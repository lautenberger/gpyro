MODULE GA_IO

USE GPYRO_VARS
USE GPYRO_FUNCS
USE GA_VARS
USE GA_SUBS

IMPLICIT NONE

! Logical units for i/o:	
INTEGER, PARAMETER :: LU = 800, LUDUMP = 110, LUGASUM = 120, LUCUMPROB = 810, LUCOPIES = 820

CONTAINS

!******************************************************************************
SUBROUTINE RESTART(IRESTART)
!******************************************************************************
INTEGER, INTENT(IN) :: IRESTART
CHARACTER(30) :: FN
CHARACTER(300) :: MESSAGE
INTEGER :: IOS

FN = 'restart.restart'
IF (IRESTART .EQ. 1) THEN !Dump restart file
   OPEN (LU,FILE=FN,FORM='UNFORMATTED',STATUS='REPLACE')
   WRITE(LU) GA%INDIV
   CLOSE(LU)
ELSE !Read restart file
   OPEN (LU,FILE=FN,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=IOS)
   IF (IOS .GT. 0) THEN
      MESSAGE='Problem opening restart file.'
      CALL SHUTDOWN_GRACEFULLY(MESSAGE)
   ENDIF
   READ (LU) GA%INDIV
   CLOSE(LU)
ENDIF

!******************************************************************************
END SUBROUTINE RESTART
!******************************************************************************

!******************************************************************************
SUBROUTINE GET_EXPERIMENTAL_DATA
!******************************************************************************

INTEGER  :: IPHI, IPHI2, ICASE, IT, NDT, IOS
CHARACTER(2) :: TWO1,TWO2
CHARACTER(30) :: FN
CHARACTER(300) :: MESSAGE

IF (GA%SIMULATED_EXPERIMENTAL_DATA) THEN

   DO IPHI = 1, GA%NPHI
      ICASE = GA%PHI(IPHI)%ICASE
      ALLOCATE(GA%EXPDATA(IPHI)%T     (1:GA%NDT(ICASE)))
      ALLOCATE(GA%EXPDATA(IPHI)%Y     (1:GA%NDT(ICASE)))
   ENDDO
         
   DO ICASE = 1, GPG%NCASES
      CALL GPYRO_DRIVER(1,ICASE,1)
      DO IPHI = 1, GA%NPHI
         IF (GA%PHI(IPHI)%ICASE .NE. ICASE) CYCLE
         GA%EXPDATA(IPHI)%NDT = GA%NDT(ICASE)
      ENDDO
   ENDDO
         
ELSE !Read experimental data files
   DO IPHI = 1, GA%NPHI
      IF ( (.NOT. GA%SPECIFY_CML_DIRECTLY) .AND. GA%PHI(IPHI)%CTYPE .EQ. 'CML') CYCLE

      ICASE = GA%PHI(IPHI)%ICASE
      WRITE (TWO1,98) IPHI
      WRITE (TWO2,98) ICASE

      FN = "exp_" // TWO1 // "_" // TWO2 // ".csv"
      
      OPEN(LU,FILE=FN,FORM='FORMATTED',STATUS='OLD',IOSTAT=IOS)
      IF (IOS .GT. 0) THEN
         MESSAGE = 'Problem opening experimental data file ' // TRIM(FN)
         CALL SHUTDOWN_GRACEFULLY(MESSAGE)
      ENDIF
      READ(LU,*) NDT
      GA%EXPDATA(IPHI)%NDT = NDT

      ALLOCATE(GA%EXPDATA(IPHI)%T(1:NDT))
      ALLOCATE(GA%EXPDATA(IPHI)%Y(1:NDT))
        
      DO IT = 1, NDT
         READ(LU,*) GA%EXPDATA(IPHI)%T(IT), GA%EXPDATA(IPHI)%Y(IT)
      ENDDO
         
      CLOSE (LU)
   ENDDO
         
! Calculate cumulative mass loss if necessary
   DO IPHI = 1, GA%NPHI
      IF (GA%PHI(IPHI)%CTYPE .NE. 'CML' .OR. GA%SPECIFY_CML_DIRECTLY) CYCLE
      ICASE = GA%PHI(IPHI)%ICASE

      DO IPHI2 = 1, GA%NPHI
         IF (GA%PHI(IPHI2)%ICASE .NE. ICASE) CYCLE
         IF (GA%PHI(IPHI2)%CTYPE .NE. 'MLR') CYCLE

         NDT = GA%EXPDATA(IPHI2)%NDT
         GA%EXPDATA(IPHI)%NDT = NDT
         ALLOCATE(GA%EXPDATA(IPHI)%T(1:NDT))
         ALLOCATE(GA%EXPDATA(IPHI)%Y(1:NDT))
              
         GA%EXPDATA(IPHI)%T(:) = 0D0
         GA%EXPDATA(IPHI)%Y(:) = 0D0

         GA%EXPDATA(IPHI)%T(1) = GA%EXPDATA(IPHI2)%T(1)
         DO IT = 2, NDT
            GA%EXPDATA(IPHI)%T(IT) = GA%EXPDATA(IPHI2)%T(IT)
            IF (GPG%ZEROD(ICASE)) THEN
               GA%EXPDATA(IPHI)%Y(IT) = GA%EXPDATA(IPHI)%Y(IT-1) + (GA%EXPDATA(IPHI2)%T(IT)-GA%EXPDATA(IPHI2)%T(IT-1)) * GA%EXPDATA(IPHI2)%Y(IT) * 1000.
            ELSE
               GA%EXPDATA(IPHI)%Y(IT) = GA%EXPDATA(IPHI)%Y(IT-1) + (GA%EXPDATA(IPHI2)%T(IT)-GA%EXPDATA(IPHI2)%T(IT-1)) * GA%EXPDATA(IPHI2)%Y(IT)
            ENDIF
         ENDDO !IT
      ENDDO !IPHI2
   ENDDO !IPHI

   CONTINUE
                  
ENDIF
            
 98  FORMAT(I2.2)

!******************************************************************************
END SUBROUTINE GET_EXPERIMENTAL_DATA
!******************************************************************************

!******************************************************************************
SUBROUTINE READ_GA_INPUT_FILES
!******************************************************************************

INTEGER :: I
CHARACTER(300) :: MESSAGE

! Variables read in from namelist group GA_GENINPUT:
INTEGER :: NINDIV,NGEN,MAXCOPIES
LOGICAL :: SIMULATED_EXPERIMENTAL_DATA,RESTART,BRUTE_FORCE,KILL_NONCONVERGED_SOLNS,ISOTROPIC_THERMAL_CONDUCTIVITY, &
           ISOTROPIC_PERMEABILITY, DUMP_INTERMEDIATE_TRIALS, DUMP_ALL_RESULTS_BEST, SPECIFY_CML_DIRECTLY
REAL(EB) :: FITMIN,FITCLIP,FITEXPONENT,WHOLEGENEFRAC,ASA,BSA
CHARACTER(30) :: OPTIMIZATION_TYPE
INTEGER :: ISEED, MINGS, NSPL, NOPT, NPS, NPG, NGS, KSTOP, MAXN
REAL(EB) :: PCENTO
REAL(EB) :: BL(1:500), BU(1:500) 

! Variables read in from namelist group GA_VARS:
INTEGER :: IOS, NGENE
INTEGER, DIMENSION (1:100) :: I1,I2,IVARPAIR
CHARACTER(60), DIMENSION(1:100) :: SHEET_NAME,VAR_TYPE
REAL(EB), DIMENSION(1:100) :: MINVAL,MAXVAL,PMUT,VMUTMAX
LOGICAL, DIMENSION (1:100) :: USE_LOG

! Variables read in from namelist group GA_PHI:
INTEGER :: NPHI
INTEGER, DIMENSION (1:100) :: ICASE
CHARACTER(3), DIMENSION (1:100) :: CTYPE
REAL(EB), DIMENSION (1:100) :: XT, YT, ZT, TSTOP, PHI, EPS

NAMELIST /GA_GENINPUT/ NINDIV,NGEN,MAXCOPIES,SIMULATED_EXPERIMENTAL_DATA,RESTART, &
                       FITCLIP,FITMIN,FITEXPONENT,WHOLEGENEFRAC,BRUTE_FORCE,KILL_NONCONVERGED_SOLNS, &
                       ASA,BSA,OPTIMIZATION_TYPE, &
                       MAXN, KSTOP, PCENTO, NGS, ISEED, NPG, &
                       NPS, NSPL, MINGS, BL, BU, NOPT, ISOTROPIC_THERMAL_CONDUCTIVITY, ISOTROPIC_PERMEABILITY, &
                       DUMP_INTERMEDIATE_TRIALS, DUMP_ALL_RESULTS_BEST, SPECIFY_CML_DIRECTLY

NAMELIST /GA_VARS/ NGENE,SHEET_NAME,I1,I2,VAR_TYPE,MINVAL,MAXVAL,USE_LOG,PMUT,VMUTMAX,IVARPAIR

NAMELIST /GA_PHI/ NPHI,ICASE,CTYPE,XT,YT,ZT,TSTOP,PHI,EPS

CALL SET_GPYRO_PROPEST_DEFAULTS
REWIND (LUINPUT)
READ (LUINPUT,NML=GA_GENINPUT,IOSTAT=IOS)
IF (IOS .GT. 0) THEN 
   MESSAGE='Error: Problem with namelist group &GA_GENINPUT.'
   CALL SHUTDOWN_GRACEFULLY(MESSAGE)
ENDIF

GA%NINDIV                         = NINDIV
GA%NGEN                           = NGEN
GA%MAXCOPIES                      = MAXCOPIES
GA%SIMULATED_EXPERIMENTAL_DATA    = SIMULATED_EXPERIMENTAL_DATA
GA%RESTART                        = RESTART
GA%FITMIN                         = FITMIN
GA%FITCLIP                        = FITCLIP
GA%FITEXPONENT                    = FITEXPONENT
GA%WHOLEGENEFRAC                  = WHOLEGENEFRAC
GA%BRUTE_FORCE                    = BRUTE_FORCE
GA%KILL_NONCONVERGED_SOLNS        = KILL_NONCONVERGED_SOLNS
GA%ASA                            = ASA
GA%BSA                            = BSA
GA%OPTIMIZATION_TYPE              = OPTIMIZATION_TYPE
GA%ISOTROPIC_THERMAL_CONDUCTIVITY = ISOTROPIC_THERMAL_CONDUCTIVITY
GA%ISOTROPIC_PERMEABILITY         = ISOTROPIC_PERMEABILITY
GA%DUMP_INTERMEDIATE_TRIALS       = DUMP_INTERMEDIATE_TRIALS
GA%DUMP_ALL_RESULTS_BEST          = DUMP_ALL_RESULTS_BEST
GA%SPECIFY_CML_DIRECTLY           = SPECIFY_CML_DIRECTLY

ALLOCATE ( SCE%BL     (1:NOPT)        ); SCE%BL(1:NOPT) = BL(1:NOPT)
ALLOCATE ( SCE%BU     (1:NOPT)        ); SCE%BU(1:NOPT) = BU(1:NOPT)

SCE%MAXN       = MAXN
SCE%KSTOP      = KSTOP
SCE%PCENTO     = PCENTO
SCE%NGS        = NGS
SCE%ISEED      = ISEED
SCE%NPG        = NPG
SCE%NPS        = NPS
SCE%NSPL       = NSPL
SCE%MINGS      = MINGS
SCE%BL(1:NOPT) = BL(1:NOPT)    
SCE%BU(1:NOPT) = BU(1:NOPT)
SCE%NOPT       = NOPT  

! Read data from namelist group GA_VARS 
REWIND(LUINPUT); READ (LUINPUT,NML=GA_VARS,IOSTAT=IOS)
IF (IOS .GT.  0) THEN
   MESSAGE='Error: Problem with namelist group &GA_VARS.'
   CALL SHUTDOWN_GRACEFULLY(MESSAGE)
ENDIF

GA%NGENE = NGENE

ALLOCATE(GA%VAR(1:GA%NGENE))

DO I = 1,GA%NGENE
   GA%VAR(I)%SHEET_NAME  = SHEET_NAME(I)
   GA%VAR(I)%I1          = I1(I)
   GA%VAR(I)%I2          = I2(I)
   GA%VAR(I)%CTYPE       = VAR_TYPE(I)
   GA%VAR(I)%MINVAL      = MINVAL(I)
   GA%VAR(I)%MAXVAL      = MAXVAL(I)
   GA%VAR(I)%USE_LOG     = USE_LOG(I)
   GA%VAR(I)%PMUT        = PMUT(I)
   GA%VAR(I)%VMUTMAX     = VMUTMAX(I)
   GA%VAR(I)%IVARPAIR    = IVARPAIR(I)
ENDDO

! Read namelist group GA_PHI
XT(:) = 0.
YT(:) = 0.
REWIND(LUINPUT); READ (LUINPUT,NML=GA_PHI,IOSTAT=IOS)
IF (IOS .GT. 0) THEN
   MESSAGE='Error: Problem with namelist group &GA_PHI.'
   CALL SHUTDOWN_GRACEFULLY(MESSAGE)
ENDIF

GA%NPHI = NPHI

ALLOCATE(GA%PHI   (1:GA%NPHI))

DO I = 1, GA%NPHI
   GA%PHI(I)%ICASE = ICASE(I)
   GA%PHI(I)%CTYPE = CTYPE(I)
   GA%PHI(I)%XT    = XT(I)
   GA%PHI(I)%YT    = YT(I)
   GA%PHI(I)%ZT    = ZT(I)
   GA%PHI(I)%TSTOP = TSTOP(I)
   GA%PHI(I)%PHI   = PHI(I)
   GA%PHI(I)%EPS   = EPS(I)
ENDDO

CLOSE(LUINPUT)
      
CONTAINS 
      
!******************************************************************************
SUBROUTINE SET_GPYRO_PROPEST_DEFAULTS
!******************************************************************************

NINDIV                         = 250
NGEN                           = 100
MAXCOPIES                      =   6               
SIMULATED_EXPERIMENTAL_DATA    = .FALSE.
RESTART                        = .FALSE.
FITMIN                         = 0D0
FITCLIP                        = 0D0
FITEXPONENT                    = 2D0
WHOLEGENEFRAC                  = 0.5
BRUTE_FORCE                    = .FALSE.
KILL_NONCONVERGED_SOLNS        = .TRUE. 
ASA                            = 1.0
BSA                            = 100.0
NGENE                          = 1
SHEET_NAME(:)                  = 'null'
I1(:)                          = 1
I2(:)                          = 0
VAR_TYPE(:)                    = 'null'
MINVAL(:)                      = 0D0
MAXVAL(:)                      = 0D0
USE_LOG(:)                     = .FALSE.
PMUT(:)                        = 0.10
VMUTMAX(:)                     = 0.50
IVARPAIR(:)                    = 0
NPHI                           = 1
ICASE(:)                       = 1
CTYPE(:)                       = 'null'
ZT(:)                          = 0D0
TSTOP(:)                       = 600D0
PHI(:)                         = 1D0
EPS(:)                         = 0.01
ISOTROPIC_THERMAL_CONDUCTIVITY = .FALSE.
ISOTROPIC_PERMEABILITY         = .FALSE. 
DUMP_INTERMEDIATE_TRIALS       = .FALSE.
DUMP_ALL_RESULTS_BEST          = .FALSE. 
SPECIFY_CML_DIRECTLY           = .FALSE. 

! Set default values of SCE input parameters:
MAXN     = 100000      ! max # of trials before optimization is terminated
KSTOP    =     10      ! # of shuffling loops in which criterion value must change by given %percentage before terminatino
PCENTO   =      0.0001 ! percentage by which criterion value must change in given number of shuffling loops
NGS      =      2      ! # complexes in initial population
ISEED    =     -3      ! initial random seed
NPG      =      5      ! # of points in each complex
NPS      =      3      ! # of points in sub-complex
NSPL     =      5      ! # of evolution steps allowed for each complex before complex shuffling
MINGS    =      2      ! min # of complexes required (if # of complexes is allowed to reduce)
BL(:)    =     -2.     ! lower bound on parameters
BU(:)    =      2.     ! upper bound on parameters
NOPT     =      2      ! # of parameters to be optimized

! *****************************************************************************
END SUBROUTINE SET_GPYRO_PROPEST_DEFAULTS
! *****************************************************************************

! *****************************************************************************
END SUBROUTINE READ_GA_INPUT_FILES
! *****************************************************************************

! *****************************************************************************
SUBROUTINE DUMP_GA(IGEN)
! *****************************************************************************

INTEGER, INTENT(IN) :: IGEN
INTEGER :: IGENE,IPHI,ICASE,IT,ILOC,IINDIV,ISKIP
CHARACTER(2) :: TWO
CHARACTER(4) :: FOUR
CHARACTER(300) :: FN
CHARACTER(10000) :: WRITESTR
CHARACTER(100) :: TEMPSTR
CHARACTER(300) :: SHELLSTR
LOGICAL :: LOPEN
REAL(EB) :: FITAVG
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: YEXP,YTRY

ALLOCATE(YEXP(1:GA%NPHI,1:GA%NDTMAX)); YEXP(:,:) = -1D0
ALLOCATE(YTRY(1:GA%NPHI,1:GA%NDTMAX)); YTRY(:,:) = -1D0
      
! Write summary.csv file
INQUIRE(UNIT=LUGASUM,OPENED=LOPEN)
IF (.NOT. LOPEN) THEN
   FN = TRIM(GPG%CASENAME) // '_summary.csv'
   OPEN(LUGASUM,FILE=FN,FORM='FORMATTED',STATUS='REPLACE')

   WRITESTR = 'Gen,Best fit,Avg fit,nrepop,nnotconverged,'
   DO IGENE = 1, GA%NGENE
      WRITE(TWO,'(I2.2)') IGENE
      TEMPSTR = 'var ' // TWO // ','
      WRITESTR = TRIM(WRITESTR) // TRIM(TEMPSTR)
   ENDDO
   WRITE(LUGASUM,'(A)') TRIM(WRITESTR)
ENDIF

! Calculate average fitness
FITAVG = 0D0
DO IINDIV = 1, GA%NINDIV
   FITAVG = FITAVG + GA%FIT(IINDIV,0)
ENDDO
FITAVG = FITAVG / REAL(GA%NINDIV,EB)

WRITE(WRITESTR,94) IGEN,GA%FIT(1,0),FITAVG,GA%NREPOPULATED,GA%NNOTCONVERGED
DO IGENE = 1, GA%NGENE
   WRITE(TEMPSTR,95) GA%INDIV(1,IGENE)
   WRITESTR = TRIM(WRITESTR) // TRIM(TEMPSTR)
ENDDO
WRITE(LUGASUM,'(A)') TRIM(WRITESTR)

!WRITE(LUGASUM,94) IGEN,GA%FIT(1,0),FITAVG,GA%NREPOPULATED,GA%NNOTCONVERGED
!WRITE(LUGASUM,95) (GA%INDIV(1,IGENE),IGENE=1,GA%NGENE)
!WRITE(LUGASUM,*)
  
IF (GA%FIT(1,0) .GT. GA%MAXFITEVER .OR. GA%DUMP_ALL_RESULTS_BEST) THEN

   WRITE(FOUR,'(I4.4)') IGEN

   FN = TRIM(GPG%CASENAME) // '_results_best_' // FOUR // '.csv'
   OPEN(LUDUMP,FILE=FN,FORM='FORMATTED',STATUS='REPLACE')
   WRITE(LUDUMP,90) (GA%INDIV(1,IGENE),IGENE=1,GA%NGENE)
         
   WRITESTR = 't,'
   DO IPHI = 1, GA%NPHI
      WRITE(TWO,'(I2.2)') IPHI
      WRITESTR = TRIM(WRITESTR) // TWO // "_act," // TWO // "_try,"
   ENDDO
   WRITE(LUDUMP,'(A)') TRIM(WRITESTR)

   CALL GET_GENES(GA%INDIV(1,:),GA%NGENE)

   DO ICASE = 1, GPG%NCASES
      IF (GPG%FDSMODE) CALL DUMP_FDS_INPUT_FILE(ICASE)
            
      CALL GPYRO_DRIVER(2,ICASE,1)

      DO IPHI = 1, GA%NPHI
         IF (GA%PHI(IPHI)%ICASE .NE. ICASE) CYCLE
            
         DO IT = 1, GA%NDT(ICASE)
            YTRY(IPHI,IT) = GA%TRYDATA(IPHI)%Y(IT)

            IF (GA%TIGO(IT) .GT. GA%PHI(IPHI)%TSTOP) CYCLE
                  
            CALL LOCATE(GA%EXPDATA(IPHI)%T,GA%EXPDATA(IPHI)%NDT,GA%TIGO(IT),ILOC)
            ILOC = MIN(ILOC,GA%EXPDATA(IPHI)%NDT)

            YEXP(IPHI,IT) = LINTERP(GA%TIGO(IT),GA%EXPDATA(IPHI)%T(ILOC),GA%EXPDATA(IPHI)%T(ILOC+1), &
            GA%EXPDATA(IPHI)%Y(ILOC),GA%EXPDATA(IPHI)%Y(ILOC+1) )
         ENDDO !IT
                         
      ENDDO !IPHI
   ENDDO !ICASE

   ISKIP = NINT(GPG%DTDUMP_GA/GPG%DT0)
   DO IT = 1, GA%NDTMAX, ISKIP
      WRITE(WRITESTR,92) GA%TIGO(IT)
      DO IPHI = 1, GA%NPHI
         WRITE(TEMPSTR,93) YEXP(IPHI,IT), YTRY(IPHI,IT)
         WRITESTR = TRIM(WRITESTR) // TRIM(TEMPSTR)
      ENDDO
      WRITE(LUDUMP,'(A)') TRIM(WRITESTR)   
   ENDDO
         
   CLOSE(LUDUMP) 

   IF (GA%FIT(1,0) .GT. GA%MAXFITEVER) THEN
      SHELLSTR = TRIM(COPYCOMMAND) // ' ' // TRIM(FN) // ' ' // TRIM(GPG%CASENAME) // '_results_best.csv'
      CALL SYSTEM(TRIM(SHELLSTR))
   ENDIF

   GA%MAXFITEVER = GA%FIT(1,0)

ENDIF !GA%FIT(1,0) .GT. GA%MAXFITEVER
      
DEALLOCATE(YTRY,YEXP)
      
 90   FORMAT(200(E16.9,','))
 92   FORMAT(F18.9,',')
 93   FORMAT(2(F16.9,','))
 94   FORMAT(I8,',',E16.9,',',E16.9,',',I8,',',I8,',')
 95   FORMAT(E16.9,',')

! *****************************************************************************
END SUBROUTINE DUMP_GA
! *****************************************************************************

! *****************************************************************************
SUBROUTINE DUMP_COPIES(IGEN)
! *****************************************************************************

INTEGER, INTENT(IN) :: IGEN
INTEGER :: J, IINDIV
LOGICAL :: LOPEN
CHARACTER(300) :: FN
CHARACTER(10000) :: WRITESTR
CHARACTER(1000) :: TEMPSTR

! Write number of copies of each individual
INQUIRE(UNIT=LUCOPIES,OPENED=LOPEN)
IF (.NOT. LOPEN) THEN
   FN = TRIM(GPG%CASENAME) // '_ncopies.csv'
   OPEN(LUCOPIES,FILE=FN,FORM='FORMATTED',STATUS='REPLACE')
   WRITESTR = ','
   DO J = 1, GA%NINDIV
      WRITE(TEMPSTR,97) J
      WRITESTR = TRIM(WRITESTR) // TRIM(TEMPSTR) 
   ENDDO
   WRITE(LUCOPIES,'(A)') TRIM(WRITESTR)
ENDIF

WRITE(WRITESTR,97) IGEN
DO IINDIV = 1, GA%NINDIV
   WRITE(TEMPSTR,97) GA%NCOPIES(IINDIV)
   WRITESTR = TRIM(WRITESTR) // TRIM(TEMPSTR)
ENDDO
WRITE(LUCOPIES,'(A)') TRIM(WRITESTR)

! Write cumulative probability file
INQUIRE(UNIT=LUCUMPROB,OPENED=LOPEN)
IF (.NOT. LOPEN) THEN
   FN = TRIM(GPG%CASENAME) // "_cumprob.csv"
   OPEN(UNIT=LUCUMPROB, FILE=FN, STATUS = 'REPLACE')
   WRITESTR = ','
   DO IINDIV = 1, GA%NINDIV
      WRITE(TEMPSTR,91) IINDIV
      WRITESTR = TRIM(WRITESTR) // TRIM(TEMPSTR)
   ENDDO
   WRITE(LUCUMPROB,'(A)') TRIM(WRITESTR)   
ENDIF

WRITE(WRITESTR,91) IGEN
DO IINDIV = 1, GA%NINDIV
   WRITE(TEMPSTR,92) GA%CUMPROB(IINDIV)
   WRITESTR = TRIM(WRITESTR) // TRIM(TEMPSTR)
ENDDO
WRITE(LUCUMPROB,'(A)') TRIM(WRITESTR)   

 91 FORMAT(I5,',')
 92 FORMAT(E16.9,',')
 97 FORMAT(I9,',')

!******************************************************************************
END SUBROUTINE DUMP_COPIES
! *****************************************************************************

!******************************************************************************
END MODULE GA_IO
! *****************************************************************************
