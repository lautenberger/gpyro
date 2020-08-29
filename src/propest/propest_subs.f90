MODULE GA_SUBS

USE GPYRO_VARS 
USE GPYRO_FUNCS
USE GPYRO_INIT
USE GPYRO_PYRO, ONLY: GPYRO_PYROLYSIS, TG_DRIVER
USE SCE_FUNCS, ONLY: RANDOM_ARRAY

USE PREC
USE GA_VARS
USE GA_MPI

USE SCE_VARS

IMPLICIT NONE

CONTAINS

! *****************************************************************************
SUBROUTINE GET_OPERATING_SYSTEM
! *****************************************************************************

CHARACTER(2000) :: PATH

CALL GET_ENVIRONMENT_VARIABLE('PATH',PATH)

IF (PATH(1:1) .EQ. '/') THEN 
   OPERATING_SYSTEM = 'linux  '
   PATH_SEPARATOR   = '/'
   DELETECOMMAND    = 'rm -f  '
   COPYCOMMAND      = 'cp -f '
ELSE
   OPERATING_SYSTEM = 'windows'
   PATH_SEPARATOR   = '\'
   DELETECOMMAND    = 'del /F '
   COPYCOMMAND      = 'copy /Y '
ENDIF

! *****************************************************************************
END SUBROUTINE GET_OPERATING_SYSTEM
! *****************************************************************************
      
!******************************************************************************
SUBROUTINE GA_INIT
!******************************************************************************

REAL(EB) :: TSTOPMAX,TI,MINDIFF,DIFF,EXPONENT
INTEGER  :: ICASE,IT,IPHI,IX,IY,IZ
      
! Set misc variables
GA%NREPOPULATED    = 0
GA%NNOTCONVERGED   = 0
GA%MAXFITEVER      = -9D9

! Initialize random number generator 
SCE%iseed1 = -abs(SCE%iseed)

! Allocate main GA arrays
ALLOCATE(GA%INDIV     (1:GA%NINDIV,1:GA%NGENE));GA%INDIV       = 0D0
ALLOCATE(GA%NCOPIES   (1:GA%NINDIV           ));GA%NCOPIES     = 0
ALLOCATE(GA%CUMPROB   (1:GA%NINDIV           ));GA%CUMPROB     = 0D0 

ALLOCATE(GA%OLDPOP    (1:GA%NINDIV,1:GA%NGENE));GA%OLDPOP      = 0D0
ALLOCATE(GA%OLDFIT    (1:GA%NINDIV           ));GA%OLDFIT      = 0D0 
ALLOCATE(GA%IPARENT   (1:GA%NINDIV           ));GA%IPARENT     = 0

! Array to track whether an individual didn't converge
ALLOCATE(GA%NOTCONVERGED(1:GA%NINDIV       ));GA%NOTCONVERGED(:) = .FALSE.
ALLOCATE(GA%REPOPULATED (1:GA%NINDIV       ));GA%REPOPULATED (:) = .FALSE.

! Allocate fitness array:
ALLOCATE(GA%FIT(1:GA%NINDIV,0:GA%NPHI)); GA%FIT(:,:)=0D0
            
! Set up array TIGO (time go)
ALLOCATE(GA%NDT(1:GPG%NCASES))
TSTOPMAX = 0D0
DO ICASE = 1, GPG%NCASES
   GA%NDT(ICASE) = GPG%TSTOP(ICASE) / GPG%DT0
   IF (GPG%TSTOP(ICASE) .GT. TSTOPMAX) THEN
      TSTOPMAX = GPG%TSTOP(ICASE)
      GA%ICASEDTMAX = ICASE
   ENDIF
ENDDO
GA%NDTMAX = TSTOPMAX / GPG%DT0 + 1
      
ALLOCATE(GA%TIGO(1:GA%NDTMAX))

TI = 0D0
DO IT = 1, GA%NDTMAX
   TI = TI + GPG%DT0
   GA%TIGO(IT) = TI
ENDDO

! Allocate arrays for storing trial data:   
ALLOCATE(GA%EXPDATA(1:GA%NPHI))
ALLOCATE(GA%TRYDATA(1:GA%NPHI))
       
DO ICASE = 1, GPG%NCASES
   DO IPHI = 1, GA%NPHI
      IF (GA%PHI(IPHI)%ICASE .NE. ICASE) CYCLE
      ALLOCATE(GA%TRYDATA(IPHI)%T(1:GA%NDT(ICASE)))
      ALLOCATE(GA%TRYDATA(IPHI)%Y(1:GA%NDT(ICASE)))
   ENDDO
ENDDO

! Determine IX, IY, IZ correspondig to each xT, yT, and zT
DO IPHI = 1, GA%NPHI   
   IF (GA%PHI(IPHI)%CTYPE .NE. 'TMP') CYCLE
   ICASE = GA%PHI(IPHI)%ICASE
   G=>GPM(GPG%IMESH(ICASE))

! Start with z:
   MINDIFF = 9D9
   DO IZ = 1, G%NCELLZ 
      DIFF = ABS(G%Z(IZ) - GA%PHI(IPHI)%ZT)
      IF (DIFF .LT. MINDIFF) THEN
         MINDIFF = DIFF
         GA%PHI(IPHI)%IZ_TMP = MAX(MIN(IZ,G%NCELLZ),1)
      ENDIF 
   ENDDO !IZ

! Now on to x:
   IF (G%NCELLX .EQ. 1) THEN
      GA%PHI(IPHI)%IX_TMP = 1
   ELSE
      MINDIFF = 9D9
      DO IX = 1, G%NCELLX 
         DIFF = ABS(G%X(IX) - GA%PHI(IPHI)%XT)
         IF (DIFF .LT. MINDIFF) THEN
            MINDIFF = DIFF
            GA%PHI(IPHI)%IX_TMP = MAX(MIN(IX,G%NCELLX),1)
         ENDIF 
      ENDDO !IX
   ENDIF

! And on to y:
   IF (G%NCELLY .EQ. 1) THEN
      GA%PHI(IPHI)%IY_TMP = 1
   ELSE
      MINDIFF = 9D9
      DO IY = 1, G%NCELLY 
         DIFF = ABS(G%Y(IY) - GA%PHI(IPHI)%YT)
         IF (DIFF .LT. MINDIFF) THEN
            MINDIFF = DIFF
            GA%PHI(IPHI)%IY_TMP = MAX(MIN(IY,G%NCELLY),1)
         ENDIF 
      ENDDO !IY
   ENDIF
   
ENDDO !PHI
      
GA%MPIFITSIZE = 3 + GA%NPHI

GA%MAXPOSSIBLEFIT = 0D0
DO IPHI = 1, GA%NPHI
   IF (GA%PHI(IPHI)%PHI .GT. 1D-9) THEN
      EXPONENT = MAX(GA%PHI(IPHI)%EPS**GA%FITEXPONENT, 1D-10)
      GA%MAXPOSSIBLEFIT = GA%MAXPOSSIBLEFIT + ( GA%PHI(IPHI)%PHI**(1D0-GA%FITEXPONENT) ) / EXPONENT
   ENDIF
ENDDO

GA%TSA = GA%BSA * GA%MAXPOSSIBLEFIT / GA%ASA

END SUBROUTINE GA_INIT
!******************************************************************************

!******************************************************************************	
SUBROUTINE CHECK_GPYRO_GA
!******************************************************************************	

INTEGER :: I,IPHI,IMESH
CHARACTER(2) :: TWO
CHARACTER(300):: MESSAGE

DO I = 1, GA%NGENE
   IF (GA%VAR(I)%MINVAL .GT. GA%VAR(I)%MINVAL) THEN
      WRITE(TWO,'(I2.2)') I
      MESSAGE='Error for Gene # ' // TWO // 'MinVal cannot be greater than MaxVal'
      CALL SHUTDOWN_GRACEFULLY(MESSAGE)
   ENDIF
ENDDO

IF (TRIM(GA%OPTIMIZATION_TYPE) .NE. 'GA' .AND. TRIM(GA%OPTIMIZATION_TYPE) .NE. 'SCE' .AND. &
    TRIM(GA%OPTIMIZATION_TYPE) .NE. 'SHC') THEN
   MESSAGE='Set OPTIMIZATION_TYPE to GA or SCE'
   CALL SHUTDOWN_GRACEFULLY(MESSAGE)
ENDIF


DO IPHI = 1, GA%NPHI
   WRITE(TWO,'(I2.2)') IPHI
   IMESH = GPG%IMESH(GA%PHI(IPHI)%ICASE)  

   IF (GA%PHI(IPHI)%ICASE .GT. GPG%NCASES) THEN
      MESSAGE='Error:  for IPHI= ' // TWO // ' ICASE is greater than NCASES.'
      CALL SHUTDOWN_GRACEFULLY(MESSAGE)             
   ENDIF

   IF (GA%PHI(IPHI)%ICASE .LT. 1) THEN
      MESSAGE='Error:  for IPHI= ' // TWO // ' ICASE is less than 1.'
      CALL SHUTDOWN_GRACEFULLY(MESSAGE)             
   ENDIF

   SELECT CASE (GA%PHI(IPHI)%CTYPE)
      CASE('TMP')

          IF (GA%PHI(IPHI)%ZT .LT. 0.) THEN
             MESSAGE='Error:  for IPHI= ' // TWO // ' z_tmp is less than zero.'
             CALL SHUTDOWN_GRACEFULLY(MESSAGE)             
          ENDIF
         
          IF (GA%PHI(IPHI)%ZT .GT. GPM(IMESH)%ZDIM) THEN
             MESSAGE='Error:  for IPHI= ' // TWO // ' z_tmp is greater than zdim.'
             CALL SHUTDOWN_GRACEFULLY(MESSAGE)             
          ENDIF

          IF (GA%PHI(IPHI)%XT .LT. 0.) THEN
             MESSAGE='Error:  for IPHI= ' // TWO // ' x_tmp is less than zero.'
             CALL SHUTDOWN_GRACEFULLY(MESSAGE)             
          ENDIF
         
          IF (GA%PHI(IPHI)%XT .GT. GPM(IMESH)%XDIM) THEN
             MESSAGE='Error:  for IPHI= ' // TWO // ' x_tmp is greater than xdim.'
             CALL SHUTDOWN_GRACEFULLY(MESSAGE)             
          ENDIF

          IF (GA%PHI(IPHI)%YT .LT. 0.) THEN
             MESSAGE='Error:  for IPHI= ' // TWO // ' y_tmp is less than zero.'
             CALL SHUTDOWN_GRACEFULLY(MESSAGE)             
          ENDIF
         
          IF (GA%PHI(IPHI)%YT .GT. GPM(IMESH)%YDIM) THEN
             MESSAGE='Error:  for IPHI= ' // TWO // ' y_tmp is greater than ydim.'
             CALL SHUTDOWN_GRACEFULLY(MESSAGE)             
          ENDIF

      CASE('MLR')
         CONTINUE
!         IF (GPM(IMESH)%NCELLX .GT. 1 .OR. GPM(IMESH)%NCELLY .GT. 1) THEN
!             MESSAGE='Error:  MLR fitness metric only valid for 0D or 1D simulation.'
!             CALL SHUTDOWN_GRACEFULLY(MESSAGE)                      
!         ENDIF

      CASE('CML')
         CONTINUE
!         IF (GPM(IMESH)%NCELLX .GT. 1 .OR. GPM(IMESH)%NCELLY .GT. 1) THEN
!             MESSAGE='Error:  CML fitness metric only valid for 1D simulation.'
!             CALL SHUTDOWN_GRACEFULLY(MESSAGE)                      
!         ENDIF

      CASE('DLT')
         IF (GPM(IMESH)%NCELLX .GT. 1 .OR. GPM(IMESH)%NCELLY .GT. 1) THEN
             MESSAGE='Error:  DLT fitness metric only valid for 1D simulation.'
             CALL SHUTDOWN_GRACEFULLY(MESSAGE)                      
         ENDIF
      
      CASE DEFAULT
         MESSAGE='For IPHI = ' // TWO // ' set type to one of TMP, MLR, CML, or DLT'
         CALL SHUTDOWN_GRACEFULLY(MESSAGE)

   END SELECT
ENDDO

!******************************************************************************	
END SUBROUTINE CHECK_GPYRO_GA
!******************************************************************************	

!******************************************************************************
SUBROUTINE GPYRO_DRIVER(ITRYACT,ICASE,IINDIV)
!******************************************************************************

! ITRYACT = 1: Simulated experimental data
! ITRYACT = 2: Trial 
INTEGER, INTENT(IN) :: ITRYACT,ICASE,IINDIV
INTEGER  :: ITIME,IPHI,IZ,IX,IY,IPHIFORTMP,ILOC
REAL(EB) :: TI,TMP
LOGICAL :: FAILED_TIMESTEP
TYPE(EXP_DATA_TYPE), POINTER, DIMENSION(:) :: METRIC

CALL INIT_GPYRO_VARS(ICASE,GPG%IMESH(ICASE))

IF (ITRYACT .EQ. 1) THEN !Simulated experimental data
   METRIC => GA%EXPDATA
ELSE !Trial data
   METRIC => GA%TRYDATA
ENDIF

IPHIFORTMP = -1
DO IPHI = 1, GA%NPHI
   METRIC(IPHI)%T     (:) = -1D0
   METRIC(IPHI)%Y     (:) =  0D0
   IF (GA%PHI(IPHI)%ICASE .EQ. ICASE .AND. IPHIFORTMP .LT. 0) IPHIFORTMP=IPHI
ENDDO
      
TI        = GPG%DT0
ITIME     = 1

DO WHILE (TI .LE. GPG%TSTOP(ICASE) ) !+ GPG%DT0)
   IF (.NOT. GPG%ZEROD(ICASE)) THEN
      FAILED_TIMESTEP = .TRUE.

      DO WHILE (FAILED_TIMESTEP)
         GPG%DT = GPG%DTNEXT
         CALL GPYRO_PYROLYSIS(GPG%IMESH(ICASE),ICASE,G%NCELLZ,G%NCELLX,G%NCELLY,TI,FAILED_TIMESTEP)
!         CALL GPYRO_PYROLYSIS(ICASE,TI,FAILED_TIMESTEP,G%NCELLX,1,G%NCELLX,G%NCELLY,1,G%NCELLY,GPG%IMESH(ICASE),G%NCELLZ)
         IF (FAILED_TIMESTEP.AND.GA%KILL_NONCONVERGED_SOLNS) THEN
            GA%FIT(IINDIV,:) = GA%FITMIN
            RETURN    
         ENDIF
                         
         IF (G%DELTA(1,1) .LT. 0.0001*G%ZDIM) RETURN

      ENDDO !FAILED_TIMESTEP
   ELSE
      GPG%DTNEXT = GPG%DT0
      FAILED_TIMESTEP = .TRUE.
      DO WHILE (FAILED_TIMESTEP)
         CALL LOCATE(GA%EXPDATA(IPHIFORTMP)%T,GA%EXPDATA(IPHIFORTMP)%NDT,TI,ILOC)
         ILOC = MIN(ILOC,GA%EXPDATA(IPHIFORTMP)%NDT-1)
         TMP  = LINTERP(TI, GA%EXPDATA(IPHIFORTMP)%T(ILOC),GA%EXPDATA(IPHIFORTMP)%T(ILOC+1), GA%EXPDATA(IPHIFORTMP)%Y(ILOC),GA%EXPDATA(IPHIFORTMP)%Y(ILOC+1) )
         CALL TG_DRIVER(ICASE,TI,FAILED_TIMESTEP,TMP)
      ENDDO
   ENDIF

! Determine if it's time to store trial solution
   IF (TI .GE. GA%TIGO(ITIME) ) THEN

      DO IPHI = 1, GA%NPHI
         IF (GA%PHI(IPHI)%ICASE .NE. ICASE) CYCLE

         METRIC(IPHI)%T(ITIME) = TI

         SELECT CASE (GA%PHI(IPHI)%CTYPE)
            CASE ('TMP') !Temperature
               IF (GPG%ZEROD(ICASE)) THEN
                  METRIC(IPHI)%Y(ITIME) = G%TPN(1,1,1) - 273.15D0
               ELSE
                  IZ = GA%PHI(IPHI)%IZ_TMP
                  IX = GA%PHI(IPHI)%IX_TMP
                  IY = GA%PHI(IPHI)%IY_TMP
                  METRIC(IPHI)%Y(ITIME) = G%TPN(IZ,IX,IY) - 273.15D0
               ENDIF
            CASE ('MLR') !Mass loss rate
               METRIC(IPHI)%Y(ITIME) = INTEGRATED_MASS_LOSS_RATE(0)
            CASE ('CML') !Cumulative mass loss
               METRIC(IPHI)%Y(ITIME) = G%INITIAL_MASS - TOTAL_MASS(0)
               IF (GPG%ZEROD(ICASE)) METRIC(IPHI)%Y(ITIME) = METRIC(IPHI)%Y(ITIME) * 1000. 
            CASE ('DLT') !Thickness (delta)
               IX = GA%PHI(IPHI)%IX_TMP
               IY = GA%PHI(IPHI)%IY_TMP
               METRIC(IPHI)%Y(ITIME) = G%DELTAN(IX,IY)
         END SELECT
            
      ENDDO
            
      ITIME = ITIME + 1
   ENDIF

   GPG%DT     = GPG%DTNEXT
   TI        = TI + GPG%DT
   GPG%DTNEXT = MIN(1.001D0*GPG%DTNEXT, GPG%DT0)

ENDDO !TIME LOOP

! *****************************************************************************      
END SUBROUTINE GPYRO_DRIVER
! *****************************************************************************

! *****************************************************************************
SUBROUTINE POPULATE(IINDIV1,IINDIV2,IGENE1,IGENE2) 
! *****************************************************************************
! Generate initial population (or repopulate)
INTEGER, INTENT(IN) :: IINDIV1,IINDIV2,IGENE1,IGENE2
INTEGER IINDIV,IGENE
REAL(EB) :: RAND(1), SPAN

DO IINDIV = IINDIV1, IINDIV2
   DO IGENE = IGENE1, IGENE2
      SPAN = GA%VAR(IGENE)%MAXVAL - GA%VAR(IGENE)%MINVAL 
      CALL RANDOM_ARRAY(RAND,1)
      GA%INDIV(IINDIV,IGENE) = GA%VAR(IGENE)%MINVAL + RAND(1)*SPAN
      CONTINUE
   ENDDO
ENDDO

END SUBROUTINE POPULATE

! *****************************************************************************
SUBROUTINE QUALITY(IGEN)
! *****************************************************************************

INTEGER  :: ICASE,IT,IINDIV,IHIGH,ISTART,ISTOP,IDEST,IRUN,IPHI,ILOC,ICOUNT
LOGICAL  :: RECVDALL,RECVD(1:GA%NINDIV),DONE
REAL(EB) :: NUMER, DENOM, TTRY 
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: YTRY, YEXP

!Variables for dump:
INTEGER, INTENT(IN) :: IGEN
INTEGER :: LU,IGENE,ISKIP
CHARACTER(2) :: TWO
CHARACTER(4) :: CGEN, CINDIV
CHARACTER(300) :: FN, WRITESTR, TEMPSTR
!REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: YEXPARR,YTRYARR

ALLOCATE(YEXP(1:GA%NPHI,1:GA%NDTMAX)); YEXP(:,:) = -1D0
ALLOCATE(YTRY(1:GA%NPHI,1:GA%NDTMAX)); YTRY(:,:) = -1D0

GA%FIT(:,:)         = 0D0
GA%NOTCONVERGED(:) = .FALSE.
      
IHIGH    = 0
IRUN     = 1
RECVD(:) = .FALSE.
RECVDALL = .FALSE. 

IF (GA%MPI .AND. GA%IRANK .EQ. 0) THEN
   DO WHILE ( (.NOT. RECVDALL) )
      IF (IHIGH .EQ. 0) THEN !Start by dishing out jobs to each node
         DO IDEST = 1, GA%NPROC-1
            IHIGH = IHIGH + 1
            CALL SEND_IINDIV_TO_RUN(IDEST,IHIGH)
         ENDDO
      ENDIF

      CALL RECV_FITNESS_ANYNODE(IDEST,IINDIV)
      RECVD(IINDIV) = .TRUE.
      IHIGH = IHIGH + 1
            
      IF (IHIGH.LE.GA%NINDIV) THEN 
         CALL SEND_IINDIV_TO_RUN(IDEST,IHIGH)
      ENDIF

! Check to see if fitnesses from all individuals have been received
      DO IINDIV = 1, GA%NINDIV
         IF ( (.NOT. RECVD(IINDIV)) ) CYCLE
         IF (IINDIV .EQ. GA%NINDIV) RECVDALL = .TRUE.
      ENDDO
   ENDDO
! Now, all individuals have been calc'd, so send slave nodes signal
! that they're done with this generation
   IHIGH = 0
   DO IDEST = 1, GA%NPROC-1
      CALL SEND_IINDIV_TO_RUN(IDEST,IHIGH)
   ENDDO

ELSE

   DONE = .FALSE.
   DO WHILE (.NOT. DONE) 
      IF (GA%MPI) THEN
         CALL RECV_IINDIV_TO_RUN(IRUN)
         ISTART = IRUN
         ISTOP  = IRUN
         IF (IRUN .EQ. 0) THEN
            DONE = .TRUE.
            ISTART = 1
            ISTOP  = 0
         ENDIF
      ELSE
         ISTART = 1
         ISTOP  = GA%NINDIV
         DONE = .TRUE.
      ENDIF

      IF (IRUN .NE. 0) THEN 
         DO IINDIV=ISTART,ISTOP

            CALL GET_GENES(GA%INDIV(IINDIV,:),GA%NGENE)
                      
            DO ICASE = 1, GPG%NCASES
               CALL GPYRO_DRIVER(2,ICASE,IINDIV)
               
               DO IPHI = 1, GA%NPHI
                  IF (GA%PHI(IPHI)%ICASE .NE. ICASE) CYCLE

                  GA%FIT(IINDIV,IPHI) = 0D0

                  ICOUNT = 0            
                  DO IT = 1, GA%NDT(ICASE)
                     TTRY = GA%TRYDATA(IPHI)%T(IT)
                     YTRY(IPHI,IT) = GA%TRYDATA(IPHI)%Y(IT) 
                  
                     IF (TTRY .GT. GA%PHI(IPHI)%TSTOP) CYCLE
                     IF (TTRY .LT. 0D0               ) CYCLE
                     ICOUNT = ICOUNT + 1
                     CALL LOCATE(GA%EXPDATA(IPHI)%T,GA%EXPDATA(IPHI)%NDT,TTRY,ILOC)
                     ILOC = MIN(ILOC,GA%EXPDATA(IPHI)%NDT)

                     IF (ILOC .GE. GA%EXPDATA(IPHI)%NDT) THEN
                        YEXP(IPHI,IT) = GA%EXPDATA(IPHI)%T(GA%EXPDATA(IPHI)%NDT)
                     ELSE
                        YEXP(IPHI,IT) = LINTERP(TTRY, GA%EXPDATA(IPHI)%T(ILOC),GA%EXPDATA(IPHI)%T(ILOC+1), &
                               GA%EXPDATA(IPHI)%Y(ILOC),GA%EXPDATA(IPHI)%Y(ILOC+1) )
                     ENDIF

                     IF (GA%PHI(IPHI)%CTYPE .EQ. 'TMP') THEN
                        YTRY(IPHI,IT) = YTRY(IPHI,IT) - (GPG%INITIAL_CONDITIONS(1)%TMP_INITIAL - 273.15D0)
                        YEXP(IPHI,IT) = YEXP(IPHI,IT) - (GPG%INITIAL_CONDITIONS(1)%TMP_INITIAL - 273.15D0)
                     ENDIF

                     NUMER = ABS(YEXP(IPHI,IT))
                     DENOM = ABS(YTRY(IPHI,IT)-YEXP(IPHI,IT)) + GA%PHI(IPHI)%EPS*NUMER
                     IF (NUMER .EQ. 0D0 .AND. DENOM .EQ. 0D0) THEN
                        GA%FIT(IINDIV,IPHI) = GA%FIT(IINDIV,IPHI) + GA%PHI(IPHI)%PHI * (1D0/(GA%PHI(IPHI)%PHI*GA%PHI(IPHI)%EPS))**GA%FITEXPONENT
                     ELSE
                        GA%FIT(IINDIV,IPHI) = GA%FIT(IINDIV,IPHI) + GA%PHI(IPHI)%PHI * (NUMER/DENOM)**GA%FITEXPONENT
                     ENDIF
                        IF (GA%FIT(IINDIV,IPHI) .NE. GA%FIT(IINDIV,IPHI)) GA%FIT(IINDIV,IPHI) = -9E9 
                  ENDDO !IT 
                 
                  GA%FIT(IINDIV,0)= GA%FIT(IINDIV,0) + GA%FIT(IINDIV,IPHI)/REAL(ICOUNT,EB)
                  CONTINUE
               ENDDO !IPHI

            ENDDO !ICASE

! Begin addition to dump all results
! These lines could be consolidated with the "results_best" dump routine executed by IRANK=0
         IF (GA%DUMP_INTERMEDIATE_TRIALS) THEN

            LU=2000+IINDIV
            WRITE(CGEN  ,'(I4.4)') IGEN
            WRITE(CINDIV,'(I4.4)') IINDIV
   
            FN = TRIM(GPG%CASENAME) // '_results_' // CGEN // '_' // CINDIV // '.csv'

            OPEN(LU,FILE=FN,FORM='FORMATTED',STATUS='REPLACE')
            WRITE(LU,90) (GA%INDIV(IINDIV,IGENE),IGENE=1,GA%NGENE)
         
            WRITESTR = 't,'
            DO IPHI = 1, GA%NPHI
               WRITE(TWO,'(I2.2)') IPHI
               WRITESTR = TRIM(WRITESTR) // TWO // "_act," // TWO // "_try,"
            ENDDO
            WRITE(LU,'(A)') TRIM(WRITESTR)

            ISKIP = NINT(GPG%DTDUMP_GA/GPG%DT0)
            DO IT = 1, GA%NDTMAX, ISKIP
               WRITE(WRITESTR,92) GA%TIGO(IT)
               DO IPHI = 1, GA%NPHI
                  WRITE(TEMPSTR,93) YEXP(IPHI,IT), YTRY(IPHI,IT)
                  WRITESTR = TRIM(WRITESTR) // TRIM(TEMPSTR)
               ENDDO
               WRITE(LU,'(A)') TRIM(WRITESTR)   
            ENDDO
         
            CLOSE(LU) 

         ENDIF
         
 90   FORMAT(200(E16.9,','))
 92   FORMAT(F16.9,',')
 93   FORMAT(2(F16.9,','))

! End addition to dump all results

            IF (GA%MPI .AND. GA%IRANK .NE. 0) CALL SEND_FITNESS_ANYNODE(GA%IRANK,IINDIV)

         ENDDO !IINDIV
      ENDIF !IRUN
   ENDDO !.NOT. DONE
ENDIF !MPI ENDIF

DEALLOCATE(YEXP) 
DEALLOCATE(YTRY)

! *****************************************************************************
END SUBROUTINE QUALITY
! *****************************************************************************

! *****************************************************************************
  SUBROUTINE GET_GENES(GENE_ARR,NGENE)
! *****************************************************************************

INTEGER, INTENT(IN) :: NGENE
REAL(EB), INTENT(IN) :: GENE_ARR(1:NGENE)
INTEGER :: IGENE,I1,I2
REAL(EB) :: VAL
LOGICAL :: USE_LOG
CHARACTER(60) :: SHEET_NAME
CHARACTER(300) :: MESSAGE
        
DO IGENE = 1, GA%NGENE
   USE_LOG    = GA%VAR(IGENE)%USE_LOG
   I1         = GA%VAR(IGENE)%I1
   I2         = GA%VAR(IGENE)%I2
   SHEET_NAME = GA%VAR(IGENE)%SHEET_NAME

   VAL = GENE_ARR(IGENE) !GA%INDIV(IINDIV,IGENE)

   IF (USE_LOG) VAL = 10D0**GENE_ARR(IGENE)
   SELECT CASE(TRIM(SHEET_NAME))
         
   CASE('sprops')
      SELECT CASE(TRIM(GA%VAR(IGENE)%CTYPE))

      CASE ('K0Z')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%K0Z(I1) = -SPROP%K0Z(I2) * VAL
         ELSE
            SPROP%K0Z(I1) = VAL
         ENDIF

      CASE ('NKZ')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%NKZ(I1) = -SPROP%NKZ(I2) * VAL
         ELSE
            SPROP%NKZ(I1) = VAL
         ENDIF

      CASE ('K0X')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%K0X(I1) = -SPROP%K0X(I2) * VAL
         ELSE
            SPROP%K0X(I1)= VAL
         ENDIF

         IF (GA%ISOTROPIC_THERMAL_CONDUCTIVITY) THEN
            SPROP%K0X(I1) = SPROP%K0Z(I1)
         ENDIF

      CASE ('NKX')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%NKX(I1) = -SPROP%K0X(I2) * VAL
         ELSE
            SPROP%NKX(I1)= VAL
         ENDIF

         IF (GA%ISOTROPIC_THERMAL_CONDUCTIVITY) THEN
            SPROP%NKX(I1) = SPROP%NKZ(I1)
         ENDIF

      CASE ('K0Y')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%K0Y(I1) = -SPROP%K0Y(I2) * VAL
         ELSE
            SPROP%K0Y(I1)= VAL
         ENDIF

         IF (GA%ISOTROPIC_THERMAL_CONDUCTIVITY) THEN
            SPROP%K0Y(I1) = SPROP%K0Z(I1)
         ENDIF

      CASE ('NKY')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%NKY(I1) = -SPROP%K0Y(I2) * VAL
         ELSE
            SPROP%NKY(I1)= VAL
         ENDIF

         IF (GA%ISOTROPIC_THERMAL_CONDUCTIVITY) THEN
            SPROP%NKY(I1) = SPROP%NKZ(I1)
         ENDIF

      CASE ('R0')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%R0(I1) = -SPROP%R0(I2) * VAL
         ELSE
            SPROP%R0(I1) = VAL
         ENDIF

      CASE ('NR')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%NR(I1) = -SPROP%NR(I2) * VAL
         ELSE
            SPROP%NR(I1) = VAL
         ENDIF

      CASE ('C0')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%C0(I1) = -SPROP%C0(I2) * VAL
         ELSE
            SPROP%C0(I1) = VAL
         ENDIF

      CASE ('NC')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%NC(I1) = -SPROP%NC(I2) * VAL
         ELSE
            SPROP%NC(I1) = VAL
         ENDIF
               
      CASE ('EMIS')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%EMIS(I1) = -SPROP%EMIS(I2) * VAL
         ELSE
            SPROP%EMIS(I1) = VAL
         ENDIF

      CASE ('KAPPA')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%KAPPA(I1) = -SPROP%KAPPA(I2) * VAL
         ELSE
            SPROP%KAPPA(I1) = VAL
         ENDIF

      CASE ('TMELT')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%TMELT(I1) = -SPROP%TMELT(I2) * VAL
         ELSE
            SPROP%TMELT(I1) = VAL
         ENDIF

      CASE ('DHMELT')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%DHMELT(I1) = -SPROP%DHMELT(I2) * VAL
         ELSE
            SPROP%DHMELT(I1) = VAL
         ENDIF

      CASE ('SIGMA2MELT')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%SIGMA2MELT(I1) = -SPROP%SIGMA2MELT(I2) * VAL
         ELSE
            SPROP%SIGMA2MELT(I1) = VAL
         ENDIF

      CASE ('GAMMA')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%GAMMA(I1) = -SPROP%GAMMA(I2) * VAL
         ELSE
            SPROP%GAMMA(I1)= VAL
         ENDIF

      CASE ('PERMZ')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%PERMZ(I1) = -SPROP%PERMZ(I2) * VAL
         ELSE
            SPROP%PERMZ(I1)= VAL
         ENDIF

      CASE ('PERMX')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%PERMX(I1) = -SPROP%PERMX(I2) * VAL
         ELSE
            SPROP%PERMX(I1)= VAL
         ENDIF

         IF (GA%ISOTROPIC_PERMEABILITY) THEN
            SPROP%PERMX(I1) = SPROP%PERMZ(I1)
         ENDIF

      CASE ('PERMY')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%PERMY(I1) = -SPROP%PERMY(I2) * VAL
         ELSE
            SPROP%PERMY(I1)= VAL
         ENDIF

         IF (GA%ISOTROPIC_PERMEABILITY) THEN
            SPROP%PERMY(I1) = SPROP%PERMZ(I1)
         ENDIF

      CASE ('RS0')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%RS0(I1) = -SPROP%RS0(I2) * VAL
         ELSE
            SPROP%RS0(I1)= VAL
         ENDIF

      CASE ('PORE_DIAMETER')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            SPROP%PORE_DIAMETER(I1) = -SPROP%PORE_DIAMETER(I2) * VAL
         ELSE
            SPROP%PORE_DIAMETER(I1)= VAL
         ENDIF

      CASE DEFAULT
         MESSAGE = 'Invalid sprops variable name ' // TRIM(GA%VAR(IGENE)%CTYPE)
         CALL SHUTDOWN_GRACEFULLY(MESSAGE)

      END SELECT

   CASE('rxns')

      SELECT CASE(TRIM(GA%VAR(IGENE)%CTYPE))

      CASE ('Z')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            RXN(I1)%Z = -RXN(I2)%Z * VAL
         ELSE
            RXN(I1)%Z = VAL
         ENDIF
               
      CASE ('E')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            RXN(I1)%E = -RXN(I2)%E * VAL
         ELSE
            RXN(I1)%E = VAL
         ENDIF

      CASE ('DHS')
         IF (I2 .NE. 0) THEN !Could this create problems?
            RXN(I1)%DHS = -RXN(I2)%DHS * VAL
         ELSE
            RXN(I1)%DHS = VAL
         ENDIF

      CASE ('DHV')
         IF (I2 .NE. 0) THEN !Could this create problems?
            RXN(I1)%DHV = -RXN(I2)%DHV * VAL
         ELSE
            RXN(I1)%DHV = VAL
         ENDIF

      CASE ('CHI')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            RXN(I1)%CHI = -RXN(I2)%CHI * VAL
         ELSE
            RXN(I1)%CHI = VAL
         ENDIF

      CASE ('ORDER')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            RXN(I1)%ORDER = -RXN(I2)%ORDER * VAL
         ELSE
            RXN(I1)%ORDER = VAL
         ENDIF

      CASE ('O2ORDER')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            RXN(I1)%ORDERO2 = -RXN(I2)%ORDERO2 * VAL
         ELSE
            RXN(I1)%ORDERO2 = VAL
         ENDIF

      CASE ('M')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            RXN(I1)%M = -RXN(I2)%M*VAL
         ELSE
            RXN(I1)%M = VAL
         ENDIF

      CASE ('KCAT')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            RXN(I1)%KCAT = -RXN(I2)%KCAT*VAL
         ELSE
            RXN(I1)%KCAT = VAL
         ENDIF

      CASE ('TCRIT')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            RXN(I1)%TCRIT = -RXN(I2)%TCRIT*VAL
         ELSE
            RXN(I1)%TCRIT = VAL
         ENDIF

      CASE DEFAULT
         MESSAGE='Invalid rxns variable name ' // TRIM(GA%VAR(IGENE)%CTYPE)
         CALL SHUTDOWN_GRACEFULLY(MESSAGE)

      END SELECT

   CASE('hgrxns')

      SELECT CASE(TRIM(GA%VAR(IGENE)%CTYPE))

      CASE ('P')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            HGRXN(I1)%P = -HGRXN(I2)%P * VAL
         ELSE
            HGRXN(I1)%P = VAL
         ENDIF
               
      CASE ('Q')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            HGRXN(I1)%Q = -HGRXN(I2)%Q * VAL
         ELSE
            HGRXN(I1)%Q = VAL
         ENDIF

      CASE ('B')
         IF (I2 .NE. 0) THEN 
            HGRXN(I1)%B = -HGRXN(I2)%B * VAL
         ELSE
            HGRXN(I1)%B = VAL
         ENDIF

      CASE ('Z')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            HGRXN(I1)%Z = -HGRXN(I2)%Z * VAL
         ELSE
            HGRXN(I1)%Z = VAL
         ENDIF

      CASE ('E')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            HGRXN(I1)%E = -HGRXN(I2)%E * VAL
         ELSE
            HGRXN(I1)%E = VAL
         ENDIF

      CASE ('DH')
         IF (VAL .LT. 0D0 .AND. I2 .NE. 0) THEN
            HGRXN(I1)%DH = -HGRXN(I2)%DH * VAL
         ELSE
            HGRXN(I1)%DH = VAL
         ENDIF

      CASE DEFAULT
         MESSAGE='Invalid hgrxns variable name ' // TRIM(GA%VAR(IGENE)%CTYPE)
         CALL SHUTDOWN_GRACEFULLY(MESSAGE)

      END SELECT

   CASE('gprops')

      SELECT CASE(TRIM(GA%VAR(IGENE)%CTYPE))

      CASE ('MW')
         GPROP%M(I1) = VAL

      CASE ('SIGMA')
         GPROP%SIGMA(I1) = VAL

      CASE ('EPSOK')
         GPROP%EPSOK(I1) = VAL

      CASE ('CPG')
         GPROP%CPG = VAL

      CASE DEFAULT
         MESSAGE='Invalid gprops variable name ' // TRIM(GA%VAR(IGENE)%CTYPE)
         CALL SHUTDOWN_GRACEFULLY(MESSAGE)

      END SELECT

   CASE('gyields')
      GPROP%YIELDS(I1,I2) = VAL
            
   CASE('hgyields')
      GPROP%HGYIELDS(I1,I2) = VAL
      
   CASE DEFAULT
      MESSAGE = 'Invalid sheet name: ' // TRIM(SHEET_NAME)
      CALL SHUTDOWN_GRACEFULLY(MESSAGE)
   END SELECT

ENDDO

END SUBROUTINE GET_GENES
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SIMULATED_ANNEALING(IGEN)
! *****************************************************************************
INTEGER, INTENT(IN) :: IGEN
REAL(EB) :: FITNEW,FITOLD,PACCEPTANCE,EXPTERM,RVEC(1),AVERAGE
INTEGER :: IINDIV,IPARENT1,IPARENT2,NLOWERSOLNS
LOGICAL :: EVEN
IF (IGEN .EQ. 1) RETURN

GA%TSA = GA%ASA * GA%TSA !Set simulated annealing temperature

AVERAGE     = 0D0
NLOWERSOLNS = 0 
DO IINDIV = 1 , GA%NINDIV
   EVEN = .FALSE.; IF (MOD(IINDIV,2) .EQ. 0) EVEN = .TRUE. 

   IF (EVEN) THEN
      IPARENT1 = GA%IPARENT(IINDIV - 1)
      IPARENT2 = GA%IPARENT(IINDIV    )
   ELSE
      IPARENT1 = GA%IPARENT(IINDIV    )
      IPARENT2 = GA%IPARENT(IINDIV + 1)
   ENDIF
   
   FITNEW = GA%FIT(IINDIV,0)
   
   IF (FITNEW.NE.FITNEW .OR. FITNEW.EQ.GPG%POSINF .OR. FITNEW.EQ.GPG%NEGINF) THEN
      FITNEW = 0D0
   ENDIF 

   IF (EVEN) THEN
      FITOLD = GA%OLDFIT(IPARENT2)
   ELSE
      FITOLD = GA%OLDFIT(IPARENT1)   
   ENDIF
   
   IF (FITOLD.NE.FITOLD .OR. FITOLD.EQ.GPG%POSINF .OR. FITOLD.EQ.GPG%NEGINF) THEN
      FITOLD = 0D0
   ENDIF 

   IF (FITNEW .GE. FITOLD) THEN
      PACCEPTANCE = 1.0
   ELSE
      IF (GA%BSA .GT. 10.0) THEN
         PACCEPTANCE = 1.0
      ELSE
         NLOWERSOLNS = NLOWERSOLNS + 1
         EXPTERM = (FITNEW - FITOLD) / GA%TSA            
         PACCEPTANCE = EXP(EXPTERM)
         AVERAGE = AVERAGE + PACCEPTANCE
      ENDIF
   ENDIF
   
   CALL RANDOM_ARRAY(RVEC,1)
   
   IF (RVEC(1) .GT. PACCEPTANCE) THEN !Revert to old solution
      GA%FIT(IINDIV,0) = FITOLD
      IF (EVEN) THEN
         GA%INDIV (IINDIV,1:GA%NGENE) = GA%OLDPOP(IPARENT2,1:GA%NGENE)
       ELSE
         GA%INDIV (IINDIV,1:GA%NGENE) = GA%OLDPOP(IPARENT1,1:GA%NGENE)      
      ENDIF
   ENDIF
   
ENDDO

!WRITE(*,*) 'Number of lower solutions: ', NLOWERSOLNS
!IF (NLOWERSOLNS .GT. 0) WRITE(*,*) 'Average Pacceptance: ', AVERAGE / REAL(NLOWERSOLNS,EB)

! *****************************************************************************
END SUBROUTINE SIMULATED_ANNEALING
! *****************************************************************************


! *****************************************************************************
SUBROUTINE RANK_INDIVIDUALS
! *****************************************************************************

REAL(EB) :: MAXFIT
INTEGER :: I,J,IINDIV
LOGICAL :: SELECTED(1:GA%NINDIV)
REAL(EB), DIMENSION(1:GA%NINDIV,1:GA%NGENE ) :: INDIVTEMP
REAL(EB), DIMENSION(1:GA%NINDIV,0:GA%NPHI  ) :: FITTEMP

! Cycle through individuals and rank them from highest to lowest
SELECTED(:) = .FALSE.
IINDIV = 0
DO I = 1, GA%NINDIV 
   MAXFIT = -9D9
   DO J = 1, GA%NINDIV
      IF (SELECTED(J)) CYCLE
      IF (GA%FIT(J,0) .GE. MAXFIT) THEN
         MAXFIT = GA%FIT(J,0)
         IINDIV = J
      ENDIF
   ENDDO

! If all fitnesses are very bad, no individual will be selected,
! so do the following to prevent exceeding array bounds
   J = 1
   DO WHILE (IINDIV .EQ. 0)
      IF ( (.NOT. SELECTED(J)) ) IINDIV = J
      J = J + 1
   ENDDO

   SELECTED(IINDIV) = .TRUE.

! Copy to temporary arrays:	    
   INDIVTEMP(I,1:GA%NGENE) = GA%INDIV (IINDIV,1:GA%NGENE)
   FITTEMP  (I,0:GA%NPHI ) = GA%FIT   (IINDIV,0:GA%NPHI )
ENDDO

! Now copy ranked population and fitnesses back 
! into original data structures
GA%INDIV(1:GA%NINDIV,1:GA%NGENE)=INDIVTEMP(1:GA%NINDIV,1:GA%NGENE)
GA%FIT  (1:GA%NINDIV,0:GA%NPHI )=FITTEMP  (1:GA%NINDIV,0:GA%NPHI )
     
! Determine which ones didn't converge and which ones need to be repopulated
DO IINDIV = 1, GA%NINDIV
   IF (GA%FIT(IINDIV,0) .LE. GA%FITMIN) THEN
      GA%REPOPULATED(IINDIV) = .TRUE.
      GA%NREPOPULATED        = GA%NREPOPULATED + 1
   ENDIF         
ENDDO 
    
END SUBROUTINE RANK_INDIVIDUALS
! *****************************************************************************
      
! *****************************************************************************
SUBROUTINE REPRODUCTION
! *****************************************************************************
USE GPYRO_FUNCS

INTEGER :: I1,I2,ILOC,IPAR,IINDIV,IGENE,NPAR,ISPLIT
REAL(EB) :: MINVAL,MAXVAL,SUMFIT,PMUT,SPAN,SUBTRACTOFF,RAND(1)
REAL(EB), DIMENSION(1:GA%NINDIV,1:GA%NGENE) :: INDIVTEMP
LOGICAL , DIMENSION(1:GA%NINDIV) :: REPRODUCED
REAL(EB), DIMENSION(1:GA%NINDIV) :: PROB,TEMPFIT
REAL(EB), DIMENSION(1:MAX(GA%NINDIV,GA%NGENE)) :: RVEC
      
GA%NCOPIES(:)   = 0
GA%CUMPROB(:)   = 0D0
GA%IPARENT(:)   = 1
INDIVTEMP (:,:) = 0D0

! If fitness is very bad, just populate a new individual
IF (GA%BRUTE_FORCE) THEN
   GA%REPOPULATED(:) = .TRUE. 
   GA%NREPOPULATED   = GA%NINDIV
ENDIF
      
DO IINDIV = 1, GA%NINDIV
   IF (GA%REPOPULATED(IINDIV)) THEN
      CALL POPULATE(IINDIV,IINDIV,1,GA%NGENE)
      GA%FIT(IINDIV,0) = GA%FITMIN
   ENDIF
ENDDO

! Calculate each individual's probability of reproduction
SUBTRACTOFF = 0D0
DO IINDIV = 1, GA%NINDIV
   IF (GA%FIT(IINDIV,0) .LT. SUBTRACTOFF) THEN
      SUBTRACTOFF =  MAX(GA%FIT(IINDIV,0),GA%FITMIN)
   ENDIF 
ENDDO
TEMPFIT(:)        = GA%FIT(:,0) - SUBTRACTOFF

SUMFIT            = SUM(TEMPFIT(1:GA%NINDIV))
PROB(1:GA%NINDIV) = TEMPFIT(1:GA%NINDIV) / SUMFIT

DO IINDIV = GA%NINDIV-1,1,-1
   GA%CUMPROB(IINDIV) = GA%CUMPROB(IINDIV+1) + PROB(IINDIV)
ENDDO

! Determine parents for offspring
! Only let parents with fitness better than GA%FITMIN reproduce
NPAR = GA%NINDIV - GA%NREPOPULATED
IF (MOD(NPAR,2) .EQ. 1) THEN !Round down to get even number
   CALL POPULATE(NPAR,NPAR,1,GA%NGENE)
   NPAR = NPAR - 1 
ENDIF
         
! Copy current population's genes to intermediate population.
! The genes for low IINDIV will be replaced with combinations of parents
! where the fitness is better than fitmin
INDIVTEMP(1:GA%NINDIV,1:GA%NGENE)=GA%INDIV(1:GA%NINDIV,1:GA%NGENE)

GA%OLDPOP(1:GA%NINDIV,1:GA%NGENE)=GA%INDIV(1:GA%NINDIV,1:GA%NGENE) !Only for simulated annealing
GA%OLDFIT(1:GA%NINDIV)           =GA%FIT  (1:GA%NINDIV,0)          !Only for simulated annealing
      
REPRODUCED(:) = .FALSE.
DO IPAR = 1, NPAR
   CALL RANDOM_ARRAY(RVEC,1)
   CALL LOCATE(GA%CUMPROB(1:NPAR),NPAR,     RVEC(1),ILOC) !IS THIS RIGHT?
   ILOC = MAX(ILOC-1,1)

   IF (GA%NCOPIES(ILOC) + 1 .LE. GA%MAXCOPIES) THEN
      GA%IPARENT(IPAR) = ILOC
   ELSE
      CALL RANDOM_ARRAY(RVEC,1)
      ILOC = NINT(RVEC(1)*REAL(NPAR,EB))
      ILOC = MAX(MIN(ILOC,NPAR),1)
      GA%IPARENT(IPAR) = ILOC
   ENDIF

   GA%NCOPIES(ILOC) = GA%NCOPIES(ILOC) + 1
   REPRODUCED(ILOC) = .TRUE.
ENDDO

! Now combine parents to make offspring
DO IINDIV = 1, NPAR, 2
   I1 = GA%IPARENT(IINDIV    )
   I2 = GA%IPARENT(IINDIV + 1)
   
   CALL RANDOM_ARRAY(RAND,1)
      
   IF (RAND(1) .LT. GA%WHOLEGENEFRAC) THEN !Crossover
      CALL RANDOM_ARRAY(RAND,1)
      ISPLIT = NINT(RAND(1)*GA%NGENE)
      ISPLIT = MAX(1       ,ISPLIT)
      ISPLIT = MIN(GA%NGENE,ISPLIT)
      RVEC (1      : ISPLIT  ) = 1D0
      RVEC (ISPLIT : GA%NGENE) = 0D0
   ELSE !Linear combination
      CALL RANDOM_ARRAY(RVEC,GA%NGENE)
   ENDIF

! If a gene is "paired" to another, set its rvec value accordingly         
   DO IGENE = 1, GA%NGENE
      IF (GA%VAR(IGENE)%IVARPAIR .NE. 0) THEN
         RVEC(IGENE) = RVEC(GA%VAR(IGENE)%IVARPAIR)
      ENDIF
   ENDDO
                  
   INDIVTEMP(IINDIV,  1:GA%NGENE) = RVEC(1:GA%NGENE) * GA%INDIV(I1,1:GA%NGENE) + & 
   (1D0-RVEC(1:GA%NGENE)) * GA%INDIV(I2,1:GA%NGENE)

   INDIVTEMP(IINDIV+1,1:GA%NGENE) = RVEC(1:GA%NGENE) * GA%INDIV(I2,1:GA%NGENE) + &
   (1D0-RVEC(1:GA%NGENE)) * GA%INDIV(I1,1:GA%NGENE)
ENDDO

! Mutate
DO IINDIV = 1, GA%NINDIV
   IF (GA%REPOPULATED(IINDIV)) CYCLE
   DO IGENE=1,GA%NGENE
      PMUT = GA%VAR(IGENE)%PMUT
      CALL RANDOM_ARRAY(RVEC,3)
      IF (RVEC(1) .LT. PMUT) THEN
         SPAN = GA%VAR(IGENE)%MAXVAL-GA%VAR(IGENE)%MINVAL
         INDIVTEMP(IINDIV,IGENE)=INDIVTEMP(IINDIV,IGENE) + (RVEC(3)-0.5D0)*GA%VAR(IGENE)%VMUTMAX*SPAN
      ENDIF
   ENDDO
ENDDO

! Check bounds
DO IGENE = 1, GA%NGENE
   MINVAL = GA%VAR(IGENE)%MINVAL 
   MAXVAL = GA%VAR(IGENE)%MAXVAL
   WHERE(INDIVTEMP(:,IGENE) .LT. MINVAL) INDIVTEMP(:,IGENE) = MINVAL
   WHERE(INDIVTEMP(:,IGENE) .GT. MAXVAL) INDIVTEMP(:,IGENE) = MAXVAL
ENDDO

! Replace parents with offspring
GA%INDIV (1:GA%NINDIV,1:GA%NGENE) = INDIVTEMP(1:GA%NINDIV,1:GA%NGENE)

! Normalize mass fractions
CALL NORMALIZE_MASS_FRACTIONS
      
! Reset counters and misc variables	
GA%REPOPULATED (:) = .FALSE.
GA%NREPOPULATED    = 0
GA%NOTCONVERGED(:) = .FALSE.
GA%NNOTCONVERGED   = 0

END SUBROUTINE REPRODUCTION
! *****************************************************************************	
      
! *****************************************************************************	
SUBROUTINE NORMALIZE_MASS_FRACTIONS
! *****************************************************************************	
INTEGER :: IINDIV,IGENE,J,I1,I2
CHARACTER(60) :: SHEET_NAME
REAL(EB) :: SUMGYIELDS(1:20)

! If we're guessing initial mass fractions, start by normalizing
! them to 1.  Assume USE_LOG = .FALSE. and that all mass fractions
! being guessed by GA add to 1.  This isn't necessarily true for the case
! of say a composite with a known glass mass fraction, so this will
! have to be modified to be able to do that. 
!
!     
! First, go through and add up mass fractions in each layer. This
! happens when J=1.  Next (when J = 2), normalize each mass fraction
! by the sum in that layer.
!DO IINDIV = 1, GA%NINDIV
!   GA%SUMYI0(:) = 0D0
!   DO J = 1, 2
!      DO IGENE = 1, GA%NGENE
!         SHEET_NAME = GA%VAR(IGENE)%SHEET_NAME
!
!         IF (TRIM(SHEET_NAME) .NE. 'layers') CYCLE
!         IF (TRIM(GA%VAR(IGENE)%CTYPE) .NE. 'YI0') CYCLE
!         I1  = GA%VAR(IGENE)%I1
!         IF (J .EQ. 1) THEN 
!            GA%SUMYI0(I1) = GA%SUMYI0(I1) + GA%INDIV(IINDIV,IGENE)
!         ELSE
!            GA%INDIV(IINDIV,IGENE) = GA%INDIV(IINDIV,IGENE) / GA%SUMYI0(I1)
!         ENDIF
!      ENDDO
!   ENDDO
!ENDDO
      
! Now do the same for the condensed-phase species yields      
DO IINDIV = 1, GA%NINDIV
   SUMGYIELDS(:) = 0D0
   DO J = 1, 2
      DO IGENE = 1, GA%NGENE
         SHEET_NAME = GA%VAR(IGENE)%SHEET_NAME

         IF (TRIM(SHEET_NAME) .NE. 'gyields') CYCLE
         I2  = GA%VAR(IGENE)%I2
         IF (J .EQ. 1) THEN 
            SUMGYIELDS(I2)= SUMGYIELDS(I2)+GA%INDIV(IINDIV,IGENE)
         ELSE
            GA%INDIV(IINDIV,IGENE) = GA%INDIV(IINDIV,IGENE) / SUMGYIELDS(I2)
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Now do the same for the homogeneous gas-phase species yields
DO IINDIV = 1, GA%NINDIV
   SUMGYIELDS(:) = 0D0
   DO J = 1, 2
      DO IGENE = 1, GA%NGENE
         SHEET_NAME = GA%VAR(IGENE)%SHEET_NAME

         IF (TRIM(SHEET_NAME) .NE. 'hgyields') CYCLE

         I1  = GA%VAR(IGENE)%I1 !IGSPEC
         I2  = GA%VAR(IGENE)%I2 !IHGRXN
               
         IF (I1 .EQ. HGRXN(I2)%IREACTANT1) THEN 
            GA%INDIV(IINDIV,IGENE) = -1D0
            CYCLE
         ENDIF
               
         IF (J .EQ. 1) THEN 
            SUMGYIELDS(I2)= SUMGYIELDS(I2)+GA%INDIV(IINDIV,IGENE)
         ELSE
            GA%INDIV(IINDIV,IGENE) = GA%INDIV(IINDIV,IGENE) / SUMGYIELDS(I2)
         ENDIF
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE NORMALIZE_MASS_FRACTIONS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SHUTDOWN_GRACEFULLY(MESSAGE)
! *****************************************************************************

CHARACTER(300), INTENT(IN) :: MESSAGE
INTEGER :: IERR,LU
LOGICAL :: LOPEN

WRITE(*,*) TRIM(MESSAGE)

IF ((GA%MPI .AND. GA%IRANK .EQ. 0) .OR. (.NOT. GA%MPI)) THEN
   DO LU = 10, 900
      INQUIRE(UNIT=LU,OPENED=LOPEN)
      IF (LOPEN) CLOSE(LU) 
   ENDDO
ENDIF

IF (GA%MPI) CALL MPI_FINALIZE(IERR)

STOP

! *****************************************************************************
END SUBROUTINE SHUTDOWN_GRACEFULLY
! *****************************************************************************

! *****************************************************************************
END MODULE GA_SUBS
! *****************************************************************************
