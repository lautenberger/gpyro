MODULE SCE_IO

USE PREC

IMPLICIT NONE

CHARACTER*4 :: XNAME(1:50)

CONTAINS

! *****************************************************************************
SUBROUTINE SCE_INPUT
! *****************************************************************************
!   THIS SUBROUTINE READS AND PRINTS THE INPUT VARIABLES FOR
!   SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
!     -- Version 2.1
!
!   WRITTEN BY QINGYUN DUAN - UNIVERSITY OF ARIZONA, APRIL 1992

USE SCE_VARS

IMPLICIT NONE

!     maxn = max no. of trials allowed before optimization is terminated
!     kstop = number of shuffling loops in which the criterion value must
!         chang by the given percentage before optimization is terminated
!     pcento = percentage by which the criterion value must change in
!         given number of shuffling loops
!     ipcnvg = flag indicating whether parameter convergence is reached
!         (i.e., check if gnrng is less than 0.001)
!         = 0, parameter convergence not satisfied
!         = 1, parameter convergence satisfied

! MAXN   = max # of trials before optimization is terminated
! KSTOP  = # of shuffling loops in which criterion value must change by given %percentage before terminatino
! PCENTO = percentage by which criterion value must change in given number of shuffling loops
! NGS    = # complexes in initial population
! ISEED  = initial random seed
! NPG    = # of points in each complex
! NPS    = # of points in sub-complex
! NSPL   = # of evolution steps allowed for each complex before complex shuffling
! MINGS  = min # of complexes required (if # of complexes is allowed to reduce)
! A      = initial parameter set
! BL     = lower bound on parameters
! BU     = upper bound on parameters
! NOPT   = # of parameters to be optimized

INTEGER :: IOS

!Local Namelist group entries:
INTEGER :: ISEED, MINGS, NSPL, NOPT, NPS, NPG, NGS, KSTOP, MAXN
REAL(EB) :: PCENTO
REAL(EB) :: BL(1:500), BU(1:500) 

NAMELIST /SCE_NAMELIST/ MAXN, KSTOP, PCENTO, NGS, ISEED, NPG, &
                        NPS, NSPL, MINGS, BL, BU, NOPT

! Set default values of input parameters:
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

!  Open and read input file and namelist group &SCE_INPUT
open(unit=SCE%LUIN, file='sce_input.data', status='OLD')

READ(SCE%LUIN,NML=SCE_NAMELIST,END=100,IOSTAT=IOS)
 100  IF (IOS > 0) WRITE(*,*) 'Error: Problem with namelist group &SCE_INPUT.'
CLOSE(SCE%LUIN)

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

! *****************************************************************************
end subroutine SCE_INPUT
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SCE_DUMP
! *****************************************************************************

USE SCE_VARS

IMPLICIT NONE

INTEGER :: J
LOGICAL, SAVE :: FIRSTCALL = .TRUE.
LOGICAL :: LOPEN

IF (FIRSTCALL) THEN

   XNAME( 1) = '  X1'
   XNAME( 2) = '  X2'
   XNAME( 3) = '  X3'
   XNAME( 4) = '  X4'
   XNAME( 5) = '  X5'
   XNAME( 6) = '  X6'
   XNAME( 7) = '  X7'
   XNAME( 8) = '  X8'
   XNAME( 9) = '  X9'
   XNAME(10) = ' X10'
   XNAME(11) = ' X11'
   XNAME(12) = ' X12'
   XNAME(13) = ' X13'
   XNAME(14) = ' X14'
   XNAME(15) = ' X15'
   XNAME(16) = ' X16'
   XNAME(17) = ' X17'
   XNAME(18) = ' X18'
   XNAME(19) = ' X19'
   XNAME(20) = ' X20'

!   open (unit=SCE%LUOUT,file='sce.out',status='REPLACE')
!
!   write(SCE%LUOUT,*) '          SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION'
!   write(SCE%LUOUT,*) '          =============================================='
!   WRITE(SCE%LUOUT,*)
!   WRITE(SCE%LUOUT,*)
!
!!  PRINT SHUFFLED COMPLEX EVOLUTION OPTIMIZATION OPTIONS
!   104 write(SCE%LUOUT,910)
!   910 format(//,2x,'SCE CONTROL',5x,'MAX TRIALS',5x, &
!        'REQUIRED IMPROVEMENT',5x,'RANDOM',/,3x,'PARAMETER',8x, &
!        'ALLOWED',6x,'PERCENT',4x,'NO. LOOPS',6x,'SEED',/, &
!        2x,11(1h-),5x,10(1H-),5x,7(1h-),4x,9(1h-),5x,6(1h-))
!
!!   write(SCE%LUOUT,912) 'USER SPEC.',SCE%MAXN,SCE%pcenta,SCE%kstop,SCE%iseed
!!   912 format(3x,a10,7x,i5,10x,f3.1,9x,i2,9x,i5)
!
!   write(SCE%LUOUT,914) SCE%ngs,SCE%npg,SCE%npt,SCE%nps,SCE%nspl
!   914 format(//,18x,'SCE ALGORITHM CONTROL PARAMETERS',/,18x,32(1H=), &
!        //,2x,'NUMBER OF',5x,'POINTS PER',5x,'POINTS IN',6x,'POINTS PER', &
!        4x,'EVOL. STEPS',/,2x,'COMPLEXES',6X,'COMPLEX',6x,'INI. POPUL.', &
!        5x,'SUB-COMPLX',4x,'PER COMPLEX',/,2x,9(1h-),5x,10(1h-),4x, &
!        11(1h-),5x,10(1h-),4x,11(1h-),5x,/,2x,5(i5,10x))
!
!   if (SCE%mings .lt. SCE%ngs) then
!      reduc = 'YES '
!   else
!      reduc = 'NO  '
!   end if
!
!   write(SCE%LUOUT,915) reduc,SCE%mings,'NO  '
!   915 format(//,15x,'COMPLX NO.',5x,'MIN COMPLEX',5x,'INI. POINT',/, &
!        15x,'REDUCTION',6x,'NO. ALLOWED',6x,'INCLUDED',/, &
!        15x,10(1h-),5x,11(1h-),5x,10(1h-),/,18x,a4,6x,i8,13x,a4)
!
!   write(SCE%LUOUT,916)
!   916 format(//,8x,'INITIAL PARAMETER VALUES AND PARAMETER BOUNDS',/, &
!               8x,45(1h=),//,2x,'PARAMETER',5x,'INITIAL VALUE',5x, &
!               'LOWER BOUND',5x,'UPPER BOUND',/,2x,9(1h-),5x,13(1h-),5x, &
!               11(1h-),5x,11(1h-))
!
!   do i = 1, SCE%nopt
!      write(SCE%LUOUT,918) xname(i),0.,SCE%bl(i),SCE%bu(i)
!   918 format(4x,a4,4x,3(6x,f10.3))
!   ENDDO
!
!   write(SCE%LUOUT,*) 'BEGIN SHUFFLED COMPLEX EVOLUTION GLOBAL SEARCH'
!
!!   write(SCE%LUOUT,*) '*** PRINT THE INITIAL POINT AND ITS CRITERION VALUE ***'
!!   write(SCE%LUOUT,510) (xname(j),j=1,nopt)
!!   write(SCE%LUOUT,520) fa,(0,j=1,nopt)
!!   510 format(/,' CRITERION',12(6x,a4),/1x,60(1h-))
!!   520 format(g10.3,12f10.3)
!
!   530 format(10x,12(6x,a4))
!   540 format(10x,12f10.3)
!
!   CLOSE (SCE%LUOUT)
      
   FIRSTCALL = .FALSE. 

ENDIF

INQUIRE(UNIT=SCE%LUSUM,OPENED=LOPEN)


IF (.NOT. LOPEN) THEN

   open (unit=SCE%LUSUM,file='gpyro_sce_summary.txt',status='REPLACE')

!  PRINT THE RESULTS FOR THE INITIAL POPULATION
   write(SCE%LUSUM,610) (xname(j),j=1,SCE%nopt)
   610 format(/,1x,'LOOP',1x,'TRIALS',1x,'COMPLXS',2x,'BEST F',3x,'WORST F',3x,'PAR RNG',1x,100(6x,a4))

ENDIF

write(SCE%LUSUM,630) SCE%nloop,SCE%icall,SCE%ngs1,SCE%bestf,SCE%worstf,SCE%gnrng,(SCE%bestx(j),j=1,SCE%nopt)
630 format(i6,2x,i7,3x,i5,3g10.3,100(f10.3))
650 format(/,1x,'POPULATION AT LOOP ',i3,/,1x,22(1h-))
!660 format(15x,g10.3,20x,8(f10.3))

! *****************************************************************************
END SUBROUTINE SCE_DUMP
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SCE_DUMP_CONVERGED
! *****************************************************************************

USE SCE_VARS

IMPLICIT NONE

INTEGER :: J

write(SCE%LUSUM,820) SCE%gnrng*100.
820 format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE POPULATION',' HAS CONVERGED INTO ',/,4x,f5.2, &
                 ' PERCENT OF THE FEASIBLE SPACE ***')

!  PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
write(SCE%LUSUM,830)
830 format(//,'*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS',' CRITERION VALUE ***')

write(SCE%LUSUM,510) (xname(j),j=1,SCE%nopt)
write(SCE%LUSUM,520) SCE%bestf,(SCE%bestx(j),j=1,SCE%nopt)
   510 format(/,' CRITERION',100(6x,a4),/1x,60(1h-))
   520 format(g10.3,100f10.3)

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

END MODULE SCE_IO