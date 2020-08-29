MODULE SCE_INIT

USE PREC

IMPLICIT NONE

CONTAINS

! *****************************************************************************
SUBROUTINE SCE_INITIALIZE(GAMPI,GAIRANK)
! *****************************************************************************

USE SCE_VARS
USE SCE_UA
USE GA_VARS, ONLY: GA

IMPLICIT NONE

LOGICAL, INTENT(IN) :: GAMPI
INTEGER, INTENT(IN) :: GAIRANK

INTEGER :: IERROR,I,j,LOOP

SCE%npt = SCE%ngs * SCE%npg

! Allocate variables:
ALLOCATE ( SCE%XX     (1:SCE%NOPT)        )
ALLOCATE ( SCE%BESTX  (1:SCE%NOPT)        )
ALLOCATE ( SCE%WORSTX (1:SCE%NOPT)        )
ALLOCATE ( SCE%XF     (1: SCE%NPT)        )
ALLOCATE ( SCE%SF     (1: SCE%NPS)        )
ALLOCATE ( SCE%XNSTD  (1:SCE%NOPT)        )
ALLOCATE ( SCE%BOUND  (1:SCE%NOPT)        )
ALLOCATE ( SCE%UNIT   (1:SCE%NOPT)        )
ALLOCATE ( SCE%X      (1:SCE%NPT, 1:SCE%NOPT) )
ALLOCATE ( SCE%S      (1:SCE%NPS, 1:SCE%NOPT) )
ALLOCATE ( SCE%LCS    (1:SCE%NPS)         )

ALLOCATE ( SCE%CF     (1:SCE%NPT )        )    !SHOULD BE NPG?
ALLOCATE ( SCE%CX     (1:SCE%NPT, 1:SCE%NOPT) )

SCE%pcenta=SCE%pcento*100.

if (SCE%iseed .eq. 0) SCE%iseed = 299989

IF (GAMPI) SCE%ISEED = SCE%ISEED - 29*GAIRANK

ierror = 0

if (SCE%ngs .lt. 1) then
   write(*,*) '**ERROR** NUMBER OF COMPLEXES IN INITIAL POPULATION ', SCE%NGS, ' IS NOT A VALID CHOICE'
   ierror = ierror + 1
end if

if (SCE%mings .lt. 1 .or. SCE%mings .gt. SCE%ngs) then
   write(*,*) '**ERROR** THE MINIMUM NUMBER OF COMPLEXES ', SCE%MINGS
   write(*,*) ' IS NOT A VALID CHOICE. SET IT TO DEFAULT'
   IERROR = IERROR + 1
end if

if (SCE%npg .lt. 2) then
   write(*,*) '**ERROR** THE NUMBER OF POINTS IN A COMPLEX ', SCE%NPG
   write(*,*) 'IS NOT A VALID CHOICE'
   IERROR = IERROR + 1
end if

if (SCE%nps.lt.2 .or. SCE%nps.gt.SCE%npg) then
   write(*,*) '**ERROR** THE NUMBER OF POINTS IN A SUB-COMPLEX ', SCE%NPS 
   write(*,*) ' IS NOT A VALID CHOICE'
   IERROR = IERROR + 1
end if

if (SCE%nspl .lt. 1) then
   write(*,*) '**ERROR** THE NUMBER OF EVOLUTION STEPS TAKEN IN EACH COMPLEX BEFORE SHUFFLING ', SCE%NSPL
   write(*,*) 'IS NOT A VALID CHOICE'
   IERROR = IERROR + 1
end if

if (GA%NPROC .NE. SCE%ngs + 1) then
   write(*,*) '**ERROR** For SCE optimization, the number of processors must equal' 
   write(*,*) 'the number of complexes (NGS) plus 1.'
   write(*,*) 'For example, if NGS = 8 then pass flag -np 9 to mpirun.'   
!   IERROR = IERROR + 1 !change this back
end if

if (ierror .ge. 1) then
   write(*,*) '*** TOTAL NUMBER OF ERROR MESSAGES IS ', IERROR
   write(*,*) '*** THE OPTIMIZATION SEARCH IS NOT CONDUCTED BECAUSE OF INPUT DATA ERROR'
   stop
end if

!Set defaults (commented for now because they're in the namelist group) :
!npg = 2*nopt + 1
!nps = nopt + 1
!nspl = npg
!mings = ngs
!iniflg = 0
!iprint = 0

!  INITIALIZE VARIABLES
SCE%nloop = 0
loop  = 0

!  INITIALIZE RANDOM SEED TO A NEGATIVE INTEGER
SCE%iseed1 = -abs(SCE%iseed)

!  COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPUALTION
SCE%ngs1 = SCE%ngs
SCE%npt1 = SCE%npt

SCE%icall = 0

!  COMPUTE THE BOUND FOR PARAMETERS BEING OPTIMIZED
do j = 1, SCE%nopt
   SCE%bound(j) = SCE%bu(j) - SCE%bl(j)
   SCE%unit(j) = 1.0
end do

!  GENERATE npt RANDOM POINTS DISTRIBUTED UNIFORMLY IN THE PARAMETER
!  SPACE, AND COMPUTE THE CORRESPONDING FUNCTION VALUES

SCE%X(:,:) = 0.
IF ( (.NOT. GAMPI) .OR. (GAMPI .AND. GAIRANK .EQ. 0) ) THEN
   do i = 1, SCE%NPT1
      call SCE_getpnt(SCE%nopt,1,SCE%iseed1,SCE%xx,SCE%bl,SCE%bu,SCE%unit,SCE%bl,SCE%NOPT)
      do j = 1, SCE%nopt
         SCE%x(i,j) = SCE%xx(j)
      end do
   end do
ENDIF

! *****************************************************************************
END SUBROUTINE SCE_INITIALIZE
! *****************************************************************************

END MODULE SCE_INIT

