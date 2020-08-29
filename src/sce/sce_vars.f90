! *****************************************************************************
MODULE SCE_VARS
! *****************************************************************************

USE PREC

TYPE :: SCE_TYPE
!Logical units:
   INTEGER :: LUIN  =  2
   INTEGER :: LUOUT = 10
   INTEGER :: LUSUM = 11

!Namelist group entries:
   INTEGER :: ISEED, MINGS, NSPL, NOPT, NPS, NPG, NGS, KSTOP, MAXN
   REAL(EB) :: PCENTO
   REAL(EB), POINTER, DIMENSION(:) :: BL, BU

!Miscellaneous globals:
   INTEGER :: NPT
   REAL(EB) :: PCENTA 
   REAL(EB) :: FA
   REAL(EB) :: BESTF
   REAL(EB) :: WORSTF
   REAL(EB) :: GNRNG 

   INTEGER :: NLOOP
   INTEGER :: ICALL
   INTEGER :: NGS1
   INTEGER :: ISEED1
   INTEGER :: NPT1
   INTEGER :: IPCNVG

!Allocatable stuff:
   REAL(EB), POINTER, DIMENSION (:)   :: XX, BESTX, WORSTX, XF, SF, CF, XNSTD, BOUND, UNIT
   REAL(EB), POINTER, DIMENSION (:,:) :: X, S, CX
   INTEGER,  POINTER, DIMENSION (:) :: LCS

!Fixed length arrays:
   INTEGER :: IRANK_STATUS(512)

END TYPE

TYPE (SCE_TYPE), SAVE  :: SCE

! *****************************************************************************
END MODULE SCE_VARS
! *****************************************************************************
