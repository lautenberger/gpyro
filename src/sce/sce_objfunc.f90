MODULE SCE_OBJFUNC

USE PREC
            
IMPLICIT NONE
      
CONTAINS

! *****************************************************************************
REAL(EB) function SCE_OBJECTIVE_FUNCTION(nopt,x)
! *****************************************************************************

USE GPYRO_VARS
USE GA_VARS
USE GPYRO_FUNCS
USE GA_SUBS

IMPLICIT NONE

integer, intent(in) :: nopt
REAL(EB), dimension(1:nopt), intent(in) :: x

INTEGER :: ICASE,IPHI,ICOUNT,IT,ILOC
REAL(EB) :: TTRY,YTRY,YEXP,NUMER,DENOM
REAL(EB), ALLOCATABLE, DIMENSION(:) :: SCE_FIT

ALLOCATE(SCE_FIT(0:GA%NPHI)); SCE_FIT(:) = 0D0

CALL GET_GENES(X(:),NOPT)

DO ICASE = 1, GPG%NCASES
   CALL GPYRO_DRIVER(2,ICASE,1)

   DO IPHI = 1, GA%NPHI
      IF (GA%PHI(IPHI)%ICASE .NE. ICASE) CYCLE

      ICOUNT = 0            
      DO IT = 1, GA%NDT(ICASE)
         TTRY = GA%TRYDATA(IPHI)%T(IT)
         YTRY = GA%TRYDATA(IPHI)%Y(IT) 
                  
         IF (TTRY .GT. GA%PHI(IPHI)%TSTOP) CYCLE
         IF (TTRY .LT. 0D0               ) CYCLE
         ICOUNT = ICOUNT + 1
         CALL LOCATE(GA%EXPDATA(IPHI)%T,GA%EXPDATA(IPHI)%NDT,TTRY,ILOC)
         ILOC = MIN(ILOC,GA%EXPDATA(IPHI)%NDT)

         YEXP = LINTERP(TTRY, GA%EXPDATA(IPHI)%T(ILOC),GA%EXPDATA(IPHI)%T(ILOC+1), &
         GA%EXPDATA(IPHI)%Y(ILOC),GA%EXPDATA(IPHI)%Y(ILOC+1) )

         IF (GA%PHI(IPHI)%CTYPE .EQ. 'TMP') THEN
            YTRY = YTRY - (GPG%INITIAL_CONDITIONS(1)%TMP_INITIAL - 273.15D0)
            YEXP = YEXP - (GPG%INITIAL_CONDITIONS(1)%TMP_INITIAL - 273.15D0)
         ENDIF

         NUMER = ABS(YEXP)
         DENOM = ABS(YTRY-YEXP) + GA%PHI(IPHI)%EPS*NUMER
         IF (NUMER .EQ. 0D0 .AND. DENOM .EQ. 0D0) THEN
            SCE_FIT(IPHI) = SCE_FIT(IPHI) + & 
            GA%PHI(IPHI)%PHI * (1D0/(GA%PHI(IPHI)%PHI*GA%PHI(IPHI)%EPS))**GA%FITEXPONENT
         ELSE
            SCE_FIT(IPHI) = SCE_FIT(IPHI) + & 
            GA%PHI(IPHI)%PHI * (NUMER/DENOM)**GA%FITEXPONENT
         ENDIF
      ENDDO !IT 
                 
      IF (ICOUNT .EQ. 0) THEN
         SCE_FIT(0) = 0D0
      ELSE
         SCE_FIT(0)= SCE_FIT(0) + SCE_FIT(IPHI)/REAL(ICOUNT,EB)
      ENDIF

   ENDDO !IPHI

IF (SCE_FIT(0) .NE. SCE_FIT(0)) THEN
   CONTINUE
ENDIF

ENDDO !ICASE

SCE_OBJECTIVE_FUNCTION = -SCE_FIT(0)

DEALLOCATE(SCE_FIT)

RETURN

! *****************************************************************************
END FUNCTION
! *****************************************************************************

! *****************************************************************************
subroutine SCE_chkcst(nopt,x,bl,bu,ibound)
! *****************************************************************************
      
INTEGER :: NOPT, IBOUND, II
REAL(EB) :: X(NOPT), BL(NOPT), BU(NOPT)  

!     This subroutine check if the trial point satisfies all
!     constraints.
!
!     ibound - violation indicator
!            = -1 initial value
!            = 0  no violation
!            = 1  violation
!     nopt = number of optimizing variables
!     ii = the ii'th variable of the arrays x, bl, and bu

ibound = -1

!     Check if explicit constraints are violated

do ii=1, nopt
   if (x(ii) .lt. bl(ii) .or. x(ii) .gt. bu(ii)) go to 10
end do

if (nopt .eq. 1) go to 9

!     Check if implicit constraints are violated
!     (no implicit constraints for this function)
!
!     No constraints are violated
      
9 ibound = 0
return

!     At least one of the constraints are violated
      
10 ibound = 1
return

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

END MODULE SCE_OBJFUNC