! *****************************************************************************
MODULE GPYRO_BC
! *****************************************************************************

USE PREC
USE GPYRO_VARS
USE GPYRO_FUNCS

IMPLICIT NONE

CONTAINS

! *****************************************************************************
SUBROUTINE GET_ALL_BOUNDARY_CONDITIONS(IMESH,T)
! *****************************************************************************

USE GPYRO_VARS

INTEGER, INTENT(IN) :: IMESH
REAL(EB), INTENT(IN) :: T
INTEGER :: IL, IH, IZ, IX, IY
REAL(EB) :: TSTART,TEND,F

CALL GET_CPU_TIME(TSTART)

G=>GPM(IMESH)

IF (IGPYRO_TYPE .NE. 2) THEN !If this isn't an FDS simulation, get user-specified boundary conditions
   DO IY = 1, G%NCELLY
   DO IX = 1, G%NCELLX
   DO IZ = 1, G%NCELLZ

      IF (G%NEEDSBCT(IZ,IX,IY)) THEN ! Do -z boundary condition
         CALL GET_BC_F(G%SURF_IDX_BCT(IZ,IX,IY),T,IL,IH,F)
         CALL INTERPOLATE_BCS(IL,IH,F,IZ,IX,IY,-3)
      ENDIF

      IF (G%NEEDSBCB(IZ,IX,IY)) THEN ! Do +z boundary condition
         CALL GET_BC_F(G%SURF_IDX_BCB(IZ,IX,IY),T,IL,IH,F)
         CALL INTERPOLATE_BCS(IL,IH,F,IZ,IX,IY, 3)
      ENDIF

! Do x boundary conditions:
      IF (G%NCELLX .GT. 1) THEN
         IF (G%NEEDSBCW(IZ,IX,IY)) THEN ! Do -x boundary condition
            CALL GET_BC_F(G%SURF_IDX_BCW(IZ,IX,IY),T,IL,IH,F)
            CALL INTERPOLATE_BCS(IL,IH,F,IZ,IX,IY,-1)
         ENDIF
        
         IF (G%NEEDSBCE(IZ,IX,IY)) THEN ! Do +x boundary condition
            CALL GET_BC_F(G%SURF_IDX_BCE(IZ,IX,IY),T,IL,IH,F)
            CALL INTERPOLATE_BCS(IL,IH,F,IZ,IX,IY, 1)
         ENDIF
      ENDIF 

! Do y boundary conditions:
      IF (G%NCELLY .GT. 1) THEN
         IF (G%NEEDSBCS(IZ,IX,IY)) THEN ! Do -y boundary condition
            CALL GET_BC_F(G%SURF_IDX_BCS(IZ,IX,IY),T,IL,IH,F)
            CALL INTERPOLATE_BCS(IL,IH,F,IZ,IX,IY,-2)
         ENDIF

         IF (G%NEEDSBCN(IZ,IX,IY)) THEN ! Do +y boundary condition
            CALL GET_BC_F(G%SURF_IDX_BCN(IZ,IX,IY),T,IL,IH,F)
            CALL INTERPOLATE_BCS(IL,IH,F,IZ,IX,IY, 2)
         ENDIF      
      ENDIF 

   ENDDO
   ENDDO
   ENDDO

ENDIF

GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)

CALL GET_CPU_TIME(TEND)
GPG%TUSED(19) = GPG%TUSED(19) + TEND - TSTART 

CONTAINS

! *****************************************************************************
SUBROUTINE GET_BC_F(SURF_IDX_TARG,T,IL,IH,F)
! *****************************************************************************

INTEGER, INTENT(IN) :: SURF_IDX_TARG
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(OUT) :: IL, IH
REAL(EB), INTENT(OUT) :: F

INTEGER :: I 

IL = -1
IH = -1
DO I = 1, GPG%NSURF_IDX
   IF (GPG%ALLBC(I)%SURF_IDX .EQ. SURF_IDX_TARG) THEN

      IF (GPG%ALLBC(I)%T .LE. T) THEN 
         IL = I
         IH = IL
         IF (I .LT. GPG%NSURF_IDX) THEN
            IF (GPG%ALLBC(I+1)%SURF_IDX .EQ. SURF_IDX_TARG) IH = IL + 1
         ENDIF
      ENDIF
      
   ENDIF
ENDDO

IF (IH .EQ. IL) THEN
   F  = 0.
ELSE
   F = (T - GPG%ALLBC(IL)%T) / (GPG%ALLBC(IH)%T - GPG%ALLBC(IL)%T)
ENDIF

! *****************************************************************************
END SUBROUTINE GET_BC_F
! *****************************************************************************

! *****************************************************************************
SUBROUTINE INTERPOLATE_BCS(IL,IH,F,IZ,IX,IY,IOR)
! *****************************************************************************

INTEGER, INTENT(IN) :: IL,IH,IZ,IX,IY,IOR
REAL(EB), INTENT(IN) :: F
INTEGER :: ISPEC
REAL(EB) :: HFIXEDLO, HFIXEDHI

! Linearly interpolate real quantities:
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%QE      = GPG%ALLBC(IL)%QE      + F * (GPG%ALLBC(IH)%QE      - GPG%ALLBC(IL)%QE      )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HC0     = GPG%ALLBC(IL)%HC      + F * (GPG%ALLBC(IH)%HC      - GPG%ALLBC(IL)%HC      )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%NHC     = GPG%ALLBC(IL)%NHC     + F * (GPG%ALLBC(IH)%NHC     - GPG%ALLBC(IL)%NHC     )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%TINF    = GPG%ALLBC(IL)%TINF    + F * (GPG%ALLBC(IH)%TINF    - GPG%ALLBC(IL)%TINF    )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%TFIXED  = GPG%ALLBC(IL)%TFIXED  + F * (GPG%ALLBC(IH)%TFIXED  - GPG%ALLBC(IL)%TFIXED  )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%PRES    = GPG%ALLBC(IL)%PRES    + F * (GPG%ALLBC(IH)%PRES    - GPG%ALLBC(IL)%PRES    )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%MFLUX   = GPG%ALLBC(IL)%MDOTPP  + F * (GPG%ALLBC(IH)%MDOTPP  - GPG%ALLBC(IL)%MDOTPP  )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HM0     = GPG%ALLBC(IL)%HM      + F * (GPG%ALLBC(IH)%HM      - GPG%ALLBC(IL)%HM      )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%QEG     = GPG%ALLBC(IL)%QEG     + F * (GPG%ALLBC(IH)%QEG     - GPG%ALLBC(IL)%QEG     )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HC0G    = GPG%ALLBC(IL)%HCG     + F * (GPG%ALLBC(IH)%HCG     - GPG%ALLBC(IL)%HCG     )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%TINFG   = GPG%ALLBC(IL)%TINFG   + F * (GPG%ALLBC(IH)%TINFG   - GPG%ALLBC(IL)%TINFG   )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%TFIXEDG = GPG%ALLBC(IL)%TFIXEDG + F * (GPG%ALLBC(IH)%TFIXEDG - GPG%ALLBC(IL)%TFIXEDG )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%PRES    = GPG%ALLBC(IL)%PRES    + F * (GPG%ALLBC(IH)%PRES    - GPG%ALLBC(IL)%PRES    )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%YJINF(1:GPROP%NGSPEC) = GPG%ALLBC(IL)%YJINF(1:GPROP%NGSPEC) + F * (GPG%ALLBC(IH)%YJINF(1:GPROP%NGSPEC) - GPG%ALLBC(IL)%YJINF(1:GPROP%NGSPEC))

! Reradiation is discrete so use value at IL:
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%RERAD = GPG%ALLBC(IL)%RERADIATION

! Set solid enthalpy:
HFIXEDLO = 0.
HFIXEDHI = 0.
DO ISPEC = 1, SPROP%NSSPEC
   HFIXEDLO = HFIXEDLO + G%YIN(ISPEC,IZ,IX,IY) * HOFT(ISPEC,GPG%ALLBC(IL)%TFIXED)
   HFIXEDHI = HFIXEDHI + G%YIN(ISPEC,IZ,IX,IY) * HOFT(ISPEC,GPG%ALLBC(IH)%TFIXED)
ENDDO
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HFIXED  = HFIXEDLO  + F * (HFIXEDHI - HFIXEDLO )

! Set gas enthalpy:
HFIXEDLO = GPROP%CPG * ( GPG%ALLBC(IL)%TFIXEDG - GPG%TDATUM)
HFIXEDHI = GPROP%CPG * ( GPG%ALLBC(IH)%TFIXEDG - GPG%TDATUM)
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HFIXEDG = HFIXEDLO + F * (HFIXEDHI - HFIXEDLO )

! *****************************************************************************
END SUBROUTINE INTERPOLATE_BCS
! *****************************************************************************

! *****************************************************************************
END SUBROUTINE GET_ALL_BOUNDARY_CONDITIONS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE ADJUST_HC_FOR_BLOWING(IZ,IX,IY,NCELLZ,NCELLX,NCELLY,CTYPE,HC0,HC0BLOWING)
! *****************************************************************************

INTEGER, INTENT(IN) :: IZ, IX, IY, NCELLZ, NCELLX, NCELLY
CHARACTER(1), INTENT(IN) :: CTYPE
REAL(EB), INTENT(IN) :: HC0
REAL(EB), INTENT(OUT) :: HC0BLOWING
REAL(EB) :: MC

IF (GPG%SOLVE_PRESSURE) THEN
   SELECT CASE(CTYPE)
      CASE('X')
         IF (IX .EQ. 1     ) MC=ABS(G%MDOTPPDARCYW(IZ,IX,IY))
         IF (IX .EQ. NCELLX) MC=ABS(G%MDOTPPDARCYE(IZ,IX,IY))
      CASE('Y')
         IF (IY .EQ. 1     ) MC=ABS(G%MDOTPPDARCYS(IZ,IX,IY))
         IF (IY .EQ. NCELLY) MC=ABS(G%MDOTPPDARCYN(IZ,IX,IY))
      CASE('Z')
         IF (IZ .EQ. 1     ) MC=ABS(G%MDOTPPDARCYT(IZ,IX,IY))
         IF (IZ .EQ. NCELLZ) MC=ABS(G%MDOTPPDARCYB(IZ,IX,IY))
   END SELECT
ELSE
   SELECT CASE(CTYPE)
      CASE('X')
         WRITE(*,*) 'Cannot adjust x-direction hc for blowing with SOLVE_PRESSURE = F' 
      CASE('Y')
         WRITE(*,*) 'Cannot adjust y-direction hc for blowing with SOLVE_PRESSURE = F' 
      CASE('Z')
         MC=ABS(G%MDOTPPZ(0,IZ,IX,IY))     
   END SELECT
ENDIF

MC         = MAX(MC*GPROP%CPG,1D-10)
HC0BLOWING = MC / (EXP(MC/HC0) - 1D0)
 
! *****************************************************************************
END SUBROUTINE 
! *****************************************************************************


! *****************************************************************************
SUBROUTINE GET_BC_INFO(IX,IZ,IY,NBCPATCH,PATCH_BCLOC,PATCH_IX1,PATCH_IX2,PATCH_IZ1, &
           PATCH_IZ2,PATCH_IY1,PATCH_IY2,PATCH_IPNEXT,PATCH_T,BCLOC,TI,IP,IPNEXT,F)
! *****************************************************************************
INTEGER, INTENT(IN) :: IX,IZ,IY,NBCPATCH
REAL(EB), INTENT(IN) :: TI
CHARACTER(1), INTENT(IN) :: BCLOC
INTEGER, INTENT(IN), DIMENSION (:) :: PATCH_IX1,PATCH_IX2,PATCH_IZ1,PATCH_IZ2,PATCH_IY1,PATCH_IY2,PATCH_IPNEXT !prevents array temporary
REAL(EB), INTENT(IN), DIMENSION (:) :: PATCH_T
CHARACTER(1), INTENT(IN), DIMENSION (:) :: PATCH_BCLOC


INTEGER, INTENT(OUT) :: IP,IPNEXT
REAL(EB), INTENT(OUT) :: F

DO IP = 1, NBCPATCH
   IF (PATCH_BCLOC(IP) .NE. BCLOC) CYCLE
   IPNEXT = PATCH_IPNEXT(IP)

   IF (IPNEXT .NE. -1) THEN
      IF (TI .GT. PATCH_T(IPNEXT)) CYCLE
   ENDIF
         
   SELECT CASE (BCLOC)
      CASE ('T','B')
         IF (IX .LT. PATCH_IX1(IP)) CYCLE
         IF (IX .GT. PATCH_IX2(IP)) CYCLE
         IF (IY .LT. PATCH_IY1(IP)) CYCLE
         IF (IY .GT. PATCH_IY2(IP)) CYCLE         
      CASE ('E','W')
         IF (IZ .LT. PATCH_IZ1(IP)) CYCLE
         IF (IZ .GT. PATCH_IZ2(IP)) CYCLE
         IF (IY .LT. PATCH_IY1(IP)) CYCLE
         IF (IY .GT. PATCH_IY2(IP)) CYCLE
      CASE ('N','S')
         IF (IZ .LT. PATCH_IZ1(IP)) CYCLE
         IF (IZ .GT. PATCH_IZ2(IP)) CYCLE
         IF (IX .LT. PATCH_IX1(IP)) CYCLE
         IF (IX .GT. PATCH_IX2(IP)) CYCLE
   END SELECT

   IF (IPNEXT .EQ. -1) THEN
      F = 1D0
   ELSE
      F = 1D0 - (TI-PATCH_T(IP)) / (PATCH_T(IPNEXT)-PATCH_T(IP))
   ENDIF
                  
   RETURN
ENDDO 

! *****************************************************************************            
END SUBROUTINE GET_BC_INFO
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SET_BC_FIXED_VALUE(AP,AT,AB,AE,AW,AN,AS,B,FIXEDVAL)
! *****************************************************************************
REAL(EB), INTENT(OUT) :: AP,AT,AB,AE,AW,AN,AS,B
REAL(EB), INTENT(IN)  :: FIXEDVAL
      
AP = 1D0
AB = 0D0
AT = 0D0
AE = 0D0
AW = 0D0
AN = 0D0
AS = 0D0
B  = FIXEDVAL

! *****************************************************************************
END SUBROUTINE SET_BC_FIXED_VALUE
! *****************************************************************************

!! *****************************************************************************
!SUBROUTINE SET_BC_SOURCE(BCLOC,DLTX,DLTZ,DLTY,FLUXIN,B,SOUT) 
!! *****************************************************************************
!CHARACTER(1), INTENT(IN) :: BCLOC
!REAL(EB), INTENT(IN) :: DLTX,DLTZ,DLTY,FLUXIN
!REAL(EB), INTENT(OUT) :: B,SOUT
!
!SELECT CASE (BCLOC)
!   CASE ('T','B')
!      SOUT = SOUT + FLUXIN / DLTZ
!      B = B + FLUXIN * DLTX * DLTY
!   CASE ('E','W')
!      SOUT = SOUT + FLUXIN / DLTX
!      B = B + FLUXIN * DLTZ * DLTY
!   CASE ('N','S')
!      SOUT = SOUT + FLUXIN / DLTY
!      B = B + FLUXIN * DLTX * DLTZ
!   CASE DEFAULT
!      CONTINUE
!END SELECT
!
!! *****************************************************************************
!END SUBROUTINE SET_BC_SOURCE
!! *****************************************************************************

! *****************************************************************************
END MODULE GPYRO_BC
! *****************************************************************************