MODULE GPYRO_FUNCS

USE PREC
USE GPYRO_VARS

IMPLICIT NONE

CONTAINS

! *****************************************************************************
SUBROUTINE TDMA_SOLVER_GENERAL(NC,PHIPN,ATDMA,BTDMA,CTDMA,DTDMA)
! *****************************************************************************

INTEGER, INTENT(IN) :: NC
REAL(EB), DIMENSION(:), INTENT(IN) :: ATDMA,BTDMA,CTDMA,DTDMA
REAL(EB), DIMENSION(:), INTENT(OUT) :: PHIPN !Deferred shape to prevent creating array temporary
REAL(EB), DIMENSION(1:NC) :: PTDMA,QTDMA
INTEGER :: IC

PTDMA (1) = BTDMA(1)/ATDMA(1)
QTDMA (1) = DTDMA(1)/ATDMA(1)

! Use recursion to zip forward to back face cell
DO IC = 2, NC
   PTDMA(IC) = BTDMA(IC) / ( ATDMA(IC) - CTDMA(IC)*PTDMA(IC-1) )
   QTDMA(IC) =  ( DTDMA(IC) + CTDMA(IC) * QTDMA(IC-1) ) / ( ATDMA(IC) - CTDMA(IC) * PTDMA(IC-1) )
ENDDO 

! Set back face phi equal to QN
PHIPN(NC) = QTDMA(NC)
            
! Back substitution
DO IC = NC-1, 1, -1
   PHIPN(IC) = PTDMA(IC) * PHIPN(IC+1) + QTDMA(IC)
ENDDO

! *****************************************************************************      
END SUBROUTINE TDMA_SOLVER_GENERAL
! *****************************************************************************

! *****************************************************************************
REAL(EB) FUNCTION INTEGRATED_MASS_LOSS_RATE(IGSPEC)
! *****************************************************************************
! Calculates total mass loss rate for current mesh for species ISPEC. 
! If ISPEC = 0, calculates total mass loss rate for all species.
!
! When SOLVE_PRESSURE=.TRUE.:  Darcian mass flux gets dotted with the surface 
! normal at boundaries so that outflows are positive and inflows are negative.
!
! When SOLVE_PRESSURE=.FALSE.: Mass loss rate is summed over all cells 
!
! In all cases, MLR gets multiplied by 1000 so units are g/m2-s in 1D, g/m-s 
! in 2D, and g/s in 3D. 

INTEGER, INTENT(IN) :: IGSPEC
REAL(EB) :: ML, V, YJ
INTEGER :: IZ, IX, IY

ML = 0.
YJ = 1.

IF (GPG%SOLVE_PRESSURE) THEN
   DO IY = 1, G%NCELLY
   DO IX = 1, G%NCELLX
   DO IZ = 1, G%NCELLZ
      IF (G%IMASK(IZ,IX,IY) .EQ. 0) CYCLE
      IF (IGSPEC .NE. 0) YJ = G%YJGN(IGSPEC,IZ,IX,IY)
      IF (G%NEEDSBCT(IZ,IX,IY)) ML = ML - G%MDOTPPDARCYT(IZ,IX,IY)*G%DXDY(IZ,IX,IY)*YJ
      IF (G%NEEDSBCB(IZ,IX,IY)) ML = ML + G%MDOTPPDARCYB(IZ,IX,IY)*G%DXDY(IZ,IX,IY)*YJ
      IF (G%NCELLX .GT. 1) THEN
         IF (G%NEEDSBCW(IZ,IX,IY)) ML = ML - G%MDOTPPDARCYW(IZ,IX,IY)*G%DYDZ(IZ,IX,IY)*YJ
         IF (G%NEEDSBCE(IZ,IX,IY)) ML = ML + G%MDOTPPDARCYE(IZ,IX,IY)*G%DYDZ(IZ,IX,IY)*YJ
      ENDIF
      IF (G%NCELLY .GT. 1) THEN
         IF (G%NEEDSBCS(IZ,IX,IY)) ML = ML - G%MDOTPPDARCYS(IZ,IX,IY)*G%DXDZ(IZ,IX,IY)*YJ
         IF (G%NEEDSBCN(IZ,IX,IY)) ML = ML + G%MDOTPPDARCYN(IZ,IX,IY)*G%DXDZ(IZ,IX,IY)*YJ
      ENDIF
   ENDDO
   ENDDO
   ENDDO
ELSE
   DO IY = 1, G%NCELLY
   DO IX = 1, G%NCELLX
   DO IZ = 1, G%NCELLZ
      IF (G%IMASK(IZ,IX,IY) .EQ. 0) CYCLE
      V = G%DXDY(IZ,IX,IY)*G%DLTZN(IZ,IX,IY) ! This works if Dz decreases
      ML = ML + G%GOMEGA(3,IGSPEC,IZ,IX,IY) * V
   ENDDO
   ENDDO
   ENDDO
ENDIF

INTEGRATED_MASS_LOSS_RATE=1000.*ML

! *****************************************************************************
END FUNCTION INTEGRATED_MASS_LOSS_RATE
! *****************************************************************************

! *****************************************************************************
REAL(EB) FUNCTION INTEGRATED_GAS_GENERATION_RATE(IGSPEC)
! *****************************************************************************

INTEGER, INTENT(IN) :: IGSPEC
REAL(EB) :: ML, V
INTEGER :: IG, IZ, IX, IY

ML = 0.

IF (IGSPEC .EQ. 0) THEN
! Sum gas generation rate for all species
   DO IY = 1, G%NCELLY
   DO IX = 1, G%NCELLX
   DO IZ = 1, G%NCELLZ
      IF (G%IMASK(IZ,IX,IY) .EQ. 0) CYCLE
      V = G%DXDY(IZ,IX,IY)*G%DLTZN(IZ,IX,IY) ! This works if Dz decreases
      DO IG = 1, GPROP%NGSPEC
         ML = ML + G%GOMEGA(3,IG,IZ,IX,IY) * V
      ENDDO
   ENDDO
   ENDDO
   ENDDO
ELSE
! Just get gas generation rate for one species
   DO IY = 1, G%NCELLY
   DO IX = 1, G%NCELLX
   DO IZ = 1, G%NCELLZ
      IF (G%IMASK(IZ,IX,IY) .EQ. 0) CYCLE
      V = G%DXDY(IZ,IX,IY)*G%DLTZN(IZ,IX,IY) ! This works if Dz decreases
      ML = ML + G%GOMEGA(3,IGSPEC,IZ,IX,IY) * V
   ENDDO
   ENDDO
   ENDDO
ENDIF

INTEGRATED_GAS_GENERATION_RATE=1000.*ML

! *****************************************************************************
END FUNCTION INTEGRATED_GAS_GENERATION_RATE
! *****************************************************************************

! *****************************************************************************
REAL(EB) FUNCTION INTEGRATED_HEAT_RELEASE_RATE()
! *****************************************************************************
! Calulates total heat release rate for current mesh. 

REAL(EB) :: HRR, V
INTEGER :: IZ, IX, IY, IRXN

HRR = 0.

DO IRXN = 1, SPROP%NRXN

   IF (RXN(IRXN)%DHS .GT. -1D-6 .AND. RXN(IRXN)%DHS .LT. 1D-6) THEN
      CONTINUE
   ELSE
      DO IZ = 1, G%NCELLZ
      DO IX = 1, G%NCELLX
      DO IY = 1, G%NCELLY
         IF (G%IMASK(IZ,IX,IY) .EQ. 0) CYCLE
         V = G%DXDY(IZ,IX,IY) * G%DLTZN(IZ,IX,IY) 
         HRR = HRR - G%OMEGASFBK(IRXN,IZ,IX,IY) * RXN(IRXN)%DHS * V
      ENDDO
      ENDDO
      ENDDO
   ENDIF

   IF (RXN(IRXN)%DHV .GT. -1D-6 .AND. RXN(IRXN)%DHV .LT. 1D-6) THEN
      CONTINUE
   ELSE
      DO IZ = 1, G%NCELLZ
      DO IX = 1, G%NCELLX
      DO IY = 1, G%NCELLY
         IF (G%IMASK(IZ,IX,IY) .EQ. 0) CYCLE
         V = G%DXDY(IZ,IX,IY) * G%DLTZN(IZ,IX,IY) 
         HRR = HRR - G%OMEGASFGK(IRXN,IZ,IX,IY) * RXN(IRXN)%DHV * V
      ENDDO
      ENDDO
      ENDDO
   ENDIF

ENDDO

! Should add homogeneous gas reactions too

INTEGRATED_HEAT_RELEASE_RATE = 0.001 * HRR

! *****************************************************************************
END FUNCTION INTEGRATED_HEAT_RELEASE_RATE
! *****************************************************************************

! *****************************************************************************
REAL(EB) FUNCTION TOTAL_MASS(ISPEC)
! *****************************************************************************
! Total mass, in g/m2 (1D), g/m (2D), or g (3D)

INTEGER, INTENT(IN) :: ISPEC
REAL(EB) :: MG,MS,M,V
INTEGER :: IZ, IX, IY

MG=0.
M =0.

DO IY = 1, G%NCELLY
DO IX = 1, G%NCELLX
DO IZ = 1, G%NCELLZ
   V = G%DLTZN(IZ,IX,IY) * G%DLTXN(IZ,IX,IY) * G%DLTYN(IZ,IX,IY)
   IF (ISPEC .EQ. 0) THEN
      MS = REAL(G%IMASK(IZ,IX,IY)) * G%RPN(IZ,IX,IY) 
      IF (GPG%SOLVE_PRESSURE) MG = REAL(G%IMASK(IZ,IX,IY))*G%POROSSN(IZ,IX,IY)*G%RGN(IZ,IX,IY)
   ELSE
      MS = REAL(G%IMASK(IZ,IX,IY)) * G%RPN(IZ,IX,IY) * G%YIN(ISPEC,IZ,IX,IY)
   ENDIF
   M = M + (MS + MG) * V
ENDDO
ENDDO
ENDDO

TOTAL_MASS = 1000.*M

! *****************************************************************************
END FUNCTION TOTAL_MASS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE CALC_MDOTPPZ(IMESH,NCELLZ,NCELLX,NCELLY)
! *****************************************************************************
! Calcuates z-direction mass flux from mass conservation. Only used in 
! 1D simulations when SOLVE_PRESSURE=.FALSE.

INTEGER, INTENT(IN) :: IMESH,NCELLZ,NCELLX,NCELLY
INTEGER :: IZ,IX,IY,IGSPEC,IOR
REAL(EB) :: TSTART,TEND

CALL GET_CPU_TIME(TSTART)

G => GPM(IMESH)
GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)

IOR=3 ! (facing +z direction)

DO IY = 1, NCELLY
DO IX = 1, NCELLX

   G%MDOTPPZ(:,:,IX,IY) = 0D0

! Account for mass flux due to internally-forced flow:
   DO IZ = G%NCELLZ, 1, -1
      IF (G%NEEDSBCT(IZ,IX,IY) .AND. GPBCP(NCELLZ,IX,IY,IOR)%PRES .LT. 0.) THEN
         DO IGSPEC = 1, GPROP%NGSPEC
            G%MDOTPPZ(IGSPEC,IZ,IX,IY) = 1D-3*GPBCP(IZ,IX,IY,IOR)%MFLUX * GPBCP(IZ,IX,IY,IOR)%YJINF(IGSPEC) 
         ENDDO
      ENDIF
   ENDDO

! Now move from back face to front face and integrate mass loss rate:
   DO IGSPEC = 1, GPROP%NGSPEC
      G%MDOTPPZ(IGSPEC,NCELLZ,IX,IY) = G%MDOTPPZ(IGSPEC,NCELLZ,IX,IY) + G%GOMEGA(3,IGSPEC,NCELLZ,IX,IY) * G%DLTZN(NCELLZ,IX,IY)
      DO IZ = NCELLZ-1, 1, -1
         G%MDOTPPZ(IGSPEC,IZ,IX,IY) = G%MDOTPPZ(IGSPEC,IZ+1,IX,IY) + G%GOMEGA(3,IGSPEC,IZ,IX,IY)*G%DLTZN(IZ,IX,IY)
      ENDDO
   ENDDO

! Now sum over species
   DO IGSPEC = 1, GPROP%NGSPEC
      DO IZ = 1, NCELLZ
         G%MDOTPPZ(0,IZ,IX,IY) = G%MDOTPPZ(0,IZ,IX,IY) + G%MDOTPPZ(IGSPEC,IZ,IX,IY) 
      ENDDO
   ENDDO

ENDDO
ENDDO

CALL GET_CPU_TIME(TEND)
GPG%TUSED(18) = GPG%TUSED(18) + TEND - TSTART

! *****************************************************************************
END SUBROUTINE CALC_MDOTPPZ
! *****************************************************************************

! *****************************************************************************
REAL(EB) FUNCTION RHOOFT(ISPEC,TMP) 
! *****************************************************************************

INTEGER,  INTENT(IN) :: ISPEC 
REAL(EB), INTENT(IN) :: TMP

IF (SPROP%NR(ISPEC) .EQ. 0D0) THEN
   RHOOFT = SPROP%R0(ISPEC)
ELSE
   RHOOFT = SPROP%R0(ISPEC) * (TMP/GPG%TREF)**SPROP%NR(ISPEC)
ENDIF

! *****************************************************************************
END FUNCTION RHOOFT
! *****************************************************************************

! *****************************************************************************
REAL(EB) FUNCTION RHOGOFT(P,M,T) 
! *****************************************************************************

REAL(EB), INTENT(IN) :: P,M,T
REAL(EB), PARAMETER :: R = 8314D0

RHOGOFT = P * M / (R * T)

! *****************************************************************************
END FUNCTION RHOGOFT
! *****************************************************************************

!******************************************************************************
REAL(EB) FUNCTION CPOFT(ISPEC,TMP)
!******************************************************************************

REAL(EB), INTENT(IN) :: TMP
INTEGER, INTENT(IN) :: ISPEC
REAL(EB), PARAMETER :: TWOPI = 2D0*PI
REAL(EB) :: TERM1,TERM2,DHMELT,SIGMA2MELT,TMELT

IF (SPROP%NC(ISPEC) .EQ. 0D0) THEN
   CPOFT  = SPROP%C0(ISPEC)
ELSE
   CPOFT  = SPROP%C0(ISPEC) * (TMP/GPG%TREF) ** SPROP%NC(ISPEC)
ENDIF

DHMELT = SPROP%DHMELT(ISPEC)
IF (DHMELT .GT. 0D0) THEN
   SIGMA2MELT = SPROP%SIGMA2MELT(ISPEC) 
   TMELT      = SPROP%TMELT(ISPEC) 
   TERM1 = DHMELT/SQRT(TWOPI*SIGMA2MELT)
   TERM2 = ((TMP-TMELT)**2D0)/(2D0*SIGMA2MELT)
   CPOFT = CPOFT + TERM1 * EXP(-TERM2)
ENDIF

END FUNCTION CPOFT

!******************************************************************************
REAL(EB) FUNCTION HOFT(ISPEC,TMP)
!******************************************************************************

REAL(EB), INTENT(IN) :: TMP
INTEGER, INTENT(IN) :: ISPEC
REAL(EB), PARAMETER :: TWOPI = 2D0*PI
REAL(EB) :: TERM1,TERM2,SIGMA2MELT,TMELT,C0,NC,TOTR,T0OTR

IF (SPROP%NC(ISPEC) .EQ. 0D0) THEN
   HOFT = SPROP%C0(ISPEC) * (TMP - GPG%TDATUM)
ELSE
   C0         = SPROP%C0(ISPEC)
   NC         = SPROP%NC(ISPEC)
   TOTR       = TMP / GPG%TREF
   T0OTR      = GPG%TDATUM / GPG%TREF
   HOFT = C0 * (TMP*TOTR**NC - GPG%TDATUM*T0OTR**NC) / (NC + 1D0)
ENDIF

IF (SPROP%DHMELT(ISPEC) .GT. 0D0) THEN
   SIGMA2MELT = SPROP%SIGMA2MELT(ISPEC) 
   TMELT      = SPROP%TMELT(ISPEC) 
   TERM1 = DERF( (TMP       - TMELT) / SQRT(2D0*SIGMA2MELT) )
   TERM2 = DERF( (GPG%TDATUM - TMELT) / SQRT(2D0*SIGMA2MELT) )
   HOFT  = HOFT + 0.5D0*SPROP%DHMELT(ISPEC)*(TERM1 - TERM2)
ENDIF

!******************************************************************************
END FUNCTION HOFT
!******************************************************************************

!******************************************************************************
REAL(EB) FUNCTION KZOFT(ISPEC,TMP)
!******************************************************************************

INTEGER, INTENT(IN) :: ISPEC
REAL(EB), INTENT(IN) :: TMP
REAL(EB) :: K0, NK

K0 = SPROP%K0Z(ISPEC)
NK = SPROP%NKZ(ISPEC)

IF (NK .EQ. 0D0) THEN
   KZOFT = K0
ELSE
   KZOFT = K0 * (TMP/GPG%TREF)**NK
ENDIF

!******************************************************************************
END FUNCTION KZOFT
!******************************************************************************

!******************************************************************************
REAL(EB) FUNCTION KXOFT(ISPEC,TMP)
!******************************************************************************

INTEGER, INTENT(IN) :: ISPEC
REAL(EB), INTENT(IN) :: TMP
REAL(EB) :: K0, NK

K0 = SPROP%K0X(ISPEC)
NK = SPROP%NKX(ISPEC)

IF (NK .EQ. 0D0) THEN
   KXOFT = K0
ELSE
   KXOFT = K0 * (TMP/GPG%TREF)**NK
ENDIF

!******************************************************************************
END FUNCTION KXOFT
!******************************************************************************

!******************************************************************************
REAL(EB) FUNCTION KYOFT(ISPEC,TMP)
!******************************************************************************

INTEGER, INTENT(IN) :: ISPEC
REAL(EB), INTENT(IN) :: TMP
REAL(EB) :: K0, NK

K0 = SPROP%K0Y(ISPEC)
NK = SPROP%NKY(ISPEC)

IF (NK .EQ. 0D0) THEN
   KYOFT = K0
ELSE
   KYOFT = K0 * (TMP/GPG%TREF)**NK
ENDIF

!******************************************************************************
END FUNCTION KYOFT
!******************************************************************************

! *****************************************************************************
REAL(EB) FUNCTION DOFT(I1,I2,TMP,PRES)
! *****************************************************************************
! Calculate binary diffusion coefficient for species I1 into species I2
! at temperature TMP and pressure PRES
INTEGER, INTENT(IN) :: I1,I2 
REAL(EB), INTENT(IN) :: TMP,PRES
REAL(EB) :: TSTAR,OMEGAD,NUMER,EPSOK,SIG2,NUMER1

EPSOK  = SQRT(GPROP%EPSOK(I1)*GPROP%EPSOK(I2))
SIG2   = (0.5D0*(GPROP%SIGMA(I1)+GPROP%SIGMA(I2)))**2D0
NUMER1 = 1D0/GPROP%M(I1) + 1D0/GPROP%M(I2)

TSTAR = TMP / EPSOK
OMEGAD = AD(1)*TSTAR**AD(2) + AD(3)*EXP(TSTAR*AD(4)) + AD(5)*EXP(AD(6)*TSTAR)
NUMER = SQRT(NUMER1 * TMP**3D0)
DOFT = 0.0188129*NUMER / (PRES*SIG2*OMEGAD)

! *****************************************************************************      
END FUNCTION DOFT
! *****************************************************************************

!******************************************************************************
REAL(EB) FUNCTION TOFH_NEWTON(HN,TGUESS,NSPEC,YN)
!******************************************************************************

REAL(EB), INTENT(IN) :: HN !Average (NEW) enthalpy
REAL(EB), INTENT(IN) :: TGUESS !Temperature guess
INTEGER,  INTENT(IN) :: NSPEC !Number of species
REAL(EB), DIMENSION(1:NSPEC), INTENT(IN) :: YN !Mass fraction array
REAL(EB), DIMENSION(1:NSPEC) :: AI,BI,DI
REAL(EB) :: RHS
INTEGER :: ISPEC
LOGICAL :: MELTING

RHS = HN
AI(:) = 0D0
BI(:) = 0D0
DI(:) = 0D0
MELTING = .FALSE.

DO ISPEC = 1, NSPEC
   BI(ISPEC) = GPG%NEWTON_B(ISPEC)
   AI(ISPEC) = YN(ISPEC) * GPG%NEWTON_A(ISPEC)
   RHS = RHS + YN(ISPEC) * GPG%NEWTON_C(ISPEC) 
   IF (SPROP%DHMELT(ISPEC) .NE. 0D0) THEN
      MELTING = .TRUE.
      RHS = RHS + YN(ISPEC) * GPG%NEWTON_G(ISPEC)
      DI(ISPEC) = YN(ISPEC) * GPG%NEWTON_D(ISPEC)
   ENDIF
ENDDO
      
TOFH_NEWTON = RTNEWT(FUNCD,290D0,3000D0,GPG%HTOL,NSPEC,TGUESS,RHS,AI,BI,DI,MELTING)

! *****************************************************************************
END FUNCTION TOFH_NEWTON
! *****************************************************************************

! *****************************************************************************
SUBROUTINE FUNCD(T,F,DF,RHS,NSPEC,AI,BI,DI,MELTING)
! *****************************************************************************
REAL(EB), INTENT(IN) :: T,RHS
INTEGER, INTENT(IN) :: NSPEC
REAL(EB), INTENT(IN) :: AI(1:NSPEC),BI(1:NSPEC),DI(1:NSPEC)
REAL(EB), INTENT(OUT) :: F,DF
LOGICAL, INTENT(IN) :: MELTING
INTEGER :: ISPEC,ILO,ITMAX,IHI
REAL(EB) :: TSTAR,MULT,EXPTRM,TCLIP
REAL(EB), PARAMETER :: SQRTRPI =  0.564189584D0 !SQRT(1D0/PI)
      
F  = -RHS
DF = 0D0

ITMAX = NINT(  10D0*(4000D0 - GPG%TDATUM) )

IF (MELTING) THEN
   DO ISPEC = 1, NSPEC
      TSTAR = (T-GPG%NEWTON_E(ISPEC)) / GPG%NEWTON_F(ISPEC)
      F =  F + AI(ISPEC)*T**BI(ISPEC) + DI(ISPEC) * DERF(TSTAR)
      MULT = SPROP%DHMELT(ISPEC) * SQRTRPI / GPG%NEWTON_F(ISPEC)
      EXPTRM = (T - GPG%NEWTON_E(ISPEC)) / GPG%NEWTON_F(ISPEC)
      DF = DF + AI(ISPEC)*   BI(ISPEC) * T**SPROP%NC(ISPEC) + MULT * EXP(-(EXPTRM*EXPTRM))
   ENDDO
ELSE
   IF (T .NE. T) THEN
      TCLIP = GPG%TAMB
   ELSE
      TCLIP = T
   ENDIF
         
   TCLIP = MAX(GPG%TDATUM,TCLIP)
   TCLIP = MIN(4D3,TCLIP)
               
   ILO = INT (10D0*(TCLIP-GPG%TDATUM) )
   ILO = MIN(ILO, ITMAX-1)
   IHI = MIN(ILO+1, ITMAX)

   DO ISPEC = 1, NSPEC
      F  =  F + AI(ISPEC)*T**BI(ISPEC)
      DF = DF + AI(ISPEC)*   BI(ISPEC) * T**SPROP%NC(ISPEC)
   ENDDO
ENDIF

! *****************************************************************************
END SUBROUTINE FUNCD
! *****************************************************************************

! *****************************************************************************
REAL(EB) FUNCTION CALC_EMIS(IZ,IX,IY)
! *****************************************************************************

INTEGER, INTENT(IN) :: IZ,IX,IY
INTEGER :: ISPEC
REAL(EB) :: E

!G=>GPM(IMESH)
!GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)

E = 0D0
DO ISPEC = 1, SPROP%NSSPEC
   E = E + G%XI(ISPEC,IZ,IX,IY) * SPROP%EMIS(ISPEC)
ENDDO
CALC_EMIS = E

! *****************************************************************************
END FUNCTION CALC_EMIS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE GET_A(NC,ISCHEME,PT,PB,BIGAT,BIGAB)
! *****************************************************************************
INTEGER, INTENT(IN) :: NC,ISCHEME
REAL(EB), DIMENSION(1:NC), INTENT(IN ) :: PT,PB
REAL(EB), DIMENSION(1:NC), INTENT(OUT) :: BIGAT,BIGAB 

SELECT CASE (ISCHEME)
   CASE(1) !Central differences
      BIGAT(:) = 1D0 - 0.5D0 * ABS(PT(:))
      BIGAB(:) = 1D0 - 0.5D0 * ABS(PB(:))
   CASE(2) !Upwind
      BIGAT(:) = 1D0
      BIGAB(:) = 1D0
   CASE(3) !Hybrid
      BIGAT(:) =  MAX (0D0, 1D0 - 0.5D0 * ABS( PT(:) ) )
      BIGAB(:) =  MAX (0D0, 1D0 - 0.5D0 * ABS( PB(:) ) )
   CASE(4) !Power law
      BIGAT(:) =  MAX (0D0,(1D0 - 0.1D0 * ABS( PT(:) ) )**5 )
      BIGAB(:) =  MAX (0D0,(1D0 - 0.1D0 * ABS( PB(:) ) )**5 )
   CASE(5) !Exponential
      BIGAT(:) = ABS(PT(:)) / (EXP(ABS(PT(:))) - 1D0)
      BIGAB(:) = ABS(PB(:)) / (EXP(ABS(PB(:))) - 1D0)
END SELECT

! *****************************************************************************      
END SUBROUTINE GET_A
! *****************************************************************************


! *****************************************************************************
SUBROUTINE GET_A_SINGLE(ISCHEME,PT,PB,BIGAT,BIGAB)
! *****************************************************************************
INTEGER,  INTENT(IN ) :: ISCHEME
REAL(EB), INTENT(IN ) :: PT,PB
REAL(EB), INTENT(OUT) :: BIGAT,BIGAB 

SELECT CASE (ISCHEME)
   CASE(1) !Central differences
      BIGAT = 1D0 - 0.5D0 * ABS(PT)
      BIGAB = 1D0 - 0.5D0 * ABS(PB)
   CASE(2) !Upwind
      BIGAT = 1D0
      BIGAB = 1D0
   CASE(3) !Hybrid
      BIGAT =  MAX (0D0, 1D0 - 0.5D0 * ABS( PT ) )
      BIGAB =  MAX (0D0, 1D0 - 0.5D0 * ABS( PB ) )
   CASE(4) !Power law
      BIGAT =  MAX (0D0,(1D0 - 0.1D0 * ABS( PT ) )**5 )
      BIGAB =  MAX (0D0,(1D0 - 0.1D0 * ABS( PB ) )**5 )
   CASE(5) !Exponential
      BIGAT = ABS(PT) / (EXP(ABS(PT)) - 1D0)
      BIGAB = ABS(PB) / (EXP(ABS(PB)) - 1D0)
END SELECT

! *****************************************************************************      
END SUBROUTINE GET_A_SINGLE
! *****************************************************************************




!******************************************************************************
REAL(EB) FUNCTION LINTERP(X,XL,XH,FXL,FXH)
!******************************************************************************

REAL(EB), INTENT(IN) :: X, XL,XH,FXL,FXH
LINTERP = FXL + (FXH-FXL)*(X-XL)/(XH-XL)

! *****************************************************************************
END FUNCTION LINTERP
!******************************************************************************

!******************************************************************************
REAL(EB) FUNCTION rtnewt(funcd,x1,x2,xacc,NSPEC,TGUESS,RHS,AI,BI,DI,MELTING)
!******************************************************************************

REAL(EB), INTENT(IN) :: TGUESS,RHS
INTEGER, INTENT(IN) :: NSPEC
LOGICAL, INTENT(IN) :: MELTING
REAL(EB) :: x1,x2,xacc,AI(1:NSPEC),BI(1:NSPEC),DI(1:NSPEC)
EXTERNAL funcd
INTEGER, PARAMETER :: NEWTMAX=50 ! Set to maximum number of iterations.
! Using the Newton-Raphson method, and the root of a function known to lie in the interval
! [x1; x2]. The root rtnewt will be retuned until its accuracy is known within .xacc. funcd
! is a user-supplied subroutine that returns both the function value and the first derivative
! of the function at the point x.
INTEGER :: j
REAL(EB) :: df,dx,f

RTNEWT=TGUESS
do j=1,NEWTMAX
   call funcd(rtnewt,f,df,RHS,NSPEC,AI,BI,DI,MELTING)
   dx=f/df
   rtnewt=rtnewt-dx
   if((x1-rtnewt)*(rtnewt-x2).lt.0.) CONTINUE
   if(abs(dx).lt.xacc) return !Convergence.
enddo 

! *****************************************************************************
END FUNCTION RTNEWT
! *****************************************************************************

! *****************************************************************************
SUBROUTINE locate(xx,n,x,j)
! *****************************************************************************
! Given an array xx(1:n), and given a value x, returns a value j such that x is between
! xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0
! or j=n is returned to indicate that x is out of range.

INTEGER j,n
REAL(EB) :: x,xx(n)
INTEGER jl,jm,ju
jl=0 !Initialize lower
ju=n+1 !and upper limits.
 10 if(ju-jl.gt.1)then !If we are not yet done,
       jm=(ju+jl)/2 !compute a midpoint,
       if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm !and replace either the lower limit
       else
          ju=jm !or the upper limit, as appropriate.
       endif
       goto 10 !Repeat until
    endif !the test condition 10 is satisfied.
      
if(x.eq.xx(1))then !Then set the output
   j=1
else if(x.eq.xx(n))then
   j=n-1
else
   j=jl
endif

return !and return.

! *****************************************************************************
END SUBROUTINE LOCATE
! *****************************************************************************
      
! *****************************************************************************
REAL(EB) FUNCTION DERF (X)
! *****************************************************************************

! *****************************************************************************
! Double precision error function so that IMSL libraries are not required. 

real(EB) :: a(0:64), b(0:64), x, w, t, y

integer :: i,k

data (a(i), i = 0, 12) /                                 &
    0.00000000005958930743d0, -0.00000000113739022964d0, &
    0.00000001466005199839d0, -0.00000016350354461960d0, &
    0.00000164610044809620d0, -0.00001492559551950604d0, &
    0.00012055331122299265d0, -0.00085483269811296660d0, &
    0.00522397762482322257d0, -0.02686617064507733420d0, &
    0.11283791670954881569d0, -0.37612638903183748117d0, &
    1.12837916709551257377d0 / 
data (a(i), i = 13, 25) /                                &
    0.00000000002372510631d0, -0.00000000045493253732d0, &
    0.00000000590362766598d0, -0.00000006642090827576d0, &
    0.00000067595634268133d0, -0.00000621188515924000d0, &
    0.00005103883009709690d0, -0.00037015410692956173d0, &
    0.00233307631218880978d0, -0.01254988477182192210d0, &
    0.05657061146827041994d0, -0.21379664776456006580d0, &
    0.84270079294971486929d0 / 
data (a(i), i = 26, 38) /                                &
    0.00000000000949905026d0, -0.00000000018310229805d0, &
    0.00000000239463074000d0, -0.00000002721444369609d0, &
    0.00000028045522331686d0, -0.00000261830022482897d0, &
    0.00002195455056768781d0, -0.00016358986921372656d0, &
    0.00107052153564110318d0, -0.00608284718113590151d0, &
    0.02986978465246258244d0, -0.13055593046562267625d0, &
    0.67493323603965504676d0 / 
 data (a(i), i = 39, 51) /                                &
    0.00000000000382722073d0, -0.00000000007421598602d0, &
    0.00000000097930574080d0, -0.00000001126008898854d0, &
    0.00000011775134830784d0, -0.00000111992758382650d0, &
    0.00000962023443095201d0, -0.00007404402135070773d0, &
    0.00050689993654144881d0, -0.00307553051439272889d0, &
    0.01668977892553165586d0, -0.08548534594781312114d0, &
    0.56909076642393639985d0 / 
data (a(i), i = 52, 64) /                                &
    0.00000000000155296588d0, -0.00000000003032205868d0, &
    0.00000000040424830707d0, -0.00000000471135111493d0, &
    0.00000005011915876293d0, -0.00000048722516178974d0, &
    0.00000430683284629395d0, -0.00003445026145385764d0, &
    0.00024879276133931664d0, -0.00162940941748079288d0, &
    0.00988786373932350462d0, -0.05962426839442303805d0, &
    0.49766113250947636708d0 / 
data (b(i), i = 0, 12) /                                  &
    -0.00000000029734388465d0, 0.00000000269776334046d0,  &
    -0.00000000640788827665d0, -0.00000001667820132100d0, &
    -0.00000021854388148686d0, 0.00000266246030457984d0,  &
    0.00001612722157047886d0, -0.00025616361025506629d0,  &
    0.00015380842432375365d0, 0.00815533022524927908d0,   &
    -0.01402283663896319337d0, -0.19746892495383021487d0, & 
    0.71511720328842845913d0 / 
data (b(i), i = 13, 25) /                                 &
    -0.00000000001951073787d0, -0.00000000032302692214d0, &
    0.00000000522461866919d0, 0.00000000342940918551d0,   &
    -0.00000035772874310272d0, 0.00000019999935792654d0,  &
    0.00002687044575042908d0, -0.00011843240273775776d0,  &
    -0.00080991728956032271d0, 0.00661062970502241174d0,  &
    0.00909530922354827295d0, -0.20160072778491013140d0,  &
    0.51169696718727644908d0 /  
data (b(i), i = 26, 38) /                                 &
    0.00000000003147682272d0, -0.00000000048465972408d0,  &
    0.00000000063675740242d0, 0.00000003377623323271d0,   &
    -0.00000015451139637086d0, -0.00000203340624738438d0, &
    0.00001947204525295057d0, 0.00002854147231653228d0,   &
    -0.00101565063152200272d0, 0.00271187003520095655d0,  &
    0.02328095035422810727d0, -0.16725021123116877197d0,  &
    0.32490054966649436974d0 / 
data (b(i), i = 39, 51) /                                 &
    0.00000000002319363370d0, -0.00000000006303206648d0,  &
    -0.00000000264888267434d0, 0.00000002050708040581d0,  &
    0.00000011371857327578d0, -0.00000211211337219663d0,  &
    0.00000368797328322935d0, 0.00009823686253424796d0,   &
    -0.00065860243990455368d0, -0.00075285814895230877d0, &
    0.02585434424202960464d0, -0.11637092784486193258d0,  &
    0.18267336775296612024d0 / 
data (b(i), i = 52, 64) /                                 &
    -0.00000000000367789363d0, 0.00000000020876046746d0,  &
    -0.00000000193319027226d0, -0.00000000435953392472d0, &
    0.00000018006992266137d0, -0.00000078441223763969d0,  &
    -0.00000675407647949153d0, 0.00008428418334440096d0,  &
    -0.00017604388937031815d0, -0.00239729611435071610d0, &
    0.02064129023876022970d0, -0.06905562880005864105d0,  &
    0.09084526782065478489d0 / 
w = abs(x)
if (w .lt. 2.2d0) then
    t = w * w
    k = int(t)
    t = t - k
    k = k * 13
    y = ((((((((((((a(k) * t + a(k + 1)) * t +             &
       a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t +    &
       a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t +    &
       a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t +   &
       a(k + 11)) * t + a(k + 12)) * w
else if (w .lt. 6.9d0) then
    k = int(w)
    t = w - k
    k = 13 * (k - 2)
    y = (((((((((((b(k) * t + b(k + 1)) * t +              &
        b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t +    &
        b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t +    &
        b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t +   &
        b(k + 11)) * t + b(k + 12)
        y = y * y
        y = y * y
        y = y * y
        y = 1 - y * y
else
    y = 1
end if

if (x .lt. 0) y = -y

derf = y

! *****************************************************************************
end function derf
! *****************************************************************************      

! *****************************************************************************      
REAL(EB) FUNCTION DERFC(X) 
! *****************************************************************************      

REAL(EB) :: pa, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11,  &
          p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, &
          x, t, u, y

parameter (                         &
    pa = 3.97886080735226000d+00,   &
    p0 = 2.75374741597376782d-01,   &
    p1 = 4.90165080585318424d-01,   &
    p2 = 7.74368199119538609d-01,   &
    p3 = 1.07925515155856677d+00,   &
    p4 = 1.31314653831023098d+00,   &
    p5 = 1.37040217682338167d+00,   &
    p6 = 1.18902982909273333d+00,   &
    p7 = 8.05276408752910567d-01,   &
    p8 = 3.57524274449531043d-01,   &
    p9 = 1.66207924969367356d-02,   &
    p10 = -1.19463959964325415d-01, & 
    p11 = -8.38864557023001992d-02)
parameter (                         &
    p12 = 2.49367200053503304d-03,  &
    p13 = 3.90976845588484035d-02,  &
    p14 = 1.61315329733252248d-02,  &
    p15 = -1.33823644533460069d-02, &
    p16 = -1.27223813782122755d-02, &
    p17 = 3.83335126264887303d-03,  &
    p18 = 7.73672528313526668d-03,  &
    p19 = -8.70779635317295828d-04, &
    p20 = -3.96385097360513500d-03, &
    p21 = 1.19314022838340944d-04,  &
    p22 = 1.27109764952614092d-03)

t = pa / (pa + abs(x))
u = t - 0.5d0
y = (((((((((p22 * u + p21) * u + p20) * u +            &
    p19) * u + p18) * u + p17) * u + p16) * u +         &
    p15) * u + p14) * u + p13) * u + p12
y = ((((((((((((y * u + p11) * u + p10) * u +           &
    p9) * u + p8) * u + p7) * u + p6) * u + p5) * u +   &
    p4) * u + p3) * u + p2) * u + p1) * u + p0) * t *   &
    exp(-x * x)

if (x .lt. 0) y = 2 - y

derfc = y

! *****************************************************************************
end function derfc
! *****************************************************************************      

! *****************************************************************************      
REAL(EB) FUNCTION WALL_CLOCK_TIME_THANKS_FDS()  ! Returns the number of seconds since January 1, 2000, including leap years
! *****************************************************************************

! Thirty days hath September,
! April, June, and November.
! February has twenty-eight alone;
! All the rest have thirty-one,
! Excepting Leap-Year, that's the time
! When February's days are twenty-nine.

INTEGER :: DATE_TIME(8),WALL_CLOCK_SECONDS
CHARACTER(10) :: BIG_BEN(3)
! X_1 = common year, X_2 = leap year
INTEGER, PARAMETER :: S_PER_YEAR_1=31536000,S_PER_YEAR_2=31622400,S_PER_DAY=86400,S_PER_HOUR=3600,S_PER_MIN=60
INTEGER, PARAMETER, DIMENSION(12) :: ACCUMULATED_DAYS_1=(/0,31,59,90,120,151,181,212,243,273,304,334/), & 
                                     ACCUMULATED_DAYS_2=(/0,31,60,91,121,152,182,213,244,274,305,335/)
INTEGER :: YEAR_COUNT

CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
WALL_CLOCK_SECONDS = 0._EB
DO YEAR_COUNT=2001,DATE_TIME(1)
   !Leap year if divisible by 4 but not 100 unless by 400 (1900 no, 1904  yes, 2000 yes)
   IF (MOD(YEAR_COUNT,4)==0 .AND. (MOD(YEAR_COUNT,100)/=0 .OR. MOD(YEAR_COUNT,400)==0)) THEN
      WALL_CLOCK_SECONDS = WALL_CLOCK_SECONDS + S_PER_YEAR_2
   ELSE
      WALL_CLOCK_SECONDS = WALL_CLOCK_SECONDS + S_PER_YEAR_1
   ENDIF
ENDDO 
IF (MOD(DATE_TIME(1),4)==0 .AND. (MOD(DATE_TIME(1),100)/=0 .OR. MOD(DATE_TIME(1),400)==0 )) THEN
   WALL_CLOCK_SECONDS = WALL_CLOCK_SECONDS + S_PER_DAY*(ACCUMULATED_DAYS_2(DATE_TIME(2))+DATE_TIME(3))
ELSE
   WALL_CLOCK_SECONDS = WALL_CLOCK_SECONDS + S_PER_DAY*(ACCUMULATED_DAYS_1(DATE_TIME(2))+DATE_TIME(3))
ENDIF
WALL_CLOCK_SECONDS = WALL_CLOCK_SECONDS +  S_PER_HOUR*DATE_TIME(5) + S_PER_MIN*DATE_TIME(6) + DATE_TIME(7)
WALL_CLOCK_TIME_THANKS_FDS    = WALL_CLOCK_SECONDS + DATE_TIME(8)*0.001_EB

! *****************************************************************************
END FUNCTION WALL_CLOCK_TIME_THANKS_FDS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE GET_CPU_TIME(T)
! *****************************************************************************
REAL(EB), INTENT(OUT) :: T
INTEGER(8) :: IT

CALL SYSTEM_CLOCK(IT)

T = REAL(IT,EB) / REAL(CLOCK_COUNT_RATE,EB)

! *****************************************************************************
END SUBROUTINE GET_CPU_TIME
! *****************************************************************************

! *****************************************************************************
INTEGER FUNCTION IJK_FROM_XYZ(COORD,NCELL,TARG,ITYPE)
! *****************************************************************************

INTEGER, INTENT(IN) :: NCELL, ITYPE
REAL(EB), INTENT(IN) :: TARG
REAL(EB), DIMENSION(:), INTENT(IN) :: COORD
REAL(EB) :: DIFF, MINDIFF
INTEGER :: I

IF (NCELL .EQ. 1) THEN
   IJK_FROM_XYZ = 1
ELSE
   MINDIFF = 9D9
   IF (ITYPE .EQ. 1) THEN
      DO I = 1, NCELL
         DIFF = ABS(COORD(I) - TARG)
         IF (DIFF .LT. MINDIFF) THEN
            MINDIFF = DIFF
            IJK_FROM_XYZ = I
         ENDIF
      ENDDO
   ELSE
      DO I = 1, NCELL
         DIFF = ABS(COORD(I) - TARG)
         IF (DIFF .LT. MINDIFF) THEN
            MINDIFF = DIFF
            IJK_FROM_XYZ = I
         ENDIF
      ENDDO
   ENDIF
ENDIF

! *****************************************************************************
END FUNCTION IJK_FROM_XYZ
! *****************************************************************************

! *****************************************************************************
END MODULE GPYRO_FUNCS
! *****************************************************************************
