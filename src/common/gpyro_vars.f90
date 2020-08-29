! *****************************************************************************
MODULE GPYRO_VARS
! *****************************************************************************

USE PREC !Following FDS, single precison is FB and double precision is EB

! For system timing:
INTEGER(8) :: CLOCK_COUNT_RATE, CLOCK_COUNT_MAX 

INTEGER :: IGPYRO_TYPE ! =1 for standalone, =2 for FDS, =3 for GA

! These are the logical units used for i/o. 
INTEGER, PARAMETER :: LUPOINT   = 200, &
                      LUSUM     = 200, & !No longer used
                      LUPROF    = 250, &
                      LUTIME    = 432, &
                      LUINPUT   = 801, &
                      LUOUTPUT  = 802, &
                      LUDOTOUT  = 901, &
                      LUSMV     = 9000, &
                      LUCSV     = 9500
                      
! Misc.
REAL(EB), PARAMETER :: PI = 3.14159265358979323846D0
      
! "alpha" coefficients for omega function (MS thesis)
REAL(EB) :: AD(1:6)

! FDS related
REAL(EB), DIMENSION(1:6) :: IGNLOC
REAL(EB) :: IGNPOW,TIGNSTART,TIGNSTOP

! Initial condition TYPE.  
TYPE :: INITIAL_CONDITION_TYPE
   REAL(EB), POINTER, DIMENSION(:) :: YI0
   REAL(EB), POINTER, DIMENSION(:) :: YJ0
   REAL(EB) :: TMP_INITIAL
   REAL(EB) :: TMPG_INITIAL
   REAL(EB) :: P_INITIAL
END TYPE 

! AllBC type: 
TYPE :: ALLBC_TYPE
   INTEGER :: SURF_IDX
   REAL(EB) :: T
   REAL(EB) :: QE
   REAL(EB) :: HC
   REAL(EB) :: NHC   
   REAL(EB) :: TINF
   LOGICAL  :: RERADIATION 
   REAL(EB) :: TFIXED
   REAL(EB) :: MDOTPP
   REAL(EB) :: PRES
   REAL(EB) :: QEG
   REAL(EB) :: HCG
   REAL(EB) :: TINFG
   REAL(EB) :: TFIXEDG
   REAL(EB) :: HM
   REAL(EB), POINTER, DIMENSION(:) :: YJINF
END TYPE

!Geometry TYPE:
TYPE :: GEOM_TYPE
   INTEGER  :: IMESH
   REAL(EB) :: Z1
   REAL(EB) :: Z2
   REAL(EB) :: X1
   REAL(EB) :: X2
   REAL(EB) :: Y1
   REAL(EB) :: Y2
   INTEGER :: ICNUM
   INTEGER :: SURF_IDX(1:6)
   REAL(EB) :: HCR(1:6)
END TYPE

! Holds info about solid to solid or solid to solid plus gas reactions. 
! Read in from rxns worksheet
TYPE :: GPYRO_REACTION_TYPE
   CHARACTER(60) :: CFROM
   CHARACTER(60) :: CTO
   INTEGER       :: IFROM !Index from
   INTEGER       :: ITOS  !Index to (solid)
   REAL(EB)      :: Z
   REAL(EB)      :: E
   REAL(EB)      :: DHS
   REAL(EB)      :: DHV
   REAL(EB)      :: CHI
   REAL(EB)      :: ORDER
   REAL(EB)      :: ORDERO2
   INTEGER       :: IKINETICMODEL
   INTEGER       :: IO2TYPE
   REAL(EB)      :: M
   REAL(EB)      :: KCAT
   INTEGER       :: ICAT
   REAL(EB)      :: TCRIT
END TYPE
TYPE (GPYRO_REACTION_TYPE),    POINTER, DIMENSION (:) :: RXN
     
! Homogeneous gaseous reactions TYPE. Holds info about homogeneous gas reactions. 
! Read in from hgrxns worksheet. 
TYPE :: HG_REACTION_TYPE
   CHARACTER(60) :: CREACTANT1
   CHARACTER(60) :: CREACTANT2
   INTEGER       :: IREACTANT1
   INTEGER       :: IREACTANT2
   REAL(EB)      :: P
   REAL(EB)      :: Q
   REAL(EB)      :: B
   REAL(EB)      :: Z
   REAL(EB)      :: E
   REAL(EB)      :: DH
END TYPE
TYPE (HG_REACTION_TYPE), POINTER, DIMENSION (:) :: HGRXN 

! Solid property TYPE.  Holds, of course, solid properties. 
! This is read in from sprops worksheet.
TYPE :: SOLID_PROPERTY_TYPE
   CHARACTER(60), POINTER, DIMENSION (:) :: NAME
   REAL(EB), POINTER, DIMENSION (:)   :: K0Z
   REAL(EB), POINTER, DIMENSION (:)   :: NKZ
   REAL(EB), POINTER, DIMENSION (:)   :: R0
   REAL(EB), POINTER, DIMENSION (:)   :: NR
   REAL(EB), POINTER, DIMENSION (:)   :: C0
   REAL(EB), POINTER, DIMENSION (:)   :: NC
   REAL(EB), POINTER, DIMENSION (:)   :: EMIS
   REAL(EB), POINTER, DIMENSION (:)   :: KAPPA
   REAL(EB), POINTER, DIMENSION (:)   :: TMELT
   REAL(EB), POINTER, DIMENSION (:)   :: DHMELT
   REAL(EB), POINTER, DIMENSION (:)   :: SIGMA2MELT
   REAL(EB), POINTER, DIMENSION (:)   :: GAMMA
   REAL(EB), POINTER, DIMENSION (:)   :: PERMZ
   REAL(EB), POINTER, DIMENSION (:)   :: RS0
   REAL(EB), POINTER, DIMENSION (:)   :: PORE_DIAMETER
   REAL(EB), POINTER, DIMENSION (:)   :: K0X
   REAL(EB), POINTER, DIMENSION (:)   :: NKX
   REAL(EB), POINTER, DIMENSION (:)   :: PERMX
   REAL(EB), POINTER, DIMENSION (:)   :: K0Y
   REAL(EB), POINTER, DIMENSION (:)   :: NKY 
   REAL(EB), POINTER, DIMENSION (:)   :: PERMY
      INTEGER  :: NSSPEC
   INTEGER  :: NRXN

END TYPE
TYPE (SOLID_PROPERTY_TYPE) :: SPROP

! Gas property type.
! This is read in from gprops worksheet.
TYPE :: GAS_PROPERTY_TYPE
   INTEGER :: NGSPEC
   INTEGER :: NHGRXN
   INTEGER :: IBG
   INTEGER :: IO2
   REAL(EB) :: CPG 
   CHARACTER(60), POINTER, DIMENSION (:)   :: NAME
   REAL(EB), POINTER, DIMENSION(:) :: M
   REAL(EB), POINTER, DIMENSION(:) :: SIGMA
   REAL(EB), POINTER, DIMENSION(:) :: EPSOK
   REAL(EB), POINTER, DIMENSION(:) :: C0
   REAL(EB), POINTER, DIMENSION(:) :: NC
   REAL(EB), POINTER, DIMENSION(:,:) :: YIELDS !Heterogeneous
   REAL(EB), POINTER, DIMENSION(:,:) :: HGYIELDS !Homogeneous
END TYPE
TYPE (GAS_PROPERTY_TYPE) :: GPROP

! Boundary condition type - holds all "calculated" quantities for each BC.
! In standalone and genetic algorithm implementations, there is only a single
! one-dimensional boundary condition, but when coupled to FDS there may be 
! hundreds or thousands of separate boundary conditions, and using this 
! "TYPE" is a convenient way to track each one. 
!
TYPE  :: GPYRO_MESH_TYPE

   REAL(EB) :: INITIAL_MASS

   REAL(EB) :: GX   !x-component of gravity vector
   REAL(EB) :: GY   !y-component of gravity vector
   REAL(EB) :: GZ   !z-component of gravity vector   

   REAL(EB) :: X0=0D0   !x absolute coordinate for FDS
   REAL(EB) :: Y0=0D0   !y absolute coordinate for FDS
   REAL(EB) :: Z0=0D0   !z absolute coordinate for FDS
   
   REAL(EB) :: ZDIM=0D0
   REAL(EB) :: XDIM=0D0
   REAL(EB) :: YDIM=0D0

   INTEGER :: NCELLZ
   INTEGER :: NCELLX
   INTEGER :: NCELLY

   REAL(EB), DIMENSION(1:3) :: MESH_CENTROIDS=0D0, MESH_EXTENTS=0D0, MESH_LL=0D0 !, MESH_G

   LOGICAL, POINTER, DIMENSION (:,:,:) :: NEEDSBCT, NEEDSBCB, NEEDSBCW, NEEDSBCE, NEEDSBCN, NEEDSBCS
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCT
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCB
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCW
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCE
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCN
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCS

!   REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRT !contact resistance (top)
!   REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRB !contact resistance (bottom)
!   REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRW !contact resistance (west)
!   REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRE !contact resistance (east)
!   REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRS !contact resistance (south)
!   REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRN !contact resistance (north)

   REAL(EB), POINTER, DIMENSION (:,:,:) :: TP  !Solid T at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: TPN !Solid T at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HP  !Solid enthalpy at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HPN !Solid enthalpy at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: TG  !Gas T at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: TGN !Gas T at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HG  !Gas enthalpy at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HGN !Gas enthalpy at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RG  !Gas density at cell P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RGN !Gas density at cell P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: P   !Pressure at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: PN  !Pressure at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RP  !Density at point P 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RPN !Density at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RSP  !Solid density at point P 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RSPN !Solid density at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RDLTZ  !Old rho*DZ
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RDLTZN !New rho*DZ
   REAL(EB), POINTER, DIMENSION (:,:,:) :: POROSS   !Porosity at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: POROSSN  !Porosity at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RS0  !Solid density @ t=0
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HCV  !Volumetric heat transfer coefficient
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RE  !Reynolds number
   REAL(EB), POINTER, DIMENSION (:,:,:) :: NU  !Nusselt number

   REAL(EB), POINTER, DIMENSION (:,:,:) :: CPS !Specific heat of the solid
            
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHP !Positive part of solid h source
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHM !Negative part of solid h source
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHGP !Positive part of gas h source
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHGM !Negative part of gas h source
            
!   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGPSIDZ  !gas density * porosity * grid size
!   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGPSIDZN !Same, but new

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: D12 !Diffusion coefficient

   INTEGER, POINTER, DIMENSION (:,:,:)   :: IMASK !Geometry mask (=0 for blocked off, =1 for not blocked off)
   INTEGER, POINTER, DIMENSION (:,:,:)   :: BCDIR !Boundary condition direction for use with IMASK (same convention as IOR in FDS)

!   INTEGER, POINTER, DIMENSION (:,:) :: IZBCT, IZBCB !, IXBCW, IXBCE, IYBCS, IYBCN

   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: YI  !Condensed-species mass fractions
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: YIN !Same, but "new"
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: HI  !Enthalpy of solid species i (hi)
         
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RRT !Temperature part of reaction rate
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RRY !Mass fraction part of reaction rate
         
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: YJG  !Gas-phase mass fractions
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: YJGN !Same, but "new"

   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZ0  !rho*Yi*Dz (at t = 0)
         
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZP  !rho*Yi*Dz at point P
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZPN !same, but "new"
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZSIGMA  !rho*Yi*Dz summation
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZSIGMAN !rho*Yi*Dz0 summation 

   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: UNREACTEDNESS ! 1 minus conversion
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: XI  !Volume fraction of solid species i
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: XIN !Same, but next 
         
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: OMSOLIDFRAC !One minus solid fraction
!   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: SOLIDFRAC   !Solid fraction
         
! OMEGASFBK is Omega (solid) - formation   of condensed species B from reaction k
! OMEGASDAK is Omega (solid) - destruction of condensed species A from reaction k
! OMEGASFGK is Omega (solid) - formation   of all gases           from reaction k
! OMEGASFG  is Omega (solid) - total formation rate of all gases (summed over all k)
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: OMEGASFBK
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: OMEGASDAK
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: OMEGASFGK
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: OMEGASFG

! OMEGASFJK is Omega (solid) - formation   of gaseous species J from reaction k
! OMEGASDJK is Omega (solid) - destruction of gaseous species J from reaction k
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: OMEGASFJK
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: OMEGASDJK
         
! SOMEGA stores the total formaton/destruction rates of solid (condensed) 
! species from condensed-phase reactions
! GOMEGA stores the total formation/destruction rates of gas-phase
! species from condensed-phase reaction
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: SOMEGA
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: GOMEGA

! HGRR is homogeneous gaseous reaction rate
! OMEGAGDJL is Omega (gaseous) - destruction of gaseous species J from reaction l (L)	   
! HGOMEGA stores the formation/destruction of gaseous species from 
! homogeneous gaseous reactions
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: HGRR 
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: OMEGAGDJL
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: HGOMEGA

   REAL(EB), POINTER, DIMENSION (:,:,:)     :: QSG !Heat transfer from solid to gas
   LOGICAL , POINTER, DIMENSION (:,:,:)     :: CONSUMED !Is a cell completely consumed?

   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: MDOTPPZ !Mass flux of each species in the z-direction
                 
   REAL(EB), POINTER, DIMENSION(:,:) :: DELTA    !Thickness
   REAL(EB), POINTER, DIMENSION(:,:) :: DELTAN   !Thickness (new)

!New section for 3D solvers:
   REAL(EB), POINTER, DIMENSION (:) :: Z ! Z coordinate
   REAL(EB), POINTER, DIMENSION (:) :: X ! X coordinate
   REAL(EB), POINTER, DIMENSION (:) :: Y ! Y coordinate

   REAL(EB), POINTER, DIMENSION (:,:,:) :: FT 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FB 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FE 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FW 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FN 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FS 

   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTZ  ! DLTZ  - Delta z 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTX  ! DLTX  - Delta x 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTY  ! DLTY  - Delta y 
   
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTZN ! DLTZN - Delta z (new), i.e. at t + Dt
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTXN ! DLTXN - Delta x (new), i.e. at t + Dt
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTYN ! DLTYN - Delta y (new), i.e. at t + Dt

   REAL(EB), POINTER, DIMENSION (:,:,:) :: DZT  !Dz top
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DZB  !Dz bottom	   
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DXW  !Dx west
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DXE  !Dx east	   
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DYN  !Dy north
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DYS  !Dy south	   

   REAL(EB), POINTER, DIMENSION (:,:,:) :: DXDYDZ !Dx*Dy*Dz
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DXDY   !Dx*Dy
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DXDZ   !Dx*Dz
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DYDZ   !Dy*Dz

   REAL(EB), POINTER, DIMENSION (:,:,:) :: KZ   !Thermal conductivity in z direction
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KX   !Thermal conductivity in x direction
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KY   !Thermal conductivity in y direction

   REAL(EB), POINTER, DIMENSION (:,:,:) :: KRZ  !Radiative conductivity in z direction
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KRX  !Radiative conductivity in x direction
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KRY  !Radiative conductivity in y direction

   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCT !(k/c) at top interface
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCB !(k/c) at bottom interface
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCE !(k/c) at east interface 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCW !(k/c) at west interface
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCN !(k/c) at north interface 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCS !(k/c) at south interface

   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIT !(k/c)*hi at top interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIB !(k/c)*hi at bottom interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIE !(k/c)*hi at east interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIW !(k/c)*hi at west interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIN !(k/c)*hi at north interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIS !(k/c)*hi at south interface

   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHYIDTDMAZ         
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHYIDTDMAX
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHYIDTDMAY

   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHYIAPZ
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHYIAPX
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHYIAPY

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMZ    !Permeability in z-direction at point P
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMX    !Permeability in x-direction at point P
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMY    !Permeability in y-direction at point P

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUT !Permeability/viscosity at top interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUB !Permeability/viscosity at bottom interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUE !Permeability/viscosity at east interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUW !Permeability/viscosity at west interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUN !Permeability/viscosity at north interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUS !Permeability/viscosity at south interface

   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYT !z-direction mass flux (top)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYB !z-direction mass flux (bottom)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYE !x-direction mass flux (east)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYW !x-direction mass flux (west)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYN !y-direction mass flux (north)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYS !y-direction mass flux (south)

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDT !porosity * gas density * D (top)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDB !porosity * gas density * D (bottom)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDE !porosity * gas density * D (east)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDW !porosity * gas density * D (west)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDN !porosity * gas density * D (north)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDS !porosity * gas density * D (south)

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: BUOYANCYZ
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: BUOYANCYX
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: BUOYANCYY

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNT
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNB
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNE
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNW
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNN
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNS

!   REAL(EB), POINTER, DIMENSION(:,:,:) :: GPYRO_CONV_FLUX,GPYRO_DIFF_FLUX,GPYRO_MASS_FLUX
   
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RESIDUAL_TMP, RESIDUAL_P, RESIDUAL_HG 
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: RESIDUAL_YIS, RESIDUAL_YJG
   
!   REAL(EB) :: TLAST !Last time at which GPYRO_PYROLYSIS was called

   REAL(EB) :: TLAST_TMP
   REAL(EB) :: TLAST_YIS
   REAL(EB) :: TLAST_YJG
   REAL(EB) :: TLAST_P
   REAL(EB) :: TLAST_HG
   
   TYPE (GPYRO_BC_TYPE), POINTER, DIMENSION(:,:,:,:) :: GPYRO_BOUNDARY_CONDITION

! Variables for orientation: 
   LOGICAL :: ORIENTATION_FILE_EXISTS
   CHARACTER(200) :: ORIENTATION_FILE
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: ORI, K_TENSOR

!WORK ARRAYS:
   REAL(EB), POINTER, DIMENSION(:,:,:) :: RWORK01, RWORK02, RWORK03, RWORK04, RWORK05, RWORK06, RWORK07, RWORK08, RWORK09, &
                                          RWORK10, RWORK11, RWORK12, RWORK13, RWORK14, RWORK15, RWORK16, RWORK17, RWORK18, &
                                          RWORK19, RWORK20, RWORK21, RWORK22, RWORK23, RWORK24, RWORK25, RWORK26, RWORK27, &
                                          RWORK28, RWORK29, RWORK30, RWORK31, RWORK32, RWORK33, RWORK34, RWORK35, RWORK36, &
                                          RWORK37, RWORK38, RWORK39, RWORK40, RWORK41 , RWORK42, RWORK43, RWORK44, RWORK45, & 
                                          RWORK46, RWORK47, RWORK48, RWORK49, RWORK50, RWORK51, RWORK52

   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: RWORK100, RWORK101, RWORK102, RWORK103, RWORK104, RWORK105, RWORK106, RWORK110, &
                                             RWORK111, RWORK112, RWORK120, RWORK121

   LOGICAL, POINTER, DIMENSION(:,:,:) :: LWORK01

   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:)   :: CONVERGED,CONVERGED_TMP,CONVERGED_HG,CONVERGED_P
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: CONVERGED_YIS,CONVERGED_YJG
   
   LOGICAL :: SMOKEVIEW_FILE_OPENED_ALREADY = .FALSE. 

END TYPE

TYPE (GPYRO_MESH_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: GPM
TYPE (GPYRO_MESH_TYPE), POINTER :: G

TYPE :: GPYRO_BC_TYPE
   REAL(EB) :: QE
   REAL(EB) :: QENET  
   REAL(EB) :: HC0
   REAL(EB) :: NHC   
   REAL(EB) :: TINF
   REAL(EB) :: TFIXED
   REAL(EB) :: HFIXED
   REAL(EB) :: PRES
   REAL(EB) :: MFLUX
   REAL(EB) :: HM0
   REAL(EB) :: QEG
   REAL(EB) :: HC0G
   REAL(EB) :: TINFG
   REAL(EB) :: TFIXEDG
   REAL(EB) :: HFIXEDG
   REAL(EB), POINTER, DIMENSION(:) :: YJINF
   LOGICAL  :: RERAD
! The remaining variables are for fds coupling:
   REAL(EB) :: EMISSIVITY
   REAL(EB) :: TMP_F
   REAL(EB) :: TMP_F_OLD
   REAL(EB) :: QRADOUT
   REAL(EB) :: QCONF
   REAL(EB), POINTER, DIMENSION(:) :: MASSFLUX
END TYPE

!TYPE (GPYRO_BC_TYPE), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:,:) :: GPYRO_BOUNDARY_CONDITION
TYPE (GPYRO_BC_TYPE), POINTER, DIMENSION(:,:,:,:) :: GPBCP

TYPE :: GPYRO_TO_FDS_COUPLING_TYPE
   INTEGER , POINTER, DIMENSION(:) :: IMESH_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IZ_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IX_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IY_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IOR_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IOR_FDS
   INTEGER , POINTER, DIMENSION(:) :: IMESH_FDS
   INTEGER , POINTER, DIMENSION(:) :: I_FDS
   INTEGER , POINTER, DIMENSION(:) :: J_FDS
   INTEGER , POINTER, DIMENSION(:) :: K_FDS
   INTEGER , POINTER, DIMENSION(:) :: IW_FDS
   INTEGER , POINTER, DIMENSION(:) :: RATIO
   REAL(EB), POINTER, DIMENSION(:) :: DIFF
   REAL(EB), POINTER, DIMENSION(:) :: XDIFF
   REAL(EB), POINTER, DIMENSION(:) :: YDIFF
   REAL(EB), POINTER, DIMENSION(:) :: ZDIFF   
   REAL(EB), POINTER, DIMENSION(:) :: Z_GPYRO
   REAL(EB), POINTER, DIMENSION(:) :: X_GPYRO
   REAL(EB), POINTER, DIMENSION(:) :: Y_GPYRO
   REAL(EB), POINTER, DIMENSION(:) :: X_FDS
   REAL(EB), POINTER, DIMENSION(:) :: Y_FDS
   REAL(EB), POINTER, DIMENSION(:) :: Z_FDS
   TYPE(GPYRO_BC_TYPE), POINTER, DIMENSION(:) :: GPBC
END TYPE
TYPE (GPYRO_TO_FDS_COUPLING_TYPE), ALLOCATABLE, DIMENSION (:) :: GPYRO2FDS

TYPE :: GPYRO_TO_GPYRO_INTERFACE
   INTEGER :: NINTERFACES
   INTEGER , POINTER, DIMENSION(:) :: IMATCH
   INTEGER , POINTER, DIMENSION(:) :: IMESH1
   INTEGER , POINTER, DIMENSION(:) :: IMESH2
   INTEGER , POINTER, DIMENSION(:) :: IZ1
   INTEGER , POINTER, DIMENSION(:) :: IZ2
   INTEGER , POINTER, DIMENSION(:) :: IX1
   INTEGER , POINTER, DIMENSION(:) :: IX2
   INTEGER , POINTER, DIMENSION(:) :: IY1
   INTEGER , POINTER, DIMENSION(:) :: IY2
   INTEGER , POINTER, DIMENSION(:) :: IOR1
   INTEGER , POINTER, DIMENSION(:) :: IOR2
   REAL(EB), POINTER, DIMENSION(:) :: DLT1
   REAL(EB), POINTER, DIMENSION(:) :: DLT2
   REAL(EB), POINTER, DIMENSION(:) :: TMP1
   REAL(EB), POINTER, DIMENSION(:) :: TMP2
   REAL(EB), POINTER, DIMENSION(:) :: K1
   REAL(EB), POINTER, DIMENSION(:) :: K2
   REAL(EB), POINTER, DIMENSION(:) :: Q
END TYPE
TYPE (GPYRO_TO_GPYRO_INTERFACE) :: G2GIFACE,IFACETEMP !GPYRO TO GPYRO INTERFACE

INTEGER, ALLOCATABLE, DIMENSION (:) :: NGPYRO_CELLS_NEEDING_BCS

TYPE :: GPYRO_GENERAL_TYPE
   INTEGER :: PROCESS_GPYRO(1:100) = 0, MYID
   LOGICAL :: USE_MPI

   REAL(EB), DIMENSION(1:100,1:3) :: MESH_CENTROIDS, MESH_EXTENTS, MESH_LL 

   INTEGER :: NUM_GPYRO_MESHES = 1
   INTEGER :: NIC !Number of initial conditions
   TYPE(INITIAL_CONDITION_TYPE), POINTER, DIMENSION(:) :: INITIAL_CONDITIONS

   INTEGER :: NSURF_IDX
   TYPE(ALLBC_TYPE), POINTER, DIMENSION(:) :: ALLBC

   INTEGER :: NOBST
   REAL(EB), POINTER, DIMENSION(:) :: ZDIM, XDIM, YDIM
   INTEGER, POINTER, DIMENSION(:) :: NCELLZ, NCELLX, NCELLY
   CHARACTER(100), POINTER, DIMENSION(:) :: GEOMETRY_FILE
   INTEGER, POINTER, DIMENSION(:,:) :: DEFAULT_SURF_IDX
   INTEGER, POINTER, DIMENSION(:) :: DEFAULT_IC
   REAL(EB), POINTER, DIMENSION(:) :: OFFSETZ, OFFSETX, OFFSETY
   TYPE(GEOM_TYPE), POINTER, DIMENSION(:) :: GEOM

   REAL(EB) :: GX   !x-component of gravity vector
   REAL(EB) :: GY   !y-component of gravity vector
   REAL(EB) :: GZ   !z-component of gravity vector
      
   REAL(EB) :: GPYRO_WALL_CLOCK_START

! Quantities read in from GPYRO_DUMP namelist group
   CHARACTER(60) :: CASENAME
   LOGICAL  :: DUMP_ENERGY_BALANCE
   INTEGER  :: N_POINT_QUANTITIES
   INTEGER  :: N_PROFILE_QUANTITIES
   INTEGER  :: N_SMOKEVIEW_QUANTITIES
   REAL(EB) :: DTDUMP_GA  
   REAL(EB) :: DTDUMP_POINT
   REAL(EB) :: DTDUMP_PROFILE
   REAL(EB) :: DTDUMP_SMOKEVIEW
   REAL(EB) :: REDUCED_DTDUMP
   REAL(EB) :: TMP_REDUCED_DTDUMP
   REAL(EB) :: DTMIN_KILL
   
   CHARACTER(60), POINTER, DIMENSION(:) :: POINT_QUANTITY,PROFILE_QUANTITY,SMOKEVIEW_QUANTITY
   REAL(EB), POINTER, DIMENSION(:) :: POINT_Z,POINT_X,POINT_Y,PROFILE_COORD1,PROFILE_COORD2,SMOKEVIEW_LOCATION
   INTEGER, POINTER, DIMENSION(:) :: POINT_IZ, POINT_IX, POINT_IY, POINT_IMESH, PROFILE_IMESH, SMOKEVIEW_IMESH, SMOKEVIEW_ICELL
   CHARACTER(1), POINTER, DIMENSION(:) :: PROFILE_DIRECTION
   CHARACTER(2), POINTER, DIMENSION(:) :: SMOKEVIEW_PLANE
   INTEGER, POINTER, DIMENSION(:) :: PROFILE_ISKIP,PROFILE_IX,PROFILE_IY,PROFILE_IZ,POINT_QUANTITY_INDEX,PROFILE_QUANTITY_INDEX, &
                                     SMOKEVIEW_QUANTITY_INDEX
               
   REAL(EB) :: TAMB !Ambient temperature
   
   REAL(EB) :: RHOINF
         
   INTEGER  :: NCASES  !Number of cases to run
   REAL(EB) :: DT0     !Initial timestep 
   REAL(EB) :: DTNEXT  !Next timestep
   REAL(EB) :: P0      !Background pressure
   
   REAL(EB) :: MAXTMP !Maximum temperature (gas or solid)
   
   LOGICAL :: NAN
   REAL(EB) :: POSINF
   REAL(EB) :: NEGINF
   
   REAL(EB) :: EPS
   REAL(EB) :: VHLC
         
   REAL(EB) :: HCV
   REAL(EB) :: NU_A
   REAL(EB) :: NU_B
   REAL(EB) :: NU_C
         
   LOGICAL  :: FDSMODE
   REAL(EB) :: FDS_HEAT_OF_COMBUSTION
   INTEGER  :: FDS_MATL_VER
   LOGICAL  :: CONSTANT_DHVOL
   LOGICAL  :: FULL_QSG
   LOGICAL  :: GASES_PRODUCED_AT_TSOLID
   LOGICAL  :: DUMP_INTERMEDIATE_TIMINGS
   LOGICAL  :: DUMP_DETAILED_CONVERGENCE
         
   LOGICAL  :: NOCONSUMPTION
   REAL(EB) :: TMPTOL
   REAL(EB) :: HTOL
   REAL(EB) :: YITOL
   REAL(EB) :: PTOL
   REAL(EB) :: YJTOL         
   REAL(EB) :: HGTOL
   REAL(EB) :: EPS_YIS
   REAL(EB) :: EPS_YJG

   REAL(EB) :: OBST_SCALING_FACTOR
   REAL(EB) :: ADD_TO_OBST_XB, ADD_TO_OBST_X, ADD_TO_OBST_Y, ADD_TO_OBST_Z
   LOGICAL  :: ANISOTROPIC_SPECIES(1:50)

   INTEGER,  POINTER, DIMENSION(:)  :: ICASE
   INTEGER,  POINTER, DIMENSION(:)  :: IMESH
   REAL(EB), POINTER, DIMENSION(:)  :: TSTOP

   LOGICAL,  POINTER, DIMENSION(:) :: ZEROD !Is this a 0D simulation
   REAL(EB), POINTER, DIMENSION(:) :: BETA !Heating rate for TGA
! Logicals 
   LOGICAL  :: THERMAL_EQUILIBRIUM
   LOGICAL  :: SOLVE_GAS_ENERGY
   LOGICAL  :: SOLVE_PRESSURE
   LOGICAL  :: SOLVE_GAS_YJ
   LOGICAL  :: EXPLICIT_T
   LOGICAL  :: USE_TOFH_NEWTON
   LOGICAL  :: SHYI_CORRECTION
   LOGICAL  :: BLOWING
   LOGICAL  :: CONVENTIONAL_RXN_ORDER
   LOGICAL  :: KOZENY_CARMAN
   LOGICAL  :: USE_ANISOTROPIC_SOLID_ENTHALPY_SOLVER
   LOGICAL  :: GEOMETRY_IS_UPSIDE_DOWN

         
! Globals 
   REAL(EB) :: TREF
   REAL(EB) :: TDATUM ! temperature datum for enthalpy
   INTEGER  :: NTDMA_ITERATIONS,NTDMAITS !For TDMA algorithm
   INTEGER  :: NSSPECIESITERNS,NCONTINUITYITERNS,NCOEFF_UPDATE_SKIP
   REAL(EB) :: ALPHA
   REAL(EB) :: ALPHA_YIS
   REAL(EB) :: ALPHA_YJG
   REAL(EB) :: ALPHA_H
   REAL(EB) :: ALPHA_P
   REAL(EB) :: ALPHA_HG   
   REAL(EB) :: DT
   INTEGER  :: CONV_DIFF_SCHEME
   REAL(EB) :: TORTUOSITY_FACTOR
   LOGICAL  :: USE_TORTUOSITY_FACTOR_FOR_FLUX
   LOGICAL  :: USE_CONSTANT_HM0
   REAL(EB) :: CONSTANT_HM0

! FDS related parameters
   REAL(EB) :: YY0(1:5)
   REAL(EB) :: TDUMPLAST_GA
   REAL(EB), DIMENSION(1:100) :: TDUMPLAST_SMOKEVIEW, TDUMPLAST_POINT, TDUMPLAST_PROFILE, TDUMPLAST_CSVDUMP !Supports up to 100 meshes
 
! Timing
   REAL(EB) :: TUSED(0:28)=0.

! Newton iteration stuff
   REAL(EB), DIMENSION(1:100) :: NEWTON_A,NEWTON_B,NEWTON_C,NEWTON_D,NEWTON_E,NEWTON_F,NEWTON_G

! Test for fds coupling
   REAL(EB), DIMENSION(1:100) :: DTNEXT_ARR = 9D9

! Used for multi-mesh i/o
   INTEGER, DIMENSION(1:100) :: FIRST_SF_IN_MESH, FIRST_PROF_IN_MESH

! Gpyro to Gpyro coupling
   REAL(EB) :: GPYRO_TO_GPYRO_TOLERANCE

! Data transfer
   LOGICAL  :: CSVDUMP
   REAL(EB) :: DT_CSVDUMP

! Solver control
   CHARACTER(4) :: SWEEP_DIRECTION
   LOGICAL :: USE_SOLID_ENTHALPY_SOLVER_TEST, FIX_DOMAIN_TEMPERATURE

END TYPE

! This is a terrible name for a variable that is used so much, but
! below GPG stands for "GPyro, General". For the longest time, this variable was 
! called FG for "FIST, general" reflecting that this model was originally developed 
! as part of the FIST project...

TYPE (GPYRO_GENERAL_TYPE), TARGET, SAVE :: GPG

! Air properties table
REAL(EB), DIMENSION(1:26,0:8) :: AIR_PROPS_TABLE

!LOGICAL, ALLOCATABLE, DIMENSION(:,:,:)   :: CONVERGED,CONVERGED_TMP,CONVERGED_HG,CONVERGED_P
!LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: CONVERGED_YIS,CONVERGED_YJG

LOGICAL, DIMENSION (0:100) :: YISCONVERGED=.FALSE., ALLYISCONVERGED=.FALSE.
LOGICAL, DIMENSION (0:100) :: YJGCONVERGED=.FALSE., ALLYJGCONVERGED=.FALSE.

! *****************************************************************************
END MODULE GPYRO_VARS
! *****************************************************************************

