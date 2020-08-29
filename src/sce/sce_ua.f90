MODULE SCE_UA

USE PREC

IMPLICIT NONE
      
CONTAINS

! *****************************************************************************
subroutine SCE_OPTIMIZATION(IGS)
! *****************************************************************************
!  SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
!     -- Version 2.1
!
!  by QINGYUN DUAN
!  DEPARTMENT OF HYDROLOGY & WATER RESOURCES
!  UNIVERSITY OF ARIZONA, TUCSON, AZ 85721
!  (602) 621-9360, email: duan@hwr.arizona.edu
!
!  WRITTEN IN OCTOBER 1990.
!  REVISED IN AUGUST 1991
!  REVISED IN APRIL 1992
!
!  STATEMENT BY AUTHOR:
!  --------------------
!
!     This general purpose global optimization program is developed at
!     the Department of Hydrology & Water Resources of the University
!     of Arizona.  Further information regarding the SCE-UA method can
!     be obtained from Dr. Q. Duan, Dr. S. Sorooshian or Dr. V.K. Gupta
!     at the address and phone number listed above.  We request all
!     users of this program make proper reference to the paper entitled
!     'Effective and Efficient Global Optimization for Conceptual
!     Rainfall-runoff Models' by Duan, Q., S. Sorooshian, and V.K. Gupta,
!     Water Resources Research, Vol 28(4), pp.1015-1031, 1992.
!
!
!  LIST OF INPUT ARGUEMENT VARIABLES
!
!     a(.) = initial parameter set
!     bl(.) = lower bound on parameters
!     bu(.) = upper bound on parameters
!     nopt = number of parameters to be optimized
!
!
!  LIST OF SCE ALGORITHMIC CONTROL PARAMETERS:

!     ngs = number of complexes in the initial population
!     npg = number of points in each complex
!     npt = total number of points in initial population (npt=ngs*npg)
!     nps = number of points in a sub-complex
!     nspl = number of evolution steps allowed for each complex before
!         complex shuffling
!     mings = minimum number of complexes required, if the number of
!         complexes is allowed to reduce as the optimization proceeds
!     iseed = initial random seed
!     iprint = flag for controlling print-out after each shuffling loop
!         = 0, print information on the best point of the population
!         = 1, print information on every point of the population
!
!
!  CONVERGENCE CHECK PARAMETERS
!
!     maxn = max no. of trials allowed before optimization is terminated
!     kstop = number of shuffling loops in which the criterion value must
!         chang by the given percentage before optimization is terminated
!     pcento = percentage by which the criterion value must change in
!         given number of shuffling loops
!     ipcnvg = flag indicating whether parameter convergence is reached
!         (i.e., check if gnrng is less than 0.001)
!         = 0, parameter convergence not satisfied
!         = 1, parameter convergence satisfied
!
!
!  LIST OF LOCAL VARIABLES
!     x(.,.) = coordinates of points in the population
!     xf(.) = function values of x(.,.)
!     xx(.) = coordinates of a single point in x
!     cx(.,.) = coordinates of points in a complex
!     cf(.) = function values of cx(.,.)
!     s(.,.) = coordinates of points in the current simplex
!     sf(.) = function values of s(.,.)
!     bestx(.) = best point at current shuffling loop
!     bestf = function value of bestx(.)
!     worstx(.) = worst point at current shuffling loop
!     worstf = function value of worstx(.)
!     xnstd(.) = standard deviation of parameters in the population
!     gnrng = normalized geometric mean of parameter ranges
!     lcs(.) = indices locating position of s(.,.) in x(.,.)
!     bound(.) = bound on ith variable being optimized
!     ngs1 = number of complexes in current population
!     ngs2 = number of complexes in last population
!     iseed1 = current random seed
!     criter(.) = vector containing the best criterion values of the last
!         10 shuffling loops

USE SCE_VARS
USE SCE_OBJFUNC
USE SCE_FUNCS
USE SCE_IO

IMPLICIT NONE

INTEGER, INTENT(IN) :: IGS
      
INTEGER :: J, K1, K2, K, LPOS, LOOP !, NGS2 

REAL(EB) :: RAND

CALL SCE_PARTITION_INTO_COMPLEXES(IGS)

do loop = 1, SCE%nspl !  BEGIN INNER LOOP - RANDOM SELECTION OF SUB-COMPLEXES

   if (SCE%nps .eq. SCE%npg) then !  CHOOSE SUB-COMPLEX (nps points) PER LINEAR DISTRIBUTION
      do k = 1, SCE%nps
         SCE%lcs(k) = k
      end do
   end if

   IF (SCE%NPS .NE. SCE%NPG) THEN
!      rand = ran1(SCE%iseed1)
      rand = RANDOM_REAL()
      SCE%lcs(1) = 1 + dint(SCE%npg + 0.5 - dsqrt( (SCE%npg+.5)**2 - SCE%npg * (SCE%npg+1) * rand ))
         
      do k = 2, SCE%nps
!60       rand = ran1(SCE%iseed1)
60       rand = RANDOM_REAL() 
         lpos = 1 + dint(SCE%npg + 0.5 - dsqrt((SCE%npg+.5)**2 - SCE%npg * (SCE%npg+1) * rand ))
         do k1 = 1, k-1
           if (lpos .eq. SCE%lcs(k1)) go to 60
         end do
         SCE%lcs(k) = lpos
      end do

      call SCE_sort1(SCE%nps,SCE%lcs) !  ARRANGE THE SUB-COMPLEX IN ORDER OF INCEASING FUNCTION VALUE
   ENDIF

!  CREATE THE SUB-COMPLEX ARRAYS 
   do k = 1, SCE%nps
      do j = 1, SCE%NOPT
         SCE%s(k,j) = SCE%cx(SCE%lcs(k),j)
      end do
      SCE%sf(k) = SCE%cf(SCE%lcs(k))
   end do

!  USE THE SUB-COMPLEX TO GENERATE NEW POINT(S)
   call SCE_cce(SCE%NOPT,SCE%nps,SCE%s,SCE%sf,SCE%bl,SCE%bu,SCE%xnstd,SCE%icall,SCE%iseed1)

!  IF THE SUB-COMPLEX IS ACCEPTED, REPLACE THE NEW SUB-COMPLEX
!  INTO THE COMPLEX
   do k = 1, SCE%nps
      do j = 1, SCE%NOPT
         SCE%cx(SCE%lcs(k),j) = SCE%s(k,j)
      end do
      SCE%cf(SCE%lcs(k)) = SCE%sf(k)
   end do

   call SCE_sort(SCE%npg,SCE%NOPT,SCE%cx,SCE%cf,SCE%NPT,SCE%NOPT)

ENDDO !loop = 1, nspl

CONTINUE

!  REPLACE THE NEW COMPLEX INTO ORIGINAL ARRAY x(.,.)

!IF (GA%MPI .AND. GA%IRANK .NE. 0) THEN
!   SCE%x(:,:) = 0.
!   SCE%xF(:) = 0.
!ENDIF

do k1 = 1, SCE%npg
   k2 = (k1-1) * SCE%ngs1 + igs
   do j = 1, SCE%NOPT
      SCE%x(k2,j) = SCE%cx(k1,j)
   end do
   SCE%xf(k2) = SCE%cf(k1)
end do

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SCE_PARTITION_INTO_COMPLEXES(IGS)
! *****************************************************************************

USE SCE_VARS

IMPLICIT NONE

INTEGER, INTENT(IN) :: IGS
INTEGER :: J, K1, K2

!  ASSIGN POINTS INTO COMPLEXES
do k1 = 1, SCE%npg
   k2 = (k1-1) * SCE%ngs1 + igs
   do j = 1, SCE%NOPT
      SCE%cx(k1,j) = SCE%x(k2,j)
   end do
   SCE%cf(k1) = SCE%xf(k2)
end do

CONTINUE

! *****************************************************************************
END SUBROUTINE SCE_PARTITION_INTO_COMPLEXES
! *****************************************************************************


! *****************************************************************************
SUBROUTINE SCE_BEST_WORST
! *****************************************************************************

USE SCE_VARS

IMPLICIT NONE

INTEGER :: J

!  RECORD THE BEST AND WORST POINTS
do j = 1, SCE%NOPT
   SCE%bestx(j) = SCE%x(1,j)
   SCE%worstx(j) = SCE%x(SCE%npt1,j)
end do

SCE%bestf = SCE%xf(1)
SCE%worstf = SCE%xf(SCE%npt1)

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

! *****************************************************************************
subroutine SCE_cce(NOPT,nps,s,sf,bl,bu,xnstd,icall,iseed)
! *****************************************************************************
!  ALGORITHM GENERATE A NEW POINT(S) FROM A SUB-COMPLEX
!  SUB-COMPLEX VARIABLES

USE SCE_OBJFUNC

IMPLICIT NONE

REAL(EB), PARAMETER :: C1 = 0.8, C2=0.4
INTEGER :: NOPT, NPS, ICALL, ISEED, N, M, J, I, IBOUND
REAL(EB) :: ALPHA, BETA, FW, FNEW
REAL(EB) :: S(NPS,NOPT), SF(NPS), BU(NOPT), BL(NOPT), XNSTD(NOPT)
REAL(EB) :: sw(NOPT),sb(NOPT),ce(NOPT),snew(NOPT)
!  LIST OF LOCAL VARIABLES
!    sb(.) = the best point of the simplex
!    sw(.) = the worst point of the simplex
!    w2(.) = the second worst point of the simplex
!    fw = function value of the worst point
!    ce(.) = the centroid of the simplex excluding wo
!    snew(.) = new point generated from the simplex
!    iviol = flag indicating if constraints are violated
!          = 1 , yes
!          = 0 , no
!
!
!  EQUIVALENCE OF VARIABLES FOR READABILTY OF CODE
n     = nps
m     = nopt
alpha = 1.0
beta  = 0.5

!  IDENTIFY THE WORST POINT wo OF THE SUB-COMPLEX s
!  COMPUTE THE CENTROID ce OF THE REMAINING POINTS
!  COMPUTE step, THE VECTOR BETWEEN wo AND ce
!  IDENTIFY THE WORST FUNCTION VALUE fw
do j = 1, m
   sb(j) = s(1,j)
   sw(j) = s(n,j)
   ce(j) = 0.0
   do i = 1, n-1
      ce(j) = ce(j) + s(i,j)
   end do
   ce(j) = ce(j)/dble(n-1)
end do
fw = sf(n)

!  COMPUTE THE NEW POINT snew

!  FIRST TRY A REFLECTION STEP
do j = 1, m
   snew(j) = ce(j) + alpha * (ce(j) - sw(j))
end do

!  CHECK IF snew SATISFIES ALL CONSTRAINTS
call SCE_chkcst(nopt,snew,bl,bu,ibound)

!  snew IS OUTSIDE THE BOUND,
!  CHOOSE A POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
!  A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
!  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
if (ibound .ge. 1) call SCE_getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb,NOPT)

!  COMPUTE THE FUNCTION VALUE AT snew
fnew = SCE_OBJECTIVE_FUNCTION(nopt,snew)
icall = icall + 1
continue

!  COMPARE fnew WITH THE WORST FUNCTION VALUE fw

if (fnew .GT. fw) THEN !  fnew IS GREATER THAN fw, SO TRY A CONTRACTION STEP
   do j = 1, m
      snew(j) = ce(j) - beta * (ce(j) - sw(j))
   end do

   fnew = SCE_OBJECTIVE_FUNCTION(nopt,snew) !  COMPUTE THE FUNCTION VALUE OF THE CONTRACTED POINT
   icall = icall + 1
ENDIF

!  COMPARE fnew TO THE WORST VALUE fw
!  IF fnew IS LESS THAN OR EQUAL TO fw, THEN ACCEPT THE POINT AND RETURN

if (fnew .GT. fw) THEN

!  IF BOTH REFLECTION AND CONTRACTION FAIL, CHOOSE ANOTHER POINT
!  ACCORDING TO A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
!  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
   call SCE_getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb,NOPT)

!  COMPUTE THE FUNCTION VALUE AT THE RANDOM POINT
   fnew = SCE_OBJECTIVE_FUNCTION(nopt,snew)
   icall = icall + 1
ENDIF

!  REPLACE THE WORST POINT BY THE NEW POINT

do j = 1, m
   s(n,j) = snew(j)
end do
sf(n) = fnew

return

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

! *****************************************************************************
subroutine SCE_getpnt(nopt,idist,iseed,x,bl,bu,std,xi,IDIM2)
! *****************************************************************************

!     This subroutine generates a new point within feasible region
!     x(.) = new point
!     xi(.) = focal point
!     bl(.) = lower bound
!     bu(.) = upper bound
!     std(.) = standard deviation of probability distribution
!     idist = probability flag
!           = 1 - uniform distribution
!           = 2 - Gaussian distribution

USE SCE_OBJFUNC, ONLY: SCE_CHKCST
USE SCE_FUNCS

!USE SCE_OBJFUNC !, ONLY: SCE_CHKCST

IMPLICIT NONE

INTEGER :: NOPT, IDIST, ISEED, J, IBOUND
REAL(EB) :: RAND
INTEGER, INTENT(IN) :: IDIM2
REAL(EB) :: x(IDIM2),bl(IDIM2),bu(IDIM2),std(IDIM2),xi(IDIM2)

1 do j=1, nopt
!   2   if (idist .eq. 1) rand = ran1(iseed)
   2   if (idist .eq. 1) rand = RANDOM_REAL()
   if (idist .eq. 2) rand = gasdev()
   x(j) = xi(j) + std(j) * rand * (bu(j) - bl(j))

!     Check explicit constraints
   call SCE_chkcst(1,x(j),bl(j),bu(j),ibound)
   if (ibound .ge. 1) go to 2
end do

!     Check implicit constraints
call SCE_chkcst(nopt,x,bl,bu,ibound)
if (ibound .ge. 1) go to 1

return

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

! *****************************************************************************
subroutine SCE_parstt(npt,nopt,x,xnstd,bound,gnrng,ipcnvg,IDIM1,IDIM2)
! *****************************************************************************

!  SUBROUTINE CHECKING FOR PARAMETER CONVERGENCE
IMPLICIT NONE
REAL(EB), PARAMETER :: DELTA = 1.0d-20,peps=1.0d-3
REAL(EB) :: GNRNG, GSUM, XSUM1, XSUM2
INTEGER :: NPT, NOPT, IPCNVG, K, I
INTEGER, INTENT(IN) :: IDIM1, IDIM2

REAL(EB) :: x(IDIM1,IDIM2),xmax(IDIM2),xmin(IDIM2)
REAL(EB) ::  xmean(IDIM2),xnstd(IDIM2),bound(IDIM2)

!  COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
gsum = 0.d0
do k = 1, nopt
   xmax(k) = -1.0d+20
   xmin(k) = 1.0d+20
   xsum1 = 0.d0
   xsum2 = 0.d0
   do i = 1, npt
      xmax(k) = dmax1(x(i,k), xmax(k))
      xmin(k) = dmin1(x(i,k), xmin(k))
      xsum1 = xsum1 + x(i,k)
      xsum2 = xsum2 + x(i,k)*x(i,k)
   end do
   xmean(k) = xsum1 / dble(npt)
   xnstd(k) = (xsum2 / dble(npt) - xmean(k)*xmean(k))
   if (xnstd(k) .le. delta) xnstd(k) = delta
   xnstd(k) = dsqrt(xnstd(k))
   xnstd(k) = xnstd(k) / bound(k)
   gsum = gsum + dlog( delta + (xmax(k)-xmin(k))/bound(k) )
end do
gnrng = dexp(gsum/dble(nopt))

!  CHECK IF NORMALIZED STANDARD DEVIATION OF PARAMETER IS <= eps
ipcnvg = 0
if (gnrng .le. peps) then
   ipcnvg = 1
end if

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

! *****************************************************************************
subroutine SCE_comp(n,npt,ngs1,ngs2,npg,a,af,b,bf,IDIM1,IDIM2)
! *****************************************************************************

!  THIS SUBROUTINE REDUCE INPUT MATRIX a(n,ngs2*npg) TO MATRIX
!  b(n,ngs1*npg) AND VECTOR af(ngs2*npg) TO VECTOR bf(ngs1*npg)
IMPLICIT NONE
INTEGER :: N, NPT, NGS1, NGS2, NPG, IGS, IPG, K1, K2, I, J
INTEGER, INTENT(IN) :: IDIM1, IDIM2
REAL(EB) :: a(IDIM1,IDIM2),af(IDIM1),b(IDIM1,IDIM2),bf(IDIM1)

do igs=1, ngs1
   do ipg=1, npg
      k1=(ipg-1)*ngs2 + igs
      k2=(ipg-1)*ngs1 + igs
         do i=1, n
            b(k2,i) = a(k1,i)
         end do
      bf(k2) = af(k1)
   end do
end do

do j=1, npt
   do i=1, n
      a(j,i) = b(j,i)
   end do
   af(j) = bf(j)
end do

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

END MODULE SCE_UA