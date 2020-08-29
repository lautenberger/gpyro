MODULE SCE_FUNCS

USE PREC

IMPLICIT NONE

CONTAINS

! *****************************************************************************
subroutine SCE_sort(n,m,rb,ra,IDIM1,IDIM2)
! *****************************************************************************
!  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
!  BY W.H. PRESS ET AL., pp. 233-234
!
!  LIST OF VARIABLES
!     ra(.) = array to be sorted
!     rb(.,.) = arrays ordered corresponding to rearrangement of ra(.)
!     wk(.,.), iwk(.) = local varibles

IMPLICIT NONE
INTEGER :: N, M, I, J
INTEGER, INTENT(IN) :: IDIM1, IDIM2
REAL(EB) :: ra(IDIM1),rb(IDIM1,IDIM2),wk(IDIM1,IDIM2)
INTEGER :: iwk(IDIM1)

call SCE_indexx(n, ra, iwk)
do 11 i = 1, n
   wk(i,1) = ra(i)
11 continue

do 12 i = 1, n
   ra(i) = wk(iwk(i),1)
12 continue
do 14 j = 1, m
   do 13 i = 1, n
      wk(i,j) = rb(i,j)
   13 continue
14 continue

do 16 j = 1, m
   do 15 i = 1, n
      rb(i,j) = wk(iwk(i),j)
   15 continue
16 continue

CONTINUE
! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

! *****************************************************************************
subroutine SCE_sort1(n,ra)
! *****************************************************************************
!  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
!  BY W.H. PRESS ET AL., pp. 231
!
!  LIST OF VARIABLES
!     ra(.) = integer array to be sorted

IMPLICIT NONE
INTEGER :: N, L, IR, I, J

integer ra(N), rra

l = (n / 2) + 1
ir = n
10 continue
if (l .gt. 1) then
   l = l - 1
   rra = ra(l)
else
   rra = ra(ir)
   ra(ir) = ra(1)
   ir = ir - 1
   if (ir .eq. 1) then
      ra(1) = rra
      return
   end if
end if
i = l
j = l + l
20 if (j .le. ir) then
   if (j .lt. ir) then
      if (ra(j) .lt. ra(j + 1)) j = j + 1
   end if
   if (rra .lt. ra(j)) then
      ra(i) = ra(j)
      i = j
      j = j + j
   else
      j = ir + 1
   end if
   goto 20
end if

ra(i) = rra
goto 10

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

! *****************************************************************************
subroutine SCE_indexx(n, arrin, indx)
! *****************************************************************************
!  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
!      implicit REAL(EB) (a-h,o-z)

IMPLICIT NONE
INTEGER :: N, J, L, IR, INDXT, I
REAL(EB) :: Q
REAL(EB) :: arrin(n)
INTEGER ::  indx(n)

do 11 j = 1, n
   indx(j) = j
11 continue
l = (n / 2) + 1
ir = n
10 continue
if (l .gt. 1) then
   l = l - 1
   indxt = indx(l)
   q = arrin(indxt)
else
   indxt = indx(ir)
   q = arrin(indxt)
   indx(ir) = indx(1)
   ir = ir - 1
   if (ir .eq. 1) then
      indx(1) = indxt
      return
   end if
end if
i = l
j = l + l
20 if (j .le. ir) then
   if (j .lt. ir) then
      if (arrin(indx(j)) .lt. arrin(indx(j + 1))) j = j + 1
   end if
   if (q .lt. arrin(indx(j))) then
      indx(i) = indx(j)
      i = j
      j = j + j
   else
      j = ir + 1
   end if
   goto 20
end if
indx(i) = indxt
goto 10

! *****************************************************************************
end SUBROUTINE
! *****************************************************************************

! *****************************************************************************
SUBROUTINE RANDOM_ARRAY(RVEC,N) 
! *****************************************************************************
! Wrapper for old "RANLUX" call to use ran1 random number generator

!USE SCE_VARS

IMPLICIT NONE

INTEGER, INTENT(IN) :: N
REAL(EB), INTENT(OUT), DIMENSION(1:N) :: RVEC
REAL(EB) :: R

INTEGER :: I

DO I = 1, N
!   RVEC(I) = RAN1(SCE%ISEED1)
   CALL RANDOM_NUMBER(R)
   RVEC(I) = R
ENDDO

! *****************************************************************************
END SUBROUTINE RANDOM_ARRAY
! *****************************************************************************

! *****************************************************************************
REAL(EB) FUNCTION RANDOM_REAL()
! *****************************************************************************

IMPLICIT NONE

REAL(EB) :: R

CALL RANDOM_NUMBER(R)
RANDOM_REAL = R

! *****************************************************************************
END FUNCTION RANDOM_REAL
! *****************************************************************************

!! *****************************************************************************
! function ran1(idum) 
!! *****************************************************************************
!
!        use PREC 
!
!        implicit none  !note after use statement 
!
!        real(EB) ran1 
!        integer, intent(inout), optional :: idum 
!        real(EB) r(97),rm1,rm2 
!        integer, parameter :: m1=259200,ia1=7141,ic1=54773 
!        integer, parameter :: m2=134456,ia2=8121,ic2=28411 
!        integer, parameter :: m3=243000,ia3=4561,ic3=51349 
!        integer j 
!        integer iff,ix1,ix2,ix3 
!        data iff /0/ 
!        save ! corrects a bug in the original routine 
!        if(present(idum))then 
!          if (idum<0.or.iff.eq.0)then 
!            rm1=1.0_EB/m1 
!            rm2=1.0_EB/m2 
!            iff=1 
!            ix1=mod(ic1-idum,m1) 
!            ix1=mod(ia1*ix1+ic1,m1) 
!            ix2=mod(ix1,m2) 
!            ix1=mod(ia1*ix1+ic1,m1) 
!            ix3=mod(ix1,m3) 
!            do j=1,97 
!                ix1=mod(ia1*ix1+ic1,m1) 
!                ix2=mod(ia2*ix2+ic2,m2) 
!                r(j)=(real(ix1,EB)+real(ix2,EB)*rm2)*rm1 
!            enddo 
!            idum=1 
!          endif 
!        endif 
!        ix1=mod(ia1*ix1+ic1,m1) 
!        ix2=mod(ia2*ix2+ic2,m2) 
!        ix3=mod(ia3*ix3+ic3,m3) 
!        j=1+(97*ix3)/m3 
!        if(j>97.or.j<1)then 
!            write(*,*)' error in ran1 j=',j 
!            stop 
!        endif 
!        ran1=r(j) 
!        r(j)=(real(ix1,EB)+real(ix2,EB)*rm2)*rm1 
!        return 
!
!! *****************************************************************************
!     end function ran1  
!! *****************************************************************************
!
!
! *****************************************************************************
REAL(EB) function gasdev()
! *****************************************************************************

!  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
IMPLICIT NONE
INTEGER, SAVE :: ISET = 0
REAL(EB) :: V1, V2, R, FAC
REAL(EB), SAVE :: GSET

if (iset .eq. 0) then
!   1 v1 = (2. * ran1(idum)) - 1.
   1 v1 = (2. * RANDOM_REAL()) - 1.
!   v2 = (2. * ran1(idum)) - 1.
   v2 = (2. * RANDOM_REAL()) - 1.
   r = (v1 ** 2) + (v2 ** 2)
   if (r .ge. 1.) goto 1
   fac = sqrt(- ((2. * log(r)) / r))
   gset = v1 * fac
   gasdev = v2 * fac
   iset = 1
else
   gasdev = gset
   iset = 0
end if

! *****************************************************************************
end FUNCTION
! *****************************************************************************

END MODULE SCE_FUNCS