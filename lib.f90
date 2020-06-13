function IsInRange(i, imin, imax)
  implicit none
  logical :: isInRange
  integer(kind=4), intent(in) :: i, imin, imax
  isInRange = .false.
  if(i<=imax .and. i>=imin) isInRange = .true.
end function IsInRange

function Rmsd(n, a1, a2)
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: a1(n), a2(n)
  real(kind=fp_kind) :: rmsd
  integer(kind=4) :: i
  if(n==0)then
    rmsd = 0.d0
    return
  end if
  rmsd = 0.d0
  rmsd = sqrt(sum((a1-a2)**2)/n)
end function Rmsd

subroutine Weighted_Mean_and_StandardDev(n,weight,x,mean,stdDev)
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: weight(n)
  real(kind=fp_kind), intent(in) :: x(n)
  real(kind=fp_kind), intent(out) :: mean, stdDev
  mean = sum(weight(:)*x(:))/sum(weight(:))
  stdDev = sqrt(sum(weight(:)/sum(weight(:))*(x(:)-mean)**2))
end subroutine Weighted_Mean_and_StandardDev

subroutine Mean_and_StandardDev(n,x,mean,stdDev)
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: x(n)
  real(kind=fp_kind), intent(out) :: mean, stdDev
  mean = sum(x)/n
  stdDev = sqrt(sum((x-mean)**2)/n)
end subroutine Mean_and_StandardDev

subroutine Mean_and_StandardErr(n,x,mean,stdErr)
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: x(n)
  real(kind=fp_kind), intent(out) :: mean, stdErr
  mean = sum(x)/n
  stdErr = sqrt(sum((x-mean)**2))/n
end subroutine Mean_and_StandardErr

function cross_correlation(n,x,y)
   use precision_m
   implicit none
   integer(kind=4), intent(in) :: n
   real(kind=fp_kind), intent(in) :: x(n), y(n)
   real(kind=fp_kind) :: cross_correlation
   real(kind=fp_kind) :: xMean, xStdDev
   real(kind=fp_kind) :: yMean, yStdDev
   real(kind=fp_kind) :: numerator, denominator
   call Mean_and_StandardDev(n, x, xMean, xStdDev)
   call Mean_and_StandardDev(n, y, yMean, yStdDev)
   numerator = sum((x(:)-xMean)*(y(:)-yMean))/n
   denominator = xStdDev*yStdDev
   cross_correlation = numerator / denominator
end function cross_correlation

function IsOdd( i )
  implicit none
  integer( kind = 4 ) :: i
  logical :: IsOdd
  IsOdd = .true.
  if( i / 2 * 2 == i ) IsOdd = .false.
end function IsOdd

function IsEven( i ) 
  implicit none
  integer( kind = 4 ) :: i
  logical :: IsEven
  IsEven = .false.
  if( i / 2 * 2 == i ) IsEven = .true.
end function IsEven

!interface swap
!  subroutine swap_int( i1, i2 )
!    implicit none
!    integer( kind = 4 ), intent( in out ) :: i1, i2
!  end subroutine swap_int 
!
!  subroutine swap_real( f1, f2 )
!    implicit none
!    real( kind = 4 ), intent( in out ) :: f1, f2
!  end subroutine swap_real
!
!  subroutine swap_double( d1, d2 )
!    implicit none
!    real( kind = 8 ), intent( in out ) :: d1, d2
!  end subroutine swap_double
!
!end interface swap

subroutine swap_int( i1, i2 )
  implicit none
  integer( kind = 4 ), intent( in out ) :: i1, i2
  integer( kind = 4 ) :: itmp
  itmp = i1
  i1 = i2
  i2 = itmp
end subroutine swap_int

subroutine swap_real( f1, f2 )
  implicit none
  real( kind = 4 ), intent( in out ) :: f1, f2
  real( kind = 4 ) :: ftmp
  ftmp = f1 
  f1 = f2
  f2 = ftmp
end subroutine swap_real

subroutine swap_double( d1, d2 )
  implicit none
  real( kind = 8 ), intent( in out ) :: d1, d2
  real( kind = 8 ) :: dtmp
  dtmp = d1 
  d1 = d2
  d2 = dtmp
end subroutine swap_double

subroutine uppercase(s,length)
implicit none
integer(kind=4),intent(in)::length
character(len=length),intent(in out)::s
integer(kind=4)::i,j,k
do i = 1, length
  if(ichar(s(i:i))<=ichar('z').and.ichar(s(i:i))>=ichar('a'))then
    s(i:i)=char(ichar('A')+ichar(s(i:i))-ichar('a'))
  end if
end do
end subroutine uppercase

subroutine cross_product3(v1,v2,v3)
  use precision_m
  implicit none
  real(kind=fp_kind),intent(in)::v1(3),v2(3)
  real(kind=fp_kind),intent(out)::v3(3)
  integer(kind=4)::i,j,k,m,n
  v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
  v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
  v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
end subroutine cross_product3

function triple_product3(v1,v2,v3)
  use precision_m
  implicit none
  real(kind=fp_kind),intent(in)::v1(3),v2(3),v3(3)
  real(kind=fp_kind)::triple_product3
  real(kind=fp_kind)::v4(3)
  call cross_product3(v2,v3,v4)
  triple_product3=dot_product(v1,v4)
end function triple_product3

subroutine normalize(v)
  use precision_m
  implicit none
  real(kind=fp_kind),intent(in out)::v(3)
  real(kind=fp_kind)::norm
  norm=sqrt(dot_product(v,v))
  v=v/norm
end subroutine normalize

