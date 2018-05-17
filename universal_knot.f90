subroutine universal_knot(np,k,N,max,min,t)

  use constants

  implicit none

  !.. Input
  integer, intent(in) :: np,k,N
  real(kind(1.d0)), intent(in) :: max,min
  real(kind(1.d0)), intent(inout) :: t(np)

  !.. Local
  real(kind(1.0d0)) :: step_size
  integer :: ii

  !.. Setting up knot vector
  step_size = (max-min)/(N-1)
 
  do ii = 1, np
     if(ii .le. k) then
        t(ii) = min
     elseif(ii .le. np-k+1) then
        t(ii) = t(ii-1)+step_size
     else
        t(ii) = t(ii-1)
     end if
   !  print *, ii, t(ii)
  end do

  return

end subroutine universal_knot

