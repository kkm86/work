subroutine alpha_knot(np2,k,NN2,alpha_1,ta)

  use constants

  implicit none

  !.. Input
  integer, intent(in) :: np2,k,NN2
  real(kind(1.d0)), intent(in) :: alpha_1
  real(kind(1.d0)), intent(inout) :: ta(np2)

  !.. Local
  real(kind(1.0d0)) :: alpha_max, alpha_min, step_size, small_step
  integer :: ii

  !.. Setting up knot vector
  alpha_max = (PI/2.d0)-Abs((PI/2.d0)-(PI/3.d0)-alpha_1)
  alpha_min = Abs((PI/3.d0)-alpha_1)
  step_size = (alpha_max-alpha_min)/(NN2-k+3)
  small_step = 1.e-6

  do ii = 1, np2
     if(ii .le. k) then
        ta(ii) = alpha_min
     elseif(ii .le. np2-k+1) then
        ta(ii) = ta(ii-1)+step_size
     else
        ta(ii) = ta(ii-1)
     end if
    ! print *, ii, ta(ii)
  end do

  return

end subroutine alpha_knot
