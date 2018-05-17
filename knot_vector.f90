subroutine knot_vector(np,k,NN,rho,t)

  use constants

  implicit none

  !.. Input
  integer, intent(in) :: np,k,NN
  real(kind(1.d0)), intent(in) :: rho
  real(kind(1.d0)), intent(inout) :: t(np)

  !.. Local
  real(kind(1.0d0)) :: alpha_max, alpha_min, step_size, small_step
  integer :: ii

  !.. Setting up knot vector
  alpha_max = PI/2.d0
  alpha_min = 0.d0
  step_size = (alpha_max-alpha_min)/(NN-k+3)
  small_step = 1.e-6

  do ii = 1, np
     if(ii .le. k) then
        t(ii) = 0.d0
     elseif(ii .le. k+2) then
        t(ii) = t(ii-1)+small_step
     elseif(ii .le. np-k) then
        t(ii) = t(ii-1)+step_size
     else
        t(ii) = alpha_max
     end if
    ! print *, ii, t(ii)
  end do

  return
end subroutine knot_vector
