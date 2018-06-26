subroutine model_potential(S,r0,r,potential)
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in)    :: S,r0,r(3,3)
  real(kind(1.d0)), intent(inout) :: potential(3,3)

  !.. Local
  integer :: ii

  potential(1,:) = S*(cosh(r(1,:)/r0))**(-2.d0)
  potential(2,:) = S*(cosh(r(2,:)/r0))**(-2.d0)
  potential(3,:) = S*(cosh(r(3,:)/r0))**(-2.d0)
 

  
  return

end subroutine model_potential
