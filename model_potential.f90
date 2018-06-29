subroutine model_potential(S,r0,r,potential)
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in)    :: S,r0,r(3)
  real(kind(1.d0)), intent(inout) :: potential(3)

  potential = S*(cosh(r/r0))**(-2.d0)
 
  return

end subroutine model_potential
