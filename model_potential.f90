subroutine model_potential(r,angconst,potential,potder)

  use constants
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in)    :: r(3),angconst(3)
  real(kind(1.d0)), intent(inout) :: potential(3),potder(3)

  potential = Potential_depth*(cosh(r/r0))**(-2.d0)

  potder = -(2.d0*Potential_depth*angconst*tanh(angconst*r/r0)*(cosh(angconst*r/r0))**(-2.d0))/r0
 
  return

end subroutine model_potential
