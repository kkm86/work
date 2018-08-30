subroutine model_potential(S,r0,r,angconst,potential,potder)
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in)    :: S,r0,r(3),angconst(3)
  real(kind(1.d0)), intent(inout) :: potential(3),potder(3)

  potential = S*(cosh(r/r0))**(-2.d0)

  potder = -(2.d0*S*angconst*tanh(angconst*r/r0)*(cosh(angconst*r/r0))**(-2.d0))/r0
 
  return

end subroutine model_potential
