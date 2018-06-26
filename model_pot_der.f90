subroutine model_pot_der(S,r0,r,angconst,potential)
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in)    :: S,r0,r(3),angconst(3)
  real(kind(1.d0)), intent(inout) :: potential(3)

  !.. Local
  integer :: ii

  potential = -(2.d0*S*angconst*tanh(angconst*r/r0)*(cosh(angconst*r/r0))**(-2.d0))/r0
  

  
  return

end subroutine model_pot_der

