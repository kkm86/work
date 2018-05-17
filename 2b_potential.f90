subroutine twobody_potential(S,b,m,rho,alpha,V)
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in) :: S,b,alpha,rho,m
  real(kind(1.d0)), intent(inout) :: V

  V = 2.d0*m*S*exp(-(sqrt(2.d0)*rho*sin(alpha)**2.d0)/(b**2.d0))

  return

end subroutine twobody_potential
  
