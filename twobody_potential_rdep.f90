subroutine twobody_potential_rdep(S,b,mass,rho,r_1,V)
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in) :: S,b,mass,rho,r_1
  real(kind(1.d0)), intent(inout) :: V

  !.. Local
  real(kind(1.0d0)) :: r

  V = -2.d0*mass*S*exp(-((r/b)**2.d0))
  
  return

end subroutine twobody_potential_rdep

  
