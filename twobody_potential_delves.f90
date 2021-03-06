subroutine twobody_potential_delves(S,r0,rho,theta,phi,V)

  use constants
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in)    :: S,r0,rho,theta,phi
  real(kind(1.d0)), intent(inout) :: V

  !.. Local
  !.. r(1)=r23, r(2)=r31, r(3)=r12
  real(kind(1.0d0)) :: r(3)

  !.. External functions/variables
  real(kind(1.d0)) :: potential(3)
  
  r(3) = sqrt(2.d0)*rho*sin(phi)
  r(1) = rho*(0.5d0*(sin(phi)**2.d0) + 1.5d0*(cos(phi)**2.d0) - sqrt(3.d0)*sin(phi)*cos(phi)*cos(theta))**0.5d0
  r(2) = rho*(0.5d0*(sin(phi)**2.d0) + 1.5d0*(cos(phi)**2.d0) + sqrt(3.d0)*sin(phi)*cos(phi)*cos(theta))**0.5d0
 
  call model_potential(S,r0,r,potential)
 

  
  V = sum(potential)
  
  return

end subroutine twobody_potential_delves

  
