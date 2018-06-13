subroutine twobody_potential(S,r0,mass,rho,theta,phi,V,points)

  use constants
  
  implicit none

  !.. Input
  integer         , intent(in)    :: points
  real(kind(1.d0)), intent(in)    :: S,r0,mass(3),rho(3,points),theta,phi
  real(kind(1.d0)), intent(inout) :: V(3,points)

  !.. Local
  !.. r(1)=r23, r(2)=r31, r(3)=r12
  real(kind(1.0d0)) :: r(3,3), d(3), redmass, totmass, channelangle(3)
  integer :: ii

  !.. External functions/variables
  real(kind(1.d0)) :: potential(3,3)

  totmass = mass(1)+mass(2)+mass(3)
  redmass = (mass(1)*mass(2)*mass(3)/totmass)**0.5d0
  d = ((mass/redmass)*(1.d0-(mass/totmass)))**0.5d0

  channelangle(3) = 0.d0
  channelangle(1) = -2.d0*atan(mass(2)/redmass)
  channelangle(2) = 2.d0*atan(mass(1)/redmass)

  do ii = 1, points
     r(1,:) = d*rho(1,ii)*(((2.d0)**(-0.5d0)))*(1.d0+sin(theta)*cos(phi+channelangle))**0.5d0
     r(2,:) = d*rho(2,ii)*(((2.d0)**(-0.5d0)))*(1.d0+sin(theta)*cos(phi+channelangle))**0.5d0
     r(3,:) = d*rho(3,ii)*(((2.d0)**(-0.5d0)))*(1.d0+sin(theta)*cos(phi+channelangle))**0.5d0
     call model_potential(S,r0,r,potential)
     V(1,ii) = sum(potential(1,:))
     V(2,ii) = sum(potential(2,:))
     V(3,ii) = sum(potential(3,:))
  end do


  return

end subroutine twobody_potential

  
