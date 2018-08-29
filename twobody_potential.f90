subroutine twobody_potential(rho,theta,phi,V,Vder,points)

  use constants
  
  implicit none

  !.. Input
  integer         , intent(in)    :: points
  real(kind(1.d0)), intent(in)    :: rho(points),theta,phi
  real(kind(1.d0)), intent(inout) :: V(points),Vder(points)

  !.. Local
  !.. r(1)=r23, r(2)=r31, r(3)=r12
  real(kind(1.0d0)) :: r(3), d(3), redmass, totmass, channelangle(3),angconst(3)
  integer :: ii

  !.. External functions/variables
  real(kind(1.d0)) :: potential(3),potder(3)

  totmass = mass(1)+mass(2)+mass(3)
  redmass = (mass(1)*mass(2)*mass(3)/totmass)**0.5d0
  d = ((mass/redmass)*(1.d0-(mass/totmass)))**0.5d0

  channelangle(3) = 0.d0
  channelangle(1) = -2.d0*atan(mass(2)/redmass)
  channelangle(2) = 2.d0*atan(mass(1)/redmass)

  angconst = d*(((2.d0)**(-0.5d0)))*(1.d0+sin(theta)*cos(phi+channelangle))**0.5d0

  do ii = 1, points
     r = d*rho(ii)*(((2.d0)**(-0.5d0)))*(1.d0+sin(theta)*cos(phi+channelangle))**0.5d0
     call model_potential(r,angconst,potential,potder)
     V(ii) = sum(potential)
     Vder(ii) = sum(potder)
  end do


  return

end subroutine twobody_potential

  
