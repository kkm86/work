subroutine twobody_potential(rho,theta,phi,V,points)

  use constants
  
  implicit none

  COMMON /twobp/  d(3), d2(3),redmass, totmass, channelangle(3)

  !.. Input
  integer         , intent(in)    :: points
  real(kind(1.d0)), intent(in)    :: rho(points),theta,phi
  real(kind(1.d0)), intent(inout) :: V(points)

  !.. Local
  !.. r(1)=r23, r(2)=r31, r(3)=r12
  real(kind(1.0d0)) :: r(3), d,d2, redmass, totmass, channelangle,angconst(3)
  integer :: ii

  !.. External functions/variables
  real(kind(1.d0)) :: potential(3),potder(3)


  angconst =d2*(1.d0+sin(theta)*cos(phi+channelangle))**0.5d0

  do ii = 1, points
     r = d2*rho(ii)*(1.d0+sin(theta)*cos(phi+channelangle))**0.5d0
     call model_potential(r,angconst,potential)
     V(ii) = sum(potential)
  end do


  return

end subroutine twobody_potential

  
