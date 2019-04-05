subroutine setup

  use constants
  
  implicit none

  COMMON /twobp/  d(3),d2(3), redmass, totmass, channelangle(3)

  !.. Local
  !.. r(1)=r23, r(2)=r31, r(3)=r12
  integer::ii
  real(kind(1.0d0)) :: d,d2, redmass, totmass, channelangle

  !.. External functions/variables
  real(kind(1.d0)) :: potential(3),potder(3)

  totmass = mass(1)+mass(2)+mass(3)
  redmass = (mass(1)*mass(2)*mass(3)/totmass)**0.5d0
  d = ((mass/redmass)*(1.d0-(mass/totmass)))**0.5d0
  d2=d*(((2.d0)**(-0.5d0)))

  channelangle(3) = 0.d0
  channelangle(1) = -2.d0*atan(mass(2)/redmass)
  channelangle(2) = 2.d0*atan(mass(1)/redmass)


  return

end subroutine setup

  
