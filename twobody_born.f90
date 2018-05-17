subroutine twobody_born(S,b,mass,r_1,V)
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in) :: S,b,mass,r_1
  real(kind(1.d0)), intent(inout) :: V

  !.. Local
  real(kind(1.0d0)) :: r

  r = r_1

  V = -2.d0*mass*S*exp(-((r/b)**2.d0))

  return

end subroutine twobody_born




  
