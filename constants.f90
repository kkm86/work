module constants

  implicit none

  real(kind(1.d0)) :: PI = 4.d0*atan(1.d0)

  !.. Parameters for the 3-body system
  real(kind(1.d0)) :: mass(3)=87.d0*1836.15d0
 
  !.. Parameters for the 2-body potential, harmonic trapping potential, and scattering lengths
  real(kind(1.d0)) :: Potential_depth = -3.086d0*10**(-8.d0)
  real(kind(1.d0)) :: r0 = 55.d0
  real(kind(1.d0)) :: osc = 731.d0
  real(kind(1.d0)) :: scatl = 228.004d0

end module constants

