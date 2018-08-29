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

contains

  function potent_omp(rho,theta,phi,points,k)
  
  implicit none

  !.. Input
  integer, intent(in) :: points,k
  real(kind(1.d0)), intent(in) :: phi(k),theta,rho(points)

  !.. Output
  real(kind(1.d0)) :: potent_omp(k,points)

  !.. Local
  real(kind(1.d0)) :: V(points),Vder(points)
  integer i
  
  potent_omp = 0.d0
  V = 0.d0
  Vder = 0.d0
  do i=1,k
  call twobody_potential(rho,theta,phi(i),V,Vder,points)
  potent_omp(i,:) = V
  end do
  return
end function potent_omp

function potentder_omp(rho,theta,phi,points,k)
  
  implicit none

  !.. Input
  integer, intent(in) :: points,k
  real(kind(1.d0)), intent(in) :: phi(k),theta,rho(points)

  !.. Output
  real(kind(1.d0)) :: potentder_omp(k,points)

  !.. Local
  real(kind(1.d0)) :: V(points),Vder(points)
  integer i

  potentder_omp = 0.d0
  V = 0.d0
  Vder = 0.d0
  do i=1,k
     call twobody_potential(rho,theta,phi(i),V,Vder,points)
     potentder_omp(i,:) = Vder
  end do
  return
end function potentder_omp



end module constants

