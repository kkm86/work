program faddeev

  use constants
  
  implicit none
  
  !.. Input
  !.. Parameters for the B-splines used in the generalized eigenvale equation
  integer, parameter :: np = 20 !.. Number of knot-points
  integer, parameter :: k = 5    !..  
  integer, parameter :: NN = 13 !.. Number of B-splines (np-k-2)

  !.. Parameters for plotting
  integer, parameter :: pp = 300
  real(kind(1.d0)) :: x(pp)
  real(kind(1.d0)) :: base(pp,NN)

  !.. Parameters for the 2-body potential
  real(kind(1.d0)) :: S
  real(kind(1.d0)) :: b
  real(kind(1.d0)) :: mass
  real(kind(1.d0)) :: V

  !.. Parameters for the energy curve
  integer, parameter :: points = 200
  real(kind(1.d0))   :: rho_vector(points), energy_curve(points), rho

  !.. Parameters for the wave function
  real(kind(1.d0))   :: wave(pp), wfn(NN)

  !.. Parameters for the knot-point grid
  real(kind(1.d0)) :: t(np), r_max, r_min, step_size, small_step

  !.. Local
  integer :: ii, jj

  !.. External
  real(kind(1.d0)) :: energy

!  write(*,*)'Give the range (b) of the potential: '
!  read(*,*)b
!  write(*,*)'Give the depth (S) of the potential: '
!  read(*,*)S
!  write(*,*)'Give the mass: '
!  read(*,*)mass

  b = 1.d0
  S = 1.d0
  mass = 1.d0

  !.. Setting up vector for plotting
  x(1) = t(k)
  x(pp) = t(np)

  step_size = (x(pp)-x(1))/pp

  do ii = 2, pp-1
     x(ii) = x(ii-1)+step_size
  end do

  rho_vector(1) = 0.d0
  do ii = 2, points
     rho_vector(ii) = ii*0.1d0
  end do

  do ii = 1, points
     rho = rho_vector(ii)

     !.. Setting up knot vector
     r_max = rho*sqrt(2.0)
     r_min = 0.d0
     step_size = (r_max-r_min)/(NN-4)
     small_step = 1.e-6

     do jj = 1, np
        if(jj .le. k) then
           t(jj) = 0.d0
        elseif(jj .le. k+2) then
           t(jj) = t(jj-1)+small_step
        elseif(jj .le. np-k) then
           t(jj) = t(jj-1)+step_size
        else
           t(jj) = r_max
        end if
        !print *, jj, t(jj)
     end do

     call Hamiltonian_rdep(np,k,NN,t,rho,energy,wfn,S,b,mass)
     energy_curve(ii) = energy
     print *, rho_vector(ii), energy_curve(ii)
  end do

  open(10,file='result_rdep.dat',status='replace')
  do ii = 1, points
     write(10,10)ii, rho_vector(ii), energy_curve(ii) 
10   format(I3,'  ',16f20.8)
  end do
  close(10)
  
  
  

end program faddeev
