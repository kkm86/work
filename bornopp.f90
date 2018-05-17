program bornopp

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
  real(kind(1.d0)) :: r_1

   !.. Parameters for the energy curve
  integer, parameter :: points = 20
  real(kind(1.d0))   :: rho_vector(points), energy_curve(points), rho

  !.. Parameters for the wave function
  real(kind(1.d0))   :: wave(pp), wfn(NN)

  !.. Parameters for the knot-point grid
  real(kind(1.d0)) :: t(np), r_max, r_min, step_size, small_step

  !.. Local
  integer :: ii, jj

  !.. External
  real(kind(1.d0)) :: energy


  b = 12.d0
  S = 1.d0
  mass = 1.d0

  !.. Setting up vector for plotting
  x(1) = t(k)
  x(pp) = t(np)

  step_size = (x(pp)-x(1))/pp

  do ii = 2, pp-1
     x(ii) = x(ii-1)+step_size
  end do

  rho_vector(1) = 0.01d0
  do ii = 2, points
     rho_vector(ii) = ii*1.d0
     print *, ii, rho_vector(ii)
  end do

  open(10,file='result_wave_born.dat',status='replace')
  do ii = 1, points
     rho = rho_vector(ii)
     call knot_vector(np,k,NN,rho,t)
     call Hamiltonian_born(np,k,NN,t,energy,wfn,S,b,mass)
     energy_curve(ii) = energy
     r_1 = t(np)
     call twobody_born(S,b,mass,r_1,V)
     print *, rho_vector(ii), energy_curve(ii), V
     write(10,10)ii, rho_vector(ii), energy_curve(ii), V
10   format(I4,'  ',16f20.8)
  end do
  close(10)
    
 
  
 ! call Hamiltonian_born(np,k,NN,t,energy,wfn,S,b,mass)

 ! print *, 'lambda_tb(0) = '
 ! print *, energy

 ! call B_spline_base(np,k,NN,t,pp,x,base)

 ! wave = matmul(base,wfn)
 
 ! open(20,file='result_wave_born.dat',status='replace')
 ! do ii = 1, pp
 !    write(20,10)ii,x(ii), wave(ii)**2.d0
!10   format(I3,'  ',16f20.8)
!  end do
!  close(20)

end program bornopp

