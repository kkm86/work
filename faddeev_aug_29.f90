program faddeev

  use constants
  
  implicit none
  
  !.. Input
  !.. Parameters for the B-splines used in the generalized eigenvale equation
  integer, parameter :: np = 190 !.. Number of knot-points
  integer, parameter :: k = 7    !..  
  integer, parameter :: NN = 181 !.. Number of B-splines (np-k-2)

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
  integer, parameter :: points = 300
  real(kind(1.d0))   :: rho_vector(points), energy_curve(points), rho

  !.. Parameters for the wave function
  real(kind(1.d0))   :: wave(pp), wfn(NN)

  !.. Parameters for the knot-point grid
  real(kind(1.d0)) :: t(np), alpha_max, alpha_min, step_size

  !.. Local
  integer :: ii, jj

  !.. External
  real(kind(1.d0)) :: energy

  write(*,*)'Give the range (b) of the potential: '
  read(*,*)b
  write(*,*)'Give the depth (S) of the potential: '
  read(*,*)S
  write(*,*)'Give the mass: '
  read(*,*)mass

  !.. Setting up knot vector
  alpha_max = PI/2.d0
  alpha_min = 0.d0
  step_size = (alpha_max-alpha_min)/(NN-4)

  do ii = 1, np
     if(ii .le. k) then
        t(ii) = 0.d0
     elseif(ii .le. np-k+1) then
        t(ii) = t(ii-1)+step_size
     else
        t(ii) = t(ii-1)
     end if
     print *, ii, t(ii)
  end do

  !.. Setting up vector for plotting
  x(1) = t(k)
  x(pp) = t(np)

  step_size = (x(pp)-x(1))/pp

  do ii = 2, pp-1
     x(ii) = x(ii-1)+step_size
  end do

  do ii = 1, points
     rho_vector(ii) = ii*0.1d0
  end do

  rho = 0.d0
  call Hamiltonian(np,k,NN,t,rho,energy,wfn)
 
  print *, 'lambda_tb(0) = '
  print *, energy

  call B_spline_base(np,k,NN,t,pp,x,base)

  wave = matmul(base,wfn)
 
  open(20,file='result_wave.dat',status='replace')
  do ii = 1, pp
     write(20,10)ii,x(ii), wave(ii)**2.d0, sin(2*x(ii))**2.d0
  end do
  close(20)

  !.. Debugging
  stop

  do ii = 1, points
     rho = rho_vector(ii)
     call Hamiltonian(np,k,NN,t,rho,energy,wfn)
     energy_curve(ii) = energy
     print *, rho_vector(ii), energy_curve(ii)
  end do

  open(10,file='result.dat',status='replace')
  do ii = 1, points
     write(10,10)ii, rho_vector(ii), energy_curve(ii) !/(rho_vector(ii)**2.d0)
10   format(I3,'  ',16f20.8)
  end do
  close(10)
  
  
  

end program faddeev
