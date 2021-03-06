program efimov2

  use constants
  
  implicit none
  
  !.. Input
  !.. Parameters for the B-splines used in the generalized eigenvalue equation
  integer, parameter :: N1 = 5   !.. Number of mesh-points in coordinate 1
  integer, parameter :: N2 = 5   !.. Number of mesh-points in coordinate 2
  integer, parameter :: k = 6    !.. B-spline order
  integer, parameter :: L = 7    !.. Number of B-splines in coordinate 1(N+k-2-cond)
  integer, parameter :: M = 7    !.. Number of B-splines in coordinate 2
  integer, parameter :: LM = 49 !.. Matrix dimension
  integer, parameter :: npl = 15  !.. Number of knot-points  N1+2(k-1)
  integer, parameter :: npm = 15  !.. Number of knot-points  N2+2(k-1)
  integer, parameter :: cc = 3    !.. Number of couplings to be calculated, for full coupling matrices cc = LM 
 

  !.. Parameters for the knot-point grids tl and tm
  real(kind(1.d0)) :: tl(npl), tm(npm), tl_max, tm_max, tl_min, tm_min

  !.. Parameters for the 3-body system
  real(kind(1.d0)) :: my

  !.. Parameters for the harmonic trapping potential
   real(kind(1.d0)) :: angfreq
 
  !.. Parameters for effective potentials and coupling matrices
  integer, parameter :: points = 10
  real(kind(1.d0)), allocatable, dimension(:) :: rho_vector
  real(kind(1.d0)), allocatable, dimension(:,:) :: energy
  real(kind(1.d0)), allocatable, dimension(:,:,:) :: H,S,Hder,Hamcoef,Pmat,P2mat,Imat

  !.. Parameters for plotting
  integer, parameter :: pp = 10
  real(kind(1.d0)) :: x(pp),y(pp),step_size,rho_min,rho_max
  real(kind(1.d0)) :: base(pp,LM),base_L(pp,L),base_M(pp,M)  

  !.. Other parameters
  real(kind(1.d0)) :: t1,t2,Vtrap(points),U(points)
  integer :: i,j,ll,mm,lj,li,mi,mj,mu,nu,n,test

  allocate(rho_vector(points))
  allocate(energy(LM,points))
  allocate(H(LM,LM,points),S(LM,LM,points),Hder(LM,LM,points),Hamcoef(LM,LM,points),Pmat(cc,cc,points),P2mat(cc,cc,points),Imat(cc,cc,points))

  !.. Declairing constants for model potential, trapping potential, and model atom
  
  my = mass(1)/sqrt(3.d0)
  angfreq = 1.d0/(mass(1)*osc**2.d0)
  

  !.. Setting up knot-vectors
  tl_min = 0.d0
  tl_max = Pi/2.d0
  tm_min = 0.d0
  tm_max = Pi/3.d0


  call universal_knot(npl,k,N1,tl_max,tl_min,tl)
  call universal_knot(npm,k,N2,tm_max,tm_min,tm)


  !.. Setting up hyperradial vector
  rho_min = 1.d0
  rho_max = 1000.d0
  step_size = (rho_max-rho_min)/(points-1)
  rho_vector(1) = rho_min
  print*, rho_vector(1)
  do i = 2,points
     rho_vector(i) = rho_vector(i-1)+step_size
  end do

  !.. Declairing rho for later use and creating the harmonic trapping potential "Vtrap"

  Vtrap = 0.5d0*my*(angfreq**2.d0)*(rho_vector**2.d0)

  !.. Calculating adiabatic potential curves and coefficients for the angular channel functions
  call CPU_TIME( t1 )
  write(6,*) 'hej5', points
     WRITE(6,*) "A",I
     call efimovham_split(npl,npm,k,L,M,LM,tl,tm,rho,my,energy,H,Hder,S,Integ,points)
     call coupling(H,Hder,Integ,energy,LM,points,cc,Pmat,P2mat,Imat)
  call CPU_TIME( t2 )
  print*, t2-t1
  write(6,*) 'hej7'

  !.. Writes adiabatic potential curves+trapping potential to file
  open(10,file='threebodypot.dat',status='replace')
  do i = 1, points
     write(10,10)i, rho_vector(i)/osc, (energy(1,i)+Vtrap(i))/angfreq,(energy(1,i)-(P2mat(1,1,i)/(2.d0*my))+Vtrap(i))/angfreq, (energy(2,i)+Vtrap(i))/angfreq,(energy(2,i)-(P2mat(2,2,i)/(2.d0*my))+Vtrap(i))/angfreq
   
10   format(I3,'  ',16f20.8)
  end do
  close(10)

  deallocate(Pmat,P2mat,Imat)

end program efimov2




   
