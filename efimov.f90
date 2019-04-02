program efimov

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
 

  !.. Parameters for the knot-point grids tl and tm
  real(kind(1.d0)) :: tl(npl), tm(npm), tl_max, tm_max, tl_min, tm_min

  !.. Parameters for the 3-body system
  real(kind(1.d0)) :: my

  !.. Parameters for the harmonic trapping potential
   real(kind(1.d0)) :: angfreq, scaling
 
  !.. Parameters for effective potentials and coupling matrices
  integer, parameter :: points = 10
  real(kind(1.d0)), allocatable, dimension(:) :: rho_vector
  real(kind(1.d0)), allocatable, dimension(:,:) :: energy
  real(kind(1.d0)), allocatable, dimension(:,:,:) :: H,S,Hder,Hamcoef,Pmat,P2mat, Imat

  !.. Parameters for plotting
  integer, parameter :: pp = 200
  real(kind(1.d0)) :: x(pp),y(pp),step_size,rho_min,rho_max
  real(kind(1.d0)) :: base(pp,LM),base_L(pp,L),base_M(pp,M)  

  !.. Other parameters
  real(kind(1.d0)) :: t1,t2,Vtrap(points),U(points)
  integer :: i,j,ll,mm,lj,li,mi,mj,mu,nu,n,test

  allocate(rho_vector(points))
  allocate(energy(LM,points))
  allocate(H(LM,LM,points),S(LM,LM,points),Hder(LM,LM,points),Hamcoef(LM,LM,points),Pmat(LM,LM,points),P2mat(LM,LM,points),Imat(LM,LM,points))

  !.. Declairing constants for model potential, trapping potential, and model atom
  
  my = mass(1)/sqrt(3.d0)
  !angfreq = 1.d0/(mass(1)*osc**2.d0)
  angfreq = 1.d0
  scaling = 10**(8.d0)
  

  !.. Setting up knot-vectors
  tl_min = 0.d0
  tl_max = Pi/2.d0
  tm_min = 0.d0
  tm_max = Pi/3.d0


  call universal_knot(npl,k,N1,tl_max,tl_min,tl)
  call universal_knot(npm,k,N2,tm_max,tm_min,tm)


  !.. Setting up hyperradial vector
  rho_min = scatl
  rho_max = 37*scatl
  step_size = (rho_max-rho_min)/(points-1)
  rho_vector(1) = rho_min
  print*, rho_vector(1)
  do i = 2,points
     rho_vector(i) = rho_vector(i-1)+step_size
  end do

  !.. Declairing rho for later use and creating the harmonic trapping potential "Vtrap"

  !Vtrap = 0.5d0*my*(angfreq**2.d0)*(rho_vector**2.d0)
  Vtrap = 0.d0

  !.. Calculating adiabatic potential curves and coefficients for the angular channel functions
  call CPU_TIME( t1 )
  write(6,*) 'hej5', points
     WRITE(6,*) "A",I
     call efimovham(npl,npm,k,L,M,LM,tl,tm,rho_vector,my,energy,H,Hder,S,Hamcoef,points,Pmat,P2mat,Imat)
  call CPU_TIME( t2 )
  print*, t2-t1
  write(6,*) 'hej7'

  !.. Writes adiabatic potential curves+trapping potential to file
  ! open(10,file='threebodypot.dat',status='replace')
!   do i = 1, points
!      write(10,10)i, rho_vector(i)/scatl, scaling*(energy(1,i)+Vtrap(i))/angfreq,scaling*(energy(1,i)-(P2mat(1,1,i)/(2.d0*my))+Vtrap(i))/angfreq, scaling*(energy(2,i)+Vtrap(i))/angfreq,scaling*(energy(2,i)-(P2mat(2,2,i)/(2.d0*my))+Vtrap(i))/angfreq

!      !, (energy(2,i)+Vtrap(i))/angfreq ,(energy(3,i)+Vtrap(i))/angfreq,(energy(4,i)+Vtrap(i))/angfreq,(energy(5,i)+Vtrap(i))/angfreq,(energy(6,i)+Vtrap(i))/angfreq, Vtrap(i)/angfreq
! 10   format(I3,'  ',16f20.8)
!   End do
!   close(10)

  !.. Writes adiabatic potential curves+trapping potential to file
  open(10,file='effective5.dat',status='replace')
  do i = 1, points
     write(10,10)i, rho_vector(i)/scatl, scaling*energy(1,i), scaling*(15.d0/4.d0)/(2.d0*my*(rho_vector(i)**2.d0)), scaling*energy(2,i), scaling*(4*(4+4.d0)+15.d0/4.d0)/(2.d0*my*(rho_vector(i)**2.d0)),scaling*energy(3,i), scaling*(6*(6+4.d0)+15.d0/4.d0)/(2.d0*my*(rho_vector(i)**2.d0)),scaling*energy(4,i) 
     10   format(I3,'  ',16f20.8)
  end do
  close(10)

  open(14,file='wave.dat',status='replace')
  do i = 1, points
     write(14,10)i,rho_vector(i)/scatl, Pmat(1,2,i), Pmat(2,1,i), P2mat(1,2,i), P2mat(1,1,i)
  end do
  close(14)

  deallocate(Pmat,P2mat,Imat)

end program efimov



   
