program efimov

  use constants
  
  implicit none
  
  !.. Input
  !.. Parameters for the B-splines used in the generalized eigenvalue equation
  integer, parameter :: N1 = 25   !.. Number of mesh-points in coordinate 1
  integer, parameter :: N2 = 25   !.. Number of mesh-points in coordinate 2
  integer, parameter :: k = 6    !.. B-spline order
  integer, parameter :: L = 28    !.. Number of B-splines in coordinate 1(N+k-2-cond)
  integer, parameter :: M = 27    !.. Number of B-splines in coordinate 2
  integer, parameter :: LM = 756  !.. Matrix dimension
  integer, parameter :: np = 35  !.. Number of knot-points  Nd+2(k-1)
 

  !.. Parameters for the knot-point grids tl and tm
  real(kind(1.d0)) :: tl(np), tm(np), tl_max, tm_max, tl_min, tm_min 
  real(kind(1.d0)) :: energy(3)

  !.. Parameters for the 2-body potential
  real(kind(1.d0)) :: d(7)
  real(kind(1.d0)) :: r0,r(3), potential(3)
  real(kind(1.d0)) :: mass(3)
  real(kind(1.d0)) :: V,theta,phi

  !.. Parameters for the energy curve
  integer, parameter :: points = 2000
  real(kind(1.d0))   :: rho_vector(points), energy_curve(7,points)

  !.. Other parameters
!  real(kind(1.d0)) :: rho, my,t1,t2,s0,lamda 
   real(kind(1.d0)) :: rho, my, H(LM,LM), S(LM,LM),t1,t2,s0,lamda
  integer :: i,j

  
  r0 = 55.d0
  mass = 87*1836.15d0
  s0 = 1.00624d0
  my = mass(1)/sqrt(3.d0)
  lamda = 4.d0

  d(1) = -6.619d0*10**(-8.d0)
  print*, d(1)
  
  tl_min = 0.d0
  tl_max = Pi/2.d0
  tm_min = 0.d0
  tm_max = 2.d0*Pi/3.d0


  call universal_knot(np,k,N1,tl_max,tl_min,tl)
  call universal_knot(np,k,N2,tm_max,tm_min,tm)



  rho_vector(1) = 10000.d0
  do i = 2, points/3
     rho_vector(i) = (i-1)*2.d0
  end do

  do i = points/3, points/2
     rho_vector(i) = (i-1)*7.d0
  end do

  do i = points/2, points
     rho_vector(i) = (i-1)*10.d0
  end do

  call CPU_TIME( t1 )
    write(6,*) 'hej5', points
  do i = 1, points
     rho = rho_vector(i)
     WRITE(6,*) "A",I
     call efimovham(np,k,L,M,LM,tl,tm,rho,my,r0,d(1),mass,energy,H,S)
     WRITE(6,*) "b",I
     energy_curve(1,i) = energy(1)
     energy_curve(2,i) = energy(2)
     energy_curve(3,i) = energy(3)
  end do
  write(6,*) 'hej6'
  call CPU_TIME( t2 )
  print*, t2-t1
  write(6,*) 'hej7'
 
  
  open(10,file='result2derny1.dat',status='replace')
  do i = 1, points
     write(10,10)i, rho_vector(i)/228.992d0, energy_curve(1,i)*(10**8.d0),energy_curve(2,i)*(10**8.d0),(lamda*(lamda+4.d0)+15*0.25d0)/(2.d0*my*(rho_vector(i)**2.d0))*(10**8.d0),(15*0.25d0)/(2.d0*my*(rho_vector(i)**2.d0))*(10**8.d0)!, -((s0**2.d0)+0.25d0)*(10**8.d0)/(2.d0*my*(rho_vector(i)**2.d0))
10   format(I3,'  ',16f20.8)
  end do
  close(10)
  write(6,*) 'hej8'
 !energy_curve(1,i)*10**(8.d0), energy_curve(2,i)*10**(8.d0)
end program efimov

   
