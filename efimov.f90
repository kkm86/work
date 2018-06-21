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
  real(kind(1.d0)) :: tld(npl), tmd(npm), tld_max, tmd_max, tld_min, tmd_min 

  !.. Parameters for the 2-body potential
  real(kind(1.d0)) :: d
  real(kind(1.d0)) :: r0
  real(kind(1.d0)) :: mass(3)
 

  !.. Parameters for the energy curve
  integer, parameter :: points = 10
  real(kind(1.d0))   :: rho_vector(3,points),energy(LM,points,3),V(3,points)

  !.. Parameters for plotting
  integer, parameter :: pp = 50
  real(kind(1.d0)) :: x(pp),y(pp),step_size,delta_rho,rho_min,rho_max
  real(kind(1.d0)) :: base(pp,LM),hh,kk

  !.. Parameters for the wave function
  real(kind(1.d0))   :: base_L(pp,L),base_M(pp,M),summa(points)
  real(kind(1.d0)), allocatable, dimension(:,:,:) :: angwfn, Pmat, twobodypot
  real(kind(1.d0)), allocatable, dimension(:,:,:,:) :: coef,angder
  real(kind(1.d0)), allocatable, dimension(:,:,:,:,:) :: Pint, Pint2 


  !.. Other parameters
  real(kind(1.d0)) :: rho(points),my,H(LM,LM,points,3),S(LM,LM,points,3),t1,t2,Vtrap(3,points),angfreq,osc,U(points),integ(points)
  integer :: i,j,ll,mm,lj,li,mi,mj,mu,nu,n

  allocate(Pint2(pp,pp,LM,LM,points), Pint(pp,pp,LM,LM,points), angwfn(pp,pp,points),Pmat(pp,pp,points),twobodypot(pp,pp,points),coef(L,M,points,LM),angder(L,M,points,LM))

  !.. Declairing constants for model potential, trapping potential, and model atom 
  r0 = 55.d0
  d = -3.086d0*10**(-8.d0)
  mass = 87.d0*1836.15d0
  my = mass(1)/sqrt(3.d0)
  osc = 731.d0
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
  delta_rho = 0.0001d0
  rho_vector(1,1) = rho_min-delta_rho
  rho_vector(2,1) = rho_min
  rho_vector(3,1) = rho_min+delta_rho
  print*, rho_vector(1,1), rho_vector(2,1), rho_vector(3,1)
  do i = 2, points
     rho_vector(1,i) = rho_vector(1,i-1)+step_size
     rho_vector(2,i) = rho_vector(2,i-1)+step_size
     rho_vector(3,i) = rho_vector(3,i-1)+step_size
  end do

  !.. Declairing rho for later use and creating the harmonic trapping potential "Vtrap"
  rho = rho_vector(2,:)
  Vtrap = 0.5d0*my*(angfreq**2.d0)*(rho_vector**2.d0)

  !.. Calculating adiabatic potential curves and coefficients for the angular channel functions
  call CPU_TIME( t1 )
  write(6,*) 'hej5', points
     WRITE(6,*) "A",I
     call efimovham(npl,npm,k,L,M,LM,tl,tm,rho_vector,my,r0,d,mass,energy,H,S,points)
  call CPU_TIME( t2 )
  print*, t2-t1
  write(6,*) 'hej7'

  !.. Writes adiabatic potential curves+trapping potential to file
  open(10,file='threebodypot.dat',status='replace')
  do i = 1, points
     write(10,10)i, rho_vector(2,i)/osc, (energy(1,i,2)+Vtrap(2,i))/angfreq, (energy(2,i,2)+Vtrap(2,i))/angfreq ,(energy(3,i,2)+Vtrap(2,i))/angfreq,(energy(4,i,2)+Vtrap(2,i))/angfreq,(energy(5,i,2)+Vtrap(2,i))/angfreq,(energy(6,i,2)+Vtrap(2,i))/angfreq, Vtrap(2,i)/angfreq
10   format(I3,'  ',16f20.8)
  end do
  close(10)

  !.. Unwinds eigenvector coefficients 
  do n = 1, LM
     do mm = 1, M
        do ll = 1, L
           i = (ll-1)*M+mm
           angder(ll,mm,:,n) = 0.5d0*(H(i,n,:,3)-H(i,n,:,1))/delta_rho
           coef(ll,mm,:,n) = H(i,n,:,2)
        end do
     end do
  end do
 
  !.. Setting up vector for plotting  
  x(1) = tl(k)
  x(pp) = tl(npl)
  y(1) = tm(k)
  y(pp) = tm(npm)
 
  hh = (x(pp)-x(1))/(pp-1)
  kk = (y(pp)-y(1))/(pp-1)
  
  do i = 2, pp-1
     x(i) = x(i-1)+hh
     y(i) = y(i-1)+kk
  end do

  call B_spline_base(npl,npm,k,L,M,tl,tm,pp,x,y,base_L,base_M)

  call CPU_TIME( t1 )
  !.. Setting up P-matrix
  do nu = 1, LM
     do mu = 1, LM
        do j = 1, pp
           do i = 1, pp
              summa = 0.0
              do lj = 1, L
                 do li = 1, L
                    do mj = 1, M
                       do mi = 1, M
                          summa = summa + sin(2.d0*x(i))*coef(lj,mj,:,nu)*base_L(i,lj)*base_M(j,mj)*angder(li,mi,:,mu)*base_L(i,li)*base_M(j,mi)
                       end do
                    end do
                 end do
              end do
              !angwfn(i,j,:) = summa
              Pint(i,j,mu,nu,:) = summa
           end do
        end do
     end do
  end do

  do nu = 1, LM
     do mu = 1, LM
        do j = 1, pp
           do i = 1, pp
              angwfn(i,j,:) = Pint(i,j,mu,nu,:)
           end do
        end do
        call T2D(points,pp,hh,kk,angwfn,integ)
        Pmat(mu,nu,:) = integ
     end do
  end do
  call CPU_TIME( t2 )
  print*, t2-t1, 'P-mat'

  print *, Pmat(1,1,1),Pmat(2,2,1), 'ok'
  print *, Pmat(1,2,1),Pmat(2,3,1), 'ok'
  print *, Pmat(2,1,1),Pmat(3,2,1), 'ok'
  

  stop

 

  open(14,file='wave.dat',status='replace')
  do i = 1, points
     write(14,10)i,rho_vector(2,i)/osc,integ(i)
  end do
  close(14)

  ! open(13,file='adia.dat',status='replace')
  ! do i = 1, points
  !    write(13,10)i, rho_vector(2,i), integ(i) !(energy(1,i,2)+Vtrap(2,i))/angfreq, (energy(1,i,2)+Vtrap(2,i)-U(i)/(2.d0*my))/angfreq, U(i)!/(2.d0*my*angfreq) 
  ! end do
  ! close(13)

  deallocate(Pint2, Pint, angwfn, Pmat, twobodypot, coef, angder)

end program efimov



   
