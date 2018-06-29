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
  integer, parameter :: points = 300
  integer, parameter :: pder = 100
  real(kind(1.d0))   :: rho_vector(points),energy(LM,points),V(points),Vder(points)

  !.. Parameters for plotting
  integer, parameter :: pp = 50
  real(kind(1.d0)) :: x(pp),y(pp),step_size,delta_rho,rho_min,rho_max
  real(kind(1.d0)) :: base(pp,LM),hh,kk

  !.. Parameters for the wave function
  real(kind(1.d0))   :: base_L(pp,L),base_M(pp,M),summa(pder)
  real(kind(1.d0)), allocatable, dimension(:,:,:) :: angwfn, Pmat,P2mat,PmatH, Imat, deriv_coef,normal_coef
  real(kind(1.d0)), allocatable, dimension(:,:,:,:) :: coef,angder
  real(kind(1.d0)), allocatable, dimension(:,:,:,:,:) :: Pint


  !.. Other parameters
  real(kind(1.d0)) :: rho(points),my,H(LM,LM,points),S(LM,LM,points),Hder(LM,LM,points),Hamcoef(LM,LM,points),t1,t2,Vtrap(points),angfreq,osc,scatl,U(points),integ(pder)
  integer :: i,j,ll,mm,lj,li,mi,mj,mu,nu,n,test

  allocate(Pint(pp,pp,LM,LM,pder),angwfn(pp,pp,pder),Pmat(LM,LM,points),P2mat(LM,LM,points),PmatH(pp,pp,pder),Imat(LM,LM,points),deriv_coef(LM,LM,pder),normal_coef(LM,LM,pder),coef(L,M,pder,LM),angder(L,M,pder,LM))

  !.. Declairing constants for model potential, trapping potential, and model atom 
  r0 = 55.d0
  d = -3.086d0*10**(-8.d0)
  scatl = 228.004
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
  delta_rho = 0.001d0
  rho_vector(1) = rho_min-delta_rho
  rho_vector(2) = rho_min
  rho_vector(3) = rho_min+delta_rho
  print*, rho_vector(1), rho_vector(2), rho_vector(3)
  do i = 4,points,3
     rho_vector(i) = rho_vector(i-3)+step_size
     rho_vector(i+1) = rho_vector(i-2)+step_size
     rho_vector(i+2) = rho_vector(i-1)+step_size
     print*, rho_vector(i), rho_vector(i+1), rho_vector(i+2)
  end do

  !.. Declairing rho for later use and creating the harmonic trapping potential "Vtrap"
  rho = rho_vector(:)
  Vtrap = 0.5d0*my*(angfreq**2.d0)*(rho_vector**2.d0)

  !.. Calculating adiabatic potential curves and coefficients for the angular channel functions
  call CPU_TIME( t1 )
  write(6,*) 'hej5', points
     WRITE(6,*) "A",I
     call efimovham(npl,npm,k,L,M,LM,tl,tm,rho,my,r0,d,mass,energy,H,Hder,S,Hamcoef,points,Pmat,P2mat,Imat)
  call CPU_TIME( t2 )
  print*, t2-t1
  write(6,*) 'hej7'

  !.. Writes adiabatic potential curves+trapping potential to file
  open(10,file='threebodypot.dat',status='replace')
  do i = 1, points
     write(10,10)i, rho_vector(i)/osc, (energy(1,i)+Vtrap(i))/angfreq,(energy(1,i)-(P2mat(1,1,i)/(2.d0*my))+Vtrap(i))/angfreq, (energy(2,i)+Vtrap(i))/angfreq,(energy(2,i)-(P2mat(2,2,i)/(2.d0*my))+Vtrap(i))/angfreq

     !, (energy(2,i)+Vtrap(i))/angfreq ,(energy(3,i)+Vtrap(i))/angfreq,(energy(4,i)+Vtrap(i))/angfreq,(energy(5,i)+Vtrap(i))/angfreq,(energy(6,i)+Vtrap(i))/angfreq, Vtrap(i)/angfreq
10   format(I3,'  ',16f20.8)
  end do
  close(10)

  !Jumping code
  test = 0
  test = 0

1 continue
  write(*,*)
  test=test+1

  if(test .gt. 2)goto 95
  if(test .eq. 1)goto 100

100 continue
  
  !.. Unwinds eigenvector coefficients for calculating derivative using finite element method

  do i = 1, 3, points
     deriv_coef(:,:,i) = 0.5d0*(H(:,:,i)-H(:,:,i+2))/delta_rho
     normal_coef(:,:,i) = H(:,:,i+1)
  end do
  
  do n = 1, LM
     do mm = 1, M
        do ll = 1, L
           i = (ll-1)*M+mm
           angder(ll,mm,:,n) = deriv_coef(i,n,:)
           coef(ll,mm,:,n) = normal_coef(i,n,:)
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
                          summa = summa+sin(2.d0*x(i))*coef(lj,mj,:,mu)*base_L(i,lj)*base_M(j,mj)*angder(li,mi,:,nu)*base_L(i,li)*base_M(j,mi)
                       end do
                    end do
                 end do
              end do
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
        call T2D(pder,pp,hh,kk,angwfn,integ)
        PmatH(mu,nu,:) = integ
     end do
  end do
  call CPU_TIME( t2 )
  print*, t2-t1, 'P-mat'

95 continue
  write(*,*) 'test = ',test,'Got to 95'


  open(14,file='wave.dat',status='replace')
  do i = 1, points
     write(14,10)i,rho_vector(i)/scatl, Pmat(1,2,i), Pmat(2,1,i), P2mat(1,2,i), P2mat(1,1,i)
  end do
  close(14)

  deallocate(Pint,angwfn,Pmat,P2mat,PmatH,Imat,coef,angder,deriv_coef,normal_coef)

end program efimov



   
