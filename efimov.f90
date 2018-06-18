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
  real(kind(1.d0)) :: d(7)
  real(kind(1.d0)) :: r0 !,r(3), potential(3)
  real(kind(1.d0)) :: mass(3)
  !real(kind(1.d0)) :: theta,phi

  !.. Parameters for the energy curve
  integer, parameter :: points = 100
  real(kind(1.d0))   :: rho_vector(3,points),energy(LM,points,3),V(3,points)

  !.. Parameters for plotting
  integer, parameter :: pp = 100
  integer, parameter :: ss = 4
  real(kind(1.d0)) :: x(pp),y(pp),step_size,delta_rho,rho_min,rho_max
  real(kind(1.d0)) :: base(pp,LM),aa,bb,cc,dd,hh,kk

  !.. Parameters for the wave function
  real(kind(1.d0))   :: wfn(pp,pp,ss),angwfn(pp,pp),f(LM,ss), c(L,M,points),cv(L,M,points),term(ss), base_L(pp,L), base_M(pp,M)


  !.. Other parameters
  real(kind(1.d0)) :: rho(points), my, H(LM,LM,points,3), S(LM,LM,points,3),t1,t2, Vtrap(3,points), angfreq, osc, U(points), integ
  integer :: i,j,jj,ii,ll,mm,n

  
  r0 = 55.d0
  mass = 87.d0*1836.15d0
  my = mass(1)/sqrt(3.d0)
  osc = 731.d0
  angfreq = 1.d0/(mass(1)*osc**2.d0)
  

  d(1) = -3.086d0*10**(-8.d0)
 
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
  rho_vector(1,1) = rho_min-delta_rho
  rho_vector(2,1) = rho_min
  rho_vector(3,1) = rho_min+delta_rho
  print*, rho_vector(1,1), rho_vector(2,1), rho_vector(3,1)
  do i = 2, points
     rho_vector(1,i) = rho_vector(1,i-1)+step_size
     rho_vector(2,i) = rho_vector(2,i-1)+step_size
     rho_vector(3,i) = rho_vector(3,i-1)+step_size
  end do

  rho = rho_vector(2,:)
  Vtrap = 0.5d0*my*(angfreq**2.d0)*(rho_vector**2.d0)

 

  
  !.. Code for plotting eigenfunctions start

  ! do n = 1, ss
!      do i = 1, LM
!         f(i,n) = H(i,n)
!      end do
!   end do

!   do ll = 1, L
!      do mm = 1, M
!         do n = 1, ss
!            i = (ll-1)*M+mm
!            c(ll,mm,n) = f(i,n)
!         end do
!      end do
!   end do
  
!   !.. Setting up vector for plotting
!   x(1) = tl(k)
!   x(pp) = tl(np)
!   y(1) = tm(k)
!   y(pp) = tm(np)

!   hh = (x(pp)-x(1))/pp
!   kk = (y(pp)-y(1))/pp

!   do ii = 2, pp-1
!      x(ii) = x(ii-1)+hh
!      y(ii) = y(ii-1)+kk
!   end do

!   call B_spline_base(np,k,L,M,tl,tm,pp,x,y,base_L,base_M)

!  do n = 1, ss
!      do j = 1, pp
!         do i = 1, pp
!            term = 0.0
!            do ll = 1, L
!               do mm = 1, M
!                  term(n) = term(n) + c(ll,mm,n)*base_L(i,ll)*base_M(j,mm)
!               end do
!            end do
!            wfn(i,j,n) = term(n)
!         end do
!      end do
!   end do

!   open(20,file='result_wave.dat',status='replace')
  
!   do j = 1, pp
!      do i = 1, pp
!         write(20,10)i, j, x(i), y(j), wfn(i,j,1)**2.d0, wfn(i,j,2)**2.d0, wfn(i,j,3)**2.d0
! 10      format(I4, I4,'  ',16f20.8)
!      end do
!   end do
!   close(20)

  !.. Code for plotting eigenfunction end
  
  
  call CPU_TIME( t1 )
  write(6,*) 'hej5', points
     WRITE(6,*) "A",I
     call efimovham(npl,npm,k,L,M,LM,tl,tm,rho_vector,my,r0,d(1),mass,energy,H,S,points)
  call CPU_TIME( t2 )
  print*, t2-t1
  write(6,*) 'hej7'


  open(10,file='threebodypot.dat',status='replace')
  do i = 1, points
     write(10,10)i, rho_vector(2,i)/osc, (energy(1,i,2)+Vtrap(2,i))/angfreq, (energy(2,i,2)+Vtrap(2,i))/angfreq ,(energy(3,i,2)+Vtrap(2,i))/angfreq,(energy(4,i,2)+Vtrap(2,i))/angfreq,(energy(5,i,2)+Vtrap(2,i))/angfreq,(energy(6,i,2)+Vtrap(2,i))/angfreq, Vtrap(2,i)/angfreq
10   format(I3,'  ',16f20.8)
  end do
  close(10)

  do j = 1, points
     do mm = 1, M
        do ll = 1, L
           i = (ll-1)*M+mm
           c(ll,mm,j) = H(i,1,j,2)
           cv(ll,mm,j) = H(i,1,j,2)
        end do
     end do
  end do

  !.. Setting up vector for plotting
  aa = tl(k)
  bb = tl(npl)
  cc = tm(k)
  dd = tm(npm)

  print*,'aabbccdd'
  print*,aa,bb,cc,dd
 
  hh = (bb-aa)/(pp-1)
  kk = (dd-cc)/(pp-1)

  x(1) = aa
  x(pp) = bb
  y(1) = cc
  y(pp) = dd
  
  do ii = 2, pp-1
     x(ii) = x(ii-1)+hh
     y(ii) = y(ii-1)+kk
  end do

  print*, 'y(n)'
  print*, y(1),y(pp),hh,kk

  call B_spline_base(npl,npm,k,L,M,tl,tm,pp,x,y,base_L,base_M)


  do j = 1, pp
     do i = 1, pp
        term = 0.0
        do ll = 1, L
           do mm = 1, M
              term(1) = term(1) + 2.d0*Pi*sin(2.d0*x(i))*c(ll,mm,1)*cv(ll,mm,1)*(base_L(i,ll)*base_M(j,mm))**2.d0
           end do
        end do
        angwfn(i,j) = term(1)
     end do
  end do

 
  call T2D(pp,hh,kk,angwfn,integ)

  stop
  ! do j = 1, 1
  !    do mm = 1, M
  !       do ll = 1, L
  !          i = (ll-1)*M+mm
  !          cv(ll,mm,j) = 0.5d0*(H(i,1,j,3)-H(i,1,j,1))/delta_rho
  !       end do
  !    end do
  ! end do

  call adiabatic(npl,npm,k,L,M,LM,tl,tm,rho,my,c,cv,U,points)

  ! open(13,file='adia.dat',status='replace')
  ! do i = 1, points
  !    write(13,10)i, rho_vector(2,i)/osc, (energy(1,i,2)+Vtrap(2,i))/angfreq, (energy(1,i,2)+Vtrap(2,i)-U(i)/(2.d0*my))/angfreq, U(i)!/(2.d0*my*angfreq) 
  ! end do
  ! close(10)

  ! open(11,file='coeffmin.dat',status='replace')
  ! do i = 1, points
  !    do j = 1, LM
  !    write(11,10) H(j,1,i,1)
  ! end do
  ! end do
  ! close(11)

  ! open(12,file='coeffmax.dat',status='replace')
  ! do i = 1, points
  !    do j = 1, LM
  !    write(12,10) H(j,1,i,3)
  ! end do
  ! end do
  ! close(12)

  

end program efimov



   
