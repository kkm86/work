program ani

  use constants
  
  implicit none
  
  !.. Input
  !.. Parameters for the B-splines used in the generalized eigenvalue equation
  integer, parameter :: N1 = 15   !.. Number of mesh-points in coordinate 1
  integer, parameter :: N2 = 15   !.. Number of mesh-points in coordinate 2
  integer, parameter :: k = 6    !.. B-spline order
  integer, parameter :: L = 18    !.. Number of B-splines in coordinate 1(N+k-2-cond)
  integer, parameter :: M = 17    !.. Number of B-splines in coordinate 2
  integer, parameter :: LM = 306  !.. Matrix dimension
  integer, parameter :: np = 25  !.. Number of knot-points  Nd+2(k-1)
 

  !.. Parameters for the knot-point grids tl and tm
  real(kind(1.d0)) :: tl(np), tm(np), tl_max, tm_max, tl_min, tm_min, energy, H(LM,LM), S(LM,LM)

  !.. Parameters for plotting eigenvectors
  integer, parameter :: pp = 100
  real(kind(1.d0)) :: f(LM,7), wfn(pp,pp,7), base_L(pp,L),base_M(pp,M), x(pp), y(pp), step_size_x, step_size_y
  real(kind(1.d0)) :: c(L,M,7), term(7)

  integer :: d, ii, i, mm, ll, n, j

  tl_min = 0.d0
  tl_max = 10.d0
  tm_min = 0.d0
  tm_max = 7.07d0

  call universal_knot(np,k,N1,tl_max,tl_min,tl)
  call universal_knot(np,k,N2,tm_max,tm_min,tm)
  
  call anisotrop(np,k,L,M,LM,tl,tm,energy,H,S)

  do i = 1, LM
     do n = 1, 7
        f(i,n) = H(i,n)
     end do
  end do

  do ll = 1, L
     do mm = 1, M
        do n = 1, 7
           i = (ll-1)*M+mm
           c(ll,mm,n) = f(i,n)
        end do
     end do
  end do    
 
  
  !.. Setting up vector for plotting
  x(1) = tl(k)
  x(pp) = 4.d0
  y(1) = tm(k)
  y(pp) = 4.d0

  ! x(1) = tl(k)
  ! x(pp) = tl(np)
  ! y(1) = tm(k)
  ! y(pp) = tm(np)

  step_size_x = (x(pp)-x(1))/pp
  step_size_y = (y(pp)-y(1))/pp

  do ii = 2, pp-1
     x(ii) = x(ii-1)+step_size_x
     y(ii) = y(ii-1)+step_size_y
  end do

  call B_spline_base(np,k,L,M,LM,tm,tl,pp,x,y,base_L,base_M)
  
  do n = 1, 7
     do j = 1, pp
        do i = 1, pp
           term = 0.0
           do ll = 1, L
              do mm = 1, M
                 term(n) = term(n) + c(ll,mm,n)*base_L(i,ll)*base_M(j,mm)
              end do
           end do
           wfn(i,j,n) = term(n)
        end do
     end do
  end do
  
  
              
 

  open(20,file='result_wave.dat',status='replace')
  
  do j = 1, pp
     do i = 1, pp
        write(20,10)i, j, x(i), y(j), wfn(i,j,1)**2.d0, wfn(i,j,6)**2.d0, wfn(i,j,7)**2.d0
10      format(I4, I4,'  ',16f20.8)
     end do
  end do
  close(20)
    

end program ani

   
