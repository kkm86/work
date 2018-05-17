program twodim

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
  real(kind(1.d0)) :: tl(np), tm(np), tl_max, tm_max, tl_min, tm_min 
  real(kind(1.d0)) :: energy2(3)

  
  tl_min = 0.d0
  tl_max = 10.d0
  tm_min = 0.d0
  tm_max = 7.07d0

  call universal_knot(np,k,N1,tl_max,tl_min,tl)
  call universal_knot(np,k,N2,tm_max,tm_min,tm)
  call ham2d(np,k,L,M,LM,tl,tm,energy2)

end program twodim
   
