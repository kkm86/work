subroutine T2D(points,pp,h,k,func,integ)
          

  use constants
  
  implicit none

  !.. This subroutine performs integration of a two variable function
  !.. Using Simpson's 2D Rule
  
  !.. Input
  integer         , intent(in) :: pp, points                !.. # subintervals in x and y
  real(kind(1.d0)), intent(in) :: func(pp,pp,points)        !.. function to integrate
  real(kind(1.d0)), intent(in) :: h,k                       !.. subinterval width h:x, k:y
  
  !.. Output
  real(kind(1.d0)), intent(inout) :: integ(points)  !.. result of integration

  !.. Local
  integer :: ii, jj, m 
  real(kind(1.d0)) :: sum4(points), sum2(points), sum16(points), sum8(points), sum42(points)
  real(kind(1.d0)) :: fsum1(points)

  sum4 = 0.d0
  sum2 = 0.d0
  sum16 = 0.d0
  sum8 = 0.d0
  sum42 = 0.d0

  m = pp/2
  
  do ii = 2, m 
     sum4 = sum4 + 9.d0*(func(2*ii-1,1,:)+func(2*ii-1,pp,:)+func(1,2*ii-1,:)+func(pp,2*ii-1,:)) 
  end do
  
  do ii = 2, m-1 
     sum2 = sum2 + 2.d0*(func(2*ii,1,:)+func(2*ii,pp,:)+func(1,2*ii,:)+func(pp,2*ii,:)) 
  end do

  do jj = 2, m
     do ii = 2, m
        sum16 = sum16 + 16.d0*func(2*ii-1,2*jj-1,:)
     end do
  end do

  do jj = 2, m-1
     do ii = 2, m
        sum8 = sum8 + 8.d0*(func(2*ii-1,2*jj,:)+func(2*jj,2*ii-1,:))
     end do
  end do

  do jj = 2, m-1
     do ii = 2, m-1
        sum42 = sum42 + 4.d0*func(2*ii,2*jj,:)
     end do
  end do

  
  fsum1 = func(1,1,:)+func(pp,1,:)+func(1,pp,:)+func(pp,pp,:)
 

  integ = h*k*(fsum1+sum4+sum2+sum16+sum8+sum42)/9.d0

  do ii = 1, points
     print *, ii, integ(ii)
  end do
  
 
  return

end subroutine T2D



  
