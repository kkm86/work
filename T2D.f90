subroutine T2D(m,h,k,func,integ)
          

  use constants
  
  implicit none

  !.. This subroutine performs integration of a two variable function
  !.. Using Trapezoidal 2D Rule
  
  !.. Input
  integer         , intent(in) :: m         !.. # subintervals in x and y
  real(kind(1.d0)), intent(in) :: func(m,m) !.. function to integrate
  real(kind(1.d0)), intent(in) :: h,k           !.. subinterval width h:x, k:y
  
  !.. Output
  real(kind(1.d0)), intent(inout) :: integ  !.. result of integration

  !.. Local
  integer :: ii, jj
  real(kind(1.d0)) :: summ1, summ2, sumn1, sumn2, summn
  real(kind(1.d0)) :: fsum1,fsum2,fsum4

  summ1 = 0.d0
  summ2 = 0.d0
  sumn1 = 0.d0
  sumn2 = 0.d0
  summn = 0.d0
  
  do ii = 2, m-1 
     summ1 = summ1 + func(ii,1)
     summ2 = summ2 + func(ii,m)
     sumn1 = sumn1 + func(1,ii)
     sumn2 = sumn2 + func(m,ii)  
  end do

  do jj = 2, m-1
     do ii = 2, m-1
        summn = summn + func(ii,jj)
     end do
  end do

  
  fsum1 = func(1,1)+func(m,1)+func(1,m)+func(m,m)
  fsum2 = 2.d0*(summ1+summ2+sumn1+sumn2)
  fsum4 = 4.d0*summn

  integ = 0.25d0*h*k*(fsum1+fsum2+fsum4)

  print *, 'T2D',integ
 
  return

end subroutine T2D



  
