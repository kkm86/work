subroutine trapezoid_uniform(xmax,x,np,NN,k,t,f,rho_f,occ_1s,occ_2s,occ_2p,occ_3s,occ_3p,occ_3d,occ_4s)

  use constants

  implicit none

  !.. Input
  integer         , intent(in) :: k,np,xmax,NN
  real(kind(1.d0)), intent(in) :: x(xmax),t(np),f(NN,3),occ_1s,occ_2s,occ_2p,occ_3s,occ_3p,occ_3d,occ_4s
  real(kind(1.d0)), intent(inout) :: rho_f(xmax)

  !.. External functions
  real(kind(1.d0)) :: bget!,PI

  !.. Local
  integer          :: ii,jj
  real(kind(1.d0)) :: Baseone(xmax,NN),res1(xmax),res2(xmax),res3(xmax),res4(xmax),res5(xmax),res6(xmax),res7(xmax),usum1,usum2,usum3,usum4,usum5,usum6,usum7,u1,u2,u3,u4,u5,u6,u7

  DO jj = 1, NN
     DO ii = 1, xmax
        Baseone(ii,jj) = bget(x(ii),t,k,np,jj+1)
     END DO
  END DO

  res1=matmul(Baseone,f(:,1))
  res2=matmul(Baseone,f(:,2))
  res3=matmul(Baseone,f(:,3))
  res4=matmul(Baseone,f(:,4))
  res5=matmul(Baseone,f(:,5))
  res6=matmul(Baseone,f(:,6))
  res7=matmul(Baseone,f(:,7))
  usum1 = (x(xmax)-x(1))*(res1(1)**2.0+res1(xmax)**2.0)/(2.0*(xmax-1))
  usum2 = (x(xmax)-x(1))*(res2(1)**2.0+res2(xmax)**2.0)/(2.0*(xmax-1))
  usum3 = (x(xmax)-x(1))*(res3(1)**2.0+res3(xmax)**2.0)/(2.0*(xmax-1))
  usum4 = (x(xmax)-x(1))*(res4(1)**2.0+res4(xmax)**2.0)/(2.0*(xmax-1))
  usum5 = (x(xmax)-x(1))*(res5(1)**2.0+res5(xmax)**2.0)/(2.0*(xmax-1))
  usum6 = (x(xmax)-x(1))*(res6(1)**2.0+res6(xmax)**2.0)/(2.0*(xmax-1))
  usum7 = (x(xmax)-x(1))*(res7(1)**2.0+res7(xmax)**2.0)/(2.0*(xmax-1))
  

  DO ii = 2, xmax-1
     u1 = ((x(xmax)-x(1))/(xmax-1))*res1(ii)**2.0
     u2 = ((x(xmax)-x(1))/(xmax-1))*res2(ii)**2.0
     u3 = ((x(xmax)-x(1))/(xmax-1))*res3(ii)**2.0
     u4 = ((x(xmax)-x(1))/(xmax-1))*res4(ii)**2.0
     u5 = ((x(xmax)-x(1))/(xmax-1))*res5(ii)**2.0
     u6 = ((x(xmax)-x(1))/(xmax-1))*res6(ii)**2.0
     u7 = ((x(xmax)-x(1))/(xmax-1))*res7(ii)**2.0
     usum1 = usum1 + u1
     usum2 = usum2 + u2
     usum3 = usum3 + u3
     usum4 = usum4 + u4
     usum5 = usum5 + u5
     usum6 = usum6 + u6
     usum7 = usum7 + u7
  END DO

  rho_f=(occ_1s*(res1**2.0)/usum1 + occ_2s*(res2**2.0)/usum2 + occ_2p*(res3**2.0)/usum3 + occ_3s*(res4**2.0)/usum4 + occ_3p*(res5**2.0)/usum5 + occ_3d*(res6**2.0)/usum6 + occ_4s*(res7**2.0)/usum7)/(4.0*(x**2.0)*PI)
   
  rho_f(1) = 0.0

  return

end subroutine trapezoid_uniform
     
      
