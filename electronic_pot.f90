subroutine electronic_pot(np,NN,k,t,xmax,kord,np1,t1,xmax_test,xabsc,weig,RHS,Vexc)

  use constants

  implicit none

  !.. Input
  integer, intent(in) :: np,NN,k,xmax,kord,np1,xmax_test
  real(kind(1.d0)),intent(in) :: xabsc(k),weig(k),t(np),t1(np1),RHS(xmax_test)
  real(kind(1.d0)),intent(inout) :: Vexc(NN,NN)

  !.. External functions
  real(kind(1.d0)) :: bget,bder,bder2,ddot

  !.. Local
  real(kind(1.d0)) :: sum1,sum2,term4,term5,base,pot,baseder,potex,potexder2,baseder2,vvv,mn(NN,k),xval,x
  integer :: i,j,m,n,lower,upper

  !.. Local for ddot
  real(kind(1.d0)) :: Bord_pot(xmax_test),Bder_pot(xmax_test),Bder2_pot(xmax_test)
  integer :: incx=1,incy=1

  Vexc = 0.0
  mn = 0.0

  do j = 1, NN
     do i = 1, NN
        sum1 = 0.d0
        sum2 = 0.d0
        term4 = 0.d0
        term5 = 0.d0
        base = 0.d0
        pot = 0.d0
        baseder = 0.d0
        potex = 0.d0
        potexder2 = 0.d0
        baseder2 = 0.d0
        vvv = 0.d0
        upper = min(i+1,j+1) + k-1 
        lower = max(i+1,j+1)
        do m = lower, upper
           term4 = 0.d0
           term5 = 0.d0
           base = 0.d0
           pot = 0.d0
           baseder = 0.d0
           potex = 0.d0
           potexder2 = 0.d0
           baseder2 = 0.d0
           vvv = 0.d0
           do n = 1, k
              if(t(m+1).gt.t(m)) then
                 mn(m,n) = 0.5d0*(t(m+1)+t(m)) + 0.5d0*(t(m+1)-t(m))*xabsc(n)
                 xval = mn(m,n)
                 call baseR(xval,xmax,kord,np1,t1,Bord_pot,Bder_pot,Bder2_pot)
                 base = ddot(xmax_test,Bord_pot,incx,RHS,incy)
                 pot=base/xval
                 term4 = term4 + weig(n)*bget(xval,t,k,np,i+1)*bget(xval,t,k,np,j+1)*pot
                ! baseder = ddot(xmax_test,Bder_pot,incx,RHS,incy)
                ! potex=baseder/(xval**2.d0)
                 baseder2 = ddot(xmax_test,Bder2_pot,incx,RHS,incy)
                 potexder2 = -baseder2/(4.d0*PI*xval)
                 vvv = (3.d0*potexder2/(8.d0*PI))**(1.d0/3.d0)
                 term5 = term5 + weig(n)*bget(xval,t,k,np,i+1)*bget(xval,t,k,np,j+1)*(-3.d0*vvv)
              end if
           end do
           sum1 = sum1 + 0.5d0*(t(m+1)-t(m))*(term4+term5)
        end do
        sum2 = sum2 + sum1
        Vexc(i,j) = sum2
     end do
  end do

  return

end subroutine electronic_pot
