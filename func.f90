 subroutine func(z,l,np,k,NN,t,Vold,Ham,Bs)

   implicit none

   !.. Input
    integer         , intent(in) :: np,k,NN          !Number of knot-points etc.
    real(kind(1.d0)), intent(in) :: z, l
    real(kind(1.d0)), intent(in) :: Vold(NN,NN), t(np)
    real(kind(1.d0)), intent(inout) :: Ham(NN,NN), Bs(NN,NN) !Hamiltonian and B-matrix

    !.. Local
    real(kind(1.d0)) :: mn(NN,k) !Abscissa-grid
    real(kind(1.d0)) :: sum1, sum2,term1,term2,term3,lower,upper,sumbsp,sumbsp2,bsp,xabsc(k),weig(k),xval
    INTEGER :: ii, jj, m, n

    !.. External functions
    real(kind(1.d0)) :: bget, bder


    CALL gauleg(k, xabsc, weig)
  
    DO jj = 1, NN
       DO ii = 1, NN
          sum1 = 0.0
          sum2 = 0.0
          sumbsp = 0.0
          sumbsp2 = 0.0
          bsp = 0.0
          term1 = 0.0
          term2 = 0.0
          term3 = 0.0
          upper = min(ii+1,jj+1) + k-1
          lower = max(ii+1,jj+1)
          do m = lower, upper
             term1 = 0.0
             term2 = 0.0
             term3 = 0.0
             bsp = 0.0
             do n = 1, k
                if(t(m+1).gt.t(m)) then
                   mn(m,n) = 0.5*(t(m+1)+t(m)) + 0.5*(t(m+1)-t(m))*xabsc(n)
                   xval = mn(m,n)
                   term1 = term1 + weig(n)*bder(xval,t,k,np,ii+1)*bder(xval,t,k,np,jj+1)/2.d0
                   term2 = term2 + weig(n)*l*(l+1.d0)*bget(xval,t,k,np,jj+1)*bget(xval,t,k,np,ii+1)/(2.d0*(xval**2.d0))
                   term3 = term3 - weig(n)*z*bget(xval,t,k,np,ii+1)*bget(xval,t,k,np,jj+1)/xval
                   bsp = bsp + weig(n)*bget(xval,t,k,np,ii+1)*bget(xval,t,k,np,jj+1)
                end if
             end do
             sum1 = sum1 + 0.5*(t(m+1)-t(m))*(term1+term2+term3)
             sumbsp = sumbsp + 0.5*(t(m+1)-t(m))*bsp
          end do
          sum2 = sum2 + sum1
          sumbsp2 = sumbsp2 + sumbsp
          Ham(ii,jj) = sum2
          Bs(ii,jj) = sumbsp2
       end do
    end do
    Ham = Ham + Vold
    
   
    return
  end subroutine func
