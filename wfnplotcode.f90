  
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
  
