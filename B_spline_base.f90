subroutine B_spline_base(np,k,L,M,tl,tm,pp,x,y,base_L,base_M)

  implicit none

  !.. Input
  integer, intent(in) :: np,k,L,M,pp
  real(kind(1.d0)), intent(in) :: tl(np),tm(np),x(pp),y(pp)
  !.. Output
  real(kind(1.d0)), intent(inout) :: base_L(pp,L), base_M(pp,M)

  !.. Local
  real(kind(1.d0)) :: B_mj,B_lj
  integer :: ii,jj,lj,mj

  !.. External subroutine
  real(kind(1.d0)) :: bget

  base_L = 0.d0
  base_M = 0.d0

  ! do jj = 1, NN
  !    do ii = 1, pp
  !       base(ii,jj)=bget(x(ii),t,k,np,jj+1)
  !    end do
  ! end do

  do lj = 1, L
     do mj = 1, M
        do ii = 1, pp
           if(mj == 1)then
              B_mj = bget(y(ii),tm,k,np,1)+bget(y(ii),tm,k,np,2)
           else if(mj == M)then
              B_mj = bget(y(ii),tm,k,np,M+1)+bget(y(ii),tm,k,np,M+2)
           else
              B_mj = bget(y(ii),tm,k,np,mj+1)
           end if
           if(lj == 1)then
              B_lj = bget(x(ii),tl,k,np,1)+bget(x(ii),tl,k,np,2)
           else if(lj == L)then
              B_lj = bget(x(ii),tl,k,np,L+1)+bget(x(ii),tl,k,np,L+2)
           else
              B_lj = bget(x(ii),tl,k,np,lj+1)     
           end if
           
           base_L(ii,lj)=B_lj
           base_M(ii,mj)=B_mj
        end do
     end do
  end do

  return
end subroutine B_spline_base

  
