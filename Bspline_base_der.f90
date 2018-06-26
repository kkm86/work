subroutine Bspline_base_der(npl,npm,k,L,M,tl,tm,pp,x,y,base_L,base_M)

  implicit none

  !.. Input
  integer, intent(in) :: npl,npm,k,L,M,pp
  real(kind(1.d0)), intent(in) :: tl(npl),tm(npm),x(pp),y(pp)
  !.. Output
  real(kind(1.d0)), intent(inout) :: base_L(pp,L), base_M(pp,M)

  !.. Local
  real(kind(1.d0)) :: B_mj,B_lj
  integer :: ii,jj,lj,mj

  !.. External subroutine
  real(kind(1.d0)) :: bder

  base_L = 0.d0
  base_M = 0.d0

 

  do lj = 1, L
     do mj = 1, M
        do ii = 1, pp
           if(mj == 1)then
              B_mj = bder(y(ii),tm,k,npm,1)+bder(y(ii),tm,k,npm,2)
           else if(mj == M)then
              B_mj = bder(y(ii),tm,k,npm,M+1)+bder(y(ii),tm,k,npm,M+2)
           else
              B_mj = bder(y(ii),tm,k,npm,mj+1)
           end if
           if(lj == 1)then
              B_lj = bder(x(ii),tl,k,npl,1)+bder(x(ii),tl,k,npl,2)
           else if(lj == L)then
              B_lj = bder(x(ii),tl,k,npl,L+1)+bder(x(ii),tl,k,npl,L+2)
           else
              B_lj = bder(x(ii),tl,k,npl,lj+1)     
           end if
           
           base_L(ii,lj)=B_lj
           base_M(ii,mj)=B_mj
        end do
     end do
  end do

  return
end subroutine Bspline_base_der


  
