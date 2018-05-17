function bder(rr,t,kord,np,index)

  implicit none

  !.. Input
  integer, intent(in) :: kord,index,np
  real(kind(1.d0)), intent(in) :: t(np),rr

  !.. Output
  real(kind(1.d0)) :: bder

  !.. Local
  integer i,j,left,it
  real(kind(1.d0)) :: Sp(kord),deri

  bder=0.d0

  if(rr.gt.t(np)) return
  if(rr.lt.t(1)) return
  do it=1,np
     if(rr.ge.t(it)) left=it
  end do

  if(abs(rr-t(np)).lt.1.d-10) then
     if(index.lt.np-kord-1) return

     if(index.eq.np-kord) then
        bder=dble(kord-1)/(t(np)-t(np-kord))
     else if(index.eq.np-kord-1) then
        bder=-dble(kord-1)/(t(np)-t(np-kord))
     end if
     return
  end if

  if(index-left+kord.lt.1.or.index-left+kord.gt.kord) return

  call bsplvb(t,kord-1,1,rr,left,Sp)
  i=index-left+kord

  if(i.eq.1) then
     deri=dble(kord-1)*(-Sp(i)/(t(index+kord)-t(index+1)))
  else if(i.eq.kord) then
     deri=dble(kord-1)*(Sp(i-1)/(t(index+kord-1)-t(index)))
  else
     deri=dble(kord-1)*(Sp(i-1)/(t(index+kord-1)-t(index))-Sp(i)/(t(index+kord)-t(index+1)))
  end if
  bder=deri
  return
end function bder
