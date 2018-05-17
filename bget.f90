function bget(rr,t,kord,np,index)

  implicit none

  !.. Input
  integer :: kord,index,np
  real(kind(1.d0)) :: t(np),rr

  !.. Output
  real(kind(1.d0)) :: bget
  
  ! .. Local
  integer i,j,left,it
  real(kind(1.d0)) :: Sp(kord)

  bget=0.d0

  if(rr.gt.t(np)) return
  if(rr.lt.t(1)) return
  if(abs(rr-t(np)).lt.1.d-10) then
     if(index.eq.np-kord) bget=1.d0
     return
  end if
  do it=1,np
     if(rr.ge.t(it)) left=it
  end do
  if(index-left+kord.lt.1.or.index-left+kord.gt.kord) return
  call bsplvb(t,kord,1,rr,left,Sp)
  bget=Sp(index-left+kord)
  return
end function bget
