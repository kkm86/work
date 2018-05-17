  function bder2(rr,t,kord,np,index)
!*     returns the value of (d/dx) Bspline(kord,index) in rr
!*     The first Bspline is called spline #1 (i.e. index=1)
!*     the first knot point is in t(1)
!*     np= number of knot points (distinct or multiple) including
!*     the ghost points: N phyical points np=N +2*(kord-1)  

      implicit none
      integer i,j,left,it,kord,index,np
      real*8 bder2,t(np),Sp(kord),rr,deri
     
      bder2=0.d0
!*     if rr=t(np) then the routine assumes that
!*     there is kord knotpoints in the last physical point and
!*     returns Bder.ne.zero if index is np-kord, or np-kord-1 
      if(rr.gt.t(np)) return
      if(rr.lt.t(1)) return
      do it=1,np
        if(rr.ge.t(it)) left=it 
      end do   

      if(abs(rr-t(np)).lt.1.d-10) then
        if(index.lt.np-kord-2) return       
        if(index.eq.np-kord) then
          bder2=dble((kord-1)*(kord-2))/(t(np-1)-t(np-kord))/(t(np-2)-t(np-kord))
        else if(index.eq.np-kord-2) then
          bder2=dble((kord-1)*(kord-2))/(t(np-1)-t(np-kord-1))/(t(np-2)-t(np-kord))
        else if(index.eq.np-kord-1) then
          bder2=-dble((kord-1)*(kord-2))/(t(np-2)-t(np-kord-1))/(t(np-2)-t(np-kord))-dble((kord-1)*(kord-2))/(t(np-1)-t(np-kord))/(t(np-2)-t(np-kord))
        end if
        return
      end if
      

      if(index-left+kord.lt.1.or.index-left+kord.gt.kord) return   
      
!*     index=left-kord+i => i=index-left+kord  for Sp with k=kord
!*     index=left-kord+i+1 => i=index-left+kord-1  for Sp with k=kord-1
!*     Sp(i-1,k-1) gives the same index as Sp(i,k)
!*     Sp(i,k-1) gives the same index as Sp(i+1.k)
      call bsplvb(t,kord-2,1,rr,left,Sp)
      i=index-left+kord
      deri=0.d0
      if(i.gt.2) deri=deri+dble(kord-1)*dble(kord-2)*Sp(i-2)/(t(index+kord-2)-t(index))/(t(index+kord-1)-t(index))
      if(i.gt.1.and.i.lt.kord) deri=deri-dble(kord-1)*dble(kord-2)*(Sp(i-1)/(t(index+kord-1)-t(index+1))/(t(index+kord-1)-t(index))+Sp(i-1)/(t(index+kord-1)-t(index+1))/(t(index+kord)-t(index+1)))
       if(i.lt.kord-1) deri=deri+dble(kord-1)*dble(kord-2)*Sp(i)/(t(index+kord)-t(index+2))/(t(index+kord)-t(index+1))     
      bder2=deri
      return
    end function bder2
