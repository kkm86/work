subroutine adiabatic(npl,npm,k,L,M,LM,tl,tm,rho,my,c,cv,U,points)
          

  use constants
  
  implicit none

  !.. Input
  integer            , intent(in) :: npl,npm,k,L,M,LM,points
  real(kind(1.d0))   , intent(in) :: tl(npl),tm(npm),rho(points),c(L,M,points),cv(L,M,points),my
  real(kind(1.d0)), intent(inout) :: U(points)
  
  !.. Local
  real(kind(1.d0))   :: coordl(L+2,k),coordm(M+2,k),xabsc(k),weig(k),theta
  real(kind(1.d0))   :: term1(points),lowerl,upperl,lowerm,upperm,t1,t2
  real(kind(1.d0))   :: sum1(points), sum2(points), sum3(points)
  real(kind(1.d0))   :: B_li, B_lj, B_mi, B_mj
  real(kind(1.d0))   :: nsize
  integer            :: li, lj, mi, mj, n, p, i, j, ll, mm

  
  !.. External functions/variables
  real(kind(1.d0)) :: bget, bder

 
  call gauleg(k, xabsc, weig)

  U = 0.d0

  call CPU_TIME( t1 )

 
  do li = 1, L
     do mj = 1, M
        sum3 = 0.0
        upperl = li+k
        lowerl = li
        upperm = mj+k
        lowerm = mj
        do ll = lowerl, upperl
           sum2 = 0.0
           do mm = lowerm, upperm
              sum1 = 0.0
              if((tl(ll+1).gt.tl(ll)) .and. (tm(mm+1).gt.tm(mm)))then
                 do n = 1, k
                    term1 = 0.0
                    do p = 1, k
                       coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
                       coordm(mm,p) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(p)
                       theta = coordl(ll,n)
              
                             
                       if(li == 1)then
                          B_li = bget(coordl(ll,n),tl,k,npl,1)+bget(coordl(ll,n),tl,k,npl,2)
                       else if(li == L)then
                          B_li = bget(coordl(ll,n),tl,k,npl,L+1)+bget(coordl(ll,n),tl,k,npl,L+2)
                       else
                          B_li = bget(coordl(ll,n),tl,k,npl,li+1)
                       end if
                                   
                             
                       if(mj == 1 .and. mj /= M)then
                          B_mj = bget(coordm(mm,p),tm,k,npm,1)+bget(coordm(mm,p),tm,k,npm,2)
                       else if(mj == M)then
                          B_mj = bget(coordm(mm,p),tm,k,npm,M+1)+bget(coordm(mm,p),tm,k,npm,M+2)
                       else
                          B_mj = bget(coordm(mm,p),tm,k,npm,mj+1)
                       end if
                             
                    

                       term1 = term1 + 2.d0*Pi*(cv(li,mj,:)*cv(li,mj,:))*B_mj*B_mj*B_li*B_li*sin(2.d0*theta)
                             
                    end do
                    sum1 = sum1 + 0.5d0*(tm(mm+1)-tm(mm))*(term1)
                 end do
                 sum2 = sum2 + 0.5d0*(tl(ll+1)-tl(ll))*sum1
              end if
           end do
           sum3 = sum3 + sum2
        end do
        U = sum3
     end do
  end do
 
   
  call CPU_TIME( t2 )
  print*,'making array', t2-t1

  do i = 1, points
     print*, 'integral value',i, rho(i),U(i)
  end do
  
  

  
  
  return

end subroutine adiabatic





