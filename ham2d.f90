subroutine ham2d(np,k,L,M,LM,tl,tm,energy2)

  use constants
  
  implicit none

  !.. Input
  integer            , intent(in) :: np,k,L,M,LM
  real(kind(1.d0))   , intent(in) :: tl(np),tm(np)
  real(kind(1.d0)), intent(inout) :: energy2(3)

  !.. Local
  real(kind(1.d0))   :: coordl(L,k),coordm(M+2,k),xabsc(k),weig(k)
  real(kind(1.d0))   :: H(LM,LM), S(LM,LM)
  real(kind(1.d0))   :: term1,term2,term3,term4,lowerl,upperl,lowerm,upperm
  real(kind(1.d0))   :: sum1, sum2, sum3, sumbsp1, sumbsp2, sumbsp3, bsp
  real(kind(1.d0))   :: B_li, B_lj, B_mi, B_mj, dB_lj, dB_mj, dB_li, dB_mi
  real(kind(1.d0))   :: nsize
  integer            :: li, lj, mi, mj, n, p, i, j, ll, mm

  !.. Paramenters for generalized eigensolver
  integer                             :: ITYPE, LDA, LDB, INFO 
  integer                 , parameter :: LWORK = 1000000
  integer                 , parameter :: LIWORK = 5000
  character*1                         :: JOBZ = 'V', UPLO = 'U'
  real(kind(1.d0))    , dimension(LM) :: W
  real(kind(1.d0)) , dimension(LWORK) :: WORK
  real(kind(1.d0)), dimension(LIWORK) :: IWORK

  !.. External functions/variables
  real(kind(1.d0)) :: bget, bder

  ITYPE = 1
  LDA = LM
  LDB = LM

  call gauleg(k, xabsc, weig)

  

  H = 0.0
  S = 0.0
 

  do lj = 1, L
     do li = 1, L
        do mj = 1, M
           do mi = 1, M
              i = (li-1)*M+mi
              j = (lj-1)*M+mj
              sum3 = 0.0
              sumbsp3 = 0.0
              upperl = min(li,lj) + k-1
              lowerl = max(li,lj)
              upperm = min(mi+1,mj+1) + k-1
              lowerm = max(mi,mj)
              do ll = lowerl, upperl
                 sum2 = 0.0
                 sumbsp2 = 0.0
                 do mm = lowerm, upperm
                    sum1 = 0.0
                    sumbsp1 = 0.0
                    if((tl(ll+1).gt.tl(ll)) .and. (tm(mm+1).gt.tm(mm)))then
                       do n = 1, k
                          term1 = 0.0
                          term2 = 0.0
                          term3 = 0.0
                          term4 = 0.0
                          bsp = 0.0
                          do p = 1, k
                             coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
                             coordm(mm,p) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(p)
                             B_li = bget(coordl(ll,n),tl,k,np,li)
                             B_lj = bget(coordl(ll,n),tl,k,np,lj)
                             dB_li = bder(coordl(ll,n),tl,k,np,li)
                             dB_lj = bder(coordl(ll,n),tl,k,np,lj)
    
                              if(mj == 1)then
                                B_mj = bget(coordm(mm,p),tm,k,np,1)+bget(coordm(mm,p),tm,k,np,2)
                                dB_mj = bder(coordm(mm,p),tm,k,np,1)+bder(coordm(mm,p),tm,k,np,2)
                             else
                                B_mj = bget(coordm(mm,p),tm,k,np,mj+1)
                                dB_mj = bder(coordm(mm,p),tm,k,np,mj+1)
                             end if
                                
                             if(mi == 1) then
                                B_mi = bget(coordm(mm,p),tm,k,np,1)+bget(coordm(mm,p),tm,k,np,2)
                                dB_mi = bder(coordm(mm,p),tm,k,np,1)+bder(coordm(mm,p),tm,k,np,2)
                             else
                                B_mi = bget(coordm(mm,p),tm,k,np,mi+1)
                                dB_mi = bder(coordm(mm,p),tm,k,np,mi+1)
                             end if
                           
                             term1 = term1 + weig(n)*weig(p)*0.5d0*dB_li*B_mi*dB_lj*B_mj*coordl(ll,n)
                             term2 = term2 + weig(n)*weig(p)*0.5d0*B_li*dB_mi*B_lj*dB_mj*coordl(ll,n)
                             term3 = term3 + weig(n)*weig(p)*0.5d0*B_li*B_mi*B_lj*B_mj*(coordl(ll,n)**3.d0)
                             term4 = term4 + weig(n)*weig(p)*2.d0*B_li*B_mi*B_lj*B_mj*(coordm(mm,p)**2.d0)*coordl(ll,n)
                             bsp = bsp + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*coordl(ll,n)
                          end do
                          sum1 = sum1 + 0.5d0*(tm(mm+1)-tm(mm))*(term1+term2+term3+term4)
                          sumbsp1 = sumbsp1 + 0.5d0*(tm(mm+1)-tm(mm))*bsp
                       end do
                       sum2 = sum2 + 0.5d0*(tl(ll+1)-tl(ll))*sum1
                       sumbsp2 = sumbsp2 + 0.5d0*(tl(ll+1)-tl(ll))*sumbsp1
                    end if
                 end do
                 sum3 = sum3 + sum2
                 sumbsp3 = sumbsp3 + sumbsp2
              end do
              H(i,j) = sum3
              S(i,j) = sumbsp3
           end do
        end do
     end do
  end do
 
  call dsygvd( ITYPE, JOBZ, UPLO, LM, H, LDA, S, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )
  print *, W(1), W(2), W(3)
  print *, 'Good job!'

  return

end subroutine ham2d
