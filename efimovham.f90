subroutine efimovham(npl,npm,k,L,M,LM,tl,tm,rho,my,r0,d,mass,energy2,H,S,ny)
          

  use constants
  
  implicit none

  !.. Input
  integer            , intent(in) :: npl,npm,k,L,M,LM
  real(kind(1.d0))   , intent(in) :: tl(npl),tm(npm),rho,my,d,r0,mass(3),ny
  real(kind(1.d0)), intent(inout) :: energy2(6)
  real(kind(1.d0)), intent(inout) :: H(LM,LM), S(LM,LM)
  
  !.. Local
  real(kind(1.d0))   :: coordl(L,k),coordm(M+2,k),xabsc(k),weig(k),theta,phi
  real(kind(1.d0))   :: term1,term2,term3,term4,term5,lowerl,upperl,lowerm,upperm
  real(kind(1.d0))   :: sum1, sum2, sum3, sumbsp1, sumbsp2, sumbsp3, bsp
  real(kind(1.d0))   :: B_li, B_lj, B_mi, B_mj, dB_lj, dB_mj, dB_li, dB_mi
  real(kind(1.d0))   :: nsize
  integer            :: li, lj, mi, mj, n, p, i, j, ll, mm

  !.. Paramenters for generalized eigensolver
  integer                             :: ITYPE, LDA, LDB, INFO 
  !integer                 , parameter :: LWORK = 4000000
  !integer                 , parameter :: LIWORK = 10000
  !character*1                         :: JOBZ = 'V', UPLO = 'U'
  integer                 , parameter :: LWORK = 8000
  integer                 , parameter :: LIWORK = 100
  character*1                         :: JOBZ = 'N', UPLO = 'U'
  real(kind(1.d0))    , dimension(LM) :: W
  real(kind(1.d0)) , dimension(LWORK) :: WORK
  real(kind(1.d0)), dimension(LIWORK) :: IWORK

  !.. External functions/variables
  real(kind(1.d0)) :: bget, bder, V

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
              upperl = min(li+1,lj+1) + k-1
              lowerl = max(li+1,lj+1)
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
                          term5 = 0.0
                          bsp = 0.0
                          do p = 1, k
                             coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
                             coordm(mm,p) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(p)
                             theta = coordl(ll,n)
                             phi = coordm(mm,p)
                             ! B_li = bget(coordl(ll,n),tl,k,npl,li+1)
                             ! B_lj = bget(coordl(ll,n),tl,k,npl,lj+1)
                             ! dB_li = bder(coordl(ll,n),tl,k,npl,li+1)
                             ! dB_lj = bder(coordl(ll,n),tl,k,npl,lj+1)
                             ! B_mi = bget(coordm(mm,p),tm,k,npm,mi+1)
                             ! B_mj = bget(coordm(mm,p),tm,k,npm,mj+1)
                             ! dB_mi = bder(coordm(mm,p),tm,k,npm,mi+1)
                             ! dB_mj = bder(coordm(mm,p),tm,k,npm,mj+1)

                             if(lj == 1)then
                                B_lj = bget(coordl(ll,n),tl,k,npl,1)+bget(coordl(ll,n),tl,k,npl,2)
                                dB_lj = bder(coordl(ll,n),tl,k,npl,1)+bder(coordl(ll,n),tl,k,npl,2)
                             else if(lj == L)then
                                B_lj = bget(coordl(ll,n),tl,k,npl,L+1)+bget(coordl(ll,n),tl,k,npl,L+2)
                                dB_lj = bder(coordl(ll,n),tl,k,npl,L+1)+bder(coordl(ll,n),tl,k,npl,L+2)
                             else
                                B_lj = bget(coordl(ll,n),tl,k,npl,lj+1)
                                dB_lj = bder(coordl(ll,n),tl,k,npl,lj+1)
                             end if

                             if(li == 1)then
                                B_li = bget(coordl(ll,n),tl,k,npl,1)+bget(coordl(ll,n),tl,k,npl,2)
                                dB_li = bder(coordl(ll,n),tl,k,npl,1)+bder(coordl(ll,n),tl,k,npl,2)
                             else if(li == L)then
                                B_li = bget(coordl(ll,n),tl,k,npl,L+1)+bget(coordl(ll,n),tl,k,npl,L+2)
                                dB_li = bder(coordl(ll,n),tl,k,npl,L+1)+bder(coordl(ll,n),tl,k,npl,L+2)
                             else
                                B_li = bget(coordl(ll,n),tl,k,npl,li+1)
                                dB_li = bder(coordl(ll,n),tl,k,npl,li+1)
                             end if

                             
                             if(mj == 1 .and. mj /= M)then
                                B_mj = bget(coordm(mm,p),tm,k,npm,1)+bget(coordm(mm,p),tm,k,npm,2)
                                dB_mj = bder(coordm(mm,p),tm,k,npm,1)+bder(coordm(mm,p),tm,k,npm,2)
                             else if(mj == M)then
                                B_mj = bget(coordm(mm,p),tm,k,npm,M+1)+bget(coordm(mm,p),tm,k,npm,M+2)
                                dB_mj = bder(coordm(mm,p),tm,k,npm,M+1)+bder(coordm(mm,p),tm,k,npm,M+2)
                             else
                                B_mj = bget(coordm(mm,p),tm,k,npm,mj+1)
                                dB_mj = bder(coordm(mm,p),tm,k,npm,mj+1)
                             end if

                             if(mi == 1 .and. mi /= M)then
                                B_mi = bget(coordm(mm,p),tm,k,npm,1)+bget(coordm(mm,p),tm,k,npm,2)
                                dB_mi = bder(coordm(mm,p),tm,k,npm,1)+bder(coordm(mm,p),tm,k,npm,2)
                             else if(mi == M)then
                                B_mi = bget(coordm(mm,p),tm,k,npm,M+1)+bget(coordm(mm,p),tm,k,npm,M+2)
                                dB_mi = bder(coordm(mm,p),tm,k,npm,M+1)+bder(coordm(mm,p),tm,k,npm,M+2)
                             else
                                B_mi = bget(coordm(mm,p),tm,k,npm,mi+1)
                                dB_mi = bder(coordm(mm,p),tm,k,npm,mi+1)
                             end if

                             call twobody_potential(d,r0,mass,rho,theta,phi,V)

                             bsp = bsp + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*sin(2.d0*theta)

                             term1 = term1 + weig(n)*weig(p)*4.d0*dB_li*B_mi*dB_lj*B_mj*sin(2.d0*theta)/(2.d0*my*rho**2.d0)
                             term2 = term2 + weig(n)*weig(p)*8.d0*B_li*dB_mi*B_lj*dB_mj*cos(theta)/(2.d0*sin(theta)*my*rho**2.d0)
                             term3 = term3 + weig(n)*weig(p)*15.d0*B_li*B_mi*B_lj*B_mj*sin(2.d0*theta)/(8.d0*my*rho**2.d0)
                             term4 = term4 + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*V*sin(2.d0*theta)

                             !term1 = term1 + weig(n)*weig(p)*4.d0*dB_li*B_mi*dB_lj*B_mj!/(2.d0*my*rho**2.d0)
                             !term2 = term2 + weig(n)*weig(p)*4.d0*B_li*dB_mi*B_lj*dB_mj/(sin(theta)**2.d0)
                             !term3 = term3 + weig(n)*weig(p)*15.d0*B_li*B_mi*B_lj*B_mj/8.d0
                             !term4 = term4 + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*V*2.d0*my*rho**2.d0
                             

                             ! term1 = term1 + weig(n)*weig(p)*4.d0*dB_li*B_mi*dB_lj*B_mj/(2.d0*my*rho**2.d0)
                             ! term2 = term2 + weig(n)*weig(p)*4.d0*B_li*dB_mi*B_lj*dB_mj/((2.d0*sin(theta)**2.d0)*my*rho**2.d0)
                             ! term3 = term3 - weig(n)*weig(p)*1.d0*B_li*B_mi*B_lj*B_mj/(8.d0*my*rho**2.d0)
                             !term4 = term4 + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*V

     
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
 
  energy2(1) = W(1)
  energy2(2) = W(2)
  energy2(3) = W(3)
  energy2(4) = W(4)
  energy2(5) = W(5)
  energy2(6) = W(6)
  

  !print*, W(1), W(2), W(3), W(4), W(5), W(6), W(7), W(8), W(9), W(10), W(11), W(12), W(13), W(14), W(15), W(16)
  !print*, 'info', INFO

 
  
  return

end subroutine efimovham

