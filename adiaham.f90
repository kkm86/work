subroutine adiaham(npl,npm,k,L,M,LM,tl,tm,rho,my,r0,d,mass,energy,H,S,points)
          

  use constants
  
  implicit none

  !.. Input
  integer            , intent(in) :: npl,npm,k,L,M,LM,points
  real(kind(1.d0))   , intent(in) :: tl(npl),tm(npm),rho(3,points),my,d,r0,mass(3)
  real(kind(1.d0)), intent(inout) :: energy(LM,points,3)
  real(kind(1.d0)), intent(inout) :: H(LM,LM,points,3), S(LM,LM,points,3)
  
  !.. Local
  real(kind(1.d0))   :: coordl(L+2,k),coordm(M+2,k),xabsc(k),weig(k),theta,phi
  real(kind(1.d0))   :: term1(3,points),term2(3,points),term3(3,points),term4(3,points),lowerl,upperl,lowerm,upperm,t1,t2
  real(kind(1.d0))   :: sum1(3,points), sum2(3,points), sum3(3,points), sumbsp1(3,points), sumbsp2(3,points), sumbsp3(3,points), bsp(3,points)
  real(kind(1.d0))   :: B_li, B_lj, B_mi, B_mj, dB_lj, dB_mj, dB_li, dB_mi
  real(kind(1.d0))   :: nsize
  real(kind(1.d0))   :: Hrez(LM,LM),Srez(LM,LM)
  integer            :: li, lj, mi, mj, n, p, i, j, ll, mm

  !.. Paramenters for generalized eigensolver
  integer                             :: ITYPE, LDA, LDB, INFO 
  integer                 , parameter :: LWORK = 168777
  integer                 , parameter :: LIWORK = 1448
  character*1                         :: JOBZ = 'V', UPLO = 'U'
  ! integer                 , parameter :: LWORK = 8000
  ! integer                 , parameter :: LIWORK = 100
  ! character*1                         :: JOBZ = 'N', UPLO = 'U'
  real(kind(1.d0))    , dimension(LM) :: W
  real(kind(1.d0)) , dimension(LWORK) :: WORK
  integer         , dimension(LIWORK) :: IWORK

  !.. External functions/variables
  real(kind(1.d0)) :: bget, bder, V(3,points)

  ITYPE = 1
  LDA = LM
  LDB = LM

  call gauleg(k, xabsc, weig)

  H = 0.0d0
  S = 0.0d0

  call CPU_TIME( t1 )
  do jj = 1, points
     do ii = 1, points
        sum1 = 0.0d0
        sum2 = 0.0d0
        sumbsp1 = 0.0d0
        sumbsp2 = 0.0d0
        upper = min(ii+1,jj+1) + k-1
        lower = max(ii+1,jj+1)
        do m = lower, upper
           bsp = 0.0d0
           term1 = 0.0d0
           term2 = 0.0d0 
           if(t(m+1).gt.t(m))then
              do n = 1, k
                 coord(ll,n) = 0.5*(t(m+1)+t(m)) + 0.5*(t(m+1)-t(m))*xabsc(n)
                             coordm(m,n) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(p)
                             theta = coordl(ll,n)
                             phi = coordm(mm,p)
                           
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
                             
                             call twobody_potential(d,r0,mass,rho,theta,phi,V,points)

                             bsp = bsp + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*sin(2.d0*theta)
                             
                             term1 = term1 + weig(n)*weig(p)*4.d0*dB_li*B_mi*dB_lj*B_mj*sin(2.d0*theta)/(2.d0*my*rho**2.d0)
                             term2 = term2 + weig(n)*weig(p)*8.d0*B_li*dB_mi*B_lj*dB_mj*cos(theta)/(2.d0*sin(theta)*my*rho**2.d0)
                             term3 = term3 + weig(n)*weig(p)*15.d0*B_li*B_mi*B_lj*B_mj*sin(2.d0*theta)/(8.d0*my*rho**2.d0)
                             term4 = term4 + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*V*sin(2.d0*theta)
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
              H(i,j,:) = sum3(:)
              S(i,j,:) = sumbsp3(:)
           end do
        end do
  call CPU_TIME( t2 )
  print*,'making array', t2-t1

  

  call CPU_TIME( t1 )
  do j = 1, 3
     do i = 1, points
        Hrez = H(:,:,i,j)
        Srez = S(:,:,i,j)
        call dsygvd( ITYPE, JOBZ, UPLO, LM, Hrez, LDA, Srez, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )
        energy(:,i,j) = W
        H(:,:,i,j) = Hrez
     end do
  end do

  call CPU_TIME( t2 )
  print*,'time2', t2-t1

  !print*, W(1), W(2), W(3), W(4), W(5), W(6), W(7), W(8), W(9), W(10), W(11), W(12), W(13), W(14), W(15), W(16)
  print*, 'info', INFO, 'LWORK', WORK(1), 'LIWORK', IWORK(1)
  
  return

end subroutine adiaham




