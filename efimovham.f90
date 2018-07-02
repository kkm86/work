subroutine efimovham(npl,npm,k,L,M,LM,tl,tm,rho,my,r0,d,mass,energy,H,Hder,S,Integ,points,Pmat,P2mat,Imat)
          

  use constants
  
  implicit none

  !.. Input
  integer            , intent(in) :: npl,npm,k,L,M,LM,points
  real(kind(1.d0))   , intent(in) :: tl(npl),tm(npm),rho(points),my,d,r0,mass(3)
  real(kind(1.d0)), intent(inout) :: energy(LM,points)
  real(kind(1.d0)), intent(inout) :: H(LM,LM,points),Hder(LM,LM,points), S(LM,LM,points),Integ(LM,LM,points),Pmat(LM,LM,points),P2mat(LM,LM,points),Imat(LM,LM,points)
  
  !.. Local
  real(kind(1.d0))   :: coordl(L+2,k),coordm(M+2,k),xabsc(k),weig(k),theta,phi
  real(kind(1.d0))   :: term(8,points)!,term2(points),term3(points),term4(points),term5(points),term6(points),term2(points),term3(points),term4(points)
  real(kind(1.d0))   :: lowerl,upperl,lowerm,upperm,t1,t2
  real(kind(1.d0))   :: sum1(points), sum2(points), sum3(points), sumbsp1(points), sumbsp2(points), sumbsp3(points), bsp(points),sumder1(points),sumder2(points),sumder3(points),sumint1(points),sumint2(points),sumint3(points)
  real(kind(1.d0))   :: B_li, B_lj, B_mi, B_mj, dB_lj, dB_mj, dB_li, dB_mi
  real(kind(1.d0))   :: nsize
  real(kind(1.d0))   :: Hrez(LM,LM),Srez(LM,LM)
  integer            :: li, lj, mi, mj, n, p, i, j, ll, mm

  !.. Paramenters for generalized eigensolver
  integer                             :: ITYPE, LDA, LDB, INFO 
  integer                 , parameter :: LWORK = 168777
  integer                 , parameter :: LIWORK = 1448
  character*1                         :: JOBZ = 'V', UPLO = 'U'
  real(kind(1.d0))    , dimension(LM) :: W
  real(kind(1.d0)) , dimension(LWORK) :: WORK
  integer         , dimension(LIWORK) :: IWORK

  !.. External functions/variables
  real(kind(1.d0)) :: bget, bder
  real(kind(1.d0)), allocatable, dimension(:) :: V, Vder

  allocate(V(points),Vder(points))

  ITYPE = 1
  LDA = LM
  LDB = LM

  call gauleg(k, xabsc, weig)

  H = 0.0d0
  Hder = 0.0d0
  S = 0.0d0


  call CPU_TIME( t1 )
  do lj = 1, L
     do li = 1, L
        do mj = 1, M
           do mi = 1, M
              i = (li-1)*M+mi
              j = (lj-1)*M+mj
              sum3 = 0.0d0
              sumder3 = 0.0d0
              sumbsp3 = 0.0d0
              upperl = min(li+1,lj+1) + k-1
              lowerl = max(li,lj)
              upperm = min(mi+1,mj+1) + k-1
              lowerm = max(mi,mj)
              do ll = lowerl, upperl
                 sum2 = 0.0d0
                 sumder2 = 0.0d0
                 sumbsp2 = 0.0d0
                 do mm = lowerm, upperm
                    sum1 = 0.0d0
                    sumder1 = 0.0d0
                    sumbsp1 = 0.0d0
                    if((tl(ll+1).gt.tl(ll)) .and. (tm(mm+1).gt.tm(mm)))then
                       do n = 1, k
                          term = 0.0d0
                          !term1 = 0.0
                          ! term2 = 0.0
                          ! term3 = 0.0
                          ! term4 = 0.0
                          bsp = 0.0d0
                          do p = 1, k
                             coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
                             coordm(mm,p) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(p)
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
                             
                             call twobody_potential(d,r0,mass,rho,theta,phi,V,Vder,points)
                             
                             bsp = bsp + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*sin(2.d0*theta)
                             
                             term(1,:) = term(1,:) + weig(n)*weig(p)*4.d0*dB_li*B_mi*dB_lj*B_mj*sin(2.d0*theta)/(2.d0*my*rho**2.d0)
                             term(2,:) = term(2,:) + weig(n)*weig(p)*8.d0*B_li*dB_mi*B_lj*dB_mj*cos(theta)/(2.d0*sin(theta)*my*rho**2.d0)
                             term(3,:) = term(3,:) + weig(n)*weig(p)*15.d0*B_li*B_mi*B_lj*B_mj*sin(2.d0*theta)/(8.d0*my*rho**2.d0)
                             term(4,:) = term(4,:) + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*V*sin(2.d0*theta)


                             term(5,:) = term(5,:) + weig(n)*weig(p)*8.d0*dB_li*B_mi*dB_lj*B_mj*sin(2.d0*theta)/(my*rho**3.d0)
                             term(6,:) = term(6,:) + weig(n)*weig(p)*8.d0*B_li*dB_mi*B_lj*dB_mj*cos(theta)/(sin(theta)*my*rho**3.d0)
                             term(7,:) = term(7,:) + weig(n)*weig(p)*15.d0*B_li*B_mi*B_lj*B_mj*sin(2.d0*theta)/(4.d0*my*rho**3.d0)
                             term(8,:) = term(8,:) - weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*Vder*sin(2.d0*theta)

                          end do
                          sum1 = sum1 + 0.5d0*(tm(mm+1)-tm(mm))*(term(1,:)+term(2,:)+term(3,:)+term(4,:))
                          sumder1 = sumder1 + 0.5d0*(tm(mm+1)-tm(mm))*(term(5,:)+term(6,:)+term(7,:)+term(8,:))
                          sumbsp1 = sumbsp1 + 0.5d0*(tm(mm+1)-tm(mm))*bsp
                       end do
                       sum2 = sum2 + 0.5d0*(tl(ll+1)-tl(ll))*sum1
                       sumder2 = sumder2 + 0.5d0*(tl(ll+1)-tl(ll))*sumder1
                       sumbsp2 = sumbsp2 + 0.5d0*(tl(ll+1)-tl(ll))*sumbsp1
                    end if
                 end do
                 sum3 = sum3 + sum2
                 sumder3 = sumder3 + sumder2
                 sumbsp3 = sumbsp3 + sumbsp2
              end do
              H(i,j,:) = sum3
              S(i,j,:) = sumbsp3
              Hder(i,j,:) = sumder3
              Integ(i,j,:) = sumbsp3
           end do
        end do
     end do
  end do
  call CPU_TIME( t2 )
  print*,'making array', t2-t1

  
  !.. Calculating Effective potentials and eigenvector coefficients
  call CPU_TIME( t1 )
  do i = 1, points
     Hrez = H(:,:,i)
     Srez = S(:,:,i)
     call dsygvd( ITYPE, JOBZ, UPLO, LM, Hrez, LDA, Srez, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )
     energy(:,i) = W
     H(:,:,i) = Hrez
  end do
  call CPU_TIME( t2 )
  print*,'DSYGVD', t2-t1

  print*, 'INFO:', INFO
  print*, 'Optimal LWORK:', WORK(1)
  print*, 'Optimal LIWORK:', IWORK(1)

  !.. Setting up coupling matrices
  call CPU_TIME( t1 )
  
  Pmat = 0.0d0
  Imat = 0.0d0

  do n = 1, LM
     do p = 1, LM
        sumder1 = 0.0d0
        sumint1 = 0.0d0
        do j = 1, LM
           do i = 1, LM
              sumder1 = sumder1 + H(i,n,:)*H(j,p,:)*Hder(i,j,:)
              sumint1 = sumint1 + H(i,n,:)*H(j,p,:)*Integ(i,j,:)
           end do
        end do
        Pmat(n,p,:) = sumder1/(energy(n,:)-energy(p,:))
        Imat(n,p,:) = sumint1
     end do
  end do

  
  do i = 1, LM
     Pmat(i,i,:) = 0.d0
  end do
  

  do p = 1, LM
     do n = 1, LM
        sum1 = 0.0d0
        do i = 1, LM
           sum1 = sum1 + Pmat(n,i,:)*Pmat(i,p,:)
        end do
        P2mat(n,p,:) = sum1
     end do
  end do

  call CPU_TIME( t2 )
  print*,'creating coupling matrices P and P²', t2-t1

  
  deallocate(V,Vder)
  return

end subroutine efimovham



