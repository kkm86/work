subroutine efimovham_split(npl,npm,k,L,M,LM,tl,tm,rho,my,energy,H,Hder,S,Integ,points)
          

  use constants
  
  implicit none

  !.. Input
  integer,             intent(in) :: npl,npm,k,L,M,LM,points
  real(kind(1.d0)),    intent(in) :: tl(npl),tm(npm),rho(points),my
  real(kind(1.d0)), intent(inout) :: energy(LM,points)
  real(kind(1.d0)), intent(inout) :: H(LM,LM,points),Hder(LM,LM,points), S(LM,LM,points),Integ(LM,LM,points)
  
  !.. Local
  real(kind(1.d0))   :: coordl(L+k,k),coordm(M+k,k),xabsc(k),weig(k),theta,phi
  real(kind(1.d0))   :: term(8,points)
  real(kind(1.d0))   :: lowerl,upperl,lowerm,upperm,t1,t2
  real(kind(1.d0))   :: sumter(3,points), sumbsp(3,points), bsp(points), sumder(3,points), sumint(3,points)
  real(kind(1.d0))   :: B_li, B_lj, B_mi, B_mj, dB_lj, dB_mj, dB_li, dB_mi
  real(kind(1.d0))   :: nsize
  real(kind(1.d0))   :: Hrez(LM,LM),Srez(LM,LM)
  integer            :: li, lj, mi, mj, n, p, i, j, ll, mm, ii, jj, map
  

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
  real(kind(1.d0)), allocatable, dimension(:) :: V, Vder, coordl_1, coordm_1, volume_element, arctan
  real(kind(1.d0)), allocatable, dimension(:,:) :: BL,dBL,BM,dBM
  real(kind(1.d0)), allocatable, dimension(:,:,:) :: V_1, Vder_1
  

  
  allocate(V(points),Vder(points))

  ITYPE = 1
  LDA = LM
  LDB = LM

  call gauleg(k, xabsc, weig)

  H = 0.0d0
  Hder = 0.0d0
  S = 0.0d0

  map = (L+k)*k
  allocate(coordl_1(map),coordm_1(map),volume_element(map),arctan(map),BL(map,L),dBL(map,L),BM(map,M),dBM(map,M),V_1(map,map,points),Vder_1(map,map,points))

  call CPU_TIME( t1 )
  do n = 1, k
     do ll = 1, L+k
        i = (ll-1)*k+n
        coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
        coordl_1(i) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
        volume_element(i) = sin(2.d0*coordl_1(i))
        arctan(i) = cos(coordl_1(i))/sin(coordl_1(i))
     end do
  end do


  do n = 1, k
     do mm = 1, M+k
        i = (mm-1)*k+n
        coordm(mm,n) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(n)
        coordm_1(i) = coordm(mm,n)
     end do
  end do

  do lj = 1, L
     do ii = 1, map
        if(lj == 1)then
           BL(ii,lj) = bget(coordl_1(ii),tl,k,npl,1)+bget(coordl_1(ii),tl,k,npl,2)
           dBL(ii,lj) = bder(coordl_1(ii),tl,k,npl,1)+bder(coordl_1(ii),tl,k,npl,2)
        else if(lj == L)then
           BL(ii,lj) = bget(coordl_1(ii),tl,k,npl,L+1)+bget(coordl_1(ii),tl,k,npl,L+2)
           dBL(ii,lj) = bder(coordl_1(ii),tl,k,npl,L+1)+bder(coordl_1(ii),tl,k,npl,L+2)
        else
           BL(ii,lj) = bget(coordl_1(ii),tl,k,npl,lj+1)
           dBL(ii,lj) = bder(coordl_1(ii),tl,k,npl,lj+1)
        end if
     end do
  end do

  do mj = 1, M
     do ii = 1, map
        if(mj == 1 .and. mj /= M)then
           BM(ii,mj) = bget(coordm_1(ii),tm,k,npm,1)+bget(coordm_1(ii),tm,k,npm,2)
           dBM(ii,mj) = bder(coordm_1(ii),tm,k,npm,1)+bder(coordm_1(ii),tm,k,npm,2)
        else if(mj == M)then
           BM(ii,mj) = bget(coordm_1(ii),tm,k,npm,M+1)+bget(coordm_1(ii),tm,k,npm,M+2)
           dBM(ii,mj) = bder(coordm_1(ii),tm,k,npm,M+1)+bder(coordm_1(ii),tm,k,npm,M+2)
        else
           BM(ii,mj) = bget(coordm_1(ii),tm,k,npm,mj+1)
           dBM(ii,mj) = bder(coordm_1(ii),tm,k,npm,mj+1)
        end if
     end do
  end do

  do jj = 1, map
     do ii = 1, map
        call twobody_potential(rho,coordl_1(jj),coordm_1(ii),V,Vder,points)
        V_1(ii,jj,:) = V
        Vder_1(ii,jj,:) = Vder
     end do
  end do

 
  do lj = 1, L
     do li = 1, L
        do mj = 1, M
           do mi = 1, M
              i = (li-1)*M+mi
              j = (lj-1)*M+mj
              sumter(3,:) = 0.0d0
              sumder(3,:) = 0.0d0
              sumbsp(3,:) = 0.0d0
              upperl = min(li+1,lj+1) + k-1
              lowerl = max(li,lj)
              upperm = min(mi+1,mj+1) + k-1
              lowerm = max(mi,mj)
              do ll = lowerl, upperl
                 sumter(2,:) = 0.0d0
                 sumder(2,:) = 0.0d0
                 sumbsp(2,:) = 0.0d0
                 do mm = lowerm, upperm
                    sumter(1,:) = 0.0d0
                    sumder(1,:) = 0.0d0
                    sumbsp(1,:) = 0.0d0
                    if((tl(ll+1).gt.tl(ll)) .and. (tm(mm+1).gt.tm(mm)))then
                       do n = 1, k
                          jj = (ll-1)*k+n
                          term = 0.0d0
                          bsp = 0.0d0
                          do p = 1, k
                             ii = (mm-1)*k+p
                       
                             bsp = bsp + weig(p)*BL(jj,li)*BM(ii,mi)*BL(jj,lj)*BM(ii,mj)
                             
                             term(1,:) = term(1,:) + weig(p)*2.d0*dBL(jj,li)*BM(ii,mi)*dBL(jj,lj)*BM(ii,mj)
                             term(2,:) = term(2,:) + weig(p)*4.d0*BL(jj,li)*dBM(ii,mi)*BL(jj,lj)*dBM(ii,mj)*arctan(jj)/volume_element(jj)
                             term(3,:) = term(3,:) + weig(p)*15.d0*BL(jj,li)*BM(ii,mi)*BL(jj,lj)*BM(ii,mj)/8.d0
                             term(4,:) = term(4,:) + weig(p)*BL(jj,li)*BM(ii,mi)*BL(jj,lj)*BM(ii,mj)*V_1(ii,jj,:)


                             term(5,:) = term(5,:) + weig(p)*8.d0*dBL(jj,li)*BM(ii,mi)*dBL(jj,lj)*BM(ii,mj)
                             term(6,:) = term(6,:) + weig(p)*8.d0*BL(jj,li)*dBM(ii,mi)*BL(jj,lj)*dBM(ii,mj)*arctan(jj)/volume_element(jj)
                             term(7,:) = term(7,:) + weig(p)*15.d0*BL(jj,li)*BM(ii,mi)*BL(jj,lj)*BM(ii,mj)/4.d0
                             term(8,:) = term(8,:) - weig(p)*BL(jj,li)*BM(ii,mi)*BL(jj,lj)*BM(ii,mj)*Vder_1(ii,jj,:)

                          end do
                          sumter(1,:) = sumter(1,:) + 0.5d0*volume_element(jj)*weig(n)*(tm(mm+1)-tm(mm))*(term(1,:)+term(2,:)+term(3,:)+term(4,:)*my*rho**2.d0)/(my*rho**2.d0)
                          sumder(1,:) = sumder(1,:) + 0.5d0*volume_element(jj)*weig(n)*(tm(mm+1)-tm(mm))*(term(5,:)+term(6,:)+term(7,:)+term(8,:)*my*rho**3.d0)/(my*rho**3.d0)
                          sumbsp(1,:) = sumbsp(1,:) + 0.5d0*volume_element(jj)*weig(n)*(tm(mm+1)-tm(mm))*bsp
                       end do
                       sumter(2,:) = sumter(2,:) + 0.5d0*(tl(ll+1)-tl(ll))*sumter(1,:)
                       sumder(2,:) = sumder(2,:) + 0.5d0*(tl(ll+1)-tl(ll))*sumder(1,:)
                       sumbsp(2,:) = sumbsp(2,:) + 0.5d0*(tl(ll+1)-tl(ll))*sumbsp(1,:)
                    end if
                 end do
                 sumter(3,:) = sumter(3,:) + sumter(2,:)
                 sumder(3,:) = sumder(3,:) + sumder(2,:)
                 sumbsp(3,:) = sumbsp(3,:) + sumbsp(2,:)
              end do
              H(i,j,:) = sumter(3,:)
              S(i,j,:) = sumbsp(3,:)
              Hder(i,j,:) = sumder(3,:)
              Integ(i,j,:) = sumbsp(3,:)
           end do
        end do
     end do
  end do
  call CPU_TIME( t2 )
  print*,'cpu_time: ', (t2-t1)
 

  
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
  
  
  
  deallocate(V,Vder)
  return

end subroutine efimovham_split




