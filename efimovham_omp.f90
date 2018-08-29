subroutine efimovham_omp(npl,npm,k,L,M,LM,tl,tm,rho,my,energy,H,Hder,S,Integ,points,Pmat,P2mat,Imat)
          
  use constants
  
  implicit none

  !.. Input
  integer            , intent(in) :: npl,npm,k,L,M,LM,points
  real(kind(1.d0))   , intent(in) :: tl(npl),tm(npm),rho(points),my
  real(kind(1.d0)), intent(inout) :: energy(LM,points)
  real(kind(1.d0)), intent(inout) :: H(LM,LM,points),Hder(LM,LM,points), S(LM,LM,points),Integ(LM,LM,points),Pmat(LM,LM,points),P2mat(LM,LM,points),Imat(LM,LM,points)
  
  !.. Local
  real(kind(1.d0))   :: coordl(L+2,k),coordm(M+2,k),xabsc(k),weig(k),theta,phi(k)
  real(kind(1.d0))   :: term(8,points)
  real(kind(1.d0))   :: lowerl,upperl,lowerm,upperm,t1,t2
  real(kind(1.d0))   :: sumter(3,points), sumbsp(3,points), bsp(points), sumder(3,points), sumint(3,points)
  real(kind(1.d0))   :: B_li, B_lj, B_mi(k), B_mj(k), dB_lj, dB_mj(k), dB_li, dB_mi(k), potent(k,points), potentder(k,points), multi(k), multider(k), firstsum(points), volume_element, arctan, multivec(points), multidervec(points)
  real(kind(1.d0))   :: nsize
  real(kind(1.d0))   :: Hrez(LM,LM),Srez(LM,LM)
  integer            :: li, lj, mi, mj, n, p, i, j, ll, mm

  !.. Parameters for generalized eigensolver
  integer                             :: ITYPE, LDA, LDB, INFO 
  integer                 , parameter :: LWORK = 168777
  integer                 , parameter :: LIWORK = 1448
  character*1                         :: JOBZ = 'V', UPLO = 'U'
  real(kind(1.d0))    , dimension(LM) :: W
  real(kind(1.d0)) , dimension(LWORK) :: WORK
  integer         , dimension(LIWORK) :: IWORK

  !.. External functions/variables
  real(kind(1.d0)) :: bget, bder
  !real(kind(1.d0)), external :: potent_omp
  !real(kind(1.d0)), allocatable, dimension(:) :: V, Vder

  !allocate(V(points),Vder(points))

  ITYPE = 1
  LDA = LM
  LDB = LM

  call gauleg(k, xabsc, weig)

  H = 0.0d0
  Hder = 0.0d0
  S = 0.0d0

  print*, 'ok1', points

  call CPU_TIME( t1 )
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
                          term = 0.0d0
                          bsp = 0.0d0
                        
                          coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
                          coordm(mm,:) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(:)
                          theta = coordl(ll,n)
                          phi(:) = coordm(mm,:)
                          volume_element = sin(2.d0*theta)
                          arctan = cos(theta)/sin(theta)
                           
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
                             B_mj(:) = bget(coordm(mm,:),tm,k,npm,1)+bget(coordm(mm,:),tm,k,npm,2)
                             dB_mj(:) = bder(coordm(mm,:),tm,k,npm,1)+bder(coordm(mm,:),tm,k,npm,2)
                          else if(mj == M)then
                             B_mj(:) = bget(coordm(mm,:),tm,k,npm,M+1)+bget(coordm(mm,:),tm,k,npm,M+2)
                             dB_mj(:) = bder(coordm(mm,:),tm,k,npm,M+1)+bder(coordm(mm,:),tm,k,npm,M+2)
                          else
                             B_mj(:) = bget(coordm(mm,:),tm,k,npm,mj+1)
                             dB_mj(:) = bder(coordm(mm,:),tm,k,npm,mj+1)
                          end if
                             
                          if(mi == 1 .and. mi /= M)then
                             B_mi(:) = bget(coordm(mm,:),tm,k,npm,1)+bget(coordm(mm,:),tm,k,npm,2)
                             dB_mi(:) = bder(coordm(mm,:),tm,k,npm,1)+bder(coordm(mm,:),tm,k,npm,2)
                          else if(mi == M)then
                             B_mi(:) = bget(coordm(mm,:),tm,k,npm,M+1)+bget(coordm(mm,:),tm,k,npm,M+2)
                             dB_mi(:) = bder(coordm(mm,:),tm,k,npm,M+1)+bder(coordm(mm,:),tm,k,npm,M+2)
                          else
                             B_mi(:) = bget(coordm(mm,:),tm,k,npm,mi+1)
                             dB_mi(:) = bder(coordm(mm,:),tm,k,npm,mi+1)
                          end if
                             
                          potent = potent_omp(rho,theta,phi,points,k)
                          potentder = potentder_omp(rho,theta,phi,points,k)

                          multi = weig*B_mi*B_mj
                          multider = weig*dB_mi*dB_mj

                          multivec = matmul(potent,multi)
                          multidervec = matmul(potentder,multi)
                         
                          bsp = weig(n)*B_li*B_lj*multi
                             
                          term(1,:) = 2.d0*dB_li*dB_lj*multi
                          term(2,:) = 4.d0*B_li*B_lj*multider*arctan/volume_element
                          term(3,:) = 15.d0*B_li*B_lj*multi/8.d0
                          term(4,:) = B_li*B_lj*multivec*my*rho*rho


                          term(5,:) = 8.d0*dB_li*dB_lj*multi
                          term(6,:) = 8.d0*B_li*B_lj*multider*arctan/volume_element
                          term(7,:) = 15.d0*B_li*B_lj*multi/4.d0
                          term(8,:) = B_li*B_lj*multidervec*my*rho*rho*rho

                        
                          sumter(1,:) = sumter(1,:) + 0.5d0*volume_element*weig(n)*(tm(mm+1)-tm(mm))*(term(1,:)+term(2,:)+term(3,:)+term(4,:))/(my*rho*rho)
                          sumder(1,:) = sumder(1,:) + 0.5d0*volume_element*weig(n)*(tm(mm+1)-tm(mm))*(term(5,:)+term(6,:)+term(7,:)+term(8,:))/(my*rho*rho*rho)
                          sumbsp(1,:) = sumbsp(1,:) + 0.5d0*volume_element*weig(n)*(tm(mm+1)-tm(mm))*bsp
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
        sumder(1,:) = 0.0d0
        sumint(1,:) = 0.0d0
        do j = 1, LM
           do i = 1, LM
              sumder(1,:) = sumder(1,:) + H(i,n,:)*H(j,p,:)*Hder(i,j,:)
              sumint(1,:) = sumint(1,:) + H(i,n,:)*H(j,p,:)*Integ(i,j,:)
           end do
        end do
        Pmat(n,p,:) = sumder(1,:)/(energy(n,:)-energy(p,:))
        Imat(n,p,:) = sumint(1,:)
     end do
  end do

  
  do i = 1, LM
     Pmat(i,i,:) = 0.d0
  end do
  

  do p = 1, LM
     do n = 1, LM
        sumter(1,:) = 0.0d0
        do i = 1, LM
           sumter(1,:) = sumter(1,:) + Pmat(n,i,:)*Pmat(i,p,:)
        end do
        P2mat(n,p,:) = sumter(1,:)
     end do
  end do

  call CPU_TIME( t2 )
  print*,'creating coupling matrices P and P²', t2-t1

  
!  deallocate(V,Vder)
  return

end subroutine efimovham_omp




