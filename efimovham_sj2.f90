!subroutine efimovham(npl,npm,k,L,M,LM,tl,tm,rho,my,energy,H,Hder,points,Pmat,P2mat)
subroutine efimovham(npl,npm,k,l,m,lm,points,tl,tm,rho,energy,S,TK,H,Hder)
          

  use constants
  
  implicit none

  !.. Input
  integer            , intent(in) :: npl,npm,points,lm,l,k,m
  real(kind(1.d0))   , intent(in) :: tl(npl),tm(npm),rho(points)
  real(kind(1.d0)), intent(inout) :: energy(LM,points)
  real(kind(1.d0)), intent(inout) :: S(LM,LM),TK(LM,LM)
  real(kind(1.d0)), intent(inout) :: H(LM,LM,points) ,Hder(LM,LM,points)
  
  !.. Local
  real(kind(1.d0))   :: coordl(L+2,k),coordm(M+2,k),xabsc(k),weig(k),theta,phi
  real(kind(1.d0))   :: term(8,points)
  real(kind(1.d0))   :: lowerl,upperl,lowerm,upperm,t1,t2
  real(kind(1.d0))   :: sumter(3,points),  sumder(3,points), sumint(3,points), volume_element, arctan
  real(kind(1.d0))   :: sumS1,sumS2,sumT11,sumT12,sumT21,sumT22
  real(kind(1.d0))   :: B_li, B_lj, B_mi, B_mj, dB_lj, dB_mj, dB_li, dB_mi
  real(kind(1.d0))   :: nsize,my
  integer            :: li, lj, mi, mj, n, p,i,j, i1, j1,i2,j2, ll, mm,mimax,limax



  !.. External functions/variables
  real(kind(1.d0)) :: bget, bder
  real(kind(1.d0)), allocatable, dimension(:) :: V, Vder

  integer :: c1,c2,cr,cm
  real(kind(1.d0)) :: rate

  !some new variables
  real(kind(1.d0))    , dimension(L,L) :: s1
  real(kind(1.d0))    , dimension(M,M) :: s2
  real(kind(1.d0))    , dimension(L,L) :: tk11,tk12
  real(kind(1.d0))    , dimension(M,M) :: tk2
  ! TK is the hyperangular momentum oprator
  !.. Initialize the system clock
  call system_clock(count_rate=cr)
  call system_clock(count_max=cm)
  rate = real(cr)
  write(*,*) 'system_clock rate', rate

  allocate(V(points),Vder(points))

  call gauleg(k, xabsc, weig)

  H = 0.0d0
  Hder = 0.0d0

  s1 = 0.d0
  s2 = 0.d0
  tk11=0.d0
  tk12=0.d0
  tk2=0.d0


  call CPU_TIME( t1 )
  call system_clock( c1 )

  my = mass(3)/sqrt(3.d0)


!! overlap and hyperangular mom. matrix in theta:
  do lj = 1, L
     limax=min(L,lj+k)
     do li = lj, limax
        upperl = min(li+1,lj+1) + k-1
        lowerl = max(li,lj)
        sumS2 = 0.0d0
        sumT12=0.d0
        sumT22=0.d0
        do ll = lowerl, upperl
           sumS1 = 0.0d0
           sumT11 = 0.0d0
           sumT21 = 0.0d0
           if(tl(ll+1).gt.tl(ll))then
              do n = 1, k
                 coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
                 theta = coordl(ll,n)
                 volume_element = sin(2.d0*theta)
                 arctan=cos(theta)/sin(theta)

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

                 sumS1 = sumS1 + volume_element*weig(n)*B_li*B_lj
                 sumT11 = sumT11 + 4.d0*volume_element*weig(n)*dB_li*dB_lj
                 sumT21 = sumT21 + 8.d0*arctan*weig(n)*B_li*B_lj
              end do ! sum over n
              sumS1 =  0.5d0*(tl(ll+1)-tl(ll))*sumS1
              sumT11 =  0.5d0*(tl(ll+1)-tl(ll))*sumT11
              sumT21 =  0.5d0*(tl(ll+1)-tl(ll))*sumT21
           end if
           sumS2 = sumS2 + sumS1
           sumT12 = sumT12 + sumT11
           sumT22 = sumT22 + sumT21
        end do ! sum over ll        
        s1(li,lj) = sumS2
        s1(lj,li)=s1(li,lj)
        tk11(li,lj) = sumT12
        tk11(lj,li)=tk11(li,lj)
        tk12(li,lj) = sumT22
        tk12(lj,li)=tk12(li,lj)
     end do
  end do

!! overlap matrix in phi:
  do mj = 1, M
     mimax=min(M,mj+k)
     do mi = mj, M
        upperm = min(mi+1,mj+1) + k-1
        lowerm = max(mi,mj)
        sumS2 = 0.0d0
        sumT12 = 0.d0
        do mm = lowerm, upperm
           sumS1 = 0.0d0
           sumT11 = 0.d0
           if(tm(mm+1).gt.tm(mm))then
              do p = 1, k
                 coordm(mm,p) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(p)
                 phi = coordm(mm,p)
                 
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


                 sumS1 = sumS1 + weig(p)*B_mi*B_mj
                 sumT11 = sumT11 + weig(p)*dB_mi*dB_mj

              end do ! sum over p
              sumS1 =  0.5d0*(tm(mm+1)-tm(mm))*sumS1
              sumT11 =  0.5d0*(tm(mm+1)-tm(mm))*sumT11

           end if
           sumS2 = sumS2 + sumS1
           sumT12 = sumT12 + sumT11

        end do ! sum over mm        
        s2(mi,mj) = sumS2
        s2(mj,mi)=s2(mi,mj)
        tk2(mi,mj) = sumT12
        tk2(mj,mi)=tk2(mi,mj)
     end do
  end do



  do lj = 1, L
     limax=min(L,lj+k)
     do li = lj, limax
        do mj = 1, M
           mimax=min(M,mj+k)
           do mi = mj, mimax
              sumter(3,:) = 0.0d0
              sumder(3,:) = 0.0d0
              upperl = min(li+1,lj+1) + k-1
              lowerl = max(li,lj)
              upperm = min(mi+1,mj+1) + k-1
              lowerm = max(mi,mj)
              do ll = lowerl, upperl
                 sumter(2,:) = 0.0d0
                 sumder(2,:) = 0.0d0
                 do mm = lowerm, upperm
                    sumter(1,:) = 0.0d0
                    sumder(1,:) = 0.0d0
                    if((tl(ll+1).gt.tl(ll)) .and. (tm(mm+1).gt.tm(mm)))then
                       do n = 1, k
                          term = 0.0d0
                          coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
                          theta = coordl(ll,n)
                          volume_element = sin(2.d0*theta)
                          arctan = 1./tan(theta)

                          if(lj == 1)then
                             B_lj = bget(coordl(ll,n),tl,k,npl,1)+bget(coordl(ll,n),tl,k,npl,2)
                          else if(lj == L)then
                             B_lj = bget(coordl(ll,n),tl,k,npl,L+1)+bget(coordl(ll,n),tl,k,npl,L+2)
                          else
                             B_lj = bget(coordl(ll,n),tl,k,npl,lj+1)
                          end if

                          if(li == 1)then
                             B_li = bget(coordl(ll,n),tl,k,npl,1)+bget(coordl(ll,n),tl,k,npl,2)
                          else if(li == L)then
                             B_li = bget(coordl(ll,n),tl,k,npl,L+1)+bget(coordl(ll,n),tl,k,npl,L+2)
                          else
                             B_li = bget(coordl(ll,n),tl,k,npl,li+1)
                          end if

                          do p = 1, k
                             coordm(mm,p) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(p)
                             phi = coordm(mm,p)
                        
                             if(mj == 1 .and. mj /= M)then
                                B_mj = bget(coordm(mm,p),tm,k,npm,1)+bget(coordm(mm,p),tm,k,npm,2)
                             else if(mj == M)then
                                B_mj = bget(coordm(mm,p),tm,k,npm,M+1)+bget(coordm(mm,p),tm,k,npm,M+2)
                             else
                                B_mj = bget(coordm(mm,p),tm,k,npm,mj+1)
                             end if
                             
                             if(mi == 1 .and. mi /= M)then
                                B_mi = bget(coordm(mm,p),tm,k,npm,1)+bget(coordm(mm,p),tm,k,npm,2)
                             else if(mi == M)then
                                B_mi = bget(coordm(mm,p),tm,k,npm,M+1)+bget(coordm(mm,p),tm,k,npm,M+2)
                             else
                                B_mi = bget(coordm(mm,p),tm,k,npm,mi+1)
                             end if
                             
                             call twobody_potential(rho,theta,phi,V,points)
                                                          
                             term(4,:) = term(4,:) + weig(p)*B_li*B_mi*B_lj*B_mj*V



                          end do
                          sumter(1,:) = sumter(1,:) + 0.5d0*volume_element*weig(n)*(tm(mm+1)-tm(mm))*term(4,:)

                       end do
                       sumter(2,:) = sumter(2,:) + 0.5d0*(tl(ll+1)-tl(ll))*sumter(1,:)
                    end if
                 end do
                 sumter(3,:) = sumter(3,:) + sumter(2,:)
              end do
              ! each (lj,li,mj,mi) corresponds to 4 (i,j) 
              i1 = (li-1)*M+mi
              j1 = (lj-1)*M+mj 
              S(i1,j1) = s1(li,lj)*s2(mi,mj)
              TK(i1,j1)=tk11(li,lj)*s2(mi,mj)+tk12(li,lj)*tk2(mi,mj)
              H(i1,j1,:) = (TK(i1,j1)+3.75d0*S(i1,j1))/(2.d0*my*rho**2.d0)+sumter(3,:)
              Hder(i1,j1,:) = sumder(3,:)
            
              S(j1,i1)=S(i1,j1)
              TK(j1,i1)=TK(i1,j1)
              H(j1,i1,:)=H(i1,j1,:)
              Hder(j1,i1,:)=Hder(i1,j1,:)
! This seems wrong for Hder, but fair enough its OK for H and S
              i2 = (li-1)*M+mj
              j2 = (lj-1)*M+mi 
              S(i2,j2) = S(i1,j1)
              TK(i2,j2)=TK(i1,j1)
              H(i2,j2,:) = H(i1,j1,:)
              Hder(i2,j2,:) = Hder(i1,j1,:)

              S(j2,i2)=S(i2,j2)
              TK(j2,i2)=TK(i2,j2)
              H(j2,i2,:)=H(i2,j2,:)
              Hder(j2,i2,:)=Hder(i2,j2,:)

           end do
        end do
     end do
  end do
  call CPU_TIME( t2 )
  call system_clock( c2 )
  print*,'making array', t2-t1
  print*,'system_clock: ', (c2-c1)/rate


  !.. Calculating Effective potentials and eigenvector coefficients
  deallocate(V,Vder)
  return

end subroutine efimovham



