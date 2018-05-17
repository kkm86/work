subroutine Hamiltonian2D(np,k,L,M,LM,tl,tm,energy)

  use constants

   implicit none

   !.. Input
   integer            , intent(in) :: np,k,L,M,LM       !..Number of knot-points etc.
   real(kind(1.d0))   , intent(in) :: tl(np), tm(np)
   real(kind(1.d0)), intent(inout) :: energy

   !.. Local
   integer, parameter :: gl = 5 !..Number of Gauss-Legendre points
   real(kind(1.d0))   :: coordl(L,gl),coordm(M,gl),Ham(LM,LM),Bs(LM,LM) !..Abscissa-grid
   real(kind(1.d0))   :: sum1,sum2,sum3,sum4,term1,term2,term3,term4,term5,lowerl,upperl,lowerm,upperm,sumbsp1,sumbsp2,sumbsp3,sumbsp4,bsp,xabsc(k),weig(k),B_li,B_lj,B_mi,B_mj,dB_lj,dB_mj,ddB_lj,ddB_mj
   real(kind(1.d0))   :: nsize
   integer            :: li, lj, mi, mj, n, p, i, j, ll, mm

   
   !.. Paramenters for generalized eigensolver
  ! integer                             :: ITYPE, LDA, LDB, INFO 
  ! integer                 , parameter :: LWORK = 1000000
  ! integer                 , parameter :: LIWORK = 5000
  ! character*1                         :: JOBZ = 'V', UPLO = 'U'
  ! real(kind(1.d0))    , dimension(LM) :: W,W2
  ! real(kind(1.d0)) , dimension(LWORK) :: WORK
  ! real(kind(1.d0)), dimension(LIWORK) :: IWORK

   !.. Paramenters for generalized eigensolver
   integer                             :: LDA, LDB, INFO 
   character*1                         :: JOBVL = 'N', JOBVR = 'N'
   integer                 , parameter :: LDVL = 100
   integer                 , parameter :: LDVR = 100
   integer                 , parameter :: LWORK = 2808
   real(kind(1.d0))    , dimension(LM) :: ALPHAR, ALPHAI, BETA,eigvr,eigvi
   real(kind(1.d0)) , dimension(LDVL,LM) :: VL
   real(kind(1.d0)) , dimension(LDVR,LM) :: VR
   real(kind(1.d0)) , dimension(max(1,LWORK)) :: WORK


   !.. External functions/variables
   real(kind(1.d0)) :: bget, bder, bder2

    !.. Initializing parameters for the generalized eigenvalue solver
   ! ITYPE = 1
    LDA = LM
    LDB = LM
   
  
    call gauleg(gl, xabsc, weig)

    Ham = 0.d0
    Bs = 0.d0 
    
    do lj = 1, L
       do li = 1, L
          do mj = 1, M
             do mi = 1, M
                sum1 = 0.0
                sum2 = 0.0
                sum3 = 0.0
                sum4 = 0.0
                sumbsp1 = 0.0
                sumbsp2 = 0.0
                sumbsp3 = 0.0
                sumbsp4 = 0.0
                bsp = 0.0
                term1 = 0.0
                term2 = 0.0
                term3 = 0.0
                term4 = 0.0
                term5 = 0.0
                upperl = min(li,lj) + k-1
                lowerl = max(li,lj)
                upperm = min(mi,mj) + k-1
                lowerm = max(mi,mj)
                do ll = lowerl, upperl
                   do mm = lowerm, upperm
                      do n = 1, gl
                         term1 = 0.0
                         term2 = 0.0
                         term3 = 0.0
                         term4 = 0.0
                         term5 = 0.0
                         bsp = 0.0
                         do p = 1, gl
                            if((tl(ll+1).gt.tl(ll)) .and. (tm(mm+1).gt.tm(mm))) then
                               coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
                               coordm(mm,p) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(p)

                               B_li = bget(coordl(ll,n),tl,k,np,li)
                               B_lj = bget(coordl(ll,n),tl,k,np,lj)
                               dB_lj = bder(coordl(ll,n),tl,k,np,lj)
                               ddB_lj = bder2(coordl(ll,n),tl,k,np,lj)
                     
                               if(mj < 3) then
                                  B_mj = bget(coordm(mm,p),tm,k,np,1)+bget(coordm(mm,p),tm,k,np,2)
                                  ddB_mj = bder2(coordm(mm,p),tm,k,np,1)+bder2(coordm(mm,p),tm,k,np,2)
                               else
                                  B_mj = bget(coordm(mm,p),tm,k,np,mj)
                                  ddB_mj = bder2(coordm(mm,p),tm,k,np,mj)
                               end if

                               if(mi < 3) then
                                  B_mi = bget(coordm(mm,p),tm,k,np,1)+bget(coordm(mm,p),tm,k,np,2)
                               else
                                  B_mi = bget(coordm(mm,p),tm,k,np,mi)
                               end if
                               
                               term1 = term1 - weig(n)*weig(p)*B_li*B_mi*dB_lj*B_mj/(2.d0*coordl(ll,n))
                               term2 = term2 - weig(n)*weig(p)*B_li*B_mi*ddB_lj*B_mj/2.d0
                               term3 = term3 - weig(n)*weig(p)*B_li*B_mi*B_lj*ddB_mj/2.d0
                               term4 = term4 + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*(coordl(ll,n)**2.d0)/2.d0
                               term5 = term5 + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj*2.d0*(coordm(ll,n)**2.d0)

                               bsp = bsp + weig(n)*weig(p)*B_li*B_mi*B_lj*B_mj
                            end if
                         end do
                         sum1 = sum1 + 0.5d0*(tm(mm+1)-tm(mm))*(term1+term2+term3+term4+term5)
                         sumbsp1 = sumbsp1 + 0.5d0*(tm(mm+1)-tm(mm))*bsp
                      end do
                      sum2 = sum2 + 0.5d0*(tl(ll+1)-tl(ll))*sum1
                      sumbsp2 = sumbsp2 + 0.5d0*(tl(ll+1)-tl(ll))*sumbsp1
                   end do
                   sum3 = sum3 + sum2
                   sumbsp3 = sumbsp3 + sumbsp2
                end do
                sum4 = sum4 + sum3
                sumbsp4 = sumbsp4 + sumbsp3
               
                i = (li-1)*M+mi
                j = (lj-1)*M+mj
                Ham(i,j) = sum4
                Bs(i,j) = sumbsp4
             end do
          end do
       end do
    end do

             
  !  do dsygvd( ITYPE, JOBZ, UPLO, LM, Ham, LDA, Bs, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )
    
    call DGGEV( JOBVL, JOBVR, LM, Ham, LDA, Bs, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )


    write(*,*) ALPHAR(1)/BETA(1), ALPHAI(1)/BETA(1)

!    stop

    open(10,file='eig.dat',status='replace')
    do i = 1, LM
       eigvr(i)=ALPHAR(i)/BETA(i)
       eigvi(i)=ALPHAI(i)/BETA(i)
    enddo

    do i = 1, LM
       write(10,*) i
       write(10,*) eigvr(i),eigvi(i)
       
       write(10,*) BETA(i)
    end do
    close(10)

    

  
    !energy = W(1)
   
    ! wfn=Ham(:,2)
   
    print *, 'Good job!',INFO,WORK(1)
    
    return
  end subroutine Hamiltonian2D
  
  
