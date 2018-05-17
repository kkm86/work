subroutine anisotrop(np,k,L,M,LM,tl,tm,energy,H,S)

  use constants

   implicit none

   !.. Input
   integer            , intent(in) :: np,k,L,M,LM       !..Number of knot-points etc.
   real(kind(1.d0))   , intent(in) :: tl(np), tm(np)
   real(kind(1.d0)), intent(inout) :: energy(3)
   real(kind(1.d0)), intent(inout) :: H(LM,LM), S(LM,LM)

   !.. Local
   real(kind(1.d0))   :: coordl(L,k),coordm(M+2,k) !..Abscissa-grid
   real(kind(1.d0))   :: Sl(L,L), Sm(M,M), V_rho(L,L), V_z(M,M), Tll(L,L), Tmm(M,M)
   real(kind(1.d0))   :: sum1,sum2,sum3
   real(kind(1.d0))   :: term1,term2,term3
   real(kind(1.d0))   :: lower,upper,xabsc(k),weig(k),B_li,B_lj,B_mi,B_mj,dB_lj,dB_mj,dB_li,dB_mi
   real(kind(1.d0))   :: nsize
   integer            :: li, lj, mi, mj, n, i, j, ll, mm

   !.. Paramenters for generalized eigensolver
   integer                             :: ITYPE, LDA, LDB, INFO 
   integer                 , parameter :: LWORK = 1000000
   integer                 , parameter :: LIWORK = 5000
   character*1                         :: JOBZ = 'V', UPLO = 'U'
   real(kind(1.d0))    , dimension(LM)  :: W
   real(kind(1.d0)) , dimension(LWORK) :: WORK
   real(kind(1.d0)), dimension(LIWORK) :: IWORK

 

   !.. External functions/variables
   real(kind(1.d0)) :: bget, bder

   ITYPE = 1
   LDA = LM
   LDB = LM
  
  
   call gauleg(k, xabsc, weig)

   Sm = 0.0
   V_z = 0.0
   Tmm = 0.0
   Sl = 0.0
   V_rho = 0.0
   Tll = 0.0
   H = 0.0
   S = 0.0
   
   do lj = 1, L
      do li = 1, L
         sum1 = 0.0
         sum2 = 0.0
         sum3 = 0.0
         upper = min(li,lj) + k-1
         lower = max(li,lj)
         do ll = lower, upper
            term1 = 0.0
            term2 = 0.0
            term3 = 0.0
            do n = 1, k
               if(tl(ll+1).gt.tl(ll)) then
                  coordl(ll,n) = 0.5*(tl(ll+1)+tl(ll)) + 0.5*(tl(ll+1)-tl(ll))*xabsc(n)
                  B_li = bget(coordl(ll,n),tl,k,np,li)
                  B_lj = bget(coordl(ll,n),tl,k,np,lj)
                  dB_li = bder(coordl(ll,n),tl,k,np,li)
                  dB_lj = bder(coordl(ll,n),tl,k,np,lj)
                  term1 = term1 + weig(n)*B_li*B_lj*coordl(ll,n)
                  term2 = term2 + weig(n)*B_li*B_lj*(coordl(ll,n)**3.d0)
                  term3 = term3 + weig(n)*dB_li*dB_lj*coordl(ll,n)
               end if
            end do
            sum1 = sum1 + 0.5d0*(tl(ll+1)-tl(ll))*term1
            sum2 = sum2 + 0.5d0*(tl(ll+1)-tl(ll))*term2
            sum3 = sum3 + 0.5d0*(tl(ll+1)-tl(ll))*term3
         end do
         Sl(li,lj) = sum1
         V_rho(li,lj) = sum2
         Tll(li,lj) = sum3
      end do
   end do

   do mj = 1, M
       do mi = 1, M
          sum1 = 0.0
          sum2 = 0.0
          sum3 = 0.0
          upper = min(mi+1,mj+1) + k-1
          lower = max(mi,mj)
          do mm = lower, upper
             term1 = 0.0
             term2 = 0.0
             term3 = 0.0
             do n = 1, k
                if(tm(mm+1).gt.tm(mm)) then
                   coordm(mm,n) = 0.5*(tm(mm+1)+tm(mm)) + 0.5*(tm(mm+1)-tm(mm))*xabsc(n)
                   if(mi == 1) then
                      B_mi = bget(coordm(mm,n),tm,k,np,1)+bget(coordm(mm,n),tm,k,np,2)
                      dB_mi = bder(coordm(mm,n),tm,k,np,1)+bder(coordm(mm,n),tm,k,np,2)
                   else
                      B_mi = bget(coordm(mm,n),tm,k,np,mi+1)
                      dB_mi = bder(coordm(mm,n),tm,k,np,mi+1)
                   end if
                   if(mj == 1) then
                      B_mj = bget(coordm(mm,n),tm,k,np,1)+bget(coordm(mm,n),tm,k,np,2)
                      dB_mj = bder(coordm(mm,n),tm,k,np,1)+bder(coordm(mm,n),tm,k,np,2)
                   else
                      B_mj = bget(coordm(mm,n),tm,k,np,mj+1)
                      dB_mj = bder(coordm(mm,n),tm,k,np,mj+1)
                   end if
                   term1 = term1 + weig(n)*B_mi*B_mj
                   term2 = term2 + weig(n)*B_mi*B_mj*2.d0*(coordm(mm,n)**2.d0)
                   term3 = term3 + weig(n)*dB_mi*dB_mj
                end if
             end do
             sum1 = sum1 + 0.5d0*(tm(mm+1)-tm(mm))*term1
             sum2 = sum2 + 0.5d0*(tm(mm+1)-tm(mm))*term2
             sum3 = sum3 + 0.5d0*(tm(mm+1)-tm(mm))*term3
          end do
          Sm(mi,mj) = sum1
          V_z(mi,mj) = sum2
          Tmm(mi,mj) = sum3
       end do
    end do

   
    do lj = 1, L
       do li = 1, L
          do mj = 1, M
             do mi = 1, M
                i = (li-1)*M+mi
                j = (lj-1)*M+mj
                H(i,j) = 0.5d0*Tll(li,lj)*Sm(mi,mj) + 0.5d0*Tmm(mi,mj)*Sl(li,lj) + V_z(mi,mj)*Sl(li,lj) + 0.5d0*V_rho(li,lj)*Sm(mi,mj)
                S(i,j) = Sl(li,lj)*Sm(mi,mj)
             end do
          end do
       end do
    end do


    call dsygvd( ITYPE, JOBZ, UPLO, LM, H, LDA, S, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )
   
    print *, W(1), W(2), W(3)
    
    print *, 'Good boy!'
    
    return
  end subroutine Anisotrop
  
  
  
