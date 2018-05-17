 subroutine Hamiltonian_born(np,k,NN,t,energy,wfn,S,b,mass)

   implicit none

   !.. Input
    integer            , intent(in) :: np,k,NN          !Number of knot-points etc.
    real(kind(1.d0))   , intent(in) :: t(np)
    real(kind(1.d0)), intent(inout) :: energy
    real(kind(1.d0)), intent(inout) :: wfn(NN)

    !.. Local
    real(kind(1.d0)) :: r(NN,k),Ham(NN,NN),Bs(NN,NN), r_1 !Abscissa-grid
    real(kind(1.d0)) :: sum1, sum2,term1,term2,term3,lower,upper,sumbsp,sumbsp2,bsp,xabsc(k),weig(k)
    integer          :: ii, jj, m, n

   
    !.. Paramenters for generalized eigensolver
    integer                             :: ITYPE, LDA, LDB, INFO 
    integer                 , parameter :: LWORK = 1000000
    integer                 , parameter :: LIWORK = 5000
    character*1                         :: JOBZ = 'V', UPLO = 'U'
    real(kind(1.d0))    , dimension(NN) :: W,W2
    real(kind(1.d0)) , dimension(LWORK) :: WORK
    real(kind(1.d0)), dimension(LIWORK) :: IWORK


    !.. External functions/variables
    real(kind(1.d0)) :: bget, bder, bder2, S, b, mass, V

    !.. Initializing parameters for the generalized eigenvalue solver
    ITYPE = 1
    LDA = NN
    LDB = NN

    call gauleg(k, xabsc, weig)

    Ham = 0.d0
    Bs = 0.d0 
    
    do jj = 1, NN
       do ii = 1, NN
          sum1 = 0.0
          sum2 = 0.0
          sumbsp = 0.0
          sumbsp2 = 0.0
          bsp = 0.0
          term1 = 0.0
          term2 = 0.0
          term3 = 0.0
          upper = min(ii+1,jj+1) + k-1
          lower = max(ii+1,jj+1)
          do m = lower, upper
             term1 = 0.0
             term2 = 0.0
             term3 = 0.0
             bsp = 0.0
             do n = 1, k
                if(t(m+1).gt.t(m)) then
                   r(m,n) = 0.5*(t(m+1)+t(m)) + 0.5*(t(m+1)-t(m))*xabsc(n)
                   r_1 = r(m,n)
                   term1 = term1 + weig(n)*bder(r(m,n),t,k,np,jj+1)*bder(r(m,n),t,k,np,ii+1)
                   term2 = term2 - weig(n)*bget(r(m,n),t,k,np,jj+1)*bder(r(m,n),t,k,np,ii+1)*2/r(m,n)
                   call twobody_born(S,b,mass,r_1,V)
                   term3 = term3 + weig(n)*bget(r(m,n),t,k,np,jj+1)*bget(r(m,n),t,k,np,ii+1)*V
                   bsp = bsp + weig(n)*bget(r(m,n),t,k,np,jj+1)*bget(r(m,n),t,k,np,ii+1)
                end if
             end do
             sum1 = sum1 + 0.5*(t(m+1)-t(m))*(term1+term2+term3)
             sumbsp = sumbsp + 0.5*(t(m+1)-t(m))*bsp
          end do
          sum2 = sum2 + sum1
          sumbsp2 = sumbsp2 + sumbsp
          Ham(ii,jj) = sum2
          Bs(ii,jj) = sumbsp2
       end do
    end do
             
    call dsygvd( ITYPE, JOBZ, UPLO, NN, Ham, LDA, Bs, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )

    !print *, W(:)
    energy = W(1)
   
    wfn=Ham(:,1)
    
    return
  end subroutine Hamiltonian_born
  
