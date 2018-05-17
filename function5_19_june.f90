PROGRAM hami


  !.. About
  !   this is 
  !
  !..

  use constants
  use mod_functions
  IMPLICIT NONE

  !.. input
  INTEGER, PARAMETER :: np = 290       !Number of knot-points
  INTEGER, PARAMETER :: np1 =286
  INTEGER, PARAMETER :: k = 6
  INTEGER, PARAMETER :: kord = 4   
  INTEGER, PARAMETER :: NN = 282       !Number of B-splines used
  INTEGER, PARAMETER :: xmax = 280                    !xmax
  integer, parameter :: xmax_test = 281
  INTEGER, PARAMETER :: xxmax = 2000

  !.. External functions
  real(kind(1.d0)) :: bget,bder,bder2,ddot

 
  !.. Local parameters
  REAL*8, DIMENSION(xxmax) :: xx
  REAL*8, DIMENSION(np) :: t          !Knot-sequence
  REAL*8, DIMENSION(np1) :: t1 !Knot-sequence, for calculating the potential
  REAL*8, ALLOCATABLE :: Mat(:,:), Vold(:,:), plot(:,:)
  REAL*8, DIMENSION(NN,NN) :: Hamil, Bmat, Hamil2, Bmat2, Ham, Bs, Ham2, Bs2
  REAL*8, DIMENSION(NN,3) :: f, P
  real(kind(1.d0)) :: Bord_pot(xmax+1),Bder_pot(xmax+1),Bder2_pot(xmax+1),Bord_norm(xmax+1),Bder_norm(xmax+1),Bder2_norm(xmax+1)
  REAL*8, ALLOCATABLE :: x(:),rho(:),RHS(:),mn(:,:),rho_plot(:),Vold_func(:,:),Vexc_func(:,:),Vexc(:,:)
  real(kind(1.d0)) :: base,baseder,baseder2,pot,potex,potexder2,xval
  real(kind(1.d0)) :: rho_f(xmax),rho_p(xxmax)
  REAL*8 :: h,hh,Rmax,rmin,z,l,l2
  real(kind(1.d0)) :: xabsc(k), weig(k),noccsum,nocc,lower, upper,sum1,sum2,term4,term5,vvv
  INTEGER :: i,j,m,n,loop

  !Setting paramenters for generalized eigensolver
  INTEGER :: ITYPE, LDA, LDB, INFO
  INTEGER, PARAMETER :: LWORK = 1000000
  INTEGER, PARAMETER :: LIWORK = 5000
  CHARACTER*1 :: JOBZ = 'V', UPLO = 'U'
  REAL*8, DIMENSION(NN) :: W, W2
  REAL*8, DIMENSION(LWORK) :: WORK
  REAL*8, DIMENSION(LIWORK) :: IWORK

  !Setting paramenters for LU
  INTEGER :: NRHS, LDA1, LDB1, INFO1
  REAL*8, ALLOCATABLE :: ipiv(:)

  !Setting parameters for dddot
  INTEGER :: incx=1,incy=1

    Rmax = 45.d0
    rmin = 0.d0
    h = (Rmax - rmin)/(NN-3)

    z = 10.d0
    l = 0.d0
    l2 = 1.d0


    ALLOCATE (Mat(xmax_test,xmax_test))
    ALLOCATE (x(xmax))
    ALLOCATE (rho(xmax))
    ALLOCATE (rho_plot(xxmax))
    ALLOCATE (RHS(xmax_test))
    ALLOCATE (ipiv(xmax_test))
    ALLOCATE (mn(NN,k))
    ALLOCATE (Vexc(NN,NN))
    ALLOCATE (Vold(NN,NN))
    ALLOCATE (Vexc_func(NN,NN))
    ALLOCATE (Vold_func(NN,NN))
    ALLOCATE (plot(xxmax,NN))
 
      
    !Setting up knot-sequence
    do i = 1, np
       if(i .le. k) then
          t(i) = 0.d0
       else if(i .le. (np-k+1)) then
          t(i) = t(i-1) + h
       else
          t(i) = t(i-1)
       end if
    end do
 
    !Setting up new knot-sequence

    do i = 1, np1
       t1(i) = t(i+2)
    end do

    do i = 1, xmax
       x(i) = t(i+5)
    end do

    hh = (10.0 - rmin)/(xxmax-1)
    do i = 2, xxmax
      xx(1) = 0.d0
      xx(i) = xx(i-1) + hh
   end do
   
   call gauleg(k, xabsc, weig)
 
   Vexc = 0.0
   Hamil = 0.0
   Bmat = 0.0
   Hamil2 = 0.0
   Bmat2 = 0.0
   RHS = 0.0
   Vold = 0.0
   !STARTING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !STARTING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !STARTING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   DO loop = 1, 1
  
      Mat = 0.d0

      call func(z,l,np,k,NN,t,Vold,Ham,Bs)
      call func(z,l2,np,k,NN,t,Vold,Ham2,Bs2)
      Hamil = Ham
      Bmat = Bs
      Hamil2 = Ham2
      Bmat2 = Bs2


  
      !Calling generalized eigenvalue problem solver
      ITYPE = 1
      LDA = NN
      LDB = NN
      CALL DSYGVD( ITYPE, JOBZ, UPLO, NN, Hamil, LDA, Bmat, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )
      CALL DSYGVD( ITYPE, JOBZ, UPLO, NN, Hamil2, LDA, Bmat2, LDB, W2, WORK, LWORK, IWORK, LIWORK, INFO )

      OPEN(33,file='eigen.dat')
      DO i = 1, NN
      WRITE(33,30)loop, W(i), W2(i)
30    FORMAT(I3,'  ',16f14.8)
      END DO
      CLOSE(33)
      PRINT *, loop, W(1), W(2), W2(1)

      DO i = 1, NN
         f(i,1) = Hamil(i,1)
         f(i,2) = Hamil(i,2)
         f(i,3) = Hamil2(i,1)
      END DO

      open(11,file='res2.dat',status='replace')
      call trapezoid_uniform(xmax,x,np,NN,k,t,f,rho_f)
 
      rho = rho_f
      do i = 1, xmax
         write(11,30)i,x(i),4.0*PI*(x(i)**2.0)*rho(i)
      end do
      close(11)

      !Testing that the integral of rho = z
      noccsum = (x(xmax)/(2.0*(xmax-1)))*(rho(1)+rho(xmax))
      do i = 2, xmax-1
         nocc = 4.d0*PI*(x(xmax)/(xmax-1))*rho(i)*x(i)**2.d0
         noccsum = noccsum + nocc
      end do

      print *,' noccsum'
      print *, noccsum

      DO i = 1,xmax+1
         IF(i .LE. xmax) THEN
            RHS(i) = -4.d0*PI*rho(i)*x(i)
         ELSE
            RHS(i) = noccsum
         END IF
      !   WRITE(*,*) i, RHS(i)
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calulating Vdir and Vexc
      DO j = 1, xmax
         DO i = 1, xmax
            Mat(i,j) = bder2(t1(i+3),t1,kord,np1,j+1)
         END DO
      END DO
      Mat(xmax,xmax+1) = bder2(t1(xmax+3),t1,kord,np1,xmax+2)
      Mat(xmax+1,xmax+1) = bget(t1(xmax+3),t1,kord,np1,xmax+2)

     ! Calling rutine for LU-factorization, returns the coefficient vector
      NRHS = 1
      LDA1 = xmax_test
      LDB1 = xmax_test
      CALL DGESV(xmax_test, NRHS, Mat, LDA1, ipiv, RHS, LDB1, INFO1)

   
      DO j = 1, NN
         DO i = 1, NN
            sum1 = 0.d0
            sum2 = 0.d0
            term4 = 0.d0
            term5 = 0.d0
            base = 0.d0
            pot = 0.d0
            baseder = 0.d0
            potex = 0.d0
            potexder2 = 0.d0
            baseder2 = 0.d0
            vvv = 0.d0
            upper = MIN(i+1,j+1) + k-1 
            lower = MAX(i+1,j+1)
            DO m = lower, upper
               term4 = 0.d0
               term5 = 0.d0
               base = 0.d0
               pot = 0.d0
               baseder = 0.d0
               potex = 0.d0
               potexder2 = 0.d0
               baseder2 = 0.d0
               vvv = 0.d0
               DO n = 1, k
                  IF(t(m+1).GT.t(m)) THEN
                     mn(m,n) = 0.5d0*(t(m+1)+t(m)) + 0.5d0*(t(m+1)-t(m))*xabsc(n)
                     xval = mn(m,n)
                     CALL baseR(xval,xmax,kord,np1,t1,Bord_pot,Bder_pot,Bder2_pot)
                     base = ddot(xmax_test,Bord_pot,incx,RHS,incy)
                     pot=base/mn(m,n)
                     term4 = term4 + weig(n)*bget(mn(m,n),t,k,np,i+1)*bget(mn(m,n),t,k,np,j+1)*pot
                     baseder = ddot(xmax_test,Bder_pot,incx,RHS,incy)
                     potex=baseder/(mn(m,n)**2.d0)
                     baseder2 = ddot(xmax_test,Bder2_pot,incx,RHS,incy)
                     potexder2 = baseder2/mn(m,n)
                     vvv = ABS((-(6.d0*pot/(mn(m,n)**2.d0))+(6.d0*potex)-(3.d0*potexder2))/(32.d0*(PI**2.d0)))**(1.d0/3.d0)
                     term5 = term5 + weig(n)*bget(mn(m,n),t,k,np,i+1)*bget(mn(m,n),t,k,np,j+1)*(-3.d0*vvv)
                  END IF 
               END DO 
               sum1 = sum1 + 0.5d0*(t(m+1)-t(m))*(term4+term5)
              ! print *, sum1,xval
            END DO
            sum2 = sum2 + sum1
            Vexc(i,j) = sum2
           ! print *,sum2
         END DO
      END DO
      !.. Debugging
     ! print *, 'shape',shape(Vexc),xabsc
     ! stop
    
      Vold = (1.d0-0.4d0)*Vexc+0.4d0*Vold

      call electronic_pot(np,NN,k,t,xmax,kord,np1,t1,xmax_test,xabsc,weig,RHS,Vexc_func)

      Vold_func = Vexc_func

       !.. Debugging
    
      print *, Vold_func(1,1),Vexc(1,1),Vold_func(NN,NN),Vexc(NN,NN)
     

      stop

   END DO

   !!. Plotting electron potential

   open(37,file='res.dat',status='replace')
   call trapezoid_uniform(xxmax,xx,np,NN,k,t,f,rho_p)
   rho_plot = rho_p
   do i = 1, xxmax
      write(37,30)i,xx(i),4.d0*PI*(xx(i)**2.d0)*rho_plot(i)
   end do
   close(37)

 end program hami
