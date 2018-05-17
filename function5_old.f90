PROGRAM hami


  !.. About
  !   this is 
  !
  !..

 
  IMPLICIT NONE

  !.. input
  INTEGER, PARAMETER :: np = 290       !Number of knot-points
  INTEGER, PARAMETER :: np1 =286
  INTEGER, PARAMETER :: k = 6
  INTEGER, PARAMETER :: kord = 4   
  INTEGER, PARAMETER :: NN = 282       !Number of B-splines used
  INTEGER, PARAMETER :: xmax = 280                    !xmax
  INTEGER, PARAMETER :: xxmax = 10000
  
  !.. Local parameters
  REAL*8, DIMENSION(xxmax) :: xx
  REAL*8, DIMENSION(np) :: t          !Knot-sequence
  REAL*8, DIMENSION(np1) :: t1 !Knot-sequence, for calculating the potential
  REAL*8, ALLOCATABLE :: Mat(:,:), Vexc(:,:), Vold(:,:), plot(:,:),Vexcsum(:)
  REAL*8, DIMENSION(NN,NN) :: Hamil, Bmat, Hamil2, Bmat2
  REAL*8, DIMENSION(NN,3) :: f, P
  REAL*8, ALLOCATABLE :: mn(:,:),Baseone(:,:),base(:),baseder(:),baseder2(:),res1(:),res2(:),res3(:),x(:),rho(:),RHS(:)
  REAL*8 :: PI = 3.14159265358979
  REAL*8 :: h,hh,Rmax,rmin,Rmax_first,rmin_first,h_first,z,l,l2
  REAL*8 :: usum1, usum2, usum3,usum1_1, usum2_1, usum3_1,usum1_2, usum2_2, usum3_2,u1,u2,u3,u1_1,u2_1,u3_1,u1_2,u2_2,u3_2,xabsc(k), weig(k),noccsum,nocc,lower, upper,sum1,sum2,term4,term5,pot,potex,potexder2,vvv
  INTEGER :: i,j,m,n,loop, count

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


    
    Rmax = 40.0
    rmin = 0.0
    h = (Rmax - rmin)/(NN-3)

    z = 10.0
    l = 0.0
    l2 = 1.0


    ALLOCATE (Mat(xmax+1,xmax+1))
    ALLOCATE (Baseone(xmax,NN))
    ALLOCATE (x(xmax))
    ALLOCATE (res1(xmax))
    ALLOCATE (res2(xmax))
    ALLOCATE (res3(xmax))
    ALLOCATE (rho(xmax))
    ALLOCATE (RHS(xmax+1))
    ALLOCATE (ipiv(xmax+1))
    ALLOCATE (mn(NN,k))
    ALLOCATE (base(xmax+1))
    ALLOCATE (baseder(xmax+1))
    ALLOCATE (baseder2(xmax+1))
    ALLOCATE (Vexc(NN,NN))
    ALLOCATE (Vold(NN,NN))
    ALLOCATE (plot(xxmax,NN))
 
      
    !Setting up knot-sequence
    DO i = 1, np
       IF(i .LE. k) THEN
          t(i) = 0.0
       ELSE IF(i .LE. (np-k+1)) THEN
          t(i) = t(i-1) + h
       ELSE
          t(i) = t(i-1)
       END IF
   !    PRINT *, i, t(i)
    END DO

    !Setting up new knot-sequence

    PRINT *,'t1'
    DO i = 1, np1
       t1(i) = t(i+2)
       PRINT *, t1
    END DO

    ! PRINT *, 'x'
    DO i = 1, xmax
       x(i) = t(i+5)
  !        PRINT *,i, x(i)
    END DO

   hh = (Rmax - rmin)/(xxmax-2)
   DO i = 2, xxmax
      xx(1) = 0.00001
      xx(i) = xx(i-1) + hh
   END DO
   
   CALL gauleg(k, xabsc, weig)

   
 !   PRINT *, xabsc
 
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
   DO loop = 1, 20
  
      Mat = 0.0

      Hamil = func(z,l,np,k,NN,t,Vold)
      Bmat =  bfunc(z,l,np,k,NN,t)
      Hamil2 = func(z,l2,np,k,NN,t,Vold)
      Bmat2 = bfunc(z,l2,np,k,NN,t)


      !..Debugg
      
    !  stop
      
     
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !TRAPEZOID INTEGRATION, !Testing that the integral of rho = z
  DO j = 1, NN
     DO i = 1, xmax
        Baseone(i,j) = bget(x(i),t,k,np,j+1)
     END DO
  END DO

  OPEN(11,file='res2.dat',status='replace')
  res1=matmul(Baseone,f(:,1))
  res2=matmul(Baseone,f(:,2))
  res3=matmul(Baseone,f(:,3))
  usum1 = Rmax*(res1(1)**2.0+res1(xmax)**2.0)/xmax
  usum2 = Rmax*(res2(1)**2.0+res2(xmax)**2.0)/xmax
  usum3 = Rmax*(res3(1)**2.0+res3(xmax)**2.0)/xmax
  DO i = 2, xmax-1
     u1 = (Rmax/(2.0*(xmax-1)))*2.0*res1(i)**2.0
     u2 = (Rmax/(2.0*(xmax-1)))*2.0*res2(i)**2.0
     u3 = (Rmax/(2.0*(xmax-1)))*2.0*res3(i)**2.0
     usum1 = usum1 + u1
     usum2 = usum2 + u2
     usum3 = usum3 + u3
  END DO
  rho=(2.0*(res1**2.0)/usum1 + 2.0*(res2**2.0)/usum2 + 6.0*(res3**2.0)/usum3)/(4.0*PI*(x**2.0))
  rho(1) = 0.0
  DO i=1, xmax
     WRITE(11,30)i, x(i), 4.0*PI*(x(i)**2.0)*rho(i)
  END DO
  WRITE(11,*)
  CLOSE(11)

  !Testing that the integral of rho = z
  noccsum = (Rmax/(xmax-1))*(rho(1)+rho(xmax))
  DO i = 2, xmax-1
     nocc = 4.0*PI*(Rmax/(2.0*(xmax-1)))*2.0*rho(i)*x(i)**2.0
     noccsum = noccsum + nocc
  END DO
  PRINT *,' noccsum'
  PRINT *, noccsum

      DO i = 1,xmax+1
         IF(i .LE. xmax) THEN
            RHS(i) = -4*PI*rho(i)*x(i)
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
      LDA1 = xmax+1
      LDB1 = xmax+1
      CALL DGESV(xmax+1, NRHS, Mat, LDA1, ipiv, RHS, LDB1, INFO1)

   
      DO j = 1, NN
         DO i = 1, NN
            sum1 = 0.0
            sum2 = 0.0
            term4 = 0.0
            term5 = 0.0
            base = 0.0
            pot = 0.0
            baseder = 0.0
            potex = 0.0
            potexder2 = 0.0
            baseder2 = 0.0
            vvv = 0.0
            upper = MIN(i+1,j+1) + k-1 
            lower = MAX(i+1,j+1)
            DO m = lower, upper
               term4 = 0.0
               term5 = 0.0
               base = 0.0
               pot = 0.0
               baseder = 0.0
               potex = 0.0
               potexder2 = 0.0
               baseder2 = 0.0
               vvv = 0.0
               DO n = 1, k
                  IF(t(m+1).GT.t(m)) THEN
                     mn(m,n) = 0.5*(t(m+1)+t(m)) + 0.5*(t(m+1)-t(m))*xabsc(n)
                     base = MATMUL(baseR(mn(m,n),xmax,kord,np1,t1),RHS)
                     pot=SUM(base)/mn(m,n)
                     term4 = term4 + weig(n)*bget(mn(m,n),t,k,np,i+1)*bget(mn(m,n),t,k,np,j+1)*pot
                     baseder = MATMUL(baseEx(mn(m,n),xmax,kord,np1,t1),RHS)
                     potex=SUM(baseder)/(mn(m,n)**2.0)
                     baseder2 = MATMUL(baseExder2(mn(m,n),xmax,kord,np1,t1),RHS)
                     potexder2 = SUM(baseder2)/mn(m,n)
                     vvv = ABS((-(6.0*pot/(mn(m,n)**2.0))+(6.0*potex)-(3.0*potexder2))/(32.0*(PI**2.0)))**(1./3.)
                     term5 = term5 + weig(n)*bget(mn(m,n),t,k,np,i+1)*bget(mn(m,n),t,k,np,j+1)*(-3.0*vvv)
                  END IF
               END DO
               sum1 = sum1 + 0.5*(t(m+1)-t(m))*(term4+term5)
            END DO
            sum2 = sum2 + sum1
            Vexc(i,j) = sum2
         END DO
      END DO

      Vold = (1-0.3)*Vexc+0.3*Vold
   END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !EXPECTATION_VALUE CALCULATION !!!!TRAPEZOID INTEGRATION of c1*P*Vee*c1*P
 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! PLOTTING ELECTRON POTENTIAL

 plot = 0.0
 
 DO j = 1, NN
    DO i = 1, xxmax
       plot(i,j) = bget(xx(i),t,k,np,j+1)
    END DO
 END DO
 
 usum1 = 0.0
 usum2 = 0.0
 usum3 = 0.0
 u1 = 0.0
 u2 = 0.0
 u3 = 0.0

 OPEN(37,file='res1.dat',status='replace')
      res1=matmul(plot,f(:,1))
      res2=matmul(plot,f(:,2))
      res3=matmul(plot,f(:,3))
      usum1 = x(xmax)*(res1(1)**2.0+res1(xxmax)**2.0)/(2.0*xxmax)
      usum2 = x(xmax)*(res2(1)**2.0+res2(xxmax)**2.0)/(2.0*xxmax)
      usum3 = x(xmax)*(res3(1)**2.0+res3(xxmax)**2.0)/(2.0*xxmax)
      DO i = 2, xxmax-1
         u1 = (x(xmax)/(xxmax-1))*res1(i)**2.0
         u2 = (x(xmax)/(xxmax-1))*res2(i)**2.0
         u3 = (x(xmax)/(xxmax-1))*res3(i)**2.0
         usum1 = usum1 + u1
         usum2 = usum2 + u2
         usum3 = usum3 + u3
      END DO
      rho=(2*(res1**2)/usum1 + 2*(res2**2)/usum2 + 6*(res3**2)/usum3)/(4*PI*(xx**2))
      rho(1) = 0.0
   
      DO i=1, xxmax
         WRITE(37,30)i, xx(i), 4*PI*(xx(i)**2)*rho(i)
      END DO
      WRITE(37,*)
      CLOSE(37)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CONTAINS

    FUNCTION baseR(x,xmax,kord,np1,t1)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: np1, xmax, kord
     REAL*8, INTENT(IN) :: x
     REAL*8,DIMENSION(np1), INTENT(IN) :: t1
     REAL*8, ALLOCATABLE :: baseR(:,:), Base2(:,:)
     INTEGER :: j,i

     ALLOCATE(Base2(1,xmax+1))
     
     DO j = 1, xmax+1
        IF (j .LE. xmax) THEN
           Base2(1,j) = bget(x,t1,kord,np1,j+1)
        ELSE
           Base2(1,j) = Base2(1,j-1)
           END IF
     END DO
     baseR=Base2
     DEALLOCATE(Base2)
     RETURN
   END FUNCTION baseR

   FUNCTION baseEx(x,xmax,kord,np1,t1)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: np1, xmax, kord
     REAL*8, INTENT(IN) :: x
     REAL*8,DIMENSION(np1), INTENT(IN) :: t1
     REAL*8, ALLOCATABLE :: baseEx(:,:), Base2(:,:)
     INTEGER :: j,i

   
     ALLOCATE(Base2(1,xmax+1))
     
     DO j = 1, xmax+1
        IF (j .LE. xmax) THEN
           Base2(1,j) = bder(x,t1,kord,np1,j+1)
        ELSE
           Base2(1,j) = Base2(1,j-1)
           END IF
     END DO
     baseEx=Base2
     DEALLOCATE(Base2)
     RETURN
   END FUNCTION baseEx

   FUNCTION baseExder2(x,xmax,kord,np1,t1)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: np1, xmax, kord
     REAL*8, INTENT(IN) :: x
     REAL*8,DIMENSION(np1), INTENT(IN) :: t1
     REAL*8, ALLOCATABLE :: baseExder2(:,:), Base2(:,:)
     INTEGER :: j,i

     ALLOCATE(Base2(1,xmax+1))
     
     DO j = 1, xmax+1
        IF (j .LE. xmax) THEN
           Base2(1,j) = bder2(x,t1,kord,np1,j+1)
        ELSE
           Base2(1,j) = Base2(1,j-1)
           END IF
     END DO
     baseExder2=Base2
     DEALLOCATE(Base2)
   END FUNCTION baseExder2

   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  FUNCTION func(z,l,np,k,NN,t,Vold)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: np,k,NN          !Number of knot-points etc.
    REAL*8, INTENT(IN) :: z, l
    REAL*8, DIMENSION(NN,NN), INTENT(IN) :: Vold
    REAL*8, DIMENSION(np), INTENT(IN) :: t
    REAL*8, ALLOCATABLE :: Ham(:,:), Bs(:,:), func(:,:), mn(:,:) !Hamiltonian and B-matrix
    REAL*8 :: sum1, sum2,term1, term2, term3,lower, upper, sumbsp,sumbsp2, bsp,  xabsc(k), weig(k)
    INTEGER :: i, j, m, n

    ALLOCATE(Ham(NN,NN))
    ALLOCATE(Bs(NN,NN))
    ALLOCATE(mn(NN,k))

   !Parameters for quantum numbers, charge and stuff
    Ham = 0.0
    Bs = 0.0

    CALL gauleg(k, xabsc, weig)
  
    DO j = 1, NN
       DO i = 1, NN
          sum1 = 0.0
          sum2 = 0.0
          term1 = 0.0
          term2 = 0.0
          term3 = 0.0
          upper = MIN(i+1,j+1) + k-1
          lower = MAX(i+1,j+1)
          DO m = lower, upper
             term1 = 0.0
             term2 = 0.0
             term3 = 0.0
             DO n = 1, k
                IF(t(m+1).GT.t(m)) THEN
                   mn(m,n) = 0.5*(t(m+1)+t(m)) + 0.5*(t(m+1)-t(m))*xabsc(n)
                   term1 = term1 + weig(n)*bder(mn(m,n),t,k,np,i+1)*bder(mn(m,n),t,k,np,j+1)/2.0
                   term2 = term2 + weig(n)*l*(l+1)*bget(mn(m,n),t,k,np,j+1)*bget(mn(m,n),t,k,np,i+1)/(2*(mn(m,n)**2))
                   term3 = term3 - weig(n)*z*bget(mn(m,n),t,k,np,i+1)*bget(mn(m,n),t,k,np,j+1)/mn(m,n)
                END IF
             END DO
             sum1 = sum1 + 0.5*(t(m+1)-t(m))*(term1+term2+term3)
          END DO
          sum2 = sum2 + sum1
          Ham(i,j) = sum2
       END DO
    END DO
    func = Ham + Vold
    
    DEALLOCATE(Ham);DEALLOCATE(Bs);DEALLOCATE(mn);
    RETURN
  END FUNCTION func

  FUNCTION bfunc(z,l,np,k,NN,t)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: np           !Number of knot-points
    INTEGER, INTENT(IN) :: k         
    INTEGER, INTENT(IN) :: NN           !Number of B-splines used
    REAL*8, DIMENSION(np), INTENT(IN) :: t
    REAL*8, ALLOCATABLE :: Bs(:,:), bfunc(:,:), mn(:,:)!B-matrix
    REAL*8, INTENT(IN) :: z, l
    REAL*8 :: lower, upper, sumbsp,sumbsp2, bsp,xabsc(k), weig(k)
    INTEGER :: i, j, m, n

    ALLOCATE(Bs(NN,NN))
    ALLOCATE(mn(NN,k))

    Bs = 0.0
    CALL gauleg(k, xabsc, weig)

    DO j = 1, NN
       DO i = 1, NN
          sumbsp = 0.0
          sumbsp2 = 0.0
          bsp = 0.0
          upper = MIN(i+1,j+1) + k-1
          lower = MAX(i+1,j+1)
          DO m = lower, upper
             bsp = 0.0
             DO n = 1, k
                IF(t(m+1).GT.t(m)) THEN
                   mn(m,n) = 0.5*(t(m+1)+t(m)) + 0.5*(t(m+1)-t(m))*xabsc(n)
                   bsp = bsp + weig(n)*bget(mn(m,n),t,k,np,i+1)*bget(mn(m,n),t,k,np,j+1)
                END IF
             END DO
             sumbsp = sumbsp + 0.5*(t(m+1)-t(m))*bsp
          END DO
          sumbsp2 = sumbsp2 + sumbsp
          Bs(i,j) = sumbsp2
       END DO
    END DO

    bfunc = Bs
    DEALLOCATE(Bs);DEALLOCATE(mn);
    RETURN
  END FUNCTION bfunc



      function bget(rr,t,kord,np,index)

      implicit none
      integer i,j,left,it,kord,index,np
      real*8 bget,t(np),Sp(kord),rr
     
      bget=0.d0

      if(rr.gt.t(np)) return
      if(rr.lt.t(1)) return
      if(abs(rr-t(np)).lt.1.d-10) then
        if(index.eq.np-kord) bget=1.d0
        return
      end if      
      do it=1,np
        if(rr.ge.t(it)) left=it 
      end do   
      if(index-left+kord.lt.1.or.index-left+kord.gt.kord) return         

      call bsplvb(t,kord,1,rr,left,Sp)
      bget=Sp(index-left+kord)
      return
    end function bget

    function bder(rr,t,kord,np,index)

      implicit none
      integer i,j,left,it,kord,index,np
      real*8 bder,t(np),Sp(kord),rr,deri
     
      bder=0.d0

      if(rr.gt.t(np)) return
      if(rr.lt.t(1)) return
      do it=1,np
        if(rr.ge.t(it)) left=it 
      end do   

      if(abs(rr-t(np)).lt.1.d-10) then
        if(index.lt.np-kord-1) return
       
        if(index.eq.np-kord) then
          bder=dble(kord-1)/(t(np)-t(np-kord))
        else if(index.eq.np-kord-1) then
          bder=-dble(kord-1)/(t(np)-t(np-kord))
        end if
        return
      end if
      

      if(index-left+kord.lt.1.or.index-left+kord.gt.kord) return   
      

      call bsplvb(t,kord-1,1,rr,left,Sp)
      i=index-left+kord


      if(i.eq.1) then
        deri=dble(kord-1)*(-Sp(i)/(t(index+kord)-t(index+1))) 
      else if(i.eq.kord) then
        deri=dble(kord-1)*(Sp(i-1)/(t(index+kord-1)-t(index)))    
      else 
        deri=dble(kord-1)*(Sp(i-1)/(t(index+kord-1)-t(index))-Sp(i)/(t(index+kord)-t(index+1)))    
      end if
      bder=deri
      return
    end function bder

      function bder2(rr,t,kord,np,index)
!*     returns the value of (d/dx) Bspline(kord,index) in rr
!*     The first Bspline is called spline #1 (i.e. index=1)
!*     the first knot point is in t(1)
!*     np= number of knot points (distinct or multiple) including
!*     the ghost points: N phyical points np=N +2*(kord-1)  

      implicit none
      integer i,j,left,it,kord,index,np
      real*8 bder2,t(np),Sp(kord),rr,deri
     
      bder2=0.d0
!*     if rr=t(np) then the routine assumes that
!*     there is kord knotpoints in the last physical point and
!*     returns Bder.ne.zero if index is np-kord, or np-kord-1 
      if(rr.gt.t(np)) return
      if(rr.lt.t(1)) return
      do it=1,np
        if(rr.ge.t(it)) left=it 
      end do   

      if(abs(rr-t(np)).lt.1.d-10) then
        if(index.lt.np-kord-2) return       
        if(index.eq.np-kord) then
          bder2=dble((kord-1)*(kord-2))/(t(np-1)-t(np-kord))/(t(np-2)-t(np-kord))
        else if(index.eq.np-kord-2) then
          bder2=dble((kord-1)*(kord-2))/(t(np-1)-t(np-kord-1))/(t(np-2)-t(np-kord))
        else if(index.eq.np-kord-1) then
          bder2=-dble((kord-1)*(kord-2))/(t(np-2)-t(np-kord-1))/(t(np-2)-t(np-kord))-dble((kord-1)*(kord-2))/(t(np-1)-t(np-kord))/(t(np-2)-t(np-kord))
        end if
        return
      end if
      

      if(index-left+kord.lt.1.or.index-left+kord.gt.kord) return   
      
!*     index=left-kord+i => i=index-left+kord  for Sp with k=kord
!*     index=left-kord+i+1 => i=index-left+kord-1  for Sp with k=kord-1
!*     Sp(i-1,k-1) gives the same index as Sp(i,k)
!*     Sp(i,k-1) gives the same index as Sp(i+1.k)
      call bsplvb(t,kord-2,1,rr,left,Sp)
      i=index-left+kord
      deri=0.d0
      if(i.gt.2) deri=deri+dble(kord-1)*dble(kord-2)*Sp(i-2)/(t(index+kord-2)-t(index))/(t(index+kord-1)-t(index))
      if(i.gt.1.and.i.lt.kord) deri=deri-dble(kord-1)*dble(kord-2)*(Sp(i-1)/(t(index+kord-1)-t(index+1))/(t(index+kord-1)-t(index))+Sp(i-1)/(t(index+kord-1)-t(index+1))/(t(index+kord)-t(index+1)))
       if(i.lt.kord-1) deri=deri+dble(kord-1)*dble(kord-2)*Sp(i)/(t(index+kord)-t(index+2))/(t(index+kord)-t(index+1))     
      bder2=deri
      return
    end function bder2

    END PROGRAM

   SUBROUTINE  gauleg(ngp, xabsc, weig)

   implicit none
   REAL*8 :: newv
   REAL*8, PARAMETER  :: EPS=3.0d-15 !EPS is the relative precision
   REAL*8,PARAMETER :: M_PI=3.141592654d0      ! Pi value
   INTEGER :: i, j, m
   REAL*8 ::  p1, p2, p3, pp, z, z1
   INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
   REAL*8, INTENT(OUT) :: xabsc(ngp), weig(ngp)
   m = (ngp + 1) / 2
   !* Roots are symmetric in the interval - so only need to find half of them  */
   do i = 1, m				! Loop over the desired roots */
      z = cos( M_PI * (i-0.25) / (ngp+0.5) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */
100     	p1 = 1.0
        	p2 = 0.0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

        	do j = 1, ngp
           	p3 = p2
           	p2 = p1
           	p1 = ((2.0*j-1.0) * z * p2 - (j-1.0)*p3) / j
        	enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
        	pp = ngp*(z*p1-p2)/(z*z-1.0)
        	z1 = z
        	z = z1 - p1/pp             ! Newton's Method  */

        	if (dabs(z-z1) .gt. EPS) GOTO  100

      	xabsc(i) =  - z    ! Roots will be bewteen -1.0 & 1.0 */
      	xabsc(ngp+1-i) =  + z           ! and symmetric about the origin  */
      	weig(i) = 2.0/((1.0-z*z)*pp*pp) ! Compute the weight and its       */
      	weig(ngp+1-i) = weig(i)         ! symmetric counterpart         */

      end do     ! i loop

    End subroutine gauleg

            SUBROUTINE bsplvb(t,jhigh,index,x,left,biatx)

      PARAMETER (JMAX=100)
      integer index,jhigh,left,i,j,jp1
      real*8 t,x,biatx,deltal,deltar,saved,term
      DIMENSION biatx(jhigh),t(left+jhigh),deltal(jmax),deltar(jmax)
      SAVE deltal,deltar
      DATA j/1/

      GO TO (10,20),index
 10   j = 1
      biatx(1) = 1.d0
      IF (j .GE. jhigh) GO TO 99

 20   CONTINUE
         jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)

         saved = 0.d0
         DO i = 1,j

             term = biatx(i)/(deltar(i) + deltal(jp1-i))
             biatx(i) = saved + deltar(i)*term
             saved = deltal(jp1-i)*term
         END DO
         biatx(jp1) = saved
         j = jp1
         IF (j .LT. jhigh) GO TO 20
 99   RETURN
       END SUBROUTINE bsplvb
  
