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
  REAL*8, DIMENSION(NN,NN) :: Hamil_0,Bmat_0,Hamil_1,Bmat_1,Hamil_2,Bmat_2,Ham_0,Bs_0,Ham_1,Bs_1,Ham_2,Bs_2
  REAL*8, DIMENSION(NN,7) :: f, P
  REAL*8, ALLOCATABLE :: x(:),rho(:),RHS(:),rho_plot(:),Vexc(:,:)
  real(kind(1.d0)) :: rho_f(xmax),rho_p(xxmax),f_new(NN,7),f_new2(1,7)
  real(kind(1.d0)) :: h,hh,Rmax,rmin,z,l_0,l_1,l_2,occ_1s,occ_2s,occ_2p,occ_3s,occ_3p,occ_3d,occ_4s
  real(kind(1.d0)) :: xabsc(k), weig(k),noccsum,nocc
  INTEGER :: i,j,m,n,loop

  !Setting paramenters for generalized eigensolver
  INTEGER :: ITYPE, LDA, LDB, INFO
  INTEGER, PARAMETER :: LWORK = 1000000
  INTEGER, PARAMETER :: LIWORK = 5000
  CHARACTER*1 :: JOBZ = 'V', UPLO = 'U'
  REAL*8, DIMENSION(NN) :: W_0,W_1,W_2
  REAL*8, DIMENSION(LWORK) :: WORK
  REAL*8, DIMENSION(LIWORK) :: IWORK

  !Setting paramenters for LU
  INTEGER :: NRHS, LDA1, LDB1, INFO1
  REAL*8, ALLOCATABLE :: ipiv(:)

  !Setting parameters for dddot
  INTEGER :: incx=1,incy=1

    Rmax = 10.d0
    rmin = 0.d0
    h = (Rmax - rmin)/(NN-3)

    write(*,*)'Number of protons = '
    read(*,*)z
    write(*,*)'Occupation number for 1s'
    read(*,*)occ_1s
    write(*,*)'Occupation number for 2s'
    read(*,*)occ_2s
    write(*,*)'Occupation number for 2p'
    read(*,*)occ_2p
    write(*,*)'Occupation number for 3s'
    read(*,*)occ_3s
    write(*,*)'Occupation number for 3p'
    read(*,*)occ_3p
    write(*,*)'Occupation number for 3d'
    read(*,*)occ_3d
    write(*,*)'Occupation number for 4s'
    read(*,*)occ_4s
    
    l_0 = 0.d0
    l_1 = 1.d0
    l_2 = 2.d0


    ALLOCATE (Mat(xmax_test,xmax_test))
    ALLOCATE (x(xmax))
    ALLOCATE (rho(xmax))
    ALLOCATE (rho_plot(xxmax))
    ALLOCATE (RHS(xmax_test))
    ALLOCATE (ipiv(xmax_test))
    ALLOCATE (Vexc(NN,NN))
    ALLOCATE (Vold(NN,NN))
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
   Hamil_0 = 0.0
   Bmat_0 = 0.0
   Hamil_1 = 0.0
   Bmat_1 = 0.0
   Hamil_2 = 0.0
   Bmat_2 = 0.0
   RHS = 0.0
   Vold = 0.0

   !.. Strating do loop
   do loop = 1, 20
  
      Mat = 0.d0

      call func(z,l_0,np,k,NN,t,Vold,Ham_0,Bs_0)
      call func(z,l_1,np,k,NN,t,Vold,Ham_1,Bs_1)
      call func(z,l_2,np,k,NN,t,Vold,Ham_2,Bs_2)
      Hamil_0 = Ham_0
      Bmat_0 = Bs_0
      Hamil_1 = Ham_1
      Bmat_1 = Bs_1
      Hamil_2 = Ham_2
      Bmat_2 = Bs_2


  
      !Calling generalized eigenvalue problem solver
      ITYPE = 1
      LDA = NN
      LDB = NN
      call dsygvd( ITYPE, JOBZ, UPLO, NN, Hamil_0, LDA, Bmat_0, LDB, W_0, WORK, LWORK, IWORK, LIWORK, INFO )
      call dsygvd( ITYPE, JOBZ, UPLO, NN, Hamil_1, LDA, Bmat_1, LDB, W_1, WORK, LWORK, IWORK, LIWORK, INFO )
      call dsygvd( ITYPE, JOBZ, UPLO, NN, Hamil_2, LDA, Bmat_2, LDB, W_2, WORK, LWORK, IWORK, LIWORK, INFO )

      open(33,file='eigen.dat')
      do i = 1, NN
      write(33,30)loop, W_0(i), W_1(i)
30    format(I3,'  ',16f14.8)
      end do
      close(33)
      print *,'loop, ','1s, ','2s, ','2p, ','3s, ','3p, ','3d, ','4s, '
      print *, loop, W_0(1), W_0(2), W_1(1), W_0(3), W_1(2), W_2(1), W_0(4)

      do i = 1, NN
         f(i,1) = Hamil_0(i,1)
         f(i,2) = Hamil_0(i,2)
         f(i,3) = Hamil_1(i,1)
         f(i,4) = Hamil_0(i,3)
         f(i,5) = Hamil_1(i,2)
         f(i,6) = Hamil_2(i,1)
         f(i,7) = Hamil_0(i,4)
      end do

      open(11,file='res.dat',status='replace')
      call trapezoid_uniform(xmax,x,np,NN,k,t,f,rho_f,occ_1s,occ_2s,occ_2p,occ_3s,occ_3p,occ_3d,occ_4s)
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

      do i = 1,xmax+1
         if(i .LE. xmax) then
            RHS(i) = -4.d0*PI*rho(i)*x(i)
         else
            RHS(i) = noccsum
         end if
      end do

      !.. Calulating Vdir and Vexc
      do j = 1, xmax
         do i = 1, xmax
            Mat(i,j) = bder2(t1(i+3),t1,kord,np1,j+1)
         end do
      end do
      Mat(xmax,xmax+1) = bder2(t1(xmax+3),t1,kord,np1,xmax+2)
      Mat(xmax+1,xmax+1) = bget(t1(xmax+3),t1,kord,np1,xmax+2)

     ! Calling rutine for LU-factorization, returns the coefficient vector
      NRHS = 1
      LDA1 = xmax_test
      LDB1 = xmax_test
      call dgesv(xmax_test, NRHS, Mat, LDA1, ipiv, RHS, LDB1, INFO1)
    
      call electronic_pot(np,NN,k,t,xmax,kord,np1,t1,xmax_test,xabsc,weig,RHS,Vexc)

      Vold = (1.d0-0.4d0)*Vexc+0.4d0*Vold

   end do

   !.. Calculating the total energy

   f_new(:,1)=matmul(Vold,f(:,1))
   f_new(:,2)=matmul(Vold,f(:,2))
   f_new(:,3)=matmul(Vold,f(:,3))
   f_new(:,4)=matmul(Vold,f(:,4))
   f_new(:,5)=matmul(Vold,f(:,5))
   f_new(:,6)=matmul(Vold,f(:,6))
   f_new(:,7)=matmul(Vold,f(:,7))
   f_new2(:,1)=ddot(NN,f(:,1),incx,f_new(:,1),incy)
   f_new2(:,2)=ddot(NN,f(:,2),incx,f_new(:,2),incy)
   f_new2(:,3)=ddot(NN,f(:,3),incx,f_new(:,3),incy)
   f_new2(:,4)=ddot(NN,f(:,4),incx,f_new(:,4),incy)
   f_new2(:,5)=ddot(NN,f(:,5),incx,f_new(:,5),incy)
   f_new2(:,6)=ddot(NN,f(:,6),incx,f_new(:,6),incy)
   f_new2(:,7)=ddot(NN,f(:,7),incx,f_new(:,7),incy)

   print *,'expectation value/2 1s',f_new2(:,1)/2.d0
   print *,'expectation value/2 2s',f_new2(:,2)/2.d0
   print *,'expectation value/2 2p',f_new2(:,3)/2.d0
   print *,'expectation value/2 3s',f_new2(:,4)/2.d0
   print *,'expectation value/2 3p',f_new2(:,5)/2.d0
   print *,'expectation value/2 3d',f_new2(:,6)/2.d0
   print *,'expectation value/2 4s',f_new2(:,7)/2.d0
 
   print *,'total energy 1s',occ_1s*W_0(1)-occ_1s*f_new2(:,1)/2.d0
   print *,'total energy 2s',occ_2s*W_0(2)-occ_2s*f_new2(:,2)/2.d0
   print *,'total energy 2p',occ_2p*W_1(1)-occ_2p*f_new2(:,3)/2.d0
   print *,'total energy 3s',occ_3s*W_0(3)-occ_3s*f_new2(:,4)/2.d0
   print *,'total energy 3p',occ_3p*W_1(2)-occ_3p*f_new2(:,5)/2.d0
   print *,'total energy 3d',occ_3d*W_2(1)-occ_3d*f_new2(:,6)/2.d0
   print *,'total energy 4s',occ_4s*W_0(4)-occ_4s*f_new2(:,7)/2.d0
   print *,'sum',occ_1s*W_0(1)-occ_1s*f_new2(:,1)/2.d0 + occ_2s*W_0(2)-occ_2s*f_new2(:,2)/2.d0 +occ_2p*W_1(1)-occ_2p*f_new2(:,3)/2.d0 + occ_3s*W_0(3)-occ_3s*f_new2(:,4)/2.d0 + occ_3p*W_1(2)-occ_3p*f_new2(:,5)/2.d0 + occ_3d*W_2(1)-occ_3d*f_new2(:,6)/2.d0 + occ_4s*W_0(4)-occ_4s*f_new2(:,7)/2.d0 

   !.. Plotting electron potential

   open(37,file='res1.dat',status='replace')
   call trapezoid_uniform(xxmax,xx,np,NN,k,t,f,rho_p,occ_1s,occ_2s,occ_2p,occ_3s,occ_3p,occ_3d,occ_4s)
   rho_plot = rho_p
   do i = 1, xxmax
      write(37,30)i,xx(i),4.d0*PI*(xx(i)**2.d0)*rho_plot(i)
   end do
   close(37)

 end program hami
