program efimov

  use constants
  
  implicit none
  INTERFACE
     subroutine efimovham(npl,npm,k,l,m,lm,points,tl,tm,rho,energy,S,TK,H,Hder)
       use constants
       implicit none
    !.. Input
       integer            , intent(in) :: npl,npm,points,lm,l,k,m
       real(kind(1.d0))   , intent(in) :: tl(npl),tm(npm),rho(points)
       real(kind(1.d0)), intent(inout) :: energy(LM,points)
       real(kind(1.d0)), intent(inout) :: S(LM,LM),TK(LM,LM)
       real(kind(1.d0)),intent(inout)  :: H(LM,LM,points),Hder(LM,LM,points)
     END subroutine efimovham
     
  END INTERFACE
  
  !.. Input
  !.. Parameters for the B-splines used in the generalized eigenvalue equation
  integer, parameter :: N1 = 10   !.. Number of mesh-points in coordinate 1
  integer, parameter :: N2 = 10   !.. Number of mesh-points in coordinate 2
  integer, parameter :: k = 6    !.. B-spline order
  integer, parameter :: L = N1+k-4    !.. Number of B-splines in coordinate 1(N+k-2-cond)
  integer, parameter :: M = N2+k-4    !.. Number of B-splines in coordinate 2
  integer, parameter :: LM = L*M !.. Matrix dimension
  integer, parameter :: npl = N1+2*(k-1)  !.. Number of knot-points  N1+2(k-1)
  integer, parameter :: npm = n2+2*(k-1)  !.. Number of knot-points  N2+2(k-1)
 

  !.. Parameters for the knot-point grids tl and tm
  real(kind(1.d0)) :: tl(npl), tm(npm), tl_max, tm_max, tl_min, tm_min


  !.. Parameters for the harmonic trapping potential
   real(kind(1.d0)) :: angfreq, scaling
 
  !.. Parameters for effective potentials and coupling matrices
  integer, parameter :: points = 300
  real(kind(1.d0)), allocatable, dimension(:) :: rho_vector
  real(kind(1.d0)), allocatable, dimension(:,:) :: energy,S,TK,Srez,Hrez
  real(kind(1.d0)), allocatable, dimension(:,:,:) :: H,Hder,Pmat,P2mat

  !.. Parameters for plotting
  integer, parameter :: pp = 200
  real(kind(1.d0)) :: x(pp),y(pp),step_size,rho_min,rho_max
  real(kind(1.d0)) :: base(pp,LM),base_L(pp,L),base_M(pp,M)  

  !.. Other parameters
  real(kind(1.d0)) :: t1,t2,my,Vtrap(points),U(points)
  integer :: i,j,ll,mm,lj,li,mi,mj,mu,nu,n,test

  !.. Paramenters for generalized eigensolver
  integer                             :: ITYPE, LDA, LDB, INFO ,lwork,liwork
  character*1                         :: JOBZ = 'V', UPLO = 'U'
  real(kind(1.d0))    , dimension(LM) :: W
  real(kind(1.d0)) ,allocatable, dimension(:) :: WORK
  integer         , allocatable,dimension(:) :: IWORK

  ITYPE = 1
  LDA = LM
  LDB = LM
  LWORK=1+6*LM+2*LM*LM
  LIWORK = 3+5*LM

  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))


  allocate(rho_vector(points))
  allocate(energy(LM,points))
  allocate(H(LM,LM,points),Hder(LM,LM,points),Pmat(LM,LM,points),P2mat(LM,LM,points))
  allocate(S(LM,LM),TK(LM,LM),Srez(LM,LM),Hrez(LM,LM))

  !.. Declairing constants for model potential, trapping potential, and model atom
  
  my = mass(3)/sqrt(3.d0)
  !angfreq = 1.d0/(mass(1)*osc**2.d0)
  angfreq = 1.d0
  scaling = 10**(8.d0)
  

  !.. Setting up knot-vectors
  tl_min = 0.d0
  tl_max = Pi/2.d0
  tm_min = 0.d0
  tm_max = Pi/3.d0


  call universal_knot(npl,k,N1,tl_max,tl_min,tl)
  call universal_knot(npm,k,N2,tm_max,tm_min,tm)


  !.. Setting up hyperradial vector
  rho_min = 1.d0
  rho_max = 3000.d0
  step_size = (rho_max-rho_min)/(points-1)
  rho_vector(1) = rho_min
  print*, rho_vector(1)
  do i = 2,points
     rho_vector(i) = rho_vector(i-1)+step_size
  end do

  !.. Declairing rho for later use and creating the harmonic trapping potential "Vtrap"

  !Vtrap = 0.5d0*my*(angfreq**2.d0)*(rho_vector**2.d0)
  Vtrap = 0.d0

  !.. Calculating adiabatic potential curves and coefficients for the angular channel functions



  ! calcualte hamiltonian

  call efimovham(npl,npm,k,l,m,lm,points,tl,tm,rho_vector,energy,S,TK,H,Hder)
  call CPU_TIME( t1 )
  do i = 1, points
     Hrez = H(:,:,i)
     Srez = S
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

!  do n = 1, 3
!     do p = 1, 3
!        sumder(1,:) = 0.0d0
!        sumint(1,:) = 0.0d0
!        do j = 1, LM
!           do i = 1, LM
!              sumder(1,:) = sumder(1,:) + H(i,n,:)*H(j,p,:)*Hder(i,j,:)
!              sumint(1,:) = sumint(1,:) + H(i,n,:)*H(j,p,:)*S(i,j)
!           end do
!        end do
!        Pmat(n,p,:) = sumder(1,:)/(energy(n,:)-energy(p,:))
!     end do
!  end do

 
!  do i = 1, LM
!     Pmat(i,i,:) = 0.d0
!  end do
 

!  do p = 1, 3
!     do n = 1, 3
!        sumter(1,:) = 0.0d0
!        do i = 1, LM
!           sumter(1,:) = sumter(1,:) + Pmat(n,i,:)*Pmat(i,p,:)
!        end do
!        P2mat(n,p,:) = sumter(1,:)
!     end do
!  end do


  call CPU_TIME( t2 )
  print*,'creating coupling matrices P and P²', t2-t1




  !.. Writes adiabatic potential curves+trapping potential to file
  ! open(10,file='threebodypot.dat',status='replace')
!   do i = 1, points
!      write(10,10)i, rho_vector(i)/scatl, scaling*(energy(1,i)+Vtrap(i))/angfreq,scaling*(energy(1,i)-(P2mat(1,1,i)/(2.d0*my))+Vtrap(i))/angfreq, scaling*(energy(2,i)+Vtrap(i))/angfreq,scaling*(energy(2,i)-(P2mat(2,2,i)/(2.d0*my))+Vtrap(i))/angfreq

!      !, (energy(2,i)+Vtrap(i))/angfreq ,(energy(3,i)+Vtrap(i))/angfreq,(energy(4,i)+Vtrap(i))/angfreq,(energy(5,i)+Vtrap(i))/angfreq,(energy(6,i)+Vtrap(i))/angfreq, Vtrap(i)/angfreq
! 10   format(I4,'  ',16f20.8)
!   end do
!   close(10)

  !.. Writes adiabatic potential curves+trapping potential to file
  open(11,file='effectivepot.dat',status='replace')
  do i = 1, points
     write(11,10)i, rho_vector(i)/scatl, scaling*energy(1,i),scaling*energy(2,i),scaling*energy(3,i),scaling*energy(4,i),scaling*energy(5,i),scaling*energy(6,i),scaling*energy(7,i),scaling*energy(8,i)
  end do
  close(11)

  !.. Writes the lowest and second lowest lambda  as defined by Braaten & Hammer
  open(12,file='lambda.dat',status='replace')
  do i = 1, points
     write(12,10)i, rho_vector(i)/scatl, 2.d0*my*rho_vector(i)*rho_vector(i)*(energy(1,i)+0.25d0), 2.d0*my*rho_vector(i)*rho_vector(i)*(energy(2,i)+0.25d0)
  end do
  close(12)

  print*, 'this is working'

  !.. Writes adiabatic potential curves+trapping potential to file
  open(10,file='effectivepot0_413_N10n.dat',status='replace')
  do i = 1, points
     write(10,10)i, rho_vector(i), (energy(1,i)*2.d0*my*(rho_vector(i)**2.d0)+0.25d0),(energy(2,i)*2.d0*my*(rho_vector(i)**2.d0)+0.25d0), -(1.00624**2.d0)
     10   format(I3,'  ',16f20.8)
  end do
  close(10)

  print*, 'this is working'


  open(14,file='wave.dat',status='replace')
  do i = 1, points
     write(14,10)i,rho_vector(i)/scatl, Pmat(1,2,i), Pmat(2,1,i), P2mat(1,2,i), P2mat(1,1,i)
  end do
  close(14)

  deallocate(Pmat,P2mat)
  deallocate(Srez,Hrez,S,TK)
  deallocate(H,Hder)

end program efimov



   
