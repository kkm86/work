subroutine coupling(H,Hder,Integ,energy,LM,points,cc,Pmat,P2mat,Imat)

  implicit none

  !.. Input
  integer, intent(in) :: LM     !.. Size of Coefficient matrices
  integer, intent(in) :: points !.. Number of hyperradial points 
  integer, intent(in) :: cc     !.. Number of couplings to be calculated

  !.. Eigenvector and differential coefficents and eigenvalue arrays
  real(kind(1.d0)), dimension(LM,LM,points),intent(in) :: H, Integ, Hder
  real(kind(1.d0)), dimension(LM,points),   intent(in) :: energy

  !.. Output
  real(kind(1.d0)), dimension(cc,cc,points),intent(inout) :: Pmat,P2mat,Imat

  !.. Local
  real(kind(1.d0)), dimension(points) :: sumter, sumder
  real(kind(1.d0)), allocatable, dimension(:,:) :: H_i, I_i
  real(kind(1.d0)) :: t1, t2
  integer :: i,j,n,p,map
  
  
  !.. Setting up coupling matrices
  call CPU_TIME( t1 )
  
  Pmat = 0.0d0
  Imat = 0.0d0

  map = cc*cc

  allocate(H_i(map,points),I_i(map,points))

  do p = 1, cc
     do n = 1, cc
        mm = (p-1)*cc+n
        sumter(1,:) = 0.d0
        sumder(1,:) = 0.d0
        do j = 1, LM
           do i = 1, LM
              sumter = sumter + H(i,n,:)*H(j,p,:)*Hder(i,j,:)
              sumder = sumder + H(i,n,:)*H(j,p,:)*Integ(i,j,:)
           end do
        end do
        H_i(mm,:) = sumter
        I_i(mm,:) = sumder
     end do
  end do

  do mm = 1, map
     n = mod((mm-1),LM) + 1
     p = (mm-1)/LM + 1
     Pmat(n,p,:) = H_ii(mm,:)/(energy(n,:)-energy(p,:))
     Imat(n,p,:) = I_ii(mm,:)
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
end subroutine coupling

