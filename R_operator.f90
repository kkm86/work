subroutine R_operator(np2,k,NN2,alpha_1,m,ii,term)

  use constants

  implicit none

  !.. Input
  integer            , intent(in) :: np2,k,NN2          !Number of knot-points etc.
  integer            , intent(in) :: m,ii
  real(kind(1.d0))   , intent(in) :: alpha_1

  !.. Output
  real(kind(1.d0))   , intent(inout) :: term

  !.. Local
    real(kind(1.d0)) :: alphaprim(m,k) !Abscissa-grid
    real(kind(1.d0)) :: xabsc(k),weig(k), innersum
    integer          :: n

  !.. External functions/variables
    real(kind(1.d0)) :: bget, ta(np2)

    call gauleg(k, xabsc, weig)
    call alpha_knot(np2,k,NN2,alpha_1,ta)

    innersum = 0.0
    do n = 1, k
       if(ta(m+1).ne.ta(m)) then
          alphaprim(m,n) = 0.5*(ta(m+1)+ta(m)) + 0.5*(ta(m+1)-ta(m))*xabsc(n)
          term = term + weig(n)*4.d0*bget(alphaprim(m,n),ta,k,np2,ii+1)/sqrt(3.d0)
       else
          term = 0.d0
       end if
       innersum = innersum + 0.5*(ta(m+1)-ta(m))*term
    end do

    innersum = term

    return

  end subroutine R_operator
  
   
