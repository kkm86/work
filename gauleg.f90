 subroutine  gauleg(ngp, xabsc, weig)

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
    
  end subroutine gauleg
