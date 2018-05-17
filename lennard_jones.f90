subroutine lennard_jones(S,b,mass,rho,alpha,V)
  
  implicit none

  !.. Input
  real(kind(1.d0)), intent(in) :: S,b,alpha,rho,mass
  real(kind(1.d0)), intent(out) :: V

  !.. Local
  real(kind(1.0d0)) :: r

   r = sqrt(2.d0)*rho*sin(alpha)
   !r = rho

   if(rho .eq. 0) then
      V = 0.d0
   else
      V = 2.d0*mass*S*((b/r)**12.d0-2.d0*(b/r)**6.d0)*rho**2.d0
   end if
   
   return

end subroutine lennard_jones
  
