subroutine baseR(x,xmax,kord,np1,t1,Bord_pot,Bder_pot,Bder2_pot)

  implicit none
  
  !.. Input
  integer         , intent(in)    :: kord, np1, xmax
  real(kind(1.d0)), intent(in)    :: x, t1(np1)
  real(kind(1.d0)), intent(inout) :: Bord_pot(xmax+1),Bder_pot(xmax+1),Bder2_pot(xmax+1)

  !.. Local
  integer :: jj

  !.. External functions
  real(kind(1.d0)) :: bget,bder,bder2

  do jj = 1, xmax+1
     if (jj .le. xmax) then
        Bord_pot(jj) = bget(x,t1,kord,np1,jj+1)
        Bder_pot(jj) = bder(x,t1,kord,np1,jj+1)
        Bder2_pot(jj) = bder2(x,t1,kord,np1,jj+1)
     else
        Bord_pot(jj) = Bord_pot(jj-1)
        Bder_pot(jj) = Bder_pot(jj-1)
        Bder2_pot(jj) = Bder2_pot(jj-1)
     end if
  end do
  return

end subroutine baseR





  
