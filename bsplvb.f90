subroutine bsplvb(t,jhigh,index,x,left,biatx)

  parameter (JMAX=100)
  integer index,jhigh,left,i,j,jp1
  real(kind(1.d0)) t,x,biatx,deltal,deltar,saved,term
  dimension biatx(jhigh),t(left+jhigh),deltal(jmax),deltar(jmax)
  save deltal,deltar

  DATA j/1/

  GO TO (10,20),index
10 j = 1
  biatx(1) = 1.d0
  if (j .ge. jhigh) go to 99

20 continue
  jp1 = j + 1
  deltar(j) = t(left+j) - x
  deltal(j) = x - t(left+1-j)

  saved = 0.d0

  do i = 1,j
     term = biatx(i)/(deltar(i) + deltal(jp1-i))
     biatx(i) = saved + deltar(i)*term
     saved = deltal(jp1-i)*term
  end do
  biatx(jp1) = saved
  j = jp1
  if (j .lt. jhigh) go to 20
99 return
end subroutine bsplvb
