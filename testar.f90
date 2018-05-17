program dot_main
  real*4 x(1), y(1), sdot, res,sasum
  real(kind(1.d0)) ddot
integer  n, incx, incy, i
!external sdot
n = 1
incx = 1
incy = 1
do i = 1, n
  x(i) = 2.0e0
  y(i) = 1.0e0
end do
res = sdot(n, x, incx, y, incy)
print*, 'SDOT = ', res
!print*, 'SDOT = ', sdot(10, x, 1, y, 1)
print*, 'SDOT = ', n, x, incx, y, incy
print*, sasum(n,x,incx)
end program dot_main
