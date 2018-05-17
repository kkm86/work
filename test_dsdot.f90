program dot_main2
  real(kind(1.d0)) :: x(10), y(10), ddot, res
  integer  n, incx, incy, i

  n = 10
  incx = 1
  incy = 1
  do i = 1, n
     x(i) = 2.0
     y(i) = 1.0
  end do
  res = ddot(n, x, incx, y, incy)
  print*, 'DDOT = ', res
  print*, 'SDOT = ', n, x, incx, y, incy

end program dot_main2
