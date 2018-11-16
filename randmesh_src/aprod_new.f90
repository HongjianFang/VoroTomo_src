subroutine aprod(mode, m, n, x, y, leniw, lenrw, iw, rw)

  implicit none

  ! mode ==1: Compute  y = y + a*x
  ! mode ==2: Compute  x = x + a(transpose)*y
  integer*4 mode
  integer*4 m, n   ! Row and column dimensions of a
  real x(n), y(m)! Input vectors
  integer*8 :: leniw
  integer*8 lenrw
  ! iw[1]  Number of non-zero elements in a
  ! iw[2:iw[1]+1]  Row indices of non-zero elements
  ! iw[iw[1]+2:2*iw[1]+1]  Column indices
  integer*4 iw(leniw) ! integer*4 work vector containing:
  real rw(lenrw)! [1..iw[1]] Non-zero elements of a

  integer*8 i1
  integer*8 j1
  integer*8 k
  integer*8 kk

  !c	set the ranges the indices in vector iw

  kk=lenrw
  i1=0
  j1=kk

  if (mode.eq.1) then
    do k = 1,kk
      y(iw(i1+k)) = y(iw(i1+k)) + rw(k)*x(iw(j1+k))
    enddo
  else
    do k = 1,kk
      x(iw(j1+k)) = x(iw(j1+k)) + rw(k)*y(iw(i1+k))
    enddo
  endif
  end
