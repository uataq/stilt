subroutine permute(ans, nax, nay, k, nkx, nky, len, lai, loi, foot)
  implicit none

  integer, intent(in)           :: lai(len), loi(len), nax, nay, nkx, nky, len
  double precision, intent(in)  :: k(nkx,nky), foot(len)
  integer                       :: lox(nkx), loy(nky)
  integer                       :: i, j, n, xs, ys
  double precision              :: ks
  double precision, intent(out) :: ans(nax, nay)

  ! Generate location index array
  lox = (/(i, i=1,nkx)/) - ((nkx + 1) / 2)
  loy = (/(i, i=1,nky)/) - ((nky + 1) / 2)

  ! Matrix permutations corresponding to gaussian kernel size
  do i = 1, nkx
    do j = 1, nky
      xs = lox(i)
      ys = loy(j)
      ks = k(i, j)

      ! Add all footprint values to their locations in arrays
      do n = 1, len
        ans(loi(n)+xs, lai(n)+ys) = ans(loi(n)+xs, lai(n)+ys) + foot(n)*ks
      end do
    end do
  end do
end subroutine permute
