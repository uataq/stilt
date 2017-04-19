      FUNCTION x_prod(A,B,C)
!  Returns as C the vector cross product of A and B:
!  C = A x B.  Also returns as function value the absolute value of C.
! CHG(03/12/03)
!
! $Id: x_prod.f,v 1.1 2010/10/26 19:03:42 jel Exp $
!
      IMPLICIT NONE

      REAL(KIND(1d0)), INTENT(IN)  :: a(3),b(3)
      REAL(KIND(1d0)), INTENT(OUT) :: C(3)
      REAL(KIND(1d0))              :: x_prod

      INTEGER, PARAMETER :: indx(4) = (/2,3,1,2/)
      INTEGER            :: k
      REAL(KIND(1d0))    :: D

      D = 0.
      do k=1,3
        C(k) = A(indx(k))*B(indx(k+1)) - A(indx(k+1))*B(indx(k))
        D = D + C(k) * C(k)
      enddo
      x_prod = SQRT(D)

      END FUNCTION x_prod
