      real*8 function csabva(a, b)
! CHG(03/12/03)
!
! $Id: limmath.f,v 1.1 2010/10/26 19:03:40 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)


!* returns (1. - dcos(a * b) ) / a / a, or its limit as a -> 0
        term = snabva(.5 * a, b)
        csabva = .5 * term * term
        return
      END

      real*8 function lnabva(a, b)

      real*8 a,b,c,t
!* returns ln (1. + a * b) / a, or its limit as a -> 0.
        c = a * b
        if (dabs(c) .gt. 0.01) then
          lnabva = log( 1. + c ) / a
        else
!* udsing (1. + t) / (1. - t) = 1. + c when t = c / (2. + c)
          t = c/(2. + c)
          t = t * t
          lnabva =  b / (2. + c) * (2.    + t *                         &
     &                             (2./3. + t *                         &
     &                             (2./5. + t *                         &
     &                             (2./7. ))))
        endif
        return
      END

      real*8 function snabva(a, b)
!* returns dsin(a * b) / a, or its limit as a -> 0.

      real*8 a,b,c,csq
        c = a * b
        csq = c * c
        if (csq .gt. .001) then
          snabva = dsin( c )  / a
        else
          snabva = b * ( 1. - csq / 6. *                                &
     &                 ( 1. - csq / 20. *                               &
     &                 ( 1. - csq / 42. )))
        endif
        return
      END

!  A series of functions which may have problems in certain ranges
!  of one of their parameters.  These functions, if written normally,
!  may suffer round-off problems.  They are used in particular for
!  Lambert Conformal projections whose parameters make them almost
!  the same as a Mercator Projection.

      real*8 function xpabva(a, b)
!* returns (exp (a * b) - 1. ) / a, or its limit as a -> 0.
      real*8 a,b,c,csq
        c = a * b
        csq = .25 * c * c
        if ( csq  .gt. .001) then
          xpabva = (dexp( c ) - 1.) / a
        else
          xpabva = dexp(.5 * c) * b * (1. + csq / 6. *                  &
     &                               (1. + csq / 20. *                  &
     &                               (1. + csq / 42. )))
        endif
        return
      END

      real*8 function atnabv(a, b, c)
!* returns atan2(a*b,1-a*c)/a, or its limit as a -> 0
      real*8 a,b,c,xi,eta,vsq,t

        xi = a * b
        eta = a * c
        vsq = xi * xi + eta * eta
        if ( vsq - 2. * eta + 1. .le. 0.) then
          atnabv = 0.
          return
        endif
        if (dabs(xi) .gt. .01 * (1. - eta) ) then
          atnabv = datan2(xi, 1. - eta) / a
        else
          t = xi/(1. - eta)
          t = t * t
          atnabv = b / (1. - eta) * (1.    - t *                        &
     &                              (1./3. - t *                        &
     &                              (1./5. - t *                        &
     &                              (1./7. ))))
        endif
        return
      END
