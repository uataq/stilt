      subroutine xy_map(stcprm, x,y, map)
! CHG(03/12/03)
!
! $Id: xy_map.f,v 1.1 2010/10/26 19:03:42 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
      real*8 map(3)
      real*8, EXTERNAL :: lnabva, atnabv
      PARAMETER (REARTH=6367.47)
        xi0 = (x - stcprm(11)) * stcprm(15) / REARTH
        eta0 = (y - stcprm(12)) * stcprm(15) / REARTH
        xi = xi0 * stcprm(13) - eta0 * stcprm(14)
        eta = xi0 * stcprm(14) + eta0 * stcprm(13)
        vsq = stcprm(1) * (xi*xi + eta*eta) - 2. * eta
        if (stcprm(1) * vsq + 1. .le. 3.e-8) then
!/* Case of point at Projection Pole */
          map(1) = 0.
          map(2) = 0.
          TEMP=1.
          map(3) = dsign(TEMP,stcprm(1))
        else
          if ( stcprm(1) .ge. 1.) then
!/* Stereographic Case, std */
            fact = 2. / (2. + vsq)
            map(2) = xi * fact
            map(1) = (1. - eta) * fact
            map(3) = fact - 1.
          else if (stcprm(1) .le. -1.) then
!/* Stereographic Case, centered at antipodes */
            fact = -2. / (2. - vsq)
            map(2) = xi * fact
            map(1) = (- 1. - eta) * fact
            map(3) = fact + 1.
          else
!/* Mercator or Lambert Case */
            ymerc = -.5 * lnabva(stcprm(1), vsq)
            if ( ymerc .gt. 0.) then
              fact = dexp(-ymerc)
              cosphi = fact / (1. + fact * fact)
              cosphi = cosphi + cosphi
              map(3) = 1. - fact * cosphi
            else
              fact = dexp(ymerc)
              cosphi = fact / (1. + fact * fact)
              cosphi = cosphi + cosphi
              map(3) = fact * cosphi - 1.
            endif
!/*      map.v[2] = 1. - fact * cosphi; */
            if (dabs(map(3)) .lt. 1.) then
              theta = atnabv(stcprm(1), xi, eta)
              map(1) = cosphi * dcos(theta)
              map(2) = cosphi * dsin(theta)
            else
              map(1) = 0.
              map(2) = 0.
            endif
          endif
        endif
      return
      END
