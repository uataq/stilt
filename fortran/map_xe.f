      subroutine map_xe(stcprm, x_map, xi, eta,flag)
! CHG(03/12/03)
!
! $Id: map_xe.f,v 1.1 2010/10/26 19:03:40 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      parameter (pi=3.14159265358979323846,pi_2=pi/2.)
      character * 1 flag
      dimension stcprm(15),x_map(3)
        if (dabs(x_map(3)) .ge. 1.) then
!* Projection pole or its antipodes
            xi = 0.
          if (stcprm(1) * x_map(3) .gt. 0.) then
            eta = 1./stcprm(1)
          else
            if (x_map(3) .gt. 0) then
              eta = 1.0e10
            else
              eta = -1.0e10
            endif
          endif
          return
        else
          if (dabs(stcprm(1)) .eq. 1.) then
!* Stereographic Case, away from pole
            fact = 1. / (stcprm(1) + x_map(3) )
            xi = x_map(2)*fact
            eta = 1./stcprm(1) - x_map(1) * fact
          else
            ymerc = .5 * dlog( ( 1. + x_map(3) ) / (1. - x_map(3) ) )
!* This Projection has a cut.  Select theta according to the rule
!* If cutflag = 'e', -PI/2 <= theta < 3PI/2, if cutflag = 'E',
!* 0 <= theta < 2PI, if cutflag = 'w', -3Pi/2 <= theta < PI/2,
!* cutflag = 'W', -2PI <= theta < 0., else -PI<=theta<PI.
          if (flag .eq. 'E') then
            theta = datan2(- x_map(2), - x_map(1) ) + pi
          else if (flag .eq. 'W') then
            theta = datan2(- x_map(2), - x_map(1) ) - pi
          else if (flag .eq. 'e') then
            theta = datan2(- x_map(1), x_map(2) )  + pi_2
          else if (flag .eq. 'w') then
            theta = datan2(x_map(1), - x_map(2) )  - pi_2
          else
            theta = datan2(x_map(2), x_map(1) )
          endif
            rhog = xpabva( stcprm(1), - ymerc)
!* rhog = ( exp( - gamma * ymerc ) - 1. ) / gamma */
            xi = (1. + stcprm(1) * rhog) * snabva(stcprm(1), theta)
            eta = stcprm(1) * (1. + stcprm(1) * rhog) *                 &
     &                        csabva(stcprm(1), theta) - rhog
          end if
        end if
        return
      END
