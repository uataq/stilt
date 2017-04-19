!  Conversions to and from the geo (geographic) and the map
!  (Map-oriented)  coordinate systems.
!
!  The geo system is a 3-D Cartesian coordinate system, with origin
!  at the Center of the Earth, with the x_3 axis pointing to the North
!  Pole, the x_1 axis pointing to the Equator at the Greenwich Meridian,
!  and the x_2 axis by the right hand rule pointing to the Equator
!  at 90 E.
!
!  In the map system, the axis of the map projection passes through the
!  center of the Earth, and through a point on the surface that acts as the
!  "Pole" of the map (Which will coincide with the North or South Pole,
!  unless the projection is Oblique or Transverse).  The x_3 axis of the
!  map coordinate system is aligned with this map pole.  In Lambert and
!  Mercator projections, a "Cut" extends from the map pole along a great
!  circle to its antipode; the unrolled map is discontinuous there.  The
!  x_1 axis of the map coordinate system is diametrically opposite this
!  cut, and 90 degrees from the pole; the x_2 axis is selected to complete
!  a right hand rule.
!
!  The coefficients of the map coordinate system relative to the geo system
!  are given in elements 2 through 10 of the stcprm array.
!
      subroutine basg2m(stcprm, xin, xout)
!  receives the vector xin, giving the components of a point in the geo
!  system; returns xout, the components of the same point in the map system.
! CHG(03/12/03)
!
! $Id: basg2m.f,v 1.1 2010/10/26 19:03:38 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      dimension stcprm(15),xin(3),xout(3)
      ifind(l,k) = -2 + 3*l + k
        do k=1,3
          xout(k) = 0.
        enddo
        do l=1,3
          do k=1,3
              xout(l) = xout(l) + xin(k) * stcprm(ifind(l,k))
            enddo
        enddo
        return
      END
