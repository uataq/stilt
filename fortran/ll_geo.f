      subroutine ll_geo(xlat, xlong, vector)
!  Given a latitude xlat and longitude xlong, returns the unit vector
!  directed to the given point in the geo (geographic) coordinate system.
!
!  The geo system is a 3-D Cartesian coordinate system, with origin
!  at the Center of the Earth, with the x_3 axis pointing to the North
!  Pole, the x_1 axis pointing to the Equator at the Greenwich Meridian,
!  and the x_2 axis by the right hand rule pointing to the Equator
!  at 90 E.
! CHG(03/12/03)
!
! $Id: ll_geo.f,v 1.1 2010/10/26 19:03:42 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      parameter (pi=3.14159265358979323846,radpdg=pi/180.)
      parameter (dgprad=180./pi)
      dimension vector(3)
        vector(3) = dsin(RADPDG * xlat)
        clat = dcos(RADPDG * xlat)
        vector(1) = clat * dcos(RADPDG * xlong)
        vector(2) = clat * dsin(RADPDG * xlong)
      return
      END
