      subroutine geo_ll(geog,plat,plon)
!  given a vector geog in the geo (geographic) coordinate system,
!  returns as plat and plon the latitude and longitude of the
!  corresponding point.
!
!  The geo system is a 3-D Cartesian coordinate system, with origin
!  at the Center of the Earth, with the x_3 axis pointing to the North
!  Pole, the x_1 axis pointing to the Equator at the Greenwich Meridian,
!  and the x_2 axis by the right hand rule pointing to the Equator
!  at 90 E.
! CHG(03/12/03)
!
! $Id: geo_ll.f,v 1.1 2010/10/26 19:03:39 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      parameter (pi=3.14159265358979,DGPRAD=180./pi,RADPDG=pi/180.)
      dimension geog(3)
      fact = geog(1)*geog(1) + geog(2)*geog(2)
      if (fact .le. 0.) then
!        plon = 0.
        TEMP=90.
        plat = dsign(TEMP,geog(3))
        plon = 90. + plat
! Change made 02/12/02 to acommodate WMO reporting conventions.  North
! pole is longitude 180., so "North" points to the Greenwich Meridian,
! South Pole is longitude 0. so "North" again points to the Greenwich
! Meridian.
      else
        fact = dsqrt(fact)
        plat = DGPRAD * datan2(geog(3),fact)
        plon = DGPRAD * datan2(geog(2),geog(1))
      endif
      return
      END
