      subroutine cgrnll(stcprm, alat,along, enx,eny,enz)
! returns a vector aligned in the direction toward the Greenwich Meridian
! at the equator.  I.e. parallel to the vector from Earth's center to
! that point.
! CHG(03/12/03)
!
! $Id: cgrnll.f,v 1.1 2010/10/26 19:03:39 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
      real*8 map(3),pole(3),geog(3)
      call ll_geo(alat,along, geog)
      call basg2m(stcprm, geog, map)
      do k=1,3
        pole(k) = stcprm(3*k - 1)
      enddo
      call proj_3d(stcprm, map, pole, enx,eny,enz)
      return
      END
