      subroutine cpolll_stilt(stcprm, alat,along, enx,eny,enz)
! returns a vector aligned with the direction toward the North Pole.  I.e.
! parallel to the vector from the Earth's center to the North Pole.
! CHG(03/12/03)
!
! $Id: cpolll_stilt.f,v 1.1 2010/10/26 19:14:23 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL :: alat,along

      real*8 stcprm(15)
      real*8 map(3),pole(3),geog(3)
      call ll_geo(dble(alat),dble(along), geog)
      call basg2m(stcprm, geog, map)
      do k=1,3
        pole(k) = stcprm(3*k + 1)
      enddo
      call proj_3d(stcprm, map, pole, enx,eny,enz)
      return
      END
