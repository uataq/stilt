      real*8 function cgszll_stilt(stcprm, alat, along)
!
! $Id: cgszll_stilt.f,v 1.1 2010/10/26 19:14:22 jel Exp $
!
      IMPLICIT NONE

      real*8 stcprm(15), alat, along
      real*8 map(3),geog(3),ymerc

      call ll_geo(alat, along, geog)
      call basg2m(stcprm, geog, map)
      if (map(3) .ge. 1.) then
         if (stcprm(1) .ge. 1.) then
            cgszll_stilt = 2. * stcprm(15)
         else
            cgszll_stilt = 0.
         endif
      else if (map(3) .le. -1.) then
         if (stcprm(1) .le. -1.) then
            cgszll_stilt = 2. * stcprm(15)
         else
            cgszll_stilt = 0.
         endif
      else if (dabs(stcprm(1)) .ge. 1.) then
         cgszll_stilt = stcprm(15) * (1. + dsign(map(3),stcprm(1)))
      else
         ymerc = -.5 * dlog( (1. - map(3))/(1. + map(3)) )
         cgszll_stilt = stcprm(15) * dexp( - (1.-stcprm(1)) * ymerc) *         &
         (1. + map(3))
      endif

      return
      END
