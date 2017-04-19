      subroutine xe_xy(stcprm, xi, eta, x, y)
! CHG(03/12/03)
!
! $Id: xe_xy.f,v 1.1 2010/10/26 19:03:41 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (REARTH=6367.47)
      dimension stcprm(15)
        x = stcprm(11) + REARTH / stcprm(15) *                          &
     &  (stcprm(13) * xi + stcprm(14) * eta)
        y = stcprm(12) + REARTH / stcprm(15) *                          &
     &  (stcprm(13) * eta - stcprm(14) * xi)
        return
      END
