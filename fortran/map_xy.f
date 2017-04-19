      subroutine map_xy(stcprm, x_map, x, y)
! CHG(03/12/03)
!
! $Id: map_xy.f,v 1.1 2010/10/26 19:03:40 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      dimension stcprm(15),x_map(3)
        call map_xe(stcprm, x_map, xi, eta, 'c')
        call xe_xy(stcprm, xi, eta, x, y)
        return
      END
