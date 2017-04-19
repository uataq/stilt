      subroutine stcm1p_stilt(grid) 

! Replace local variables:
!       stcprm, x1, y1, xlat1, xlong1,                  &
!     & xlatrf, xlonrf, gridsz, orient, proj)
!
! by components of gset derived type:
!      GRID(KG)%GBASE, GRID(KG)%SYNC_XP,  GRID(KG)%SYNC_YP,               &
!      GRID(KG)%SYNC_LAT, GRID(KG)%SYNC_LON,              &
!      GRID(KG)%REF_LAT,  GRID(KG)%REF_LON,               &
!      GRID(KG)%SIZE,     GRID(KG)%ORIENT,                &
!      GRID(KG)%PROJ)

!
! CHG(03/12/03)
!
! $Id: stcm1p_stilt.f,v 1.1 2010/10/26 19:14:23 jel Exp $
!
      use module_defgrid, only : gset, vmiss, vmissle

      IMPLICIT REAL*8 (A-H,O-Z)

      type(gset), intent(inout) :: grid

      parameter (pi=3.14159265358979323846d0,radpdg=pi/180d0)
      parameter (dgprad=180d0/pi)
!      real*8 stcprm(15)
      real*8 enx,eny,enz,norm
      real x1a,y1a
      REAL(KIND(1d0)), EXTERNAL :: cgszll_stilt

      c_or = dcos(RADPDG * grid%orient)
      s_or = - dsin(RADPDG * grid%orient)
!  stcprm->x0 = stcprm->y0 = stcprm->srotate = 0;
      grid%gbase(11) = 0.
      grid%gbase(12) = 0.
      grid%gbase(13) = 1.
      grid%gbase(14) = 0.
      grid%gbase(15) = 1.
!  stcprm->crotate = stcprm -> gridszeq = 1.0;
      call cpolll_stilt(grid%gbase, grid%ref_lat, grid%ref_lon, enx, eny, enz)
      norm = dsqrt(enx*enx + eny*eny)
      if (norm .eq. 0.) then
        call cgrnll(grid%gbase,dble(grid%ref_lat),dble(grid%ref_lon),enx,eny,enz)
        norm = dsqrt (enx* enx + eny*eny)
      endif
      enx = enx/norm
      eny = eny/norm
!  stcprm->gridszeq *= gridsz / cgszll(stcprm, xlatrf,xlonrf);
      grid%gbase(15) = grid%gbase(15) * grid%size / cgszll_stilt(grid%gbase,dble(grid%ref_lat),dble(grid%ref_lon))

!  stcprm -> crotate = eny * c_or - enx * s_or;
      grid%gbase(13) = eny * c_or - enx * s_or
!  stcprm -> srotate = -eny * s_or - enx * c_or;
      grid%gbase(14) = -eny * s_or - enx * c_or
      call cll2xy_wps(grid%gbase, grid%sync_lat,grid%sync_lon, x1a,y1a, grid%proj)

!  stcprm->x0 += x1 - x1a;
!  stcprm->y0 += y1 - y1a;
      grid%gbase(11) = grid%gbase(11) + grid%sync_xp - x1a
      grid%gbase(12) = grid%gbase(12) + grid%sync_yp - y1a
      return
      END
