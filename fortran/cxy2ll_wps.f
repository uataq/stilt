
subroutine cxy2ll_wps (STCPRM, X,Y, XLAT,XLONG, proj)

!
! $Id: cxy2ll_wps.f,v 1.4 2010/10/26 19:10:52 jel Exp $
!
  use map_utils
  use module_defgrid, only : vmiss, vmissle

  IMPLICIT none
  type(proj_info), intent(in) :: proj
  real*8 :: stcprm(15),xlatD,xlongD
  real   :: x, y, xlat, xlong
  real*8 :: map(3),geog(3)
  if (all(stcprm .le. vmissle)) then
     ! use WPS mapping routines
     call ij_to_latlon(proj, x, y, xlat, xlong)
  else
      call xy_map(stcprm, dble(x),dble(y), map)
      call basm2g(stcprm, map, geog)
      call geo_ll(geog, xlatD,xlongD)
      xlat=xlatD
      xlong=xlongD
  end if
  return
end subroutine cxy2ll_wps
