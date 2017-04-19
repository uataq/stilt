 subroutine cll2xy_wps (STCPRM, XLAT,XLONG, X,Y, proj)
  ! $Id: cll2xy_wps.f,v 1.3 2010/10/26 19:10:52 jel Exp $
  !
  use map_utils
  use module_defgrid, only : vmiss, vmissle

  IMPLICIT none
  type(proj_info), intent(in) :: proj
  real*8   :: stcprm(15),xD,yD
  real     :: xlat, xlong,x,y
  real*8   :: geog(3),temp(3)
  if (all(stcprm .le. vmissle)) then
     ! use WPS mapping routines
     call latlon_to_ij(proj, xlat, xlong, x, y)
  else
      call ll_geo(dble(xlat),dble(xlong), geog)
      call basg2m(stcprm, geog, temp)
      call map_xy(stcprm, temp, xD,yD)
      x=xD
      y=yD
  end if
  return
end subroutine cll2xy_wps
