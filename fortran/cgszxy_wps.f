      real function CGSZXY_wps (stcprm, XS,YS, grid)
      use map_utils
      use module_defgrid, only : gset, vmiss, vmissle
      implicit none
      type(gset), intent(in) :: grid
      real*8, intent(in) :: stcprm(15)
      real, intent(in)   :: xS,yS
      
      real*8 :: map(3),ymerc,x,y

      integer :: i_mapfactor, j_mapfactor
      x=xS
      y=yS

      if (all(stcprm .le. vmissle)) then
     !     use WPS mapping routines
         i_mapfactor = nint(x)
         i_mapfactor = min(grid%nx,max(1,i_mapfactor))
         j_mapfactor = nint(y)
         j_mapfactor = min(grid%ny,max(1,j_mapfactor))
         cgszxy_wps = grid%proj%dx * grid%mapfactor(i_mapfactor,j_mapfactor) / 1000.
      else
         call xy_map(stcprm, x,y,map)
         if (map(3) .ge. 1.)  then
            if (stcprm(1) .ge. 1.) then
               cgszxy_wps = 2.*stcprm(15)
            else
               cgszxy_wps = 0.
            endif
         else if (map(3) .le. -1.) then
            if (stcprm(1) .le. -1.) then
               cgszxy_wps = 2.*stcprm(15)
            else
               cgszxy_wps = 0.
            endif
         else if (dabs(stcprm(1)) .ge. 1.) then
            cgszxy_wps = stcprm(15) * (1. + dsign(map(3),stcprm(1)))
         else
            ymerc = -.5 * dlog( (1. - map(3))/(1. + map(3)) )
            cgszxy_wps = stcprm(15) * dexp( - (1. - stcprm(1)) * ymerc) *       &
            (1. + map(3))
         endif
      end if
      return
      end function CGSZXY_wps
