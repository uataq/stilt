!###############################################################################
! GEMGRD - Creates the latitude/longitude dependent grid dimensions for
! computing the finite difference fluxes and converting cell mass to air 
! concentration. In the current global lat-lon grid configuration the 1,1 grid
! point is at the southwest corner and the maximum points are at the northwest
! corner.  Lat-lon grids can start at either 000 or 180E longitude.
!-------------------------------------------------------------------------------
! GRID COORDINATE SYSTEM
! 
!     *      *      *        where * are the meteorlogical grid points
!                            of dimensions nx,ny     
!         +-----+
!         |     |            and the box represents the concentration
!     *   |  *  |   *        grid cell whose dimensions are computed 
!         |     |            between the + verticies
!         +-----+
!
!     *      *      *
! ------------------------------------------------------------------------------
! LAST REVISED: 16 May 2008 (RRD) - initial version
!               26 Sep 2008 (MDC) - added clat1 to grid size equation (gxy)
!-------------------------------------------------------------------------------

SUBROUTINE gemgrd

  USE funits
  USE gemkon  
  USE gemvar  
  USE gemcfg  

  IMPLICIT NONE

! establish the default computational horzontal grid spacing for latitudes 70-90
! (near the poles) in terms of the number of meteorological grid points comprise
! each concentration grid cell for both the 1-deg and 2.5 deg grids

  INTEGER*4 :: nh10(181),nh25(73)
  DATA nh10 /360,2*180,4*90,2*45,5*15,6*3,141*1,6*3,5*15,2*45,4*90,2*180,360/
  DATA nh25 /144,72,36,18,6,6,2,2,57*1,2,2,6,6,18,36,72,144/                          
                          
  INTEGER*4 :: j
  REAL*4    :: clat,base,pole,abot,atop

  INTEGER*4 :: nx,ny,nz,np,kgrd
  REAL*4    :: clat1,clon1,dlat,dlon

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  IF(KINIT.LT.0)RETURN

! establish the horizontal grid aggregation factors to limit the number
! of computational points near the poles
  
  IF(ny.EQ.181)THEN
     ngp = nh10
  ELSEIF(ny.EQ.73)THEN
     ngp = nh25
  ELSE
     ngp = 1
     WRITE(KF21,*)'*ERROR* gemgrd: grid spacing aggregation not defined - ',ny
     STOP 800
  END IF

! Distance on the earth's surface is given by the earth's circumference
! converted to distance per degree times the cell spacing. The grid spacing
! is used in the finite difference equations and therefore should apply
! at the meteorological grid points.

  base = 2.0*pi*rearth/360.0 ! base distance on earth in meters/deg-latitude

  gsy  = dlat*base           ! spacing (meters) same all latitudes

! Distance in the horizontal uses a similar approach by computing the 
! circumference at each latitude. Horizontal distance = 0 at the poles.

  gsx = 0.0

  DO j=2,(ny-1)   
     clat = (j-1)*dlat+clat1
     gsx(:,j) = dlon*base*cos(clat*radpdeg)  ! the same for all longitudes
  END DO

! The area of the cell surrounding each meteorological grid point is computed
! from the difference in the sector areas as measured from the pole to the 
! the base and top of each cell. The sector area is defined as (pi P**2),
! where P is the linear distance from the pole to the latitude. Northern
! hemisphere distances are the same as southern hemisphere distances.

  DO j=2,ny/2
     clat = (j-1)*dlat-(0.5*dlat)+clat1   ! bot of box lat
     base=rearth*cos(clat*radpdeg)
     pole=rearth*(1.0-sin(clat*radpdeg))
     abot=pi*(base*base+pole*pole)          

     clat = j*dlat-(0.5*dlat)+clat1       ! top of box lat
     base=rearth*cos(clat*radpdeg)
     pole=rearth*(1.0-sin(clat*radpdeg))
     atop=pi*(base*base+pole*pole)          

     gxy(:,j)      = abs(atop-abot)/nx    ! area (sq-meters) per grid cell
     gxy(:,ny-j+1) = abs(atop-abot)/nx
  END DO

! special case at the equator

  gxy(:,1+ny/2) = gsx(:,1+ny/2) * gsy(:,1+ny/2) 

! The area of the polar grid cell (there is only one, but all are the same in
! the array) is computed from the radius using the above equations.

  clat=90.0-0.5*dlat
  base=rearth*cos(clat*radpdeg)
  pole=rearth*(1.0-sin(clat*radpdeg))
  gxy(:,1)  = pi*(base*base+pole*pole)          
  gxy(:,ny) = pi*(base*base+pole*pole)          

END SUBROUTINE gemgrd
