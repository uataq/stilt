      subroutine mpstrt(stcprm, conang, p_lat, p_long, r_lat, r_long)
! CHG(03/12/03)
!
! $Id: mpstrt.f,v 1.1 2010/10/26 19:03:40 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      parameter (pi=3.14159265358979,RADPDG=pi/180.,DGPRAD=180./pi)

!*  REARTH=6356.766  from U.S. Standard Atmosphere, 1976
!*  REARTH=6367.47   for spherical earths in GRIB grids
!*  REARTH=6371.2    original assumption for CMAPF routines.
!*                   source lost, probably old NWPC grids.

! CHG(09/15/03) adapt to saulos version
!      PARAMETER (REARTH=6367.47)
      PARAMETER (REARTH=6367.00)
      dimension stcprm(15)
!typedef struct {
! 1:     gamma,
! 2-10:  rotate(3,3),
! 11-12: x0,y0,
! 13-14: crotate,srotate,
!  15:   grdszq
!

!*
!*  General Purpose routine to Set Map Parameters.  Called by the
!*  special purpose routines for oblique stereographic, oblique and
!*  transverse Mercator, oblique Lambert Conformal, etc. Projections
!*  Inputs: p_lat,p_long - Latitude and Longitude of the Pole
!*            Point for the Projection
!*          r_lat,r_long - Latitude and Longitude of a
!*            reference point - 180degrees around the Pole Point
!*            from the cut
!*          conang - angle between the Projection Pole axis and
!*            the generators (sides) of the cone on which the Earth
!*            is projected. + or - 90 degrees indicates a plane,
!*            hence Stereographic; 0 degrees indicates a cylinder
!*            and hence Mercator.
!*  Outputs: stcprm - map parameters
!*/
      dimension temp(3)
      REAL(KIND(1d0)), EXTERNAL :: x_prod
      
      ifind(l,k) = -2 + 3*l + k
        call ll_geo(p_lat,p_long,temp)
        do k=1,3
          stcprm(ifind(3,k)) = temp(k)
        enddo
!  for (k=0;k<3;k++) stcprm->rotate[2][k] = temp.v[k];
        call ll_geo(r_lat,r_long,temp)
        do k=1,3
          stcprm(ifind(1,k)) = temp(k)
        enddo
!  for (k=0;k<3;k++) stcprm->rotate[0][k] = temp.v[k];
        tnorm = x_prod(stcprm(ifind(3,1):ifind(3,1)+2),stcprm(ifind(1,1):ifind(1,1)+2), &
     &                 stcprm(ifind(2,1):ifind(2,1)+2))
!  norm =
!  x_product (stcprm->rotate[2],stcprm->rotate[0],stcprm->rotate[1]);
        do k=1,3
          stcprm(ifind(2,k)) = stcprm(ifind(2,k)) /tnorm
        enddo
!  for (k=0;k<3;k++) stcprm->rotate[1][k] /= norm
        tnorm = x_prod(stcprm(ifind(2,1):ifind(2,1)+2),stcprm(ifind(3,1):ifind(3,1)+2), &
     &                 stcprm(ifind(1,1):ifind(1,1)+2))
!  x_product (stcprm->rotate[1],stcprm->rotate[2],stcprm->rotate[0]);
!  stcprm->x0 = stcprm->y0 = stcprm->srotate = 0;
        stcprm(11)=0.
        stcprm(12)=0.
        stcprm(14)=0.
        stcprm(13)=1.
        stcprm(15) = REARTH
        stcprm(1) = dsin(RADPDG * conang )
! stcprm->crotate = 1.;
! stcprm->gridszeq = REARTH;
! stcprm->gamma = sin(RADPDEG * cone_ang);
!*geographic triple : i = equator @ std merid,
!       j = equator @ 90E,
!       k = North Pole */
!*map base triple : i' = M-prime meridian @ M-equator,
!     j' = M-East at location i',
!     k' = M-Pole */
!
      return
      END
