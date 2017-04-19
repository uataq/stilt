      subroutine sobstr(stcprm, p_lat, p_lon)
!*
!*  Set Map Parameters for an OBlique STeReographic Projection
!*  Inputs: p_lat,P_lon - latitude and longitude of the
!*          projection pole (tangent point)
!*  Outputs: stcprm - map parameters
!*/
! CHG(03/12/03)
!
! $Id: sobstr.f,v 1.1 2010/10/26 19:03:41 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      dimension stcprm(15)
      TEMP=90.
! CHG(09/15/03) pass on pole lat and lon (temp2 and temp4) opposing ref_lat and lon
! might need adjustment for northern hemisphere
      TEMP2=-(180.0+p_lat)
      TEMP4=p_lon+180.0
!      call mpstrt(stcprm, TEMP, p_lat,p_lon, TEMP2,p_lon)  Original: NAN's
! CHG(09/15/03) pass on pole lat and lon (temp2 and temp4) opposing ref_lat and lon
!      WRITE(*,*)'mpstrargs:',TEMP,TEMP2,TEMP4, p_lat,p_lon
      call mpstrt(stcprm, TEMP, TEMP2,TEMP4, p_lat,p_lon)
!      call mpstrt(stcprm, TEMP, p_lat,p_lon,TEMP2,TEMP4) !change pole and ref => doesn't work (nan's)

! JCL(03/28/03) according to mpstrt.f:
!     mpstrt(stcprm, conang, p_lat, p_long, r_lat, r_long)
      return
      END
