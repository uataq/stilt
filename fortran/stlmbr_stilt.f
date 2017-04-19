      subroutine stlmbr_stilt(stcprm, reflat, reflon)
!*
!*  Set Map Parameters for a North Polar LaMBeRt Conic Conformal
!*    Projection
!*  Inputs: reflat - tangent latitude of the tangent latitude of
!*                   cone.
!*          reflon - midrange longitude (180 degrees from cut)
!*  Outputs: stcprm - map parameters
!*/
! CHG(03/12/03)
!
! $Id: stlmbr_stilt.f,v 1.1 2010/10/26 19:14:24 jel Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
! CHG(03/11/03) according to mpstrt.f:
!     mpstrt(stcprm, conang, p_lat, p_long, r_lat, r_long)
      TEMP=90.
      TEMP2=0.
! JCL(03/11/03)
      CALL mpstrt(stcprm,reflat,TEMP,TEMP2,reflat,reflon)
      return
      END

      subroutine stcmap(stcprm, tnglat, reflon)
!*
!*  Set Map Parameters for a North Polar LaMBeRt Conic Conformal
!*    Projection
!*  included for compatibliity with previous version
!*  Inputs: tnglat - tangent latitude of the tangent latitude of
!*                   cone.
!*          reflon - midrange longitude (180 degrees from cut)
!*  Outputs: stcprm - map parameters
!*/
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
      TEMP=90.
      TEMP2=0.
!      CALL mpstrt(stcprm,tnglat,TEMP,reflon,TEMP2,reflon)
! JCL(03/11/03)
!      XXX3=35.0
      CALL mpstrt(stcprm,tnglat,TEMP,TEMP2,tnglat,reflon)
!      CALL mpstrt(stcprm,tnglat,TEMP,TEMP2,XXX3,reflon)
      return
      END
