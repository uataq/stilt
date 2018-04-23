!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVINT           ADVection INTerpolation of 3D meteo variable
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION INTERPOLATION OF 3D METEO VARIABLES TO A POINT IN X,Y,
!   USING LINEAR METHODS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL ADVINT(S,NXS,NYS,NZM,X1,Y1,ZX,SS)
!   INPUT ARGUMENT LIST:
!     S           - real      input variable of dimensions NXS,NYS,NZM
!     X1,Y1 - real      position of interpolated value
!     ZX    - real      vertical interpolation index fraction
!   OUTPUT ARGUMENT LIST:
!     SS    - real      value of S at x1,y1,z1
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: advintZML.f,v 1.1 2009/10/26 15:36:52 jel Exp $
!
!$$$

! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
! JCL:4/28/00 Have a special interpolation routine for extracting
!      mixed-layer height info.--declare arrays as having MAX dimensions
!     Most of code is identical to ADVINTZML subroutine

      SUBROUTINE ADVINTZML (S,NXS,NYS,   X1,Y1,ZX,GLOBAL,NXP,NYP,SS)
!     SUBROUTINE ADVINTZML(S,NXS,NYS,NZM,X1,Y1,ZX,SS)

      use module_defsize
      IMPLICIT REAL*8 (A-H,O-Z)

! JCL:4/28/00
!     array size information

      CHARACTER(2) GLOBAL
      REAL*8  S(NXM,NYM,1)

!---------------------------------------------------------------------------------------------------
!     grid index values
      I1 = INT(X1)
      J1 = INT(Y1)
      K1 = INT(ZX)

!     fractional interpolation factors to x1,y1,zx
      XF = X1-I1
      YF = Y1-J1
      ZF = ZX-K1
      IF (ZF > TINY(ZF)) STOP 'advintZML: vertical interpolation not allowed. STOP.'

! JCL:(07/12/2004) added cyclic boundary condition flag
!     set upper index
      I1P = I1+1
      J1P = J1+1
!     cyclic boundary conditons
 !tk(20160317)
!      IF (GLOBAL) THEN
      IF (GLOBAL.EQ."gl".OR.GLOBAL.EQ."nh") THEN
        IF (I1P > NXP) I1P = 1
        IF (J1P > NYP) J1P = NYP
! JCL:(07/26/2004) occasionally even I1 or J1 exceeds the limit
        IF (I1 > NXP) THEN
           I1 = 1
           I1P = 2
        END IF
        IF (J1 > NYP) THEN
           J1 = NYP
           J1P = NYP-1
        END IF
      END IF
      IF (GLOBAL.EQ."sh") THEN
        IF (I1P > NXP) I1P = 1
!        IF (J1P > NYP) J1P = NYP
! JCL:(07/26/2004) occasionally even I1 or J1 exceeds the limit
        IF (I1 > NXP) THEN
           I1 = 1
           I1P = 2
        END IF
!        IF (J1 > NYP) THEN
!           J1 = NYP
!           J1P = NYP-1
!        END IF
      END IF


      !--- index range checks in x and y

      IF (I1 < 1  .OR.  I1 > NXS  .OR.  (XF > 0d0  .AND. I1P > NXS)) THEN
         PRINT *, 'Subroutine ADVINTZML: one or more array indices exceed declared dimension.'
         PRINT '(a,2(i0,1x),g12.4,1x,i0)', 'NXS, I1, XF, I1P: ', NXS, I1, XF, I1P
      END IF
      IF (J1 < 1  .OR.  J1 > NYS  .OR.  (YF > 0d0  .AND. J1P > NYS)) THEN
         PRINT *, 'Subroutine ADVINTZML: one or more array indices exceed declared dimension.'
         PRINT '(a,2(i0,1x),g12.4,1x,i0)', 'NYS, J1, YF, J1P: ', NYS, J1, YF, J1P
      END IF


!     value at x1 along bottom of lowest layer
!     SB = (S(I1+1,J1  ,K1)-S(I1,J1  ,K1))*XF+S(I1,J1,K1)
      SB = (S(I1P,J1  ,K1)-S(I1,J1  ,K1))*XF+S(I1,J1,K1)
!     value at x1 along top of lowest layer
!     ST = (S(I1+1,J1+1,K1)-S(I1,J1+1,K1))*XF+S(I1,J1+1,K1)
      ST = (S(I1P,J1P,K1)-S(I1,J1P,K1))*XF+S(I1,J1P,K1)
!     value at x1,y1 of lowest layer
      SS = (ST-SB)*YF+SB

!     no vertical interpolation required
      IF (ZF == 0.0) RETURN

!     value at x1 along bottom of upper layer
!     SB = (S(I1+1,J1  ,K1+1)-S(I1,J1  ,K1+1))*XF+S(I1,J1,K1+1)
      SB = (S(I1P,J1  ,K1+1)-S(I1,J1  ,K1+1))*XF+S(I1,J1,K1+1)
!     value at x1 along top of upper layer
!     ST = (S(I1+1,J1+1,K1+1)-S(I1,J1+1,K1+1))*XF+S(I1,J1+1,K1+1)
      ST = (S(I1P,J1P,K1+1)-S(I1,J1P,K1+1))*XF+S(I1,J1P,K1+1)
!     value at x1,y1 of upper layer
      SH = (ST-SB)*YF+SB

!     value at x1,y1,z1
      SS = (SH-SS)*ZF+SS

      RETURN
      END
