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
! $Id: advint2D.f,v 1.1 2009/10/26 15:36:51 jel Exp $
!
!$$$

! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
! JCL:(10/31/02) interpolation of TWO-DIMENSIONAL array (instead of 3D); based on ADVINT
      SUBROUTINE ADVINT2D (S,NXS,NYS,X1,Y1,GLOBAL,NXP,NYP,SS)

      IMPLICIT REAL*8 (A-H,O-Z)

      CHARACTER(2) GLOBAL
! JCL:(6/6/2000)
      REAL*8   SS
      REAL*8   S(NXS,NYS)

!---------------------------------------------------------------------------------------------------
!     grid index values
      I1 = INT(X1)
      J1 = INT(Y1)

!     fractional interpolation factors to x1,y1,zx
      XF = X1-I1
      YF = Y1-J1

! JCL:(07/12/2004) added cyclic boundary condition flag
!     set upper index
      I1P = I1+1
      J1P = J1+1
!     cyclic boundary conditons
!      IF (GLOBAL.EQ."gl") THEN
!tk(20160317)
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
         PRINT *, 'Subroutine ADVINT2D: one or more array indices exceed declared dimension.'
         PRINT '(a,2(i0,1x),g12.4,1x,i0)', 'NXS, I1, XF, I1P: ', NXS, I1, XF, I1P
      END IF
      IF (J1 < 1  .OR.  J1 > NYS  .OR.  (YF > 0d0  .AND. J1P > NYS)) THEN
         PRINT *, 'Subroutine ADVINT2D: one or more array indices exceed declared dimension.'
         PRINT '(a,2(i0,1x),g12.4,1x,i0)', 'NYS, J1, YF, J1P: ', NYS, J1, YF, J1P
      END IF


      IF (XF == 0.0) THEN
!     value at x1 along bottom of lowest layer
        SB = S(I1,J1)
!     value at x1 along top of lowest layer
!       IF(YF /= 0.0) ST = S(I1,J1+1)
        IF (YF /= 0.0) ST = S(I1,J1P)
      ELSE
!     value at x1 along bottom of lowest layer
!       SB = (S(I1+1,J1)-S(I1,J1))*XF+S(I1,J1)
        SB = (S(I1P,J1)-S(I1,J1))*XF+S(I1,J1)
!     value at x1 along top of lowest layer
!       IF(YF /= 0.0) ST = (S(I1+1,J1+1)-S(I1,J1+1))*XF+S(I1,J1+1)
        IF (YF /= 0.0) ST = (S(I1P,J1P)-S(I1,J1P))*XF+S(I1,J1P)
      END IF
!     value at x1,y1 of lowest layer
      IF (YF == 0.0) THEN
        SS = SB
      ELSE
        SS = (ST-SB)*YF+SB
      END IF

      RETURN
      END
