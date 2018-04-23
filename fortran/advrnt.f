!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVRNT           ADVection Regression iNTerpolation
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION INTERPOLATION OF 3D METEO VARIABLES TO A POINT IN X,Y,
!   USING LINEAR REGRESSION TO FIT A POLYNOMIAL.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 21 Jan 2003 (RRD) - initial version from adv3nt
!
! USAGE:  CALL ADVRNT(S,X1,Y1,ZX,SS,GLOBAL,NXP,NYP)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            none
!   OUTPUT FILES:           none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE ADVRNT(S,X1,Y1,ZX,SS,GLOBAL,NXP,NYP)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
  REAL,      INTENT(IN)    :: x1,y1         ! position of interpolated value
  REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  CHARACTER(2),   INTENT(IN)    :: global        ! cyclic boundary condition flag
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundaries

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER                  :: i,j,i1,j1,k1,ii,jj,kk 
  REAL                     :: xf,yf,zf,sb,st,sh 
  REAL                     :: val(4),b0,b1,b2

! grid index values
  I1=INT(X1)
  J1=INT(Y1)
  K1=INT(ZX)

! linear vertical interpolation factor
  ZF=ZX-FLOAT(K1)

! find value at X1 for all J about Y1
  DO J=J1-1,J1+2
!     IF(GLOBAL)THEN
 !tk(20160317)
    IF(GLOBAL .EQ. "gl")THEN
        IF(J.GT.NYP)JJ=2*NYP-J
        IF(J.LT.1  )JJ=1-J
 !tk(20160317)
     ELSE IF (GLOBAL .EQ. "nh")THEN
        IF(J.GT.NYP)JJ=2*NYP-J
      !  IF(J.LT.1  )JJ=1-J
     ELSE IF (GLOBAL .EQ. "sh")THEN
      !  IF(J.GT.NYP)JJ=2*NYP-J
        IF(J.LT.1  )JJ=1-J
     ELSE
        JJ=J
     END IF

     CALL POLYNT(0,XF,YF,B0,B1,B2)
     DO I=I1-1,I1+2
  !      IF(GLOBAL)THEN
 !tk(20160317)
      IF(GLOBAL.EQ."gl" .OR. GLOBAL.EQ."nh" .OR. GLOBAL.EQ."sh")THEN
           IF(I.GT.NXP)II=I-NXP
           IF(I.LT.1  )II=NXP+I
        ELSE
 !tk(20160317)
           II=I
        END IF
        CALL POLYNT(1,FLOAT(II),S(II,JJ,K1),B0,B1,B2)
     END DO

     CALL POLYNT(-1,XF,YF,B0,B1,B2)
     VAL(J-J1+2)=B2*X1*X1+B1*X1+B0
  END DO

! find value at Y1 from previous X1 values along J
  CALL POLYNT(0,XF,YF,B0,B1,B2)
  DO J=J1-1,J1+2
 !    IF(GLOBAL)THEN
 !tk(20160317)
      IF(GLOBAL .EQ. "gl")THEN
        IF(J.GT.NYP)JJ=2*NYP-J
        IF(J.LT.1  )JJ=1-J
      ELSE IF (GLOBAL .EQ. "nh")THEN
        IF(J.GT.NYP)JJ=2*NYP-J
       ! IF(J.LT.1  )JJ=1-J
      ELSE IF (GLOBAL .EQ. "sh")THEN
       ! IF(J.GT.NYP)JJ=2*NYP-J
        IF(J.LT.1  )JJ=1-J
     ELSE
        JJ=J
     END IF
     CALL POLYNT(1,FLOAT(JJ),VAL(J-J1+2),B0,B1,B2)
  END DO

  CALL POLYNT(-1,XF,YF,B0,B1,B2)
  SS=B2*Y1*Y1+B1*Y1+B0

! value requested only on the lower boundary
  IF(ZF.EQ.0.0)RETURN

! find value at X1 for all J about Y1
  DO J=J1-1,J1+2
 !    IF(GLOBAL)THEN
 !tk(20160317)
      IF(GLOBAL .EQ. "gl")THEN
        IF(J.GT.NYP)JJ=2*NYP-J
        IF(J.LT.1  )JJ=1-J
     ELSE IF (GLOBAL .EQ. "nh")THEN
        IF(J.GT.NYP)JJ=2*NYP-J
       ! IF(J.LT.1  )JJ=1-J
     ELSE IF (GLOBAL .EQ. "sh")THEN
       ! IF(J.GT.NYP)JJ=2*NYP-J
        IF(J.LT.1  )JJ=1-J
     ELSE
        JJ=J
     END IF

     CALL POLYNT(0,XF,YF,B0,B1,B2)
     DO I=I1-1,I1+2
   !     IF(GLOBAL)THEN
  !tk(20160317)
 !tk(20160317)
      IF(GLOBAL.EQ."gl" .OR. GLOBAL.EQ."nh" .OR. GLOBAL.EQ."sh")THEN
           IF(I.GT.NXP)II=I-NXP
           IF(I.LT.1  )II=NXP+I
        ELSE
           II=I
        END IF
        CALL POLYNT(1,FLOAT(II),S(II,JJ,K1+1),B0,B1,B2)
     END DO

     CALL POLYNT(-1,XF,YF,B0,B1,B2)
     VAL(J-J1+2)=B2*X1*X1+B1*X1+B0
  END DO

! find value at Y1 from previous X1 values along J
  CALL POLYNT(0,XF,YF,B0,B1,B2)
  DO J=J1-1,J1+2
   !  IF(GLOBAL)THEN
 !tk(20160317)
      IF(GLOBAL .EQ. "gl")THEN
        IF(J.GT.NYP)JJ=2*NYP-J
        IF(J.LT.1  )JJ=1-J
      ELSE IF (GLOBAL .EQ. "nh")THEN
        IF(J.GT.NYP)JJ=2*NYP-J
       ! IF(J.LT.1  )JJ=1-J
      ELSE IF (GLOBAL .EQ. "sh")THEN
       ! IF(J.GT.NYP)JJ=2*NYP-J
        IF(J.LT.1  )JJ=1-J
     ELSE
        JJ=J
     END IF
     CALL POLYNT(1,FLOAT(JJ),VAL(J-J1+2),B0,B1,B2)
  END DO

  CALL POLYNT(-1,XF,YF,B0,B1,B2)
  SH=B2*Y1*Y1+B1*Y1+B0

! vertical interpolation to value at x1,y1,z1
  SS=(SH-SS)*ZF+SS

END SUBROUTINE advrnt

!------------------------------------------------------------

SUBROUTINE POLYNT(K,X,Y,B0,B1,B2)

! solution of the equation: Y = B0 + B1 X + B2 X^2
! for regression solutions

  SAVE A,B,C,D,E,C1,C2,C3

  IF(K.EQ.0)THEN
!    initialize regression sums
     A=0.0
     B=0.0
     C=0.0
     D=0.0
     E=0.0
     C1=0.0
     C2=0.0
     C3=0.0

  ELSEIF(K.LT.0)THEN
!    computes solution for set of simultaneous linear equations
!    where x=B0; y=B1; z=B2 as defined in the argument list
!    Ax + By + Cz = C1
!    Bx + Cy + Dz = C2
!    Cx + Dy + Ez = C3

     B2=( (C1*B/A-C2)*(B*C/A-D)-(C1*C/A-C3)*(B*B/A-C) ) /                      &
        ( (B*C/A-D)**2 - (C*C/A-E)*(B*B/A-C) )

     B1=(C1*B/A-C2)/(B*B/A-C) - B2*(B*C/A-D)/(B*B/A-C)

     B0=C1/A - B1*B/A - B2*C/A

  ELSE
!    perform summation
     A=A+1.0
     B=B+X
     C=C+X*X
     D=D+X*X*X
     E=E+X*X*X*X
     C1=C1+Y
     C2=C2+X*Y
     C3=C3+X*X*Y

  END IF

END SUBROUTINE polynt
