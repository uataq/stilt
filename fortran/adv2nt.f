!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADV2NT           ADVection INTerpolation of 2D meteo variable
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION INTERPOLATION OF 3D METEO VARIABLES TO A POINT IN X,Y,
!   USING LINEAR METHODS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 28 Sep 2000 (RRD) - fortran90 upgrade
!                 12 Mar 2001 (RRD) - cyclic boundary conditions
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 11 May 2005 (RRD) - improved cyclic boundary
!                 15 Nov 2007 (RRD) - north pole global interpolate
!
! USAGE:  CALL ADVINT(S,X1,Y1,SS,GLOBAL,NXP,NYP)
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

SUBROUTINE ADV2NT(S,X1,Y1,SS,GLOBAL,NXP,NYP)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,      INTENT(IN)    :: s(:,:)        ! field for interpolation
  REAL,      INTENT(IN)    :: x1,y1         ! position of interpolated value
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  CHARACTER(2), INTENT(IN)    :: global        ! cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundaries  

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER                  :: i1,j1,i1p,j1p
  REAL                     :: xf,yf,sb,st

! grid index values
  I1=INT(X1)
  J1=INT(Y1)

! assign dummy interpolation factors
  XF=X1
  YF=Y1

! cyclic boundary conditions
!tk(20160317)
!   IF(GLOBAL)THEN  
  IF(GLOBAL .EQ."gl")THEN
     IF(I1.GT.NXP)THEN
        I1=1
        XF=XF-NXP
     END IF
     IF(I1.LT.1)THEN
        I1=NXP
        XF=NXP+XF
     END IF
     IF(J1.LT.1)THEN
        J1=1
        YF=2.0-YF
     END IF
     IF(J1.EQ.NYP)THEN
        YF=2*NYP-YF
        J1=INT(YF)
        IF(J1.EQ.NYP)THEN
           J1=J1-1
           YF=FLOAT(NYP)
        END IF
     END IF
  END IF

!tk(20160317)
  IF(GLOBAL .EQ."nh")THEN
! northern hemisphere
     IF(I1.GT.NXP)THEN
        I1=1
        XF=XF-NXP
     END IF
     IF(I1.LT.1)THEN
        I1=NXP
        XF=NXP+XF
     END IF
  !   IF(J1.LT.1)THEN
  !     J1=1
  !      YF=2.0-YF
  !   END IF
     IF(J1.EQ.NYP)THEN
        YF=2*NYP-YF
        J1=INT(YF)
        IF(J1.EQ.NYP)THEN
           J1=J1-1
           YF=FLOAT(NYP)
        END IF
     END IF
  END IF

!tk(20160317)
  IF(GLOBAL .EQ."sh")THEN
! southern hemisphere
     IF(I1.GT.NXP)THEN
        I1=1
        XF=XF-NXP
     END IF
     IF(I1.LT.1)THEN
        I1=NXP
        XF=NXP+XF
     END IF
     IF(J1.LT.1)THEN
       J1=1
        YF=2.0-YF
     END IF
 !    IF(J1.EQ.NYP)THEN
 !       YF=2*NYP-YF
 !       J1=INT(YF)
 !       IF(J1.EQ.NYP)THEN
 !          J1=J1-1
 !          YF=FLOAT(NYP)
 !       END IF
  !   END IF
  END IF


! set upper index
  I1P=I1+1
  J1P=J1+1

! cyclic boundary conditions
!  IF(GLOBAL.AND.I1P.GT.NXP)I1P=1
!tk(20160317)
  IF(GLOBAL.EQ.'gl'.AND.I1P.GT.NXP)I1P=1
  IF((GLOBAL.EQ.'nh' .or. GLOBAL.EQ.'sh').AND.I1P.GT.NXP)I1P=1

  if (i1p > size(s,1) .or. j1p > size(s,2)) then
     write (*,*) 'subroutine adv2nt: i1p and/or j1p exceeds allowable limit:', &
          & i1p, j1p, size(s,1), size(s,2)
!!$     I1P=min(I1+1,size(s,1))
!!$     J1P=min(J1+1,size(s,2))
  end if

! fractional interpolation factors to x1,y1,zx
  XF=XF-FLOAT(I1)
  YF=YF-FLOAT(J1)

! value at x1 bottom
  SB=(S(I1P,J1 )-S(I1,J1 ))*XF+S(I1,J1)
! value at x1 top 
  ST=(S(I1P,J1P)-S(I1,J1P))*XF+S(I1,J1P)
! value at x1,y1 
  SS=(ST-SB)*YF+SB

END SUBROUTINE adv2nt
