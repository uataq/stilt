!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GRFBDY           GRAphics map BounDarY 
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:98-12-07
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CREATES A MAP BACKGROUND IN LAT/LON COORDINATES. THIS PROGRAM
!   IS A SIMPLIFIED VERSION OF MAPBDY
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 14 May 2003 (RRD)
!                  17 Jul 2003 (RRD) - revised directory search 
!
! USAGE:  CALL GRFBDY(XL,YL)
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!     NONE
!   INPUT FILES:
!     UNIT 15 - ascii map background vector file
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE GRFBDY(XL,YL)

  IMPLICIT NONE

  REAL, INTENT(IN) :: XL(5),YL(5)  ! map domain limits

  LOGICAL       :: FTEST
  CHARACTER(80) :: FNAME
  INTEGER       :: N,K,KPTS,KRET,KLINE
  REAL          :: XX,XA(1000),PLON(1000)
  REAL          :: YY,YA(1000),PLAT(1000)

  FNAME='arlmap'
  INQUIRE(FILE=FNAME,EXIST=FTEST)
  IF(FTEST)THEN
     OPEN(15,FILE=FNAME,STATUS='OLD')
  ELSE
     FNAME='graphics/arlmap'
     INQUIRE(FILE=FNAME,EXIST=FTEST)
     IF(FTEST)THEN
        OPEN(15,FILE=FNAME,STATUS='OLD')
     ELSE
        FNAME='../graphics/arlmap'
        INQUIRE(FILE=FNAME,EXIST=FTEST)
        IF(FTEST)THEN
           OPEN(15,FILE=FNAME,STATUS='OLD')
        ELSE
           WRITE(*,*)'WARNING: map background file not found'
           RETURN
        END IF
     END IF
  END IF

  KRET=0
vector : DO WHILE (KRET.EQ.0)

     READ(15,'(2I5)',IOSTAT=KRET)KLINE,KPTS
     IF(KRET.NE.0) EXIT vector

     READ(15,'(10F6.2)')(PLAT(K),K=1,KPTS)
     READ(15,'(10F7.2)')(PLON(K),K=1,KPTS)
     IF(KPTS.LT.2) CYCLE vector

     K=0
     DO N=1,KPTS
        CALL MAP2XY(PLON(N),PLAT(N),XX,YY)
        IF(XX.LT.XL(1).OR.XX.GT.XL(3).OR.YY.LT.YL(1).OR.YY.GT.YL(3))THEN
           IF(K.GT.0) CALL SLDCRV(XA,YA,K,0.005)
           K=0
        ELSE
           K=K+1
           XA(K)=XX
           YA(K)=YY
        END IF
     END DO
     IF(K.GT.0) CALL SLDCRV(XA,YA,K,0.01)

  END DO vector

  CLOSE (15)
  END SUBROUTINE grfbdy
