!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVSRT           ADVection SoRTting to minimize data I/O
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION SORTTING TO MINIMIZE DATA I/O WHEN METEOROLOGICAL
!   SUBGRID IS DEFINED.  MAY BE CALLED BEFORE EACH ADVECTION CYCLE
!   WHEN THE SUB-GRID IS SMALLER THAN THE FULL METEOROLOGICAL GRID.
!   ALL ELEMENTS THAT ARE ON THE CURRENT SUB-GRID ARE GROUPED TOGETHER
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 23 Jul 1998 (RRD)
!                 02 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL ADVSRT(KPM,NXP,NYP,LX1,LY1,XPOS,YPOS,NSORT)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE ADVSRT(KPM,NXP,NYP,LX1,LY1,XPOS,YPOS,NSORT)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,   INTENT(IN)    :: kpm         ! number of points
  INTEGER,   INTENT(IN)    :: nxp,nyp     ! max meteo grid size
  INTEGER,   INTENT(IN)    :: lx1,ly1     ! lower left corner of subgrid
  REAL,      INTENT(IN)    :: xpos  (:)   ! puff center positions (grid units)
  REAL,      INTENT(IN)    :: ypos  (:)   ! puff center positions (grid units)
  INTEGER,   INTENT(OUT)   :: nsort (:)   ! sortted array by position

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER    :: nxs2,nys2,kk,kpt
  REAL       :: xnew,ynew,xcnt,ycnt

!-------------------------------------------------------------------------------

! need at least three elements to sort
  IF(KPM.LE.2)RETURN
! exit if subgrid not yet set
  IF(LX1.LE.0.OR.LY1.LE.0)RETURN

! index for points within subgrid (eventually equals KPM)
  KPT=0

! set scan range for subgrid
  NXS2=(NXP-1)/2-2
  NYS2=(NYP-1)/2-2

! default initial position to center of current subgrid
  XCNT=FLOAT(LX1+NXS2)
  YCNT=FLOAT(LY1+NYS2)

! zero out index array to indicate unsorted particles
  DO KK=1,KPM
     NSORT(KK)=0
  END DO

! loop through particle array until none left
  DO WHILE (KPT.LT.KPM)

!    reset the subgrid center position each pass
     XNEW=-1.0
     YNEW=-1.0

!    determine if position within range of subgrid
     DO KK=1,KPM
!       only test particles not previously identified
        IF(NSORT(KK).EQ.0)THEN
           IF(XPOS(KK).GT.XCNT-FLOAT(NXS2).AND.XPOS(KK).LT.XCNT+FLOAT(NXS2).AND. &
              YPOS(KK).GT.YCNT-FLOAT(NYS2).AND.YPOS(KK).LT.YCNT+FLOAT(NYS2))THEN
!             identify elements within current subgrid
              KPT=KPT+1
              NSORT(KPT)=KK
           ELSE
!             save position of first particle not in array
!             which will be the center of the next subgrid
              IF(XNEW.LT.0.0.AND.YNEW.LT.0.0)THEN
                 XNEW=XPOS(KK)
                 YNEW=YPOS(KK)
              END IF
           END IF
        END IF
     END DO

!    set new sort center point to last outside
     XCNT=XNEW
     YCNT=YNEW

! while loop
  END DO

END SUBROUTINE advsrt
