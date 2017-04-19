!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  MAP2XY           MAP 2 XY converts user units to map units
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:98-12-07
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ROUTINE TO CONVERT USER PLOT UNITS TO PAPER PLOTTING UNITS IN INCHES
!   EACH CALL TO MAPSET SCALES BETWEEN RELATIVE AND USER UNITS WHERE
!   RELATIVE IS ALWAYS BETWEEN 0-1 AND USER UNITS ARE THE CONFORMAL MAP
!   PROJECTION UNITS.  SUBSEQUENT CALLS TO MAP2XY CONVERT THE USER CONFORMAL
!   MAP UNITS TO ACTUAL INCHES ON THE PAPER FOR THE PSPLOT LIBRARY
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 07 Dec 1998 - RRD
!
! USAGE:  CALL MAPSET(PMODE,XR1,XR2,YR1,YR2,XU1,XU2,YU1,YU2)
!         CALL MAP2XY(XU,YU,XP,YP)
!   INPUT ARGUMENT LIST:
!     PMODE - logical that is true for portrait mode
!     XR1,XR2,YR1,YR2 - window coordinates (0-1) for user references frame
!     XU1,XU2,YU1,YU2 - user coordinate system that is mapped into window
!     XU,YU - input point in user units
!   OUTPUT ARGUMENT LIST:
!     XP,YP - input point mapped into output paper coordinates
!     MAPBND - common block data passed to MAPBDY subroutine
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE MAPSET(PMODE,XR1,XR2,YR1,YR2,XU1,XU2,YU1,YU2)

      LOGICAL PMODE
      COMMON /MAPBND/ XMIN,XMAX,YMIN,YMAX
      SAVE XW1,XW2,YW1,YW2,XYMIN,XDIR,YDIR,XREF,YREF

!     default paper size represents does not include 0.5 inch margin
!     relative scale 0-1 is assumed over the xdir,ydir distance

!     set directional paper size
      IF(PMODE)THEN
         XDIR=7.5
         YDIR=10.0
         XREF=0.0
         YREF=2.0
      ELSE
         XDIR=10.0
         YDIR=7.5
         XREF=3.0
         YREF=0.0
      END IF
      XYMIN=MIN(XDIR,YDIR)

!     save window limits
      XW1=XR1
      XW2=XR2
      YW1=YR1
      YW2=YR2

!     save user map limits in common area
      XMIN=XU1
      XMAX=XU2
      YMIN=YU1
      YMAX=YU2
      RETURN

!==>convert user units to paper units

      ENTRY MAP2XY(XU,YU,XP,YP)

!     convert user coordinates to user fractional units
      XUF=(XU-XMIN)/(XMAX-XMIN)
      YUF=(YU-YMIN)/(YMAX-YMIN)

!     convert user fraction to absolute fractional units
      XAF=XUF*(XW2-XW1)+XW1
      YAF=YUF*(YW2-YW1)+YW1

!     convert absolute fractional units to absolute paper units
      XP=XAF*XYMIN+XREF
      YP=YAF*XYMIN+YREF

      RETURN
      END
