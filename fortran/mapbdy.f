!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  MAPBDY           MAP BounDarY creates a map background
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:98-12-07
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CREATES A MAP BACKGROUND ACCORDING TO THE PROJECTION
!   SPECIFIED BY THE PARMAP ARGUMENT
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 08 Dec 1998 (RRD)
!                  18 Feb 1999 (RRD) - converted map file to ascii
!                  28 May 1999 (RRD) - changed argument list
!                  09 Jun 1999 (RRD) - lat/lon lines position fix
!                  15 Jun 1999 (RRD) - search local directory first for arlmap
!                  03 Feb 2000 (RRD) - split thickness from argument list
!                  29 Mar 2000 (RRD) - prevent lat/lon redraw
!                  20 Apr 2001 (RRD) - support lat/lon grids
!                  31 Oct 2001 (RRD) - argument list gbl routines
!                  01 Mar 2002 (RRD) - added clip region for drawing
!                  11 Apr 2002 (RRD) - simplified ptest for absoft 
!                  11 Apr 2003 (RRD) - clipping optional (absoft)
!                  17 Jul 2003 (RRD) - revised directory search
!                  19 Aug 2003 (GDR) - mod loc of lon label above 80 deg lat
!                  04 Dec 2003 (RRD) - replaced llint with aintr for <1 deg
!                  09 Feb 2004 (RRD) - correction on drawing lat/lon lines 
!                  15 Dec 2004 (RRD) - g95 compatibility
!                  10 Jan 2005 (RRD) - dynamic array allocation
!                  09 Mar 2005 (RRD) - added array size tests
!                  12 Jan 2007 (RRD) - cyclindrical equidistant
!                  18 Oct 2007 (RRD) - prime meridian issues on global plots
!
! USAGE:  CALL MAPBDY(KPROJ,LLGRID,PARMAP,AINTR,XTHICK,FNAME)
!   INPUT ARGUMENT LIST:
!     LLGRID - logical flag
!     PARMAP - real 9 element array defining the map projection
!     AINTR - real lat/lon line intervals (0 for none)
!     XTHICK - map background thickness (0.5 = half  1.0 = normal  2.0=twice)
!     GFILE - background file name
!   OUTPUT ARGUMENT LIST:
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

      SUBROUTINE MAPBDY(KPROJ,LLGRID,PARMAP,AINTR,XTHICK,FNAME)

      INTEGER           :: KPROJ
      LOGICAL           :: LLGRID, FTEST, PTEST1, PTEST2, CLIP
      CHARACTER(80)     :: FNAME
      CHARACTER(20)     :: POINTS
      CHARACTER(10)     :: COLOR
      CHARACTER(6)      :: LATLON
      REAL              :: PARMAP(9)
      REAL, ALLOCATABLE :: PLAT(:),PLON(:),XARR(:),YARR(:)

      COMMON /MAPBND/ XU1,XU2,YU1,YU2
      COMMON /MAPCLP/ XC1,XC2,YC1,YC2

!==>limits check

      PTEST1(XP,YP)=(XP.GE.XU1.AND.XP.LE.XU2.AND.YP.GE.YU1.AND.YP.LE.YU2) 
      PTEST2(XP,YP)=(XP.GE.XU2.AND.XP.LE.XU1.AND.YP.GE.YU2.AND.YP.LE.YU1)

!     check for boundary clip region
      CLIP=.TRUE.

!=>initial  allocations

      MARR=5000 ! temporary array space for drawing
      ALLOCATE (XARR(MARR),YARR(MARR))
      MPTS=500  ! map background vector array
      ALLOCATE (PLAT(MPTS),PLON(MPTS))

!==>check range on line thickness multiplier

      THICK=AMAX1(0.01, AMIN1(100.0, XTHICK))

!==>set map paper units clip region

      CALL MAP2XY(XU1,YU1,XC1,YC1)
      CALL MAP2XY(XU2,YU2,XC2,YC2)

!------------------------------------------------------------------------------
! mark lat/lon line labels
!------------------------------------------------------------------------------

!     determine position of map center
      IF(LLGRID)THEN
         CALL GBL2LL(1,1,10.5*(XU1+XU2),0.5*(YU1+YU2),CLAT,CLON)
      ELSEIF(KPROJ.EQ.4)THEN 
         CALL CYL2LL(0.5*(XU1+XU2),0.5*(YU1+YU2),CLAT,CLON)
      ELSE
         CALL CXY2LL(PARMAP,0.5*(XU1+XU2),0.5*(YU1+YU2),CLAT,CLON)
      END IF

      INTRA=AINTR*10.0
!     print out lat/lon labels
      IF(AINTR.GT.0.0)THEN 

         CLON=NINT(CLON/AINTR)*AINTR
         CLAT=NINT(CLAT/AINTR)*AINTR

!        longitude labels
         IF(CLAT.GT.80.0.OR.CLAT.LT.-80.0)THEN
           XLAT=CLAT/2.0+AINTR/2.0
         ELSE
           XLAT=CLAT+AINTR/2.0
         END IF

         DO LONX=-(1800-INTRA),(1800-INTRA),INTRA
            XLON=LONX/10.0
            IF(LLGRID)THEN
               CALL GBL2XY(1,1,XLAT,XLON,XP,YP)
            ELSEIF(KPROJ.EQ.4)THEN
               CALL CYL2XY(XLAT,XLON,XP,YP)
            ELSE
               CALL CLL2XY(PARMAP,XLAT,XLON,XP,YP)
            END IF

            IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
               CALL MAP2XY(XP,YP,XX,YY)
               IF(AINTR.LT.1.0)THEN
                  WRITE(LATLON,'(F6.1)')XLON
                  CALL KEKSYMC(XX,YY,0.10,LATLON,0.0,6,1)
               ELSE
                  WRITE(LATLON,'(I4)')INT(XLON)
                  CALL KEKSYMC(XX,YY,0.10,LATLON,0.0,4,1)
               END IF
            END IF
         END DO

!        latitude labels
         XLON=CLON+AINTR/2.0

         DO LATX=-(900-INTRA),(900-INTRA),INTRA
            XLAT=LATX/10.0 
            IF(LLGRID)THEN
               CALL GBL2XY(1,1,XLAT,XLON,XP,YP)
            ELSEIF(KPROJ.EQ.4)THEN
               CALL CYL2XY(XLAT,XLON,XP,YP)
            ELSE
               CALL CLL2XY(PARMAP,XLAT,XLON,XP,YP)
            END IF

            IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
               CALL MAP2XY(XP,YP,XX,YY)
               IF(AINTR.LT.1.0)THEN
                  WRITE(LATLON,'(F5.1)')XLAT
                  CALL KEKSYMC(XX,(YY-0.05),0.10,LATLON,0.0,5,1)
               ELSE
                  WRITE(LATLON,'(I3)')INT(XLAT)
                  CALL KEKSYMC(XX,(YY-0.05),0.10,LATLON,0.0,3,1)
               END IF
            END IF
         END DO
      END IF

      IF(AINTR.NE.0.0)THEN 

!------------------------------------------------------------------------------
! draw lines of longitude  
!------------------------------------------------------------------------------

         DO LONX=-1800,(1800-INTRA),INTRA
            XLON=LONX/10.0
            N=0
            NPTS=0
            FTEST=.FALSE.

            DO LATX=-890,890,MAX(MIN(1,INTRA/10),1)
               XLAT=LATX/10.0
               N=N+1
  
               IF(LLGRID)THEN
                  CALL GBL2XY(1,1,XLAT,XLON,XP,YP)
               ELSEIF(KPROJ.EQ.4)THEN
                  CALL CYL2XY(XLAT,XLON,XP,YP)
               ELSE
                  CALL CLL2XY(PARMAP,XLAT,XLON,XP,YP)
               END IF

               IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
                  IF(.NOT.FTEST.AND.N.GT.1.AND.CLIP)THEN
!                    add previous position just outside domain
                     CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
                     CALL MAP2XY(XNEW,YNEW,XX,YY)
                     NPTS=MIN(NPTS+1,MARR)
                     XARR(NPTS)=XX
                     YARR(NPTS)=YY
                  END IF
                  CALL MAP2XY(XP,YP,XX,YY)
                  NPTS=MIN(NPTS+1,MARR)
                  XARR(NPTS)=XX
                  YARR(NPTS)=YY
                  FTEST=.TRUE.

               ELSEIF(FTEST)THEN
                  IF(CLIP)THEN
!                    add current off domain position to line
                     CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
                     CALL MAP2XY(XNEW,YNEW,XX,YY)
                     NPTS=MIN(NPTS+1,MARR)
                     XARR(NPTS)=XX
                     YARR(NPTS)=YY
                  END IF
                  CALL DSHCRV(XARR,YARR,NPTS,20,0.005*THICK)
                  NPTS=0
                  FTEST=.FALSE.
               END IF

!              save previous value
               XOLD=XP
               YOLD=YP
            END DO
            IF(FTEST)CALL DSHCRV(XARR,YARR,NPTS,20,0.005*THICK)
            IF(NPTS.EQ.MARR)WRITE(*,*)'Exceeding (XARR,YARR) array limits'
         END DO

!------------------------------------------------------------------------------
! draw lines of latitude  
!------------------------------------------------------------------------------

         DO LATX=-(900-INTRA),(900-INTRA),INTRA
            XLAT=LATX/10.0
            N=0
            NPTS=0
            FTEST=.FALSE.

            DO LONX=-1800,1800,MAX(MIN(1,INTRA/10),1)
               XLON=LONX/10.0
               N=N+1
  
               IF(LLGRID)THEN
                  CALL GBL2XY(1,1,XLAT,XLON,XP,YP)
               ELSEIF(KPROJ.EQ.4)THEN
                  CALL CYL2XY(XLAT,XLON,XP,YP)
               ELSE
                  CALL CLL2XY(PARMAP,XLAT,XLON,XP,YP)
               END IF

               IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
                  IF(.NOT.FTEST.AND.N.GT.1.AND.CLIP)THEN
!                    add previous position just outside domain
                     CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
                     CALL MAP2XY(XNEW,YNEW,XX,YY)
                     NPTS=MIN(NPTS+1,MARR)
                     XARR(NPTS)=XX
                     YARR(NPTS)=YY
                  END IF
                  CALL MAP2XY(XP,YP,XX,YY)

                  NPTS=MIN(NPTS+1,MARR)
                  XARR(NPTS)=XX
                  YARR(NPTS)=YY
                  FTEST=.TRUE.

               ELSEIF(FTEST)THEN
                  IF(CLIP)THEN
!                    add current off domain position to line
                     CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
                     CALL MAP2XY(XNEW,YNEW,XX,YY)
                     NPTS=MIN(NPTS+1,MARR)
                     XARR(NPTS)=XX
                     YARR(NPTS)=YY
                  END IF
                  CALL DSHCRV(XARR,YARR,NPTS,20,0.005*THICK)
                  NPTS=0
                  FTEST=.FALSE.
               END IF

!              save previous value
               XOLD=XP
               YOLD=YP
            END DO
            IF(FTEST)CALL DSHCRV(XARR,YARR,NPTS,20,0.005*THICK)
            IF(NPTS.EQ.MARR)WRITE(*,*)'Exceeding (XARR,YARR) array limits'
         END DO
      END IF

!------------------------------------------------------------------------------
! open map background vector data file
! standard file names:
! arlmap-boundaries, arlmap2-roads, arlmap3-counties, arlmap4-water
!------------------------------------------------------------------------------

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
               WRITE(*,*)'WARNING: map background file not found ',FNAME
               RETURN
            END IF
         END IF
      END IF

!==>read vectors and draw lines

   50 READ(15,'(A)',END=800)POINTS
      READ(POINTS,'(2I5)')KLINE,KPTS
      COLOR=POINTS(11:)

      IF(KPTS.GT.MPTS)THEN
         DEALLOCATE (PLAT,PLON)
         MPTS=MPTS+500
         ALLOCATE (PLAT(MPTS),PLON(MPTS))
      END IF

      READ(15,'(10F6.2)')(PLAT(K),K=1,KPTS)
      READ(15,'(10F7.2)')(PLON(K),K=1,KPTS)
      IF(KPTS.LE.1)GO TO 50

      SELECT CASE (COLOR)
         CASE ('  BOUNDARY') 
            call setcolr(0.,0.,0.)
            th=0.01
         CASE ('  COUNTIES') 
            call setcolr(.8,.8,.8)
            th=0.008
         CASE ('     ROADS') 
            call setcolr(.8,.0,.0)
            th=0.008
         CASE ('    RIVERS') 
            call setcolr(.0,.0,.8)
            th=0.008
         CASE DEFAULT
            th=0.01
      END SELECT

      NPTS=0
      FTEST=.FALSE.
      dline : DO N=1,KPTS
  
         IF(LLGRID)THEN
            CALL GBL2XY(1,1,PLAT(N),PLON(N),XP,YP)
         ELSEIF(KPROJ.EQ.4)THEN
            CALL CYL2XY(PLAT(N),PLON(N),XP,YP)
         ELSE
            CALL CLL2XY(PARMAP,PLAT(N),PLON(N),XP,YP)
         END IF

!        test for too long vectors due to coordinate change (prime meridian)
!        modification introduced 10/18/2007
         IF(N.GT.1)THEN
            IF(ABS(XP-XOLD).GT.0.5*(XU2-XU1).OR.    &
               ABS(YP-YOLD).GT.0.5*(YU2-YU1))THEN
               IF(NPTS.GT.1) CALL SLDCRV(XARR,YARR,NPTS,TH*THICK)
               NPTS=0
               FTEST=.FALSE.
               CYCLE dline 
            END IF
         ELSE
            XOLD=XP
            YOLD=YP
         END IF

         IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
            IF(.NOT.FTEST.AND.N.GT.1.AND.CLIP)THEN
!              add previous position just outside domain
               CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
               CALL MAP2XY(XNEW,YNEW,XX,YY)
               NPTS=MIN(NPTS+1,MARR)
               XARR(NPTS)=XX
               YARR(NPTS)=YY
            END IF
            CALL MAP2XY(XP,YP,XX,YY)
            NPTS=MIN(NPTS+1,MARR)
            XARR(NPTS)=XX
            YARR(NPTS)=YY
            FTEST=.TRUE.

         ELSEIF(FTEST)THEN
            IF(CLIP)THEN
!              add current off domain position to line
               CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
               CALL MAP2XY(XNEW,YNEW,XX,YY)
               NPTS=MIN(NPTS+1,MARR)
               XARR(NPTS)=XX
               YARR(NPTS)=YY
            END IF
!           draw curve if anything in array
            CALL SLDCRV(XARR,YARR,NPTS,TH*THICK)
!           start new line if out of window
            NPTS=0
            FTEST=.FALSE.
         END IF

!        save previous value
         XOLD=XP
         YOLD=YP
      END DO dline

!     complete curve if anything in array
      IF(FTEST) CALL SLDCRV(XARR,YARR,NPTS,TH*THICK)
      IF(NPTS.EQ.MARR)WRITE(*,*)'Exceeding (XARR,YARR) array limits'
      GO TO 50

  800 CLOSE (15)

!==>draw map boundary box in black before returning

      CALL SETCOLR(0.,0.,0.)
      CALL MAP2XY(XU1,YU1,XARR(1),YARR(1))
      CALL MAP2XY(XU1,YU2,XARR(2),YARR(2))
      CALL MAP2XY(XU2,YU2,XARR(3),YARR(3))
      CALL MAP2XY(XU2,YU1,XARR(4),YARR(4))
      CALL MAP2XY(XU1,YU1,XARR(5),YARR(5))
      CALL SLDCRV(XARR,YARR,5,0.015*THICK)

      DEALLOCATE (PLAT,PLON,XARR,YARR)
      RETURN
      END

!==>find clip point on boundary

      SUBROUTINE BOUNDS(X1,Y1,X2,Y2,XP,YP)

      LOGICAL PTEST3, PTEST4
      COMMON /MAPBND/ XU1,XU2,YU1,YU2

      PTEST3(XD,YD)=(XD.GE.X1.AND.XD.LE.X2).OR.(XD.GE.X2.AND.XD.LE.X1)
      PTEST4(XD,YD)=(YD.GE.Y1.AND.YD.LE.Y2).OR.(YD.GE.Y2.AND.YD.LE.Y1)

      DELTA=X2-X1
      IF(DELTA.EQ.0.0)THEN
         XP=X1
         YP=YU1
         IF((YP.GE.Y1.AND.YP.LE.Y2).OR.(YP.GE.Y2.AND.YP.LE.Y1))RETURN
         YP=YU2
         RETURN  
      END IF

      SLOPE=(Y2-Y1)/DELTA   
      ZEROY=Y1-SLOPE*X1

      XP=XU1
      YP=SLOPE*XP+ZEROY
      IF(PTEST3(XP,YP).AND.PTEST4(XP,YP))RETURN

      XP=XU2
      YP=SLOPE*XP+ZEROY
      IF(PTEST3(XP,YP).AND.PTEST4(XP,YP))RETURN

      YP=YU1
      XP=(YP-ZEROY)/SLOPE
      IF(PTEST3(XP,YP).AND.PTEST4(XP,YP))RETURN

      YP=YU2
      XP=(YP-ZEROY)/SLOPE

      END SUBROUTINE bounds
