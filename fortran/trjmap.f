!$$$  SUB PROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TRJMAP           MAIN ROUTINE FOR PLOTTING
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!     basic trajectory plotting program for hysplit4 format files
!     Most  labels and customized features are set in this subroutine
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 08 Jan 1999 (RRD)
!                 28 May 1999 (RRD) - call to mapbdy argument list change
!                 11 Apr 2000 (RRD) - option to draw base map
!                 12 Dec 2000 (RRD) - fortran90 upgrade
!                 03 Jan 2001 (RRD) - vertical coordinate flag
!                 20 Feb 2001 (RRD) - starting height array added
!                 20 Apr 2001 (RRD) - mapbdy argument list
!                 02 Nov 2001 (RRD) - default label changes
!                                   - match rsmc plotting defaults
!                 11 Feb 2002 (RRD) - multiple trajectories in time
!                 11 Mar 2002 (RRD) - various label modifications
!                 20 Jun 2002 (RRD) - process ID suffix
!                 15 Jul 2002 (RRD) - vary traj color by file 
!                 08 Aug 2002 (RRD) - expanded text in maptext file
!                 29 Aug 2002 (RRD) - test for zero points in trajectory
!                 23 Oct 2002 (RRD) - improved split traj plot
!                 17 Dec 2002 (BBS) - origin point multiple files
!                 04 Feb 2003 (RRD) - multiple source location label
!                 17 Jul 2003 (RRD) - revised NOAA in default label
!                 19 Nov 2003 (RRD) - endpoint intervals <1 and diag meteo
!                 04 Dec 2003 (RRD) - converted llint to aintr
!                 01 Mar 2004 (BS)  - trajectory colors,no vertical plot
!                 30 Mar 2004 (RRD) - color allocation test
!                 19 May 2004 (RRD) - multi-traj time labels
!                 15 Dec 2005 (RRD) - g95 compatibility
!                 10 Jan 2005 (RRD) - multi color issues
!                 20 Apr 2005 (BS)  - symbols at time intervals wrt traj start
!                 31 May 2005 (RRD) - terrain height plot
!                 29 Jul 2005 (BS)  - identifier label (idlbl)
!                 28 Dec 2005 (BS)  - STATIONPLOT.CFG input file
!                 29 Sep 2006 (BS)  - enhance LABELS.CFG
!                 30 Oct 2006 (RRD) - sun compiler compatibility
!                 12 Jan 2007 (RRD) - mapbdy argument list change
!                 27 Aug 2007 (BS)  - 3 digits vertical/time plot labeling 
!                 03 Oct 2007 (RRD) - pass through initial terrain height
!                                   - horizontal map boundary test 
!                 05 Nov 2007 (BS)  - cluster labeling
!                 11 Jun 2008 (RRD) - reversed scale for backward trajectory
!                 12 Aug 2008 (RRD) - shapefile mapping option
!                 30 Oct 2008 (RRD) - added starting minutes to label (kk)
!                 06 Feb 2009 (RRD) - fix terrain display for back trajectory
!
! USAGE:  CALL TRJMAP(PMODE,KOLOR,KTHRS,BASE,NTRAJ,OLAT,OLON,OLVL,OTER,KAGL,KV,
!              TKOL,IBYR,IBMO,IBDA,IBHR,JFMO,JFDA,JFHR,MODEL,DIRCTN,MOTION,
!              LDIAG,PBOT,PTOP,TLAT,TLON,THGT,TERR,KFHR,MM,DD,HH,KK,RR,NP,
!              NXP,NYP,PARMAP,ALATB,ALONL,ALATT,ALONR,GFILE,PROCESS,KFN,NRAD,
!  
!   INPUT PARAMETERS:   see below
!   OUTPUT PARAMETERS:  see below
!   INPUT FILES:        unit 10 - STATIONPLOT.CFG
!   OUTPUT FILES:       none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE TRJMAP(PMODE,KOLOR,KTHRS,BASE,NTRAJ,OLAT,OLON,OLVL,OTER,KAGL,KV,    &
                  TKOL,IBYR,IBMO,IBDA,IBHR,JFMO,JFDA,JFHR,MODEL,DIRCTN,        &
                  MOTION,LDIAG,PBOT,PTOP,TLAT,TLON,THGT,TERR,KFHR,MM,DD,HH,    &
                  KK,RR,NP,NXP,NYP,PARMAP,ALATB,ALONL,ALATT,ALONR,GFILE,       &
                  PROCESS,KFN,NRAD,DIST,KLR,SYMBOL,KTER,IDLBL,MAPDEL,IPROJ)

  IMPLICIT NONE

  COMMON /MAPCLP/ XC1,XC2,YC1,YC2

!-------------------------------------------------------------------------------
! special simulation setup parameters read from namelist
!-------------------------------------------------------------------------------
  INTEGER :: HGT
  NAMELIST/SETUP/HGT

!-------------------------------------------------------------------------------
  LOGICAL,      INTENT(IN) :: PMODE       ! indicate portrait mode
  INTEGER,      INTENT(IN) :: KOLOR       ! flag to set B&W(0) or Color(1)
  INTEGER,      INTENT(IN) :: KTHRS       ! time interval trajectory label
  LOGICAL,      INTENT(IN) :: BASE        ! option to draw base map and labels
  INTEGER,      INTENT(IN) :: NTRAJ       ! numb of trajectories in this file
  REAL,         INTENT(IN) :: OLAT(:)     ! starting location 
  REAL,         INTENT(IN) :: OLON(:)     ! starting location 
  REAL,         INTENT(IN) :: OLVL(:)     ! starting height from endpts
  REAL,         INTENT(IN) :: OTER(:)     ! starting terrain height     
  INTEGER,      INTENT(IN) :: KAGL        ! vertical coordinate flag
  INTEGER,      INTENT(IN) :: KV          ! numb of diff levels
  REAL,         INTENT(IN) :: TKOL(:)     ! height color index
  INTEGER,      INTENT(IN) :: IBYR,IBMO   ! starting date/time 
  INTEGER,      INTENT(IN) :: IBDA,IBHR   ! starting date/time 
  INTEGER,      INTENT(IN) :: JFMO,JFDA   ! starting hour meteo 
  INTEGER,      INTENT(IN) :: JFHR        ! starting hour meteo 
  CHARACTER(8), INTENT(IN) :: MODEL       ! meteo data file used
  CHARACTER(8), INTENT(IN) :: DIRCTN      ! integration direction
  CHARACTER(8), INTENT(IN) :: MOTION      ! vertical motion method
  CHARACTER(8), INTENT(IN) :: LDIAG       ! last diagnostic label
  REAL,         INTENT(IN) :: PBOT,PTOP   ! height limit
  REAL,         INTENT(IN) :: TLAT(:,:)   ! endpoint positions
  REAL,         INTENT(IN) :: TLON(:,:)   ! endpoint positions
  REAL,         INTENT(IN) :: THGT(:,:)   ! endpoint positions
  REAL,         INTENT(IN) :: TERR(:,:)   ! endpoint positions
  INTEGER,      INTENT(IN) :: KFHR        ! maximum forecast hour
  INTEGER,      INTENT(IN) :: MM(:,:)     ! endpoint month
  INTEGER,      INTENT(IN) :: DD(:,:)     ! day
  INTEGER,      INTENT(IN) :: HH(:,:)     ! hour
  INTEGER,      INTENT(IN) :: KK(:)       ! minute
  REAL,         INTENT(IN) :: RR(:)       ! age (h)
  INTEGER,      INTENT(IN) :: NP(:,:)     ! numb of endpoints with each traj
  INTEGER,      INTENT(IN) :: NXP,NYP     ! grid for mapping endpoints
  REAL,         INTENT(IN) :: PARMAP(9)   ! defining map projection
  REAL,         INTENT(IN) :: ALATB,ALONL ! map limits
  REAL,         INTENT(IN) :: ALATT,ALONR ! map limits
  CHARACTER(80),INTENT(IN) :: GFILE       ! map background file
  CHARACTER(8), INTENT(IN) :: PROCESS     ! process ID suffix
  INTEGER,      INTENT(IN) :: KFN         ! sequential file number    
  INTEGER,      INTENT(IN) :: NRAD        ! number of radius circles
  REAL,         INTENT(IN) :: DIST        ! distance between circels (km)
  CHARACTER(1), INTENT(IN) :: KLR(:)      ! trajectory colors
  LOGICAL,      INTENT(IN) :: SYMBOL      ! symbol at traj origin 
  INTEGER,      INTENT(IN) :: KTER        ! terrain height index 
  CHARACTER(8), INTENT(IN) :: IDLBL       ! identifier label
  INTEGER,      INTENT(IN) :: MAPDEL      ! enhanced map labelling options
  INTEGER,      INTENT(IN) :: IPROJ       ! map projection flag
!-------------------------------------------------------------------------------

  CHARACTER(80) :: FIELD, TITLE, LABEL, RUNID
  CHARACTER(8)  :: IDENT
  CHARACTER(2)  :: TRJID,NCL
  CHARACTER(3)  :: MONTH
  LOGICAL       :: FTEST, LTEST, SYNOP
  REAL          :: XARR(4),YARR(4)
  REAL          :: RCOL(7),GCOL(7),BCOL(7)
  INTEGER       :: KSYM(7)
  INTEGER       :: KMAX = 0      
  INTEGER       :: KDIR = 0      
  INTEGER       :: CLN,NTPC(50),NTOT,PCT

  INTEGER, PARAMETER :: NLAB = 10
  REAL               :: AINTR(NLAB)

  INTEGER, PARAMETER :: NLIN = 8
  CHARACTER(80)      :: TXTBOXL(NLIN)

  INTEGER :: k,kp,kl,kn,kol,klen,klen1,klen2,klen3,ippp
  REAL    :: ymax,xmin,ymin,xmax,ppp,pmax,pmin,pinc,delt,ptopf,pbotf,inc
  REAL    :: xp,yp,x0,y0,x1,y1,x2,y2,xlonr,slat,slon,qsize,xc1,xc2,yc1,yc2,tsize
  integer :: i,ir,nflg,nlines,lmax
  real    :: xold,yold,xnew,ynew
  integer :: xtraj,zflg

!-------------------------------------------------------------------------------
! color definition table
!         trajectories  map  special
!            red blue green cyan magn yel olive 
  DATA RCOL/ 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.2/
  DATA GCOL/ 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.6/
  DATA BCOL/ 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.8/

! decimal character code for special multiple trajectory plot symbols
! filled: triangle, square, circle
  DATA KSYM/115,    110,    108,    115,110,108,  115/

! potential lat/lon label intervals
  DATA AINTR/45.0,30.0,20.0,15.0,10.0,5.0,2.0,1.0,0.5,0.2/

! maximum dimension along the time axis for multiple plots,direction
  SAVE kmax,kdir

!--------------------------------------------------------------------
! set defaults
  HGT=0
  NFLG=0

!--------------------------------------------------------------------

  XMAX=NXP
  YMAX=NYP
  IF(KMAX.EQ.0)KMAX=MAXVAL(NP(:,2))

!-------------------------------------------------------------------------------
! initialize custom labels if labels.cfg in root directory
!-------------------------------------------------------------------------------

! default values
! TITLE='NATIONAL OCEANIC ATMOSPHERIC ADMINISTRATION&'
  TITLE='NOAA HYSPLIT MODEL&'
  NLINES=0

! optional process ID suffix
  IF(PROCESS(1:2).EQ.'ps')THEN
     LABEL='LABELS.CFG'
  ELSE
     LABEL='LABELS.'//PROCESS
  END IF

  INQUIRE(FILE=LABEL,EXIST=FTEST)
  IF(FTEST)THEN

!    check NTXBOXL and TXBOXL lines
     OPEN(50,FILE=LABEL)
 65  READ(50,*,END=75)IDENT,FIELD
     KLEN1=INDEX(FIELD,'&')
     KLEN2=INDEX(IDENT,'&')

     IF(IDENT(1:KLEN2).EQ.'TXBOXL&')THEN 
        IF(NFLG.EQ.0)THEN
          WRITE(*,*)'ERROR - NTXBOXL line must preceed TXBOXL line(s)'
          WRITE(*,*)'        Fix ',LABEL
          STOP 10
        ELSE
          WRITE(*,*)'ERROR - Too many text box lines. Should be',NLINES
          WRITE(*,*)'        Fix ',LABEL
          STOP 10
        END IF

     ELSEIF(IDENT(1:KLEN2).EQ.'NTXBOXL&')THEN
!       check number of text lines given (NLINES) vs. compiles max (NLIN)
!       remove &
        NFLG=1
        KLEN1=KLEN1-1
        READ(FIELD(1:KLEN1),'(I2)') NLINES
        IF(NLINES.GT.NLIN)THEN
           WRITE(*,*)'ERROR - Text box lines ',NLINES,' exceed compiled ',NLIN
           WRITE(*,*)'        Fix ',LABEL
           STOP 20
        END IF

!       check for NLINES of text
!       NTXBOXL (number of text box lines) must be just before TXTBOXL lines
        DO K=1,NLINES
           READ(50,*,END=70)IDENT,FIELD
           KLEN1=INDEX(FIELD,'&')
           KLEN2=INDEX(IDENT,'&')
           IF(IDENT(1:KLEN2).NE.'TXBOXL&')THEN
              WRITE(*,*)  &
             'ERROR -',NLINES,' text box lines must follow NTXBOXL in ',LABEL
              STOP 1
           END IF        
        END DO
     END IF
     GOTO 65 
 
 70  WRITE(*,*)  &
    'ERROR - Number of textbox lines does not match specified:',NLINES
     WRITE(*,*)'        Fix ',LABEL
     STOP 2

 75  CLOSE(50)
!    read data in LABELS.CFG (correct number of text box lines)
     OPEN(50,FILE=LABEL)

100  READ(50,*,END=200)IDENT,FIELD
     KLEN1=INDEX(FIELD,'&')
     KLEN2=INDEX(IDENT,'&')

     IF(IDENT(1:KLEN2).EQ.'TITLE&')THEN
       TITLE=FIELD(1:KLEN1)

     ELSEIF(IDENT(1:KLEN2).EQ.'VUNIT&'.AND.FIELD(1:KLEN1).EQ.'FEET&')THEN
!       hgt=1 is a flag to label the vertical projection height in feet
!       assuming height in tdump file is in meters
        HGT=1

     ELSEIF(IDENT(1:KLEN2).EQ.'NTXBOXL&')THEN
!       read labels and determine max line length
        LMAX=0
        DO K=1,NLINES
           READ(50,*,END=190)IDENT,FIELD
           KLEN1=INDEX(FIELD,'&')
!          remove & symbol from text line
           TXTBOXL(K)=FIELD(1:KLEN1-1)
           LMAX=MAX(LMAX,KLEN1-1)
        END DO
     END IF
     GOTO 100

190  WRITE(*,*)'ERROR - reading text lines'
     STOP 4
200  CLOSE(50)
  END IF

!-------------------------------------------------------------------------------
! draw map background and set display
!-------------------------------------------------------------------------------

! check for dateline/prime meridian switch
  XLONR=ALONR
  IF(SIGN(1.0,ALONL).NE.SIGN(1.0,ALONR).AND.ALONR.LT.0.0) XLONR=360.0+ALONR

! set interval to have at least 4 lat/lon lines on a map
  DELT=MAX(ABS(ALATT-ALATB)/4.0,ABS(ALONL-XLONR)/4.0,AINTR(NLAB))
  KL=1
  DO WHILE (KL.LT.NLAB.AND.AINTR(KL).GT.DELT)
     KL=KL+1
  END DO

! enhanced labelling options (10 Jan 2007)
  IF(MAPDEL.EQ.0)THEN
!    labels turned off
     DELT=0.0
  ELSEIF(MAPDEL.LT.0)THEN
!    interval set on command line integer tenths
     DELT=ABS(FLOAT(MAPDEL)/10.0)
  ELSE
!    default auto label option mapdel=1
     DELT=AINTR(KL)
  END IF

! map limits in user coordinates (lat/lon conversion)
  XMIN=1.0
  YMIN=1.0

  CALL MAPSET(PMODE,.11,.89,.20,.80,XMIN,XMAX,YMIN,YMAX)
  IF(KOLOR.NE.0) CALL SETCOLR(0.4,0.6,0.8) 
  IF(BASE)THEN
!    map background
     KLEN=LEN_TRIM(GFILE)
     IF(KLEN.GT.13)THEN
        IF(GFILE(klen-13:klen).EQ.'shapefiles.txt')THEN
           CALL SHPBDY(IPROJ,.FALSE.,PARMAP,DELT,1.0,GFILE)
        ELSE
           CALL MAPBDY(IPROJ,.FALSE.,PARMAP,DELT,1.0,GFILE)
        END IF
     ELSE
        CALL MAPBDY(IPROJ,.FALSE.,PARMAP,DELT,1.0,GFILE)
     END IF
  END IF

!-------------------------------------------------------------------------------
! draw optional concentric circles (sub xcircle in conmap)
!-------------------------------------------------------------------------------

! optional concentric circles
  IF(NRAD.GT.0)THEN
     IF(IPROJ.EQ.4)THEN
        CALL CYL2XY(OLAT(1),OLON(1),XP,YP)
     ELSE
        CALL CLL2XY(PARMAP,OLAT(1),OLON(1),XP,YP)
     END IF
     CALL MAP2XY(XP,YP,X1,Y1)

     IF(IPROJ.EQ.4)THEN
        CALL CYL2XY(OLAT(1)+DIST/111.0,OLON(1),XP,YP)
     ELSE
        CALL CLL2XY(PARMAP,OLAT(1)+DIST/111.0,OLON(1),XP,YP)
     END IF
     CALL MAP2XY(XP,YP,X2,Y2)

     DO K=1,NRAD
!       CALL CIRCLE(X1,Y1,(ABS(Y2-Y1)*K),.FALSE.)
        CALL XCIRCLE(X1,Y1,(ABS(Y2-Y1)*K),LTEST)
        IF(LTEST)THEN
           WRITE(LABEL,'(I4,A)')NINT(DIST*K),' km '
           CALL KEKSYMC(X1,(Y1-ABS(Y2-Y1)*K),0.12,LABEL,0.0,8,1)
        END IF
     END DO
  END IF

!-------------------------------------------------------------------------------
! symbols along trajectory - synoptic intervals or wrt start time
!-------------------------------------------------------------------------------
  SYNOP=.TRUE.
  IF(KTHRS.LT.0) SYNOP=.FALSE.

!-------------------------------------------------------------------------------
! If the file STATIONPLOT.CFG exists in the root directory, plot character(s) 
! specified in that file or it plots a symbol if no characters are given.  
! The format is F6.2,1X,F7.2,1X,A :
!-------------------------------------------------------------------------------
  IF(PROCESS(1:2).EQ.'ps')THEN
     LABEL='STATIONPLOT.CFG'
  ELSE
     LABEL='STATIONPLOT.'//PROCESS
  END IF

  INQUIRE(FILE=LABEL,EXIST=FTEST)
  IF(FTEST)THEN
     OPEN(10,FILE='STATIONPLOT.CFG')
  95 READ(10,'(F6.2,1X,F7.2,1X,A)',END=99)SLAT,SLON,LABEL
        IF(LABEL.EQ.'')THEN
         ! octal 154 = 108 (closed circle if no LABEL specified)
           LABEL=CHAR(108)
           KLEN=1
           QSIZE=0.06 
           CALL SETFNT(35)
        ELSE
           KLEN=INDEX(LABEL,' ')-1
           QSIZE=0.10 
           CALL SETFNT(20)
        END IF
      ! write(*,*)slat,slon,label(1:klen),klen
 
!       set color back to black
        CALL SETCOLR(0.,0.,0.)
        IF(IPROJ.EQ.4)THEN
           CALL CYL2XY(SLAT,SLON,XP,YP)
        ELSE
           CALL CLL2XY(PARMAP,SLAT,SLON,XP,YP)
        END IF
        CALL MAP2XY(XP,YP,X1,Y1)
        IF(X1.GE.XC1.AND.X1.LE.XC2.AND.Y1.GE.YC1.AND.Y1.LE.YC2)  &
           CALL KEKSYMC(X1,(Y1-0.06),QSIZE,LABEL,0.0,KLEN,1)	

     GOTO 95

  99 CLOSE(10)
   ! set default font (Helvetica)
     CALL SETFNT(20)

  END IF

!-------------------------------------------------------------------------------
! for mean clusters, calculate percentage trajectories in each cluster
!   ntpc=number of trajectories per cluster
!   ntot=total number of trajectories clustered
! read number of clusters in CLUSLIST_N,
!        where N=#clusters not counting cluster 0
!-------------------------------------------------------------------------------
   IF(IDLBL.EQ.'MERGMEAN')THEN

     IF(NTRAJ.LT.10)THEN
        WRITE(NCL,'(I1)')NTRAJ
        LABEL='CLUSLIST_'//NCL(1:1)
     ELSE
        WRITE(NCL,'(I2)')NTRAJ
        LABEL='CLUSLIST_'//NCL(1:2)
     END IF
     INQUIRE(FILE=LABEL,EXIST=FTEST)
     IF(FTEST)THEN
     ! input Cmean file is for clusters 1-NTRAJ
       ZFLG=0
     ELSE
       WRITE(*,*)'NOTICE: file not found ',LABEL
       XTRAJ=NTRAJ-1
       IF(XTRAJ.LT.10)THEN
          WRITE(NCL,'(I1)')XTRAJ
          LABEL='CLUSLIST_'//NCL(1:1)
       ELSE
          WRITE(NCL,'(I2)')XTRAJ
          LABEL='CLUSLIST_'//NCL(1:2)
       END IF
       INQUIRE(FILE=LABEL,EXIST=FTEST)
       IF(FTEST)THEN
       ! input Cmean is for clusters 0-NTRAJ
         ZFLG=1
       ELSE
         WRITE(*,*)'ERROR: file not found ',LABEL
         STOP 101
       END IF

     END IF

     write(*,*)'Reading ',LABEL
     OPEN(50,FILE=LABEL)
     NTOT=0
 300 READ(50,*,END=310)CLN,NTPC(CLN)
     IF((CLN.EQ.0.AND.ZFLG.EQ.1).OR.CLN.GT.0)THEN
       !write(*,*)'Cluster:',cln,' Number traj in cluster:',ntpc(cln)
        NTOT=NTOT+NTPC(CLN)
     END IF
     DO I=1,NTPC(CLN)-1
        READ(50,*)
     END DO
     GO TO 300
 310 CLOSE(50)
    !write(*,*)'Total number of trajs:',ntot

 END IF

!-------------------------------------------------------------------------------
! draw horizontal-plane trajectory
!-------------------------------------------------------------------------------

! set special character font
  CALL SETFNT(35)

  hloop : DO KN=1,NTRAJ
 
!    plot source location position (octal 110 = 72)
     IF(SYMBOL)THEN
     IF(.NOT.(OLAT(KN).EQ.99.0.AND.OLON(KN).EQ.99.0))THEN
!       set color back to black
        CALL SETCOLR(0.,0.,0.)
        IF(IPROJ.EQ.4)THEN
           CALL CYL2XY(OLAT(KN),OLON(KN),XP,YP)
        ELSE
           CALL CLL2XY(PARMAP,OLAT(KN),OLON(KN),XP,YP)
        END IF
        CALL MAP2XY(XP,YP,X1,Y1)
        CALL KEKSYMC(X1,(Y1-0.06),0.12,CHAR(72),0.0,1,1)
     END IF
     END IF

     IF(NP(KN,1).EQ.0.OR.NP(KN,2).EQ.0) CYCLE hloop

!    set starting position      
     IF(IPROJ.EQ.4)THEN
        CALL CYL2XY(OLAT(KN),OLON(KN),XP,YP)
     ELSE
        CALL CLL2XY(PARMAP,OLAT(KN),OLON(KN),XP,YP)
     END IF
     CALL MAP2XY(XP,YP,X1,Y1)

!    set trajectory colors
     IF(KOLOR.EQ.1)THEN
        IF(KV.EQ.1)THEN
!          rotating trajectory index for colors (only 3 valid)
           KOL=MOD((KN-1)+(KFN-1),3)+1
           KOL=MOD(KOL-1,7)+1
        ELSE
!          color based upon starting height
           KOL=MOD(INT(TKOL(KN))-1,7)+1
        END IF
        CALL SETCOLR(RCOL(KOL),GCOL(KOL),BCOL(KOL))
     ELSEIF(KOLOR.EQ.2)THEN
!       itemized colors
        READ(KLR(KN),'(I1)')KOL
        KOL=MOD(KOL-1,7)+1        ! force 1 to 7
        CALL SETCOLR(RCOL(KOL),GCOL(KOL),BCOL(KOL))
     ELSE
!       set color back to black
        CALL SETCOLR(0.,0.,0.)
        KOL=1				! Stunder 6-22-05
     END IF

!    draw trajectory
     sloop : DO KP=NP(KN,1),NP(KN,2)

!       save previous position in user units
        XOLD=XP
        YOLD=YP

!       convert lat/lon to user projection units
        IF(IPROJ.EQ.4)THEN
           CALL CYL2XY(TLAT(KP,KN),TLON(KP,KN),XP,YP)
        ELSE
           CALL CLL2XY(PARMAP,TLAT(KP,KN),TLON(KP,KN),XP,YP)
        END IF

!       off map test using bounds from mapbdy routine (rrd 10/3/2007)
        IF(XP.LT.XMIN.OR.XP.GT.XMAX.OR.YP.LT.YMIN.OR.YP.GT.YMAX) THEN       
           CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
           CALL MAP2XY(XNEW,YNEW,X2,Y2)
           CALL SLDLIN(X1,Y1,X2,Y2,0.02)
!          terminate trajectory when it runs off the map
           EXIT sloop
 
        ELSE
!          continue drawing trajectory while within map domain
           CALL MAP2XY(XP,YP,X2,Y2)
           CALL SLDLIN(X1,Y1,X2,Y2,0.02)
        END IF

!       special plot symbol 
        IF(KK(KP).EQ.0.AND.IABS(KTHRS).NE.0)THEN
         IF(SYNOP)THEN
!          at synoptic intervals
           IF(HH(KP,KN).EQ.0)THEN
!             special symbol at 0000 GMT
              CALL KEKSYMC(X2,(Y2-0.05),0.09,CHAR(KSYM(KOL)),0.0,1,1)
           ELSEIF(MOD(HH(KP,KN),IABS(KTHRS)).EQ.0)THEN
!             smaller symbol at output interval
              CALL KEKSYMC(X2,(Y2-0.05),0.06,CHAR(KSYM(KOL)),0.0,1,1)
           END IF
         ELSE
!          hours since trajectory start
           IF(MOD(ABS(RR(KP)),24.0).EQ.0)THEN
!             special symbol every 24 hours
              CALL KEKSYMC(X2,(Y2-0.05),0.09,CHAR(KSYM(KOL)+1),0.0,1,1)
           ELSEIF(MOD(ABS(RR(KP)),FLOAT(IABS(KTHRS))).EQ.0)THEN
!             smaller symbol at output interval
              CALL KEKSYMC(X2,(Y2-0.05),0.06,CHAR(KSYM(KOL)),0.0,1,1)
           END IF
         END IF
        END IF

        X1=X2
        Y1=Y2

!    segments loop
     END DO sloop

!    set color back to black
     CALL SETCOLR(0.,0.,0.)

!    -----------------------------------------------------------------
!    for cluster mean plot (MERGMEAN in tdump) label cluster (traj) #
     IF(IDLBL.EQ.'MERGMEAN')THEN

     ! set default font (Helvetica)
       CALL SETFNT(20)

       LABEL='-- (---%)&'
     ! traj 1 is cluster 0, if including cluster 0; else traj 1 is cluster 1
       IF(ZFLG.EQ.0)THEN
         WRITE(LABEL(1:2),'(I2)')KN
       ! percentage of trajectories in this cluster
         PCT=NTPC(KN)*100./NTOT+0.5
       ELSE
         WRITE(LABEL(1:2),'(I2)')KN-1
         PCT=NTPC(KN-1)*100./NTOT+0.5
       END IF
       WRITE(LABEL(5:7),'(I3)')PCT
       KLEN=INDEX(LABEL,'&')-1
       CALL KEKSYMC(X1,Y1,0.10,LABEL(1:KLEN),0.0,KLEN,2)

     ! set special character font
       CALL SETFNT(35)

     END IF
!    -----------------------------------------------------------------

! trajectory loop
  END DO hloop

!-------------------------------------------------------------------------------
! determine levels for vertical axis scale
!-------------------------------------------------------------------------------
  IF(KAGL.NE.4)THEN

  IF(KAGL.EQ.0)THEN
!    pressure coordinate
     PMIN=(INT(PBOT/50.0)+1)*50.0
     PMAX=(INT(PTOP/50.0)-1)*50.0
     PINC=-MAX(50.0,INT((ABS(PMAX-PMIN)/5.0)/50.0)*50.0)
     PMAX=MIN(PMAX,PMIN+4.0*PINC)
  ELSEIF(KAGL.EQ.1)THEN
!    height coordinate 
     PMIN=(INT(PBOT/500.0)-1)*500.0
     IF(KTER.NE.0) PMIN=0.0
     PMAX=(INT(PTOP/500.0)+1)*500.0 
     PINC=MAX(500.0,INT((ABS(PMAX-PMIN)/5.0)/500.0)*500.0)
     PMIN=MAX(0.0,PMIN)
     PMAX=MAX(PMAX,PMIN+4.0*PINC)
   ! set units to feet 
     IF(HGT.EQ.1)THEN
     ! cannot simply set PTOP=PTOP/0.3048, PBOT=... because of INTENT(IN)
       PBOTf=PBOT/0.3048
       PTOPf=PTOP/0.3048
       PMIN=(INT(PBOTf/500.0)-1)*500.0
       IF(KTER.NE.0) PMIN=0.0
       PMAX=(INT(PTOPf/500.0)+1)*500.0 
       PINC=MAX(500.0,INT((ABS(PMAX-PMIN)/5.0)/500.0)*500.0)
       PMIN=MAX(0.0,PMIN)
       PMAX=MAX(PMAX,PMIN+4.0*PINC)
     END IF
  ELSEIF(KAGL.EQ.2)THEN
!    isentropic coordinate
     PMIN=(INT(PBOT/5.0)-1)*5.0
     PMAX=(INT(PTOP/5.0)+1)*5.0
     PINC=MAX(5.0,INT((ABS(PMAX-PMIN)/5.0)/5.0)*5.0)
     PMAX=MAX(PMAX,PMIN+4.0*PINC)
  ELSEIF(KAGL.EQ.3)THEN
!    meteorological data
     IF(ABS(PTOP-PBOT).GT.1000)THEN
        PPP=1000.0
     ELSEIF(ABS(PTOP-PBOT).GT.100.0)THEN
        PPP=100.0
     ELSEIF(ABS(PTOP-PBOT).GT.10.0)THEN
        PPP=10.0
     ELSE
        PPP=1.0
     END IF
     PMIN=(INT(PBOT/PPP)-1)*PPP
     PMAX=(INT(PTOP/PPP)+1)*PPP 
     PINC=MAX(PPP,INT((ABS(PMAX-PMIN)/5.0)/PPP)*PPP)
     PMIN=MAX(0.0,PMIN)
     PMAX=MAX(PMAX,PMIN+4.0*PINC)
  ELSE
     WRITE(*,*)'Invalid vertical plot coordinate: ',KAGL
     STOP 900
  END IF

  END IF

!-------------------------------------------------------------------------------
! draw vertical projection trajectory
!-------------------------------------------------------------------------------

  IF(KAGL.NE.4)THEN

  XMIN=0.0  
  XMAX=KMAX 
  CALL MAPSET(PMODE,.20,.80,.05,.20,XMIN,XMAX,PMIN,PMAX)

! special case of superimposed forward and backward trajectories
  IF(DIRCTN(1:7).EQ.'FORWARD')THEN
     IF(KDIR.EQ.-1)CALL MAPSET(PMODE,.20,.80,.05,.20,XMAX,XMIN,PMIN,PMAX)
  ELSE
     IF(KDIR.EQ. 1)CALL MAPSET(PMODE,.20,.80,.05,.20,XMAX,XMIN,PMIN,PMAX)
  END IF

  vloop : DO KN=1,NTRAJ

     IF(NP(KN,1).EQ.0.OR.NP(KN,2).EQ.0) CYCLE vloop

!    place pointer at trajectory origin
     CALL MAP2XY((NP(KN,1)-1.0),OLVL(KN),X1,Y1)
     IF(HGT.EQ.1) CALL MAP2XY((NP(KN,1)-1.0),OLVL(KN)/0.3048,X1,Y1)

!    set color back to black
     CALL SETCOLR(0.,0.,0.)
!    plot source on vertical
     IF(KTHRS.GE.0) CALL KEKSYMC(X1,(Y1-0.06),0.12,CHAR(72),0.0,1,1)

!    set color for trajectory
     IF(KOLOR.EQ.1)THEN
        IF(KV.EQ.1)THEN
!          rotating trajectory index for colors (only 3 valid)
           KOL=MOD((KN-1)+(KFN-1),3)+1
           KOL=MOD(KOL-1,7)+1
        ELSE
!          color based upon starting height
           KOL=MOD(INT(TKOL(KN))-1,7)+1
        END IF
        CALL SETCOLR(RCOL(KOL),GCOL(KOL),BCOL(KOL))
     ELSEIF(KOLOR.EQ.2)THEN
!       itemized colors
        READ(KLR(KN),'(I1)')KOL
        KOL=MOD(KOL-1,7)+1        ! force 1 to 7
        CALL SETCOLR(RCOL(KOL),GCOL(KOL),BCOL(KOL))
     ELSE
!       set color back to black
        CALL SETCOLR(0.,0.,0.)
     END IF

     DO KP=NP(KN,1),NP(KN,2)
        CALL MAP2XY(FLOAT(KP),THGT(KP,KN),X2,Y2)
        IF(HGT.EQ.1)CALL MAP2XY(FLOAT(KP),THGT(KP,KN)/0.3048,X2,Y2)
        CALL SLDLIN(X1,Y1,X2,Y2,0.02)

!       special plot symbol 
        IF(KK(KP).EQ.0.AND.IABS(KTHRS).NE.0)THEN
           IF(SYNOP)THEN
!             at synoptic intervals
              IF(HH(KP,KN).EQ.0)THEN
!                special symbol at 0000 GMT
                 CALL KEKSYMC(X2,(Y2-0.05),0.09,CHAR(KSYM(KOL)),0.0,1,1)
              ELSEIF(MOD(HH(KP,KN),IABS(KTHRS)).EQ.0)THEN
!                smaller symbol at output interval
                 CALL KEKSYMC(X2,(Y2-0.05),0.06,CHAR(KSYM(KOL)),0.0,1,1)
              END IF
           ELSE
!             hours since trajectory start
              IF(MOD(ABS(RR(KP)),24.0).EQ.0)THEN
!                special symbol every 24 hours
                 CALL KEKSYMC(X2,(Y2-0.05),0.09,CHAR(KSYM(KOL)+1),0.0,1,1)
              ELSEIF(MOD(ABS(RR(KP)),FLOAT(IABS(KTHRS))).EQ.0)THEN
!                smaller symbol at output interval
                 CALL KEKSYMC(X2,(Y2-0.05),0.06,CHAR(KSYM(KOL)),0.0,1,1)
              END IF
           END IF
        END IF

        X1=X2
        Y1=Y2
!    segments loop
     END DO

! trajectory loop
  END DO vloop

  END IF

!-------------------------------------------------------------------------------
! terrain height plot section
!-------------------------------------------------------------------------------

  IF(KAGL.EQ.1.AND.KTER.NE.0)THEN

  KOL=1
  CALL SETCOLR(0.,0.,0.)    ! set color back to black

  XMIN=0.0  
  XMAX=KMAX 
  CALL MAPSET(PMODE,.20,.80,.05,.20,XMIN,XMAX,PMIN,PMAX)

  IF(DIRCTN(1:7).EQ.'FORWARD')THEN
     IF(KDIR.EQ.-1)CALL MAPSET(PMODE,.20,.80,.05,.20,XMIN,XMAX,PMIN,PMAX)
  ELSE
     IF(KDIR.EQ. 1)CALL MAPSET(PMODE,.20,.80,.05,.20,XMAX,XMIN,PMIN,PMAX)
  END IF

  zhts : DO KN=1,NTRAJ

     IF(NP(KN,1).EQ.0.OR.NP(KN,2).EQ.0) CYCLE zhts 

!    place pointer at trajectory origin (see revised code below)
!    CALL MAP2XY((NP(KN,1)-1.0),TERR(NP(KN,1),KN),X1,Y1)
!    IF(HGT.EQ.1)CALL MAP2XY((NP(KN,1)-1.0),TERR(INT(NP(KN,1)/0.3048),KN),X1,Y1)

!    correction for initial terrain height (rrd 10/2/2007)
     CALL MAP2XY((NP(KN,1)-1.0),OTER(KN),X1,Y1)
     IF(HGT.EQ.1)CALL MAP2XY((NP(KN,1)-1.0),(OTER(KN)/0.3048),X1,Y1)

!    plot source on vertical
     IF(KTHRS.GE.0) CALL KEKSYMC(X1,(Y1-0.06),0.12,CHAR(72),0.0,1,1)

     DO KP=NP(KN,1),NP(KN,2)
        CALL MAP2XY(FLOAT(KP),TERR(KP,KN),X2,Y2)
        IF(HGT.EQ.1)CALL MAP2XY(FLOAT(KP),TERR(KP,KN)/0.3048,X2,Y2)
        CALL SLDLIN(X1,Y1,X2,Y2,0.02)

!       special plot symbol 
        IF(KK(KP).EQ.0.AND.IABS(KTHRS).NE.0)THEN
           IF(SYNOP)THEN
!             at synoptic intervals
              IF(HH(KP,KN).EQ.0)THEN
!                special symbol at 0000 GMT
                 CALL KEKSYMC(X2,(Y2-0.05),0.09,CHAR(KSYM(KOL)),0.0,1,1)
              ELSEIF(MOD(HH(KP,KN),IABS(KTHRS)).EQ.0)THEN
!                smaller symbol at output interval
                 CALL KEKSYMC(X2,(Y2-0.05),0.06,CHAR(KSYM(KOL)),0.0,1,1)
              END IF
           ELSE
!             hours since trajectory start
              IF(MOD(ABS(RR(KP)),24.0).EQ.0)THEN
!                special symbol every 24 hours
                 CALL KEKSYMC(X2,(Y2-0.05),0.09,CHAR(KSYM(KOL)+1),0.0,1,1)
              ELSEIF(MOD(ABS(RR(KP)),FLOAT(IABS(KTHRS))).EQ.0)THEN
!                smaller symbol at output interval
                 CALL KEKSYMC(X2,(Y2-0.05),0.06,CHAR(KSYM(KOL)),0.0,1,1)
              END IF
           END IF
        END IF

        X1=X2
        Y1=Y2
!    segments loop
     END DO

! trajectory loop
  END DO zhts 

  END IF

!-------------------------------------------------------------------------------

! set color back to black
  CALL SETCOLR(0.,0.,0.)
! set default font (Helvetica)
  CALL SETFNT(20)
! return if only drawing trajectories
  IF(.NOT.BASE)RETURN

!-------------------------------------------------------------------------------
! vertical projection labels (invert axis)
!-------------------------------------------------------------------------------

  IF(KAGL.NE.4)THEN

  XMIN=0.0
  XMAX=1.0
  CALL MAPSET(PMODE,.11,.89,.05,.20,XMIN,XMAX,PMIN,PMAX)

! projection box
  CALL MAP2XY(XMIN,PMAX,XARR(1),YARR(1))
  CALL MAP2XY(XMIN,PMIN,XARR(2),YARR(2))
  CALL MAP2XY(XMAX,PMIN,XARR(3),YARR(3))
  CALL MAP2XY(XMAX,PMAX,XARR(4),YARR(4))
  CALL SLDCRV(XARR,YARR,4,0.015)

! reference lines and value
  DO IPPP=INT(PMIN+PINC),INT(PMAX-PINC),INT(PINC)
     PPP=FLOAT(IPPP)
     WRITE(LABEL,'(I5)')INT(PPP)
     CALL MAP2XY(0.10,PPP,X1,Y1)
     CALL MAP2XY(0.90,PPP,X2,Y2)
     CALL KEKSYMC((X2+0.06),(Y2-0.05),0.10,LABEL,0.0,5,0)
     CALL DSHLIN(X1,Y1,X2,Y2,6,0.0)
  END DO

! starting height value
  PPP=-1.0
  zloop : DO KN=1,NTRAJ  
     IF(NP(KN,1).EQ.0.OR.NP(KN,2).EQ.0) CYCLE zloop
     IF(NP(KN,1).EQ.1)THEN
        CALL MAP2XY(0.02,OLVL(KN),XP,YP)
        IF(HGT.EQ.1)CALL MAP2XY(0.02,OLVL(KN)/0.3048,XP,YP)
        WRITE(LABEL,'(I5)')INT(OLVL(KN))
        IF(HGT.EQ.1)WRITE(LABEL,'(I5)')INT(OLVL(KN)/0.3048+0.5)
        IF(ABS(PPP-YP).GT.0.10)THEN
           CALL KEKSYMC((XP-0.04),(YP-0.05),0.10,LABEL,0.0,5,0)
           PPP=YP
        END IF
     END IF
  END DO zloop

  END IF

!-------------------------------------------------------------------------------
! time projection box
!-------------------------------------------------------------------------------

  IF(KAGL.NE.4)THEN

  PMIN=0.0
  PMAX=1.0
  CALL MAPSET(PMODE,.11,.89,0.0,.05,XMIN,XMAX,PMIN,PMAX)

! projection box
  CALL MAP2XY(XMIN,PMAX,XARR(1),YARR(1))
  CALL MAP2XY(XMIN,PMIN,XARR(2),YARR(2))
  CALL MAP2XY(XMAX,PMIN,XARR(3),YARR(3))
  CALL MAP2XY(XMAX,PMAX,XARR(4),YARR(4))
  CALL SLDCRV(XARR,YARR,4,0.015)

  END IF

!-------------------------------------------------------------------------------
! time projection labels based upon the first trajectory        
!-------------------------------------------------------------------------------

  IF(KAGL.NE.4)THEN

  KL=0
  XMIN=0.0  
  XMAX=KMAX
  CALL MAPSET(PMODE,.20,.80,0.0,.05,XMIN,XMAX,PMIN,PMAX)

  DO KN=1,NTRAJ
  DO KP=NP(KN,1),NP(KN,2)

     IF(KK(KP).EQ.0.AND.IABS(KTHRS).NE.0)THEN
      IF(SYNOP)THEN
        IF(MOD(HH(KP,KN),IABS(KTHRS)).EQ.0.AND.KP.GT.KL)THEN
           CALL MAP2XY(FLOAT(KP),0.80,X1,Y1)
           CALL MAP2XY(FLOAT(KP),0.90,X2,Y2)
           CALL MAP2XY(FLOAT(KP),PMAX,XP,YP)
           CALL SLDLIN(X2,Y2,XP,YP,0.0)

           CALL MAP2XY(FLOAT(KP),0.50,X0,Y0)
           WRITE(LABEL,'(I2.2)')HH(KP,KN)    
           CALL KEKSYMC(X0,Y0,0.09,LABEL,0.0,2,1)

           IF(HH(KP,KN).EQ.0)THEN
              CALL SLDLIN(X1,Y1,XP,YP,0.0)

              CALL MAP2XY(FLOAT(KP),0.10,X0,Y0)
              WRITE(LABEL,'(I2.2,A1,I2.2)')MM(KP,KN),'/',DD(KP,KN)    
              CALL KEKSYMC(X0,Y0,0.09,LABEL,0.0,5,1)
           END IF
           KL=MAX(KL,KP)
        END IF
      ELSE
      ! with respect to traj start 
        IR=ABS(RR(KP))+0.5
        IF(MOD(IABS(IR),IABS(KTHRS)).EQ.0)THEN
           CALL MAP2XY(FLOAT(KP),0.80,X1,Y1)
           CALL MAP2XY(FLOAT(KP),0.90,X2,Y2)
           CALL MAP2XY(FLOAT(KP),PMAX,XP,YP)
           CALL SLDLIN(X2,Y2,XP,YP,0.0)

           CALL MAP2XY(FLOAT(KP),0.50,X0,Y0)
           IF(IR.LT.100)THEN
             WRITE(LABEL,'(I2.2)')IR
             CALL KEKSYMC(X0,Y0,0.09,LABEL,0.0,2,1)
           ELSE
             WRITE(LABEL,'(I3.3)')IR
             CALL KEKSYMC(X0,Y0,0.09,LABEL,0.0,3,1)
           END IF
        END IF
      END IF
     END IF

  END DO
  END DO

  END IF

!-------------------------------------------------------------------------------
! top of panel labels
!-------------------------------------------------------------------------------

! set absolute frame for text
  CALL MAPSET(PMODE,.05,.95,.10,.90,0.,1.,0.,1.)

  LABEL=TITLE
  KLEN=INDEX(LABEL,'&')-1
  CALL MAP2XY(.50,.98,XP,YP)
  CALL KEKSYMC(XP,YP,0.12,LABEL,0.0,KLEN,1)

! trajectory cluster or merged-file labels, else regular labels
  IF(IDLBL.EQ.'MERGMEAN'.OR.IDLBL.EQ.'MERGLIST')THEN
     IF(DIRCTN(1:7).EQ.'FORWARD')THEN
        IF(IDLBL.EQ.'MERGMEAN')THEN
           LABEL='----- forward trajectories&'
           WRITE(LABEL(1:5),'(I5)')NTOT
        ELSEIF(IDLBL.EQ.'MERGLIST')THEN
           LABEL='----- forward trajectories starting at various times&'
           WRITE(LABEL(1:5),'(I5)')NTRAJ
        END IF
     ELSE
        IF(IDLBL.EQ.'MERGMEAN')THEN
           LABEL='----- backward trajectories&'
           WRITE(LABEL(1:5),'(I5)')NTOT
        ELSEIF(IDLBL.EQ.'MERGLIST')THEN
           LABEL='----- backward trajectories ending at various times&'
           WRITE(LABEL(1:5),'(I5)')NTRAJ
        END IF
     END IF
  ELSE
   ! add tests here for number of trajs for proper label 12-10-96
     IF(DIRCTN(1:7).EQ.'FORWARD')THEN
        IF(NTRAJ.GE.2) THEN
           LABEL='Forward trajectories starting at ---- UTC -- --- --&'
        ELSE
           LABEL='  Forward trajectory starting at ---- UTC -- --- --&'
        END IF
     ELSE
        IF(NTRAJ.GE.2) THEN
           LABEL=' Backward trajectories ending at ---- UTC -- --- --&'
        ELSE
           LABEL='   Backward trajectory ending at ---- UTC -- --- --&'
        END IF
     END IF

     CALL MONSET(IBMO,MONTH,1)
     WRITE(LABEL(34:37),'(2I2.2)')IBHR,KK(SIZE(KK))
     WRITE(LABEL(43:48),'(I2.2,1X,A3)')IBDA,MONTH
     WRITE(LABEL(50:51),'(I2.2)')IBYR
  END IF

  KLEN=INDEX(LABEL,'&')-1
  CALL MAP2XY(.50,.94,XP,YP)
  CALL KEKSYMC(XP,YP,0.14,LABEL,0.0,KLEN,1)

! forecast hour determines label
  CALL MONSET(JFMO,MONTH,1)
  IF(KFHR.GT.12)THEN
     RUNID='-- UTC -- ---'//MODEL//'  Forecast Initialization&'
     WRITE(RUNID(1:2),'(I2.2)')JFHR
     WRITE(RUNID(8:13),'(I2.2,1X,A3)')JFDA,MONTH
  ELSE
     RUNID=MODEL//' Meteorological Data&'
  END IF

  KLEN=INDEX(RUNID,'&')-1
  CALL MAP2XY(.50,.90,XP,YP)
  CALL KEKSYMC(XP,YP,0.12,RUNID,0.0,KLEN,1)

!-------------------------------------------------------------------------------
! left edge labels
!-------------------------------------------------------------------------------

  LABEL='Source&'
  KLEN=INDEX(LABEL,'&')-1
! CALL MAP2XY(.05,.20,XP,YP)
  CALL MAP2XY(.05,.33,XP,YP)
  CALL KEKSYMC(XP,YP,0.12,LABEL,90.0,KLEN,0)

! enhanced source symbol
  IF(KTHRS.GE.0)THEN
     CALL MAP2XY(.05,.45,XP,YP)
     CALL SETFNT(35)
     CALL KEKSYMC(XP,YP,0.12,CHAR(72),90.0,1,0)
     CALL SETFNT(20)
  END IF

  IF(MAXVAL(OLAT).EQ.MINVAL(OLAT).AND.MAXVAL(OLON).EQ.MINVAL(OLON))THEN
!    all sources at the same location
     LABEL='at    0.00 N    0.00 E&'
     WRITE(LABEL(5:10),'(F6.2)')ABS(OLAT(1))
     IF(OLAT(1).LT.0.0)WRITE(LABEL(12:12),'(A1)')'S'
     WRITE(LABEL(14:20),'(F7.2)')ABS(OLON(1))
     IF(OLON(1).LT.0.0)WRITE(LABEL(22:22),'(A1)')'W'
  ELSE
     LABEL='at multiple locations&'
  END IF
  KLEN=INDEX(LABEL,'&')-1
  CALL MAP2XY(.05,.50,XP,YP)
  CALL KEKSYMC(XP,YP,0.12,LABEL,90.0,KLEN,0)

!-------------------------------------------------------------------------------
! right edge labels (3/1/2002)
!-------------------------------------------------------------------------------

! LABEL='Vertical Motion Method - '//MOTION//'&'
! KLEN=INDEX(LABEL,'&')-1
! CALL MAP2XY(.95,.50,XP,YP)
! CALL KEKSYMC(XP,YP,0.12,LABEL,-90.0,KLEN,1)

!-------------------------------------------------------------------------------
! vertical projection label
!-------------------------------------------------------------------------------

  CALL MAP2XY(.05,.04,XP,YP)
  IF(KAGL.EQ.0)THEN
     CALL KEKSYMC(XP,YP,0.12,'hPa',90.0,3,1)
  ELSEIF(KAGL.EQ.1)THEN
     IF(KTER.EQ.0)THEN
        IF(HGT.EQ.1) THEN
           CALL KEKSYMC(XP,YP,0.12,'Feet AGL',90.0,8,1)
        ELSE
           CALL KEKSYMC(XP,YP,0.12,'Meters AGL',90.0,10,1)
        END IF
     ELSE
        IF(HGT.EQ.1) THEN
           CALL KEKSYMC(XP,YP,0.12,'Feet MSL',90.0,8,1)
        ELSE
           CALL KEKSYMC(XP,YP,0.12,'Meters MSL',90.0,10,1)
        END IF
     END IF
  ELSEIF(KAGL.EQ.2)THEN
     CALL KEKSYMC(XP,YP,0.12,'Theta',90.0,5,1)
  ELSEIF(KAGL.EQ.3)THEN
     CALL KEKSYMC(XP,YP,0.12,LDIAG,90.0,8,1)
  END IF

!-------------------------------------------------------------------------------
! optional RSMC text box with supplemental information
!-------------------------------------------------------------------------------

  PMIN=0.0
  PMAX=7.0
  XMIN=0.0
  XMAX=1.0
  CALL MAPSET(PMODE,.11,.89,-0.11,0.0,XMIN,XMAX,PMIN,PMAX)

! optional process ID suffix
  IF(PROCESS(1:2).EQ.'ps')THEN
     LABEL='MAPTEXT.CFG'
  ELSE
     LABEL='MAPTEXT.'//PROCESS
  END IF

  INQUIRE(FILE=LABEL,EXIST=FTEST)
  IF(FTEST)THEN

!    projection box
     CALL MAP2XY(XMIN,PMAX,XARR(1),YARR(1))
     CALL MAP2XY(XMIN,PMIN,XARR(2),YARR(2))
     CALL MAP2XY(XMAX,PMIN,XARR(3),YARR(3))
     CALL MAP2XY(XMAX,PMAX,XARR(4),YARR(4))
     CALL SLDCRV(XARR,YARR,4,0.015)

     OPEN(50,FILE=LABEL)
     DO K=1,15
        READ(50,'(A80)',END=500)LABEL
        IF(K.EQ.1.OR.(K.GE.3.AND.K.LE.5).OR.K.EQ.9.OR.K.EQ.15) THEN
           PMAX=PMAX-1.1
           CALL MAP2XY(0.05,PMAX,XP,YP)
           IF(K.EQ.4)LABEL=' '
           CALL KEKSYMC(XP,YP,0.09,LABEL,0.0,80,0)
        END IF
     END DO
 500 CLOSE (50)

  END IF

!-------------------------------------------------------------------------------
! optional alternate text box (text in LABELS.CFG)
!-------------------------------------------------------------------------------
  IF(NLINES.GT.0)THEN

! black
  CALL SETCOLR(0.,0.,0.)

  PMIN=0.0
  PMAX=7.0
  XMIN=0.0
  XMAX=1.0
  CALL MAPSET(PMODE,.11,.89,-0.11,0.0,XMIN,XMAX,PMIN,PMAX)

!    projection box
     CALL MAP2XY(XMIN,PMAX,XARR(1),YARR(1))
     CALL MAP2XY(XMIN,PMIN,XARR(2),YARR(2))
     CALL MAP2XY(XMAX,PMIN,XARR(3),YARR(3))
     CALL MAP2XY(XMAX,PMAX,XARR(4),YARR(4))
     CALL SLDCRV(XARR,YARR,4,0.015)

   ! vertical text spacing increment
     INC=(PMAX-PMIN)/(NLINES+0.5)

     DO K=1,NLINES
        PMAX=PMAX-INC
        CALL MAP2XY(0.05,PMAX,XP,YP)
        LABEL=TXTBOXL(K)
      ! ---
      ! text size (0.09 is good for up to 6 lines)
        TSIZE=0.09
        IF(NLINES.GT.6)THEN
           TSIZE=0.07
        ELSEIF(NLINES.LE.4)THEN
         ! Alaska volcano trajs
           IF(LMAX.LE.67)TSIZE=0.13
        END IF
      ! small font if "http" starting in first 15 chars (Alaska trajs)
      !    web runs URL near end of label
        DO I=1,10
           IF(LABEL(I:I+3).EQ.'http')TSIZE=AMAX1(TSIZE-0.04,0.06)
        END DO
      ! small font if date-time stamp starts "Created" (Alaska trajs)
        IF(LABEL(1:7).EQ.'Created')TSIZE=AMAX1(TSIZE-0.04,0.06)
       !write(*,*)'TSIZE:',tsize
      ! ---
        CALL KEKSYMC(XP,YP,TSIZE,LABEL,0.0,80,0)
     END DO

  END IF

!---------------------------------------------------------------
! set flag for special case of superimposed fwd and back traj
! only one pair permitted, multiple fwd/back mix not supported
!---------------------------------------------------------------

  IF(DIRCTN(1:7).EQ.'FORWARD')THEN
     KDIR=1
  ELSE
     KDIR=-1
  END IF

END SUBROUTINE trjmap

