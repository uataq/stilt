!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONTUR           MAP CONTOUR DRAW
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-12-11
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            DETERMINES THE CONTOURS AND INTERVAL ACCORDING TO THE
!            KONTUR FLAG VARIABLE THEN CALLS COLOR FILL ROUTINES
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 06 Feb 1998 (RRD)
!                 08 Dec 1998 (RRD) - converted NCAR graphics to psplot
!                 02 Nov 2000 (BBS) - ctmax fix
!                 21 Nov 2000 (RRD) - fortran90 upgrade
!                 06 Feb 2001 (RRD) - fixed contour option
!                 02 Nov 2001 (RRD) - symbol plot instead of contour
!                 16 Nov 2001 (RRD) - revised fixed contour method
!                 10 Dec 2001 (GDR) - added Arcview option
!                 05 Jan 2002 (RRD) - gis dimension fix
!                 21 Jun 2002 (RRD) - local min for display
!                 26 Sep 2002 (RRD) - GIS format modification
!                 05 May 2003 (RRD) - avoid conturing when below minimum
!                 10 Feb 2003 (RRD) - option to drop contour outline
!                 26 Oct 2004 (RJP) - added dbf capability
!                 23 Dec 2004 (RRD) - expanded contour format in generate
!                 03 Jun 2005 (RRD) - linear contour option
!                 06 Mar 2006 (GDR) - Google Earth option added
!                                   - added source lat/lon and heights
!                 01 Dec 2006 (GDR) - valid time to GE label field
!                 12 Jan 2007 (RRD) - cylindrical equidistant
!                 01 Feb 2007 (RRD) - added minutes field
!                 23 Feb 2007 (GDR) - move Google Earth routines to geplot.f
!                 14 Sep 2007 (GDR) - agl/msl issues
!                 09 Oct 2007 (GDR) - pass IYR to geplot
!                 22 Oct 2007 (BS)  - -k3 option
!                 02 Jan 2008 (GDR) - uncommented section for duplicate contours so all contours are checked
!                 27 Jan 2009 (RRD) - force colors from command line
!
! USAGE: CALL CONTUR( see below )
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            GISTMP, GISOUT, GISATT
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CONTUR(KPROJ,LEVEL1,LEVEL2,PTYPE,PARMAP,KONTUR,NXP,NYP,XCON,CMIN,   &
                  CMAX,TMAX,TMIN,CVAL,NVAL,KOLOR,NGPT,SYM,IYR,IMO,IDA,IHR,IMN, &
                  JYR,JMO,JDA,JHR,JMN,JFH,OLAT,OLON,OLVL1,OLVL2,GOPEN,KMAP,    &
                  KAVG,MAPN,PROCESS)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER, INTENT(IN)     :: kproj         ! projection                    
  INTEGER, INTENT(IN)     :: level1,level2 ! range of display levels
  CHARACTER(4),INTENT(IN) :: ptype         ! ident of pollutant displayed

  REAL,    INTENT(IN)    :: parmap (9)    ! conformal map projection parameters
  INTEGER, INTENT(IN)    :: kontur        ! flag auto or fixed contours
  INTEGER, INTENT(IN)    :: nxp,nyp       ! size of the concentration array
  REAL,    INTENT(IN)    :: xcon (:,:)    ! concentration array
  INTEGER, INTENT(IN)    :: nval          ! number of contours
  INTEGER, INTENT(IN)    :: kolor         ! flag for B&W=0 or color=1

  REAL,    INTENT(IN)    :: cmin,cmax     ! global min max for contours
  REAL,    INTENT(OUT)   :: tmax,tmin     ! local max min for display
  REAL,    INTENT(INOUT) :: cval(:)       ! contour array
  INTEGER, INTENT(IN)    :: ngpt          ! number of grid points per input 
  REAL,    INTENT(IN)    :: sym           ! symbol plot threshold value (>=0)  

  CHARACTER(8),INTENT(IN):: process       ! process ID suffix
  INTEGER, INTENT(IN)    :: imo,ida,ihr   ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: imn,iyr       ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: jmo,jda,jhr   ! stop date/time of the sample
  INTEGER, INTENT(IN)    :: jmn,jyr       ! stop date/time of the sample
  INTEGER, INTENT(IN)    :: jfh           ! stop time of the sample
  INTEGER, INTENT(IN)    :: mapn          ! current plot map number        
  REAL,    INTENT(IN)    :: olat(:)       ! source latitude
  REAL,    INTENT(IN)    :: olon(:)       ! source longitude
  REAL,    INTENT(IN)    :: olvl1,olvl2   ! release height range
  INTEGER, INTENT(IN)    :: kmap,kavg     ! flag to indicate concen or deposit
  LOGICAL, INTENT(INOUT) :: gopen         ! Google Earth output file open

!-------------------------------------------------------------------------------

  REAL     :: COLOR(3,6), COLOR1(3,6), COLOR2(3,6), TVAL(0:10), GRYLEV(6)
  REAL     :: COLORX(3,6), GRYLEV1(6), GRYLEV2(6), XARR(5), YARR(5)

  INTEGER  :: n,ii,jj,nexp,kret,k,kge 
  INTEGER  :: ispec,ioffp,ilegg,ilabb,lscal,ldash,nlbll,nhii,ndeccn 
  REAL     :: cint,spval,hgtlab,ctmax,ctmin 
  REAL     :: x1,y1,x2,y2,rcon,xmax,ymax,xlen,ylen 

! Rolph added for GIS
  REAL          :: xp,yp,xr,yr,clat,clon,cv,xra,yra,cvl,clvl(5)
  CHARACTER(80) :: gisout,gisatt,gistmp,gislbl
  INTEGER       :: igis,icnt,nn,nc,np,maxct,maxpt,iclvl,iprv,jlen
! number of contour points
  REAL     :: XPNT(9999), YPNT(9999)
  REAL     :: XLAT(9999), XLON(9999),AREA
  LOGICAL  :: CLSPLY
!-------------------------------------------------------------------------------

! conrec contouring common block
  COMMON /CONPAR/ ISPEC,IOFFP,SPVAL,ILEGG,ILABB,NHII,NDECCN,NLBLL,LSCAL,       &
                  LDASH,HGTLAB

! colors from command line
  COMMON /RGBVAL/ COLORX

! Rolph added for GIS
! GIS option common block
  COMMON /GIS/GISOUT,GISATT,GISTMP,IGIS

! RGB contours: red  green  blue
! order: white,cyan,green,blue,yellow,white
  DATA COLOR1/  1.0,  1.0,  1.0,                                           &
                0.0,  1.0,  1.0,                                           &
                0.0,  1.0,  0.0,                                           &
                0.0,  0.0,  1.0,                                           &
                1.0,  1.0,  0.0,                                           &
                1.0,  1.0,  1.0/

! order: white,white,yellow,orange,red,white    
  DATA COLOR2/  1.0,  1.0,  1.0,                                           &
                0.9,  0.9,  0.9,                                           &
                1.0,  1.0,  0.0,                                           &
                1.0,  0.5,  0.0,                                           &
                1.0,  0.0,  0.0,                                           &
                1.0,  1.0,  1.0/

! set equivalent gray scale values
  DATA GRYLEV1/ 1.0, 0.9, 0.7, 0.5, 0.3, 1.0/
  DATA GRYLEV2/ 1.0, 1.0, 0.7, 0.5, 0.3, 1.0/

!-------------------------------------------------------------------------------
  INTERFACE
!-------------------------------------------------------------------------------
  SUBROUTINE GEPLOT(KGE,KPROJ,LEVEL1,LEVEL2,PTYPE,PARMAP,NXP,NYP,   &
                  TMAX,TMIN,XMAX,YMAX,CVAL,NVAL,NGPT,               &
                  IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,          &
                  OLAT,OLON,OLVL1,OLVL2,GOPEN,KMAP,                 &
                  MAPN,PROCESS,XLEN,YLEN)

  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: kge           ! 0=normal, 1=initialize, 2=zero concentrations
  INTEGER, INTENT(IN)    :: kproj         ! projection
  INTEGER, INTENT(IN)    :: level1,level2 ! range of display levels
  CHARACTER(4),INTENT(IN):: ptype         ! ident of pollutant displayed
  REAL,    INTENT(IN)    :: parmap (9)    ! conformal map projection parameters
  INTEGER, INTENT(IN)    :: nxp,nyp       ! size of the concentration array
  INTEGER, INTENT(IN)    :: nval          ! number of contours
  INTEGER, INTENT(IN)    :: ngpt          ! number of grid points per input 
  REAL,    INTENT(IN)    :: xmax,ymax     ! minimum and maximum concentrations
  REAL,    INTENT(IN)    :: tmax,tmin     ! local time max for display
  REAL,    INTENT(INOUT) :: cval(:)       ! contour array
  INTEGER, INTENT(IN)    :: imo,ida,ihr   ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: imn,iyr       ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: jmo,jda,jhr   ! stop date/time of the sample
  INTEGER, INTENT(IN)    :: jmn,jyr       ! stop date/time of the sample
  INTEGER, INTENT(IN)    :: mapn          ! current plotting map number
  REAL,    INTENT(IN)    :: olat (:)      ! starting location
  REAL,    INTENT(IN)    :: olon (:)      ! starting location
  REAL,    INTENT(IN)    :: olvl1,olvl2   ! release height range
  INTEGER, INTENT(IN)    :: kmap          ! flag to indicate concen or deposit
  LOGICAL, INTENT(INOUT) :: gopen         ! Google Earth output file open
  CHARACTER(8),INTENT(IN):: process       ! process ID suffix
  REAL,    INTENT(IN)    :: xlen,ylen     ! size of plot in paper units

  END SUBROUTINE geplot
!-------------------------------------------------------------------------------
  END INTERFACE

!-------------------------------------------------------------------------------
! set contouring defaults
!-------------------------------------------------------------------------------

  IF(NVAL.EQ.0)RETURN

  ISPEC=1
  ILEGG=0
  ILABB=0
  NHII=-1

!-------------------------------------------------------------------------------
! The optional symbol plot above threshold is only used by special versions of 
! concplot (e.g. volcplot) to create blank maps or concentration regions 
! represented by a special character.  This option is not currently enabled 
! within concplot.
!-------------------------------------------------------------------------------

  IF(SYM.GT.0.0)THEN
     IF(KOLOR.EQ.1) CALL SETCOLR(1.,0.,1.)
     CALL SETFNT(35)
     DO II=1,NXP
     DO JJ=1,NYP
        RCON=XCON(II,JJ)
        IF(RCON.GE.SYM.AND.RCON.GT.0.0)THEN
           CALL MAP2XY(FLOAT(II),FLOAT(JJ),X1,Y1)
!          filled square=110;  open square=111
           CALL KEKSYMC(X1,Y1,0.025,CHAR(110),0.0,1,1)
        END IF
     END DO
     END DO
     IF(KOLOR.EQ.1) CALL SETCOLR(0.,0.,0.)
     CALL SETFNT(20)
     RETURN
  ELSEIF(SYM.EQ.0.0)THEN
!    blank map option
     RETURN
  END IF

! Open Google Earth format XML file (.kml) and write defaults
  IF(IGIS.EQ.3.AND..NOT.GOPEN)THEN
     KGE=1        
     CALL GEPLOT(KGE,KPROJ,LEVEL1,LEVEL2,PTYPE,PARMAP,NXP,NYP,   &
               TMAX,TMIN,XMAX,YMAX,CVAL,NVAL,NGPT,               &
               IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,          &
               OLAT,OLON,OLVL1,OLVL2,GOPEN,KMAP,                 &
               MAPN,PROCESS,XLEN,YLEN)
   END IF

!-------------------------------------------------------------------------------
! find max and mins for labels or to set contours
!-------------------------------------------------------------------------------

  CTMAX=0.0
  CTMIN=1.0E+25
  DO II=1,NXP
  DO JJ=1,NYP
     RCON=XCON(II,JJ)
     IF(RCON.GT.0.0)THEN
        CTMIN=AMIN1(RCON,CTMIN)
        IF(RCON.GT.CTMAX)THEN
           CTMAX=RCON
           XMAX=II
           YMAX=JJ
        END IF
     END IF
  END DO
  END DO

! return local time max for display purposes
  ! tmax tmin are max min of conc on conformal grid 
  ! and so will likely differ from max min in cdump file (lat-lon grid)
  TMAX=CTMAX
  TMIN=CTMIN

! Avoid contouring nothing, however create Google Earth label before exiting
  IF(CTMAX.EQ.0.0)THEN
   IF(IGIS.EQ.3)THEN
       KGE=2
       CALL GEPLOT(KGE,KPROJ,LEVEL1,LEVEL2,PTYPE,PARMAP,NXP,NYP,   &
                  TMAX,TMIN,XMAX,YMAX,CVAL,NVAL,NGPT,              &
                  IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,         &
                  OLAT,OLON,OLVL1,OLVL2,GOPEN,KMAP,                &
                  MAPN,PROCESS,XLEN,YLEN)
   END IF
   RETURN
  END IF

!    plot grid size box at maximum value location

! when the dynamic contour option is set, the max/min is calculated 
! above for each map and over-rides the global value determined for
! all time periods in the conndx subroutine

  IF(KONTUR.EQ.1.OR.KONTUR.EQ.3.OR.KONTUR.EQ.51)THEN
!    with the fixed option the contours are the same for all maps
     CTMAX=CMAX
     CTMIN=CMIN
  END IF

!-------------------------------------------------------------------------------
! set contour intervals: always the first time but sometimes
! each time depending upon the value of KONTUR
!-------------------------------------------------------------------------------

! test for contour values on command line 
  IF(KONTUR.EQ.4)THEN

!    check if contour order needs to be reversed (nval=max value)
     IF(CVAL(1).GT.CVAL(NVAL))THEN
        DO N=1,NVAL
           TVAL(N)=CVAL(NVAL-N+1)
        END DO
        DO N=1,NVAL
           CVAL(N)=TVAL(N)
        END DO
     END IF

!    compute contour increment
     CINT=CVAL(NVAL)-CVAL(NVAL-1)

!    for threshold runs, set minimum contour to minimum concentration
     IF(KMAP.EQ.4)CVAL(1)=CTMIN

!    fill in remaining values for offset contour
     DO N=1,NVAL
        TVAL(N-1)=CVAL(N)
     END DO
!    TVAL(NVAL)=CVAL(NVAL)+CINT
     TVAL(NVAL)=CTMAX
     TVAL(NVAL+1)=TVAL(NVAL)+CINT

!    check to see if 2 contours have the same value...if so subtract 1000th from it
     DO N=1,NVAL-1
        IF(TVAL(N+1).EQ.TVAL(N))THEN
          TVAL(N)=TVAL(N)-0.001*TVAL(N)
        END IF
        IF(CVAL(N+1).EQ.CVAL(N))THEN
          CVAL(N)=CVAL(N)-0.001*CVAL(N)
        END IF
     END DO

!    check if min contur > max conc   
     IF(KMAP.NE.4.AND.TVAL(1).GT.CTMAX)RETURN

! exponential contour interval

  ELSEIF(KONTUR.EQ.0.OR.KONTUR.EQ.1.OR.KONTUR.EQ.50.OR.KONTUR.EQ.51)THEN

!    determine if contours should be at 10's or 100's
!      force 10's if kontur=50 or 51
     CINT=10.0
     IF(.NOT.(KONTUR.EQ.50.OR.KONTUR.EQ.51))THEN
        IF(CTMAX/CTMIN.GT.1.0E+08)CINT=100.0
     END IF

!    set the basic intervals: 1-min nval=max
     NEXP=INT(ALOG10(CTMAX))
     IF(NEXP.LT.0)NEXP=NEXP-1

!    set max values (color fill contour below value)
     CVAL(NVAL)  =10.0**NEXP
     TVAL(NVAL)  =CVAL(NVAL)*CINT
     TVAL(NVAL+1)=CVAL(NVAL)*CINT*CINT

!    set remaining contours
     DO N=(NVAL-1),1,-1
        CVAL(N)=CVAL(N+1)/CINT
        TVAL(N)=TVAL(N+1)/CINT
     END DO

!    manually force min contour to min conc when kontur.ne.4
!      (cannot set contour to 0 since a contour shows bdry between 2 values)
!      (map domain is scaled based on plume size (conc>0) and conc grid span, 
!         see conndx,mapini)
   ! CVAL(1)=CTMIN

!    set the lowest contour to the minimum
     TVAL(0)=TVAL(1)/CINT

! linear contour interval
  ELSEIF(KONTUR.EQ.2.OR.KONTUR.EQ.3)THEN

     NEXP=NINT(ALOG10(CTMAX/4.0))
     IF(NEXP.LT.0)NEXP=NEXP-1
     CINT=10.0**NEXP
     IF(CTMAX/CINT.GT.6)CINT=2.0*CINT

!    fill in remaining values for offset contour
     DO N=1,NVAL
        CVAL(N)=CINT*N
     END DO
     DO N=1,NVAL
        TVAL(N-1)=CVAL(N)
     END DO
     TVAL(NVAL)=CVAL(NVAL)+CINT
     TVAL(NVAL+1)=TVAL(NVAL)+CINT

  ELSE
     WRITE(*,*)'ERROR contur: contour values not properly set'
     STOP 900
  END IF

!-------------------------------------------------------------------------------
! call contour fill routine (one time for each value)
!-------------------------------------------------------------------------------

! determine size of plot in paper units
  CALL MAP2XY(1.0,1.0,X1,Y1)
  CALL MAP2XY(FLOAT(NXP),FLOAT(NYP),X2,Y2)
  XLEN=X2-X1
  YLEN=Y2-Y1

! shift plot coordinate reference
  CALL PLOT(X1,Y1,-3)

! set the contour colors based on type of output
  IF(KMAP.EQ.4)THEN
   COLOR=COLOR2
   GRYLEV=GRYLEV2
!  check to see if some contour levels are missing
   DO N=3,5
    IF(CVAL(N-1).EQ.0.0)THEN
        COLOR(1,N)=0.9
        COLOR(2,N)=0.9
        COLOR(3,N)=0.9
        GRYLEV(N)=0.95
    END IF
   END DO
  ELSE
   COLOR=COLOR1
   GRYLEV=GRYLEV1
   IF(COLORX(1,1).GE.0.0) COLOR=COLORX
  END IF

! B&W or color fill routine
  IF(KOLOR.EQ.0.OR.KOLOR.EQ.3)THEN
     CALL CONFILL(XCON,NXP,NXP,NYP,XLEN,YLEN,TVAL,GRYLEV,NVAL+2,1,0.)
  ELSE
     CALL CONCOLR(XCON,NXP,NXP,NYP,XLEN,YLEN,TVAL,COLOR,NVAL+2,1,0.)
  END IF

! open GIS file if requested
  IF(IGIS.GE.1)OPEN(16,FILE=GISTMP,STATUS="REPLACE")

  IF(.NOT.(KOLOR.EQ.2.OR.KOLOR.EQ.3))THEN
!   contouring routine for black outline on fill pattern
    CALL SETLW(0.002)
    CALL CONREC(XCON,NXP,NXP,NYP,XLEN,YLEN,CVAL,NVAL)
    CALL SETLW(0.0)
  END IF

! shift reference frame back
  CALL PLOT(-X1,-Y1,-3)

!-------------------------------------------------------------------------------
! special color at max conc point
!-------------------------------------------------------------------------------

  CALL MAP2XY((XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(1),YARR(1))
  CALL MAP2XY((XMAX+0.5*NGPT),(YMAX-0.5*NGPT),XARR(2),YARR(2))
  CALL MAP2XY((XMAX+0.5*NGPT),(YMAX+0.5*NGPT),XARR(3),YARR(3))
  CALL MAP2XY((XMAX-0.5*NGPT),(YMAX+0.5*NGPT),XARR(4),YARR(4))
  CALL MAP2XY((XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(5),YARR(5))

  IF(KOLOR.EQ.0.OR.KOLOR.EQ.3)THEN
     CALL COLBOX(XARR,YARR,5,0.,0.,0.)
  ELSE IF(KMAP.EQ.4)THEN
!    make max box black
     CALL COLBOX(XARR,YARR,5,0.,0.,0.)
  ELSE
!    make max box red
     CALL COLBOX(XARR,YARR,5,1.,0.,0.)
  END IF

! convert GIS contours from paper units to lat/lon 
  IF(IGIS.GT.0.AND.IGIS.LE.2)THEN
     REWIND(16)
     OPEN(17,FILE=GISOUT)

     OPEN(18,FILE=GISATT)
     WRITE(18,'(A)')'#CONC,NAME,DATE,TIME,LLEVEL,HLEVEL'

     ICNT=0
     KRET=0
     loop : DO WHILE (kret.EQ.0)
        READ(16,'(E14.7,I6)',IOSTAT=kret) CV,NP
        IF(kret.NE.0) EXIT loop
        IF(IGIS.EQ.1) THEN
           CVL=ALOG10(CV)
        ELSE
           CVL=CV
        END IF
        DO NN=1,NP
           READ(16,'(2F10.7)') XRA,YRA
           XR=XRA/XLEN
           YR=YRA/YLEN
           XP=XR*(NXP-1.0)+1.0
           YP=YR*(NYP-1.0)+1.0
           IF(KPROJ.EQ.4)THEN
              CALL CYL2LL(XP,YP,CLAT,CLON)
           ELSE
              CALL CXY2LL(PARMAP,XP,YP,CLAT,CLON)
           END IF
           IF(NN.EQ.1)THEN
              IF(IGIS.EQ.1)THEN
                 WRITE(17,'(F10.5,2(A1,F10.5))')CVL,',',CLON,',',CLAT 
              ELSE
                 WRITE(17,'(E10.3,2(A1,F10.5))')CVL,',',CLON,',',CLAT 
              END IF
           ELSE
              WRITE(17,'(F10.5,A1,F10.5)')CLON,',',CLAT 
           END IF 
        END DO
        WRITE(17,'(A)')'END'

!       Write output to attributes file (rjp - dec 2004)
        IF(JYR.LT.40)THEN
           WRITE(18,'(E10.3,A1,A,A1,I4,2I2.2,A1,I4.4,A1,I5.5,A1,I5.5)') &
           CVL,',',PTYPE,',',JYR+2000,JMO,JDA,',',JHR*100,',',LEVEL1,',',LEVEL2
        ELSE
           WRITE(18,'(E10.3,A1,A,A1,I4,2I2.2,A1,I4.4,A1,I5.5,A1,I5.5)') &
           CVL,',',PTYPE,',',JYR+1900,JMO,JDA,',',JHR*100,',',LEVEL1,',',LEVEL2
        END IF
     END DO loop

     WRITE(17,'(A)')'END'
     CLOSE (16)
     CLOSE (17)
     CLOSE (18)
  ELSE IF(IGIS.EQ.3)THEN
     KGE=0
!    Avoid contouring nothing, however create Google Earth label before exiting
!    write(6,*)'ctmax,tval(1),tval(2),tval(3),tval(4)',ctmax,tval(1),tval(2),tval(3),tval(4)
     xloop: DO N=0,4
       IF(TVAL(N).GT.0.0.AND.CTMAX.LE.TVAL(N))THEN
         KGE=2
         EXIT xloop
       ELSE
         EXIT xloop
       END IF
     END DO xloop
     CALL GEPLOT(KGE,KPROJ,LEVEL1,LEVEL2,PTYPE,PARMAP,NXP,NYP,   &
                  TMAX,TMIN,XMAX,YMAX,CVAL,NVAL,NGPT,            &
                  IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,       &
                  OLAT,OLON,OLVL1,OLVL2,GOPEN,KMAP,              &
                  MAPN,PROCESS,XLEN,YLEN)
  END IF

END SUBROUTINE contur
