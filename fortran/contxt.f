!$$$ SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONTXT           ADD CUSTOMIZED LABELS TO BASIC MAP
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            ADDS LABELS AND OTHER TEXT INFORMATION TO THE FINISHED
!            MAP.  EACH GRAPHICAL APPLICATION PROGRAM WILL HAVE A
!            SLIGHTLY DIFFERENT VERSION OF THIS SUBROUTINE ATTACHED
!            TO THE MAIN PROGRAM FOR COMPILATION
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 (RRD)
!                 09 Dec 1998 (RRD) - converted NCAR graphics to psplot
!                 10 Mar 2000 (RRD) - label change for exposure plot
!                                   - use labels.cfg file for custom plot
!                 21 Nov 2000 (RRD) - fortran90 upgrade      
!                 05 Dec 2000 (RRD) - enhanced labels
!                 02 Nov 2001 (RRD) - default label change
!                                   - incorporate RSMC plotting standards
!                 21 Nov 2001 (BS)  - add extra digit to source height label
!                 11 Jan 2002 (RRD) - source plotting symbol variable
!                 28 Feb 2002 (BS)  - KAVG=3 forces layer-avg label
!                 04 Mar 2002 (RRD) - various label modifications
!                 20 Jun 2002 (RRD) - process ID for supplemental files
!                 08 Aug 2002 (RRD) - expanded text in maptext file
!                 23 Dec 2002 (RRD) - inconsistent forecast time with no data 
!                 04 Feb 2003 (RRD) - multiple source location label
!                 17 Jul 2003 (RRD) - revised default NOAA label
!                 20 Nov 2003 (RRD) - added minutes to emission label
!                 26 Jan 2007 (RRD) - added minutes to averaging period label
!                 30 Jul 2007 (RRD) - turn off left edge source if symbol=blank
!                 29 Oct 2007 (BS)  - KMAP=5 volcanic ash (VA)
!                 30 Nov 2007 (GDR) - Passed IYR through from CONMAP
!                                   - output integration time Google Earth label
!                 21 Mar 2008 (RRD) - foward/backward date test
!                 27 Jan 2009 (RRD) - command line color mapping
!
! ***************  NOTE: volcplot.f needs to be edited for its version of CONTXT
! ***************        if changes are made here.
!
!
! USAGE:  CALL CONTXT(LEVEL1,LEVEL2,KMAP,KAVG,MODEL,PTYPE,NVAL,CMAX,CMIN,CVAL,
!              KOLOR,PMODE,IBYR,IBMO,IBDA,IBHR,IBMN,OLAT,OLON,OLVL1,OLVL2,
!              IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,JFH,MAPN,QCHAR,PROCESS,CLBL)
!
!   INPUT PARAMETERS:     see below
!   OUTPUT PARAMETERS:    see below
!   INPUT FILES:          none
!   OUTPUT FILES:         none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CONTXT(LEVEL1,LEVEL2,KMAP,KAVG,MODEL,PTYPE,NVAL,CMAX,CMIN,CVAL,     &
           KOLOR,PMODE,QUNIT,IBYR,IBMO,IBDA,IBHR,IBMN,OLAT,OLON,OLVL1,OLVL2,   &
           IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,JFH,MAPN,QCHAR,PROCESS,CLBL)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: level1,level2     ! range of display levels 
  INTEGER, INTENT(IN)  :: kmap,kavg         ! map type (concen/deposit)
  INTEGER, INTENT(IN)  :: nval              ! number of contour values
  REAL,    INTENT(IN)  :: cmax,cmin         ! max and min concentrations   
  REAL,    INTENT(IN)  :: cval(:)           ! concentration contour values
  INTEGER, INTENT(IN)  :: kolor             ! flag to set B&W=0 or color=1
  LOGICAL, INTENT(IN)  :: pmode             ! true indicates portrait display
  INTEGER, INTENT(IN)  :: ibyr (:)          ! emission start year 
  INTEGER, INTENT(IN)  :: ibmo (:)          ! emission start month
  INTEGER, INTENT(IN)  :: ibda (:)          ! emission start day
  INTEGER, INTENT(IN)  :: ibhr (:)          ! emission start hour
  INTEGER, INTENT(IN)  :: ibmn (:)          ! emission start minute
  REAL,    INTENT(IN)  :: olat (:)          ! starting location
  REAL,    INTENT(IN)  :: olon (:)          ! starting location 
  REAL,    INTENT(IN)  :: olvl1,olvl2       ! release height range
  INTEGER, INTENT(IN)  :: imn,jmn           ! minutes field    
  INTEGER, INTENT(IN)  :: imo,ida,ihr       ! starting date/time of the sample
  INTEGER, INTENT(IN)  :: iyr               ! starting date/time of the sample
  INTEGER, INTENT(IN)  :: jmo,jda,jhr       ! stop date/time of the sample
  INTEGER, INTENT(IN)  :: jyr,jfh           ! stop date/time of the sample
  INTEGER, INTENT(IN)  :: mapn              ! current plotting map number

  CHARACTER(4),  INTENT(IN) :: model        ! ident of meteorological data
  CHARACTER(4),  INTENT(IN) :: ptype        ! ident of pollutant displayed
  CHARACTER(80), INTENT(IN) :: qunit        ! mass units for pollutant
  CHARACTER(1),  INTENT(IN) :: qchar        ! source plotting symbol
  CHARACTER(8),  INTENT(IN) :: process      ! process ID suffix
  CHARACTER(8),  INTENT(IN) :: clbl(:)      ! concentration contour labels

!-------------------------------------------------------------------------------

  LOGICAL       :: FTEST
  INTEGER       :: PREFOR = 0               ! previous forecast time minutes
  INTEGER       :: CURFOR = 0               ! current forecast time minutes
  REAL          :: XARR(5),YARR(5)
  REAL          :: COLORX(3,6),COLOR(3,5),COLOR1(3,5),COLOR2(3,5)
  CHARACTER(80) :: LABEL, MAPID, LAYER, UNITS, FIELD, TITLE
  CHARACTER(3)  :: MONTH
  CHARACTER(8)  :: IDENT, VOLUM
  CHARACTER(80) :: gisout,gisatt,gistmp,gislbl

  REAL      :: xp,yp,xp2,yp2,xpt,xtab,ypt,delx,dely,dely2,gray 
  INTEGER   :: n,k,kk,klen,klen1,klen2,klen6,klen7,klen8,klen9,kyr,kmo,kda,khr,igis
  INTEGER   :: iyear,jyear

! colors from command line
  COMMON /RGBVAL/ COLORX

! GIS option common block
  COMMON /GIS/GISOUT,GISATT,GISTMP,IGIS

!-------------------------------------------------------------------------------

! RGB contours: red  green  blue
! order: cyan,green,blue,yellow,white
  DATA COLOR1/  0.0,  1.0,  1.0,                                           &
                0.0,  1.0,  0.0,                                           &
                0.0,  0.0,  1.0,                                           &
                1.0,  1.0,  0.0,                                           &
                1.0,  1.0,  1.0/

! order: gray,yellow,orange,red,white    
  DATA COLOR2/  0.9,  0.9,  0.9,                                           &
                1.0,  1.0,  0.0,                                           &
                1.0,  0.5,  0.0,                                           &
                1.0,  0.0,  0.0,                                           &
                1.0,  1.0,  1.0/

!-------------------------------------------------------------------------------

  SAVE PREFOR

! set window to relative for final label information
  CALL MAPSET(PMODE,.01,.95,.10,.90,0.,1.,0.,1.)

!-------------------------------------------------------------------------------
! initialize custom labels if labels.cfg in root directory
!-------------------------------------------------------------------------------

! default values
  TITLE='NOAA HYSPLIT MODEL&'

! default map/unit/volume labels 
  UNITS=QUNIT
  IF(KMAP.EQ.1.AND.UNITS(1:3).EQ.'ppm')THEN
     MAPID='Concentration &'
     VOLUM='&'
  ELSEIF(KMAP.EQ.1.OR.KMAP.EQ.5)THEN
     MAPID='Concentration &'
     VOLUM='/m3&'
  ELSEIF(KMAP.EQ.2)THEN
     MAPID='Exposure &'
     VOLUM='-s/m3&'
  ELSEIF(KMAP.EQ.3)THEN
     MAPID='Deposition &'
     VOLUM='/m2&'
  ELSEIF(KMAP.EQ.4)THEN
     MAPID='Concentration &'
     VOLUM='&'
  END IF

! default averaging layer labels 
  LAYER=' &'
  IF(KAVG.EQ.1)THEN
     LAYER=' at level &'
  ELSEIF(KAVG.EQ.2.OR.KAVG.EQ.3)THEN
     LAYER=' averaged between &'
  END IF

! optional process ID suffix
  IF(PROCESS(1:2).EQ.'ps')THEN
     LABEL='LABELS.CFG' 
  ELSE
     LABEL='LABELS.'//PROCESS
  END IF

  INQUIRE(FILE=LABEL,EXIST=FTEST)
  IF(FTEST)THEN
     OPEN(50,FILE=LABEL)
100  READ(50,*,END=200)IDENT,FIELD
        KLEN1=INDEX(FIELD,'&')
        KLEN2=INDEX(IDENT,'&')
        IF(IDENT(1:KLEN2).EQ.'TITLE&')TITLE=FIELD(1:KLEN1)
        IF(IDENT(1:KLEN2).EQ.'UNITS&')UNITS=FIELD(1:KLEN1)
        IF(IDENT(1:KLEN2).EQ.'VOLUM&')VOLUM=FIELD(1:KLEN1)
        IF(IDENT(1:KLEN2).EQ.'LAYER&')LAYER=FIELD(1:KLEN1)
        IF(IDENT(1:KLEN2).EQ.'MAPID&')MAPID=FIELD(1:KLEN1)
     GOTO 100
200  CLOSE(50)
  END IF
  KLEN6=INDEX(MAPID,'&')-1
  KLEN7=INDEX(LAYER,'&')-1
  KLEN8=INDEX(UNITS,'&')-1
  KLEN9=INDEX(VOLUM,'&')-1

!-------------------------------------------------------------------------------
! top of map identification
!-------------------------------------------------------------------------------

  LABEL=TITLE
  KLEN=INDEX(LABEL,'&')-1
  CALL MAP2XY(0.50,0.99,XP,YP)
  CALL KEKSYMC(XP,YP,0.12,LABEL,0.0,KLEN,1)

!---------------------

  IF(KMAP.LE.2.OR.KMAP.EQ.4.OR.KMAP.EQ.5)THEN
!    concentration and exposure
     IF(KAVG.EQ.1)THEN  
        WRITE(LABEL,'(A,I5,A)')                            &
              MAPID(:KLEN6)//' ('//                        &
              UNITS(:KLEN8)//VOLUM(:KLEN9)//')'//          &
              LAYER(:KLEN7), LEVEL2,' m&'
     ELSEIF(KAVG.EQ.2.OR.KAVG.EQ.3)THEN
        WRITE(LABEL,'(A,2(I5,A))')                        &
              MAPID(:KLEN6)//' ('//                        &
              UNITS(:KLEN8)//VOLUM(:KLEN9)//')'//          &
              LAYER(:KLEN7),LEVEL1,' m and ',LEVEL2,' m&'
     END IF

  ELSEIF(KMAP.EQ.3)THEN
!    deposition
     WRITE(LABEL,'(A)')                                    &
           MAPID(:KLEN6)//' ('//                           &
           UNITS(:KLEN8)//VOLUM(:KLEN9)//')'//             &
           ' at ground-level&'
  END IF

  KLEN=INDEX(LABEL,'&')-1
  CALL MAP2XY(0.50,0.955,XP,YP)
  CALL KEKSYMC(XP,YP,0.10,LABEL,0.0,KLEN,1)

!---------------------

! check integration direction (rrd: 2/28/2002)
  CALL TM2MIN(IYR,IMO,IDA,IHR,IMN,K)
  CALL TM2MIN(JYR,JMO,JDA,JHR,JMN,KK)

! base integration time label 
  IF(KK.GE.K)THEN
     LABEL='Integrated from ---- -- --- to ---- -- --- -- (UTC)&'
  ELSE
     LABEL='Integrated from ---- -- --- to ---- -- --- -- (UTC) [backward]&'
  END IF

! write sample start times into record
  CALL MONSET(IMO,MONTH,1)
  WRITE(LABEL(17:20),'(2I2.2)')IHR,IMN
  WRITE(LABEL(22:27),'(I2.2,1X,A3)')IDA,MONTH

! write sample start times into record
  CALL MONSET(JMO,MONTH,1)
  WRITE(LABEL(32:35),'(2I2.2)')JHR,JMN
  WRITE(LABEL(37:45),'(I2.2,1X,A3,1X,I2.2)')JDA,MONTH,JYR

! write record
  KLEN=INDEX(LABEL,'&')-1
  CALL MAP2XY(0.50,0.925,XP,YP)
  CALL KEKSYMC(XP,YP,0.10,LABEL,0.0,KLEN,1)

! write contour levels to text file for Google Earth overlay label
! modified gdr 11/30/2007
  IYEAR=IYR+2000
  JYEAR=JYR+2000
  IF(IYR.GT.40)IYEAR=IYR+1900
  IF(JYR.GT.40)JYEAR=JYR+1900
  IF(IGIS.EQ.3)THEN   
     WRITE(18,'(I1)')KMAP
     WRITE(18,'(A)')UNITS(:KLEN8)//VOLUM(:KLEN9)//'&'
     CALL MONSET(IMO,MONTH,2)
     WRITE(18,'(A,2I2.2,A,A3,A,I2.2,A,I4,A)')'Integrated:  ',IHR,IMN,' UTC  ',MONTH,'  ',IDA,'   ',IYEAR,'&'
     CALL MONSET(JMO,MONTH,2)
     WRITE(18,'(A,2I2.2,A,A3,A,I2.2,A,I4,A)')'        to:  ',JHR,JMN,' UTC  ',MONTH,'  ',JDA,'   ',JYEAR,'&'
     WRITE(18,'(2((1P,E7.1),1X))')CMAX,CMIN
     WRITE(18,'(4((1P,E7.1),1X))')(CVAL(K),K=4,1,-1)
     WRITE(18,'(4(A6,1X))')(CLBL(K),K=1,4)
  END IF

!---------------------

! release starting information only for source #1
! added minutes field through argument list (20 Nov 2003)

  LABEL=PTYPE//' Release started at ---- -- --- -- (UTC)&'
  CALL MONSET(IBMO(1),MONTH,1)
  WRITE(LABEL(25:28),'(2I2.2)')IBHR(1),IBMN(1)
  WRITE(LABEL(30:38),'(I2.2,1X,A3,1X,I2.2)')IBDA(1),MONTH,IBYR(1)
  KLEN=INDEX(LABEL,'&')-1
  CALL MAP2XY(0.50,0.895,XP,YP)
  CALL KEKSYMC(XP,YP,0.10,LABEL,0.0,KLEN,1)

!-------------------------------------------------------------------------------
! right panel labels
!-------------------------------------------------------------------------------

  XP =0.68
  XPT=1.00
  YP =0.125
  YPT=0.875

  CALL MAP2XY(XP ,YPT,XARR(1),YARR(1))
  CALL MAP2XY(XP ,YP ,XARR(2),YARR(2))
  CALL MAP2XY(XPT,YP ,XARR(3),YARR(3))
  CALL MAP2XY(XPT,YPT,XARR(4),YARR(4))
  CALL MAP2XY(XP ,YPT,XARR(5),YARR(5))
  CALL SLDCRV(XARR,YARR,5,0.015)

! lower left corner position and box size     
  XPT=0.74
  YPT=0.85
  DELX=0.04
  DELY=0.02

! disclaimer (optional)
  IF(KMAP.EQ.4.OR.KMAP.EQ.5)THEN
     LABEL='Not for Public Dissemination&'
     IF(KMAP.EQ.5)LABEL='*** Hypothetical eruption ***&'
     KLEN=INDEX(LABEL,'&')-1
     CALL MAP2XY(XPT+0.01,YPT,XP,YP)
     CALL SETCOLR(1.0,0.0,0.0)
     CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
     CALL SETCOLR(0.0,0.0,0.0)
     YPT=YPT-2.0*DELY
  END IF

! set the contour colors based on type of output
  IF(KMAP.EQ.4)THEN
     COLOR=COLOR2
  ELSE
     COLOR=COLOR1
     IF(COLORX(1,1).GE.0.0)THEN
        DO K=1,5
           COLOR(:,K)=COLORX(:,K+1)
        END DO
     END IF
  END IF

! color fill box size and left shift from numbers

  DO K=1,NVAL
     KK=NVAL-K+1

     CALL MAP2XY(XPT-DELX,YPT-DELY,XARR(1),YARR(1))
     CALL MAP2XY(XPT+DELX,YPT-DELY,XARR(2),YARR(2))
     CALL MAP2XY(XPT+DELX,YPT+DELY,XARR(3),YARR(3))
     CALL MAP2XY(XPT-DELX,YPT+DELY,XARR(4),YARR(4))
     CALL MAP2XY(XPT-DELX,YPT-DELY,XARR(5),YARR(5))

     IF(KOLOR.EQ.0.OR.KOLOR.EQ.3)THEN
!       cval kk=1 at highest value
!       set gray scale high=dark=0  low=white=1
        GRAY=0.2*K+0.1
        IF(KMAP.EQ.4.AND.K.EQ.4)GRAY=1.0
        CALL FILLBOX(XARR,YARR,5,GRAY)
     ELSE
        CALL COLBOX(XARR,YARR,5,COLOR(1,KK),COLOR(2,KK),COLOR(3,KK))
     END IF

     YPT=YPT-2.0*DELY
  END DO

!---------------------

  CALL SETCOLR(0.,0.,0.)

  XPT=0.80
  YPT=0.84
  IF(KMAP.EQ.4.OR.KMAP.EQ.5)YPT=0.80

  DO K=1,NVAL
     KK=NVAL-K+1
!    Special labels for Emergency Activation Levels
     IF(KMAP.EQ.4.OR.KMAP.EQ.5)THEN
       IF(KK.LT.6)THEN
         IF(KK.EQ.1.OR.CVAL(KK).NE.0.0)THEN
            WRITE(LABEL(1:18), '(1A,1P,E7.1,A)')'>',CVAL(KK),' '//UNITS(:KLEN8)//VOLUM(:KLEN9)//'&'
         ELSE
            WRITE(LABEL(1:18), '(A3)')'NR&'
         END IF
         KLEN=INDEX(LABEL,'&')-1
         CALL MAP2XY(XPT,YPT,XP,YP)
         CALL KEKSYMC(XP,YP,0.10,LABEL,0.0,KLEN,0)
         CALL MAP2XY(XPT-0.1,YPT,XP,YP)
         KLEN=INDEX(CLBL(K),' ')-1
         IF(KK.NE.1)CALL KEKSYMC(XP,YP,0.085,CLBL(K),0.0,KLEN,0)
       END IF
     ELSE
       LABEL='> 0.0E 00&'
       WRITE(LABEL(1:19), '(1A,1P,E7.1,A)')'>',CVAL(KK),' '//UNITS(:KLEN8)//VOLUM(:KLEN9)//'&'
       KLEN=INDEX(LABEL,'&')-1
       CALL MAP2XY(XPT-0.02,YPT,XP,YP)
       CALL KEKSYMC(XP,YP,0.10,LABEL,0.0,KLEN,0)
     END IF
     YPT=YPT-2.0*DELY
  END DO

  IF(CMAX.GT.0.0)THEN
    XPT=0.70
    WRITE(LABEL,'(A,1P,E7.1,A)') 'Maximum: ',CMAX,'&'
    KLEN=INDEX(LABEL,'&')-1
    CALL MAP2XY(XPT,YPT,XP,YP)
    CALL KEKSYMC(XP,YP,0.09,LABEL,0.0,KLEN,0)
    WRITE(LABEL,'(A)') '(identified as a square)&'
    KLEN=INDEX(LABEL,'&')-1
    YPT=YPT-1.4*DELY
    CALL MAP2XY(XPT+0.09,YPT,XP,YP)
    CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)

    WRITE(LABEL,'(A,1P,E7.1,A)') 'Minimum: ',CMIN,'&'
    KLEN=INDEX(LABEL,'&')-1
    YPT=YPT-1.8*DELY
    CALL MAP2XY(XPT,YPT,XP,YP)
    CALL KEKSYMC(XP,YP,0.09,LABEL,0.0,KLEN,0)
  END IF

! Explanation of Emergency Activation Levels
! here
  IF(KMAP.EQ.4)THEN
    XPT=0.69
    XTAB=0.78
    YPT=0.54
    DELY2=0.013
    IF(CLBL(1)(1:4).EQ.'AEGL'.OR.CLBL(1)(1:4).EQ.'ERPG'        &
      .OR.CLBL(1)(1:4).EQ.'TEEL')THEN
      LABEL='ACUTE (SHORT-TERM) EFFECTS&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-2.0*DELY
        CALL MAP2XY(XPT+0.01,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
        YPT=YPT-0.5*DELY
        CALL MAP2XY(XPT,YPT,XP,YP)
        CALL MAP2XY(XPT+0.3,YPT,XP2,YP2)
        CALL SLDLIN(XP,YP,XP2,YP2,0.)
      WRITE(LABEL,'(A1,A,A)')'>',CLBL(1),': Death or irreversible health&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY
        CALL MAP2XY(XPT,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      LABEL='effects possible.&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY2
        CALL MAP2XY(XTAB,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      WRITE(LABEL,'(A1,A,A)')'>',CLBL(2),': Serious health effects or&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-2.0*DELY
        CALL MAP2XY(XPT,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      LABEL='or impaired ability to take&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY2
        CALL MAP2XY(XTAB,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      LABEL='protective action.&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY2
        CALL MAP2XY(XTAB,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      WRITE(LABEL,'(A1,A,A)')'>',CLBL(3),': Minor reversible health&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-2.0*DELY
        CALL MAP2XY(XPT,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      LABEL='effects. Possible odor.&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY2
        CALL MAP2XY(XTAB,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
        YPT=YPT-0.5*DELY
        CALL MAP2XY(XPT,YPT,XP,YP)
        CALL MAP2XY(XPT+0.3,YPT,XP2,YP2)
        CALL SLDLIN(XP,YP,XP2,YP2,0.)
    END IF
  ELSEIF(KMAP.EQ.5)THEN
    XPT=0.69
    XTAB=0.74
    YPT=0.54
    DELY2=0.013
        YPT=YPT-0.5*DELY
        CALL MAP2XY(XPT,YPT,XP,YP)
        CALL MAP2XY(XPT+0.3,YPT,XP2,YP2)
        CALL SLDLIN(XP,YP,XP2,YP2,0.)
        CALL SETCOLR(1.0,0.0,0.0)
      LABEL='Initial ash mass is 1.0 unit&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY2
        CALL MAP2XY(XPT,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      LABEL='Multiply by estimated mass&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY2
        CALL MAP2XY(XTAB,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      LABEL='or choose contour.&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY2
        CALL MAP2XY(XTAB,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      LABEL='For real eruption, see&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY
        CALL MAP2XY(XPT,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
      LABEL='SIGMET and VAAC products&'
        KLEN=INDEX(LABEL,'&')-1
        YPT=YPT-1.5*DELY2
        CALL MAP2XY(XTAB,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,KLEN,0)
        CALL SETCOLR(0.0,0.0,0.0)
        YPT=YPT-0.5*DELY
        CALL MAP2XY(XPT,YPT,XP,YP)
        CALL MAP2XY(XPT+0.3,YPT,XP2,YP2)
        CALL SLDLIN(XP,YP,XP2,YP2,0.)
  END IF

!-------------------------------------------------------------------------------
! left edge labels
!-------------------------------------------------------------------------------
  IF(QCHAR.NE.' ')THEN
     IF(MAXVAL(OLAT).EQ.MINVAL(OLAT).AND.MAXVAL(OLON).EQ.MINVAL(OLON))THEN
!       location label only for one release point
        LABEL='Source       0.000 N    0.000 E&'
        WRITE(LABEL(12:18),'(F6.3)')ABS(OLAT(1))
        IF(OLAT(1).LT.0.0)WRITE(LABEL(20:20),'(A1)')'S'
        WRITE(LABEL(22:29),'(F7.3)')ABS(OLON(1))
        IF(OLON(1).LT.0.0)WRITE(LABEL(31:31),'(A1)')'W'
     ELSE
        LABEL='Source       at multiple locations&'
     END IF
     KLEN=INDEX(LABEL,'&')-1
     CALL MAP2XY(0.02,0.12,XP,YP)
     CALL KEKSYMC(XP,YP,0.12,LABEL,90.0,KLEN,0)

!    enhanced source point symbol
     CALL MAP2XY(0.0,0.23,XP,YP)
     CALL SETFNT(35)
     CALL KEKSYMC(XP,YP,0.12,QCHAR,0.0,1,0)
     CALL SETFNT(20)
  END IF

!---------------------

  IF(OLVL1.EQ.OLVL2)THEN
     LABEL='from 00000 m&'
     WRITE(LABEL(6:10),'(I5)') NINT(OLVL1)
  ELSE
     LABEL='00000 m to 00000 m&'
     WRITE(LABEL(1:5),'(I5)') NINT(OLVL1)
     WRITE(LABEL(12:16),'(I5)') NINT(OLVL2)
  END IF

  KLEN=INDEX(LABEL,'&')-1
  CALL MAP2XY(0.02,0.88,XP,YP)
  CALL KEKSYMC(XP,YP,0.12,LABEL,90.0,KLEN,2)


!-------------------------------------------------------------------------------
! Bottom edge label 
!-------------------------------------------------------------------------------

! determine hindcast or forecast
  KYR=JYR
  KMO=JMO
  KDA=JDA
  KHR=JHR
  CALL TMPLUS(KYR,KMO,KDA,KHR,-JFH)
  CALL TM2MIN(KYR,KMO,KDA,KHR,JMN,CURFOR)
  CALL MONSET(KMO,MONTH,1)
  IF(JFH.GT.12.AND.(CURFOR.EQ.PREFOR.OR.PREFOR.EQ.0))THEN
     LABEL='---- -- --- -- '//MODEL//' FORECAST INITIALIZATION&'
     WRITE(LABEL(1:4),'(I4.4)')(KHR*100)
     WRITE(LABEL(6:14),'(I2.2,1X,A3,1X,I2.2)')KDA,MONTH,KYR
  ELSE
     LABEL=MODEL//' METEOROLOGICAL DATA&'
  END IF
  PREFOR=CURFOR

  KLEN=INDEX(LABEL,'&')-1
  CALL MAP2XY(0.42,0.09,XP,YP)
  CALL KEKSYMC(XP,YP,0.11,LABEL,0.0,KLEN,1)

!-------------------------------------------------------------------------------
! optional bottom text box with supplemental information
!-------------------------------------------------------------------------------

  YP =1.0
  YPT=15.0
  XP =0.0
  XPT=1.0
  CALL MAPSET(PMODE,.045,0.95,0.00,0.16,XP,XPT,YP,YPT)

! optional process ID suffix
  IF(PROCESS(1:2).EQ.'ps')THEN
     LABEL='MAPTEXT.CFG' 
  ELSE
     LABEL='MAPTEXT.'//PROCESS
  END IF

  INQUIRE(FILE=LABEL,EXIST=FTEST)
  IF(FTEST)THEN

!    projection box
     CALL MAP2XY(XP ,YPT,XARR(1),YARR(1))
     CALL MAP2XY(XP ,YP ,XARR(2),YARR(2))
     CALL MAP2XY(XPT,YP ,XARR(3),YARR(3))
     CALL MAP2XY(XPT,YPT,XARR(4),YARR(4))
     CALL MAP2XY(XP ,YPT,XARR(5),YARR(5))
     CALL SLDCRV(XARR,YARR,5,0.015)

     OPEN(50,FILE=LABEL)
     LABEL='                                                                                                    '
aa : DO K=1,15
        READ(50,'(A80)',END=300)LABEL
        IF(K.LT.3.OR.K.EQ.4.OR.K.EQ.6.OR.K.EQ.7) CYCLE aa
        IF(LABEL(1:4).EQ.'Traj') CYCLE aa 
        YPT=YPT-1.45
        CALL MAP2XY(0.05,YPT,XP,YP)
        CALL KEKSYMC(XP,YP,0.09,LABEL,0.0,80,0)
     END DO aa
300  CLOSE (50)
  END IF

END SUBROUTINE contxt
