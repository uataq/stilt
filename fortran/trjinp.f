!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TRJINP           TRAJECTORY INPUT FILE
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            READS THE TRAJECTORY ENDPOINTS FILE AND PLACES POINTS
!            IN AN ARRAY
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Feb 1997 (RRD)
!                 04 May 1999 (RRD) - included direction test
!                 05 Jul 2000 (RRD) - number of diagnostic variables
!                 12 Dec 2000 (RRD) - fortran90 upgrade
!                 03 Jan 2001 (RRD) - vertical coordinate option
!                 20 Feb 2001 (RRD) - saved starting height data
!                 21 Jun 2001 (RRD) - multiple diagnostic variables
!                 18 Nov 2001 (RRD) - vertical starting position 
!                 11 Feb 2002 (RRD) - multiple trajectories in time
!                 18 Mar 2002 (RRD) - diagnostic variable format
!                 28 Aug 2002 (RRD) - splitting trajectory origins
!                 20 Nov 2003 (RRD) - minutes field & mult diagnostic
!                 11 May 2004 (RRD) - array space min alloc all times
!                 19 May 2004 (RRD) - multi-time traj labels
!                 07 Jun 2004 (BS)  - backward initialization time
!                 04 Nov 2004 (BS)  - last hour endpoint to read,output hours
!                 31 May 2005 (RRD) - updated trajectory file format
!                 06 Oct 2005 (RRD) - four digit year test
!                 02 Oct 2007 (RRD) - pass through initial terrain height
!                 30 Oct 2008 (RRD) - save initial times in date array
!
! USAGE:  CALL TRJINP(DIRCTN,TFMT,NTRAJ,KPTS,NDIAG,JFYR,JFMO,JFDA,JFHR,KAGL,
!                     OLAT,OLON,OLVL,OTER,TLAT,TLON,THGT,TERR,KFHR,YY,MM,DD,
!                     HH,KM,FF,RR,NP,LHRS,KTER)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             UNIT 10 - trajectory endpoints
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$
 
SUBROUTINE TRJINP(DIRCTN,TFMT,NTRAJ,KPTS,NDIAG,JFYR,JFMO,JFDA,JFHR,KAGL,OLAT,  &
           OLON,OLVL,OTER,TLAT,TLON,THGT,TERR,KFHR,YY,MM,DD,HH,KM,FF,RR,NP,    &
           LHRS,KTER)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
  CHARACTER(8), INTENT(IN)    :: DIRCTN       ! direction variable
  INTEGER,      INTENT(IN)    :: TFMT         ! trajectory file format 
  INTEGER,      INTENT(IN)    :: NTRAJ        ! number of trajectories in file
  INTEGER,      INTENT(INOUT) :: KPTS         ! number of points each traj
  INTEGER,      INTENT(IN)    :: NDIAG        ! number of diagnostic variables
  INTEGER,      INTENT(IN)    :: KAGL         ! vertical coordinate flag
  INTEGER,      INTENT(OUT)   :: JFYR,JFMO    ! meteo initialization time
  INTEGER,      INTENT(OUT)   :: JFDA,JFHR    ! meteo initialization time
  REAL,         INTENT(INOUT) :: OLAT (:)     ! starting latitude
  REAL,         INTENT(INOUT) :: OLON (:)     ! starting longitude
  REAL,         INTENT(OUT)   :: OLVL (:)     ! starting heights from endpoints
  REAL,         INTENT(OUT)   :: OTER (:)     ! starting terrain height        
  REAL,         INTENT(OUT)   :: TLAT (:,:)   ! array of endpoint positions
  REAL,         INTENT(OUT)   :: TLON (:,:)   ! array of endpoint positions
  REAL,         INTENT(OUT)   :: THGT (:,:)   ! array of endpoint height    
  REAL,         INTENT(OUT)   :: TERR (:,:)   ! array of endpoint terrain 
  INTEGER,      INTENT(OUT)   :: KFHR         ! maximum forecast hour 
  INTEGER,      INTENT(OUT)   :: YY(:)        ! endpoint year
  INTEGER,      INTENT(OUT)   :: MM(:,:)      ! endpoint month
  INTEGER,      INTENT(OUT)   :: DD(:,:)      ! endpoint day
  INTEGER,      INTENT(OUT)   :: HH(:,:)      ! endpoint hour
  INTEGER,      INTENT(OUT)   :: KM(:)        ! endpoint minute
  INTEGER,      INTENT(OUT)   :: FF(:)        ! endpoint forecast hour
  REAL,         INTENT(OUT)   :: RR(:)        ! endpoint duration
  INTEGER,      INTENT(OUT)   :: NP(:,:)      ! number of endpoints each traj
  INTEGER,      INTENT(IN)    :: LHRS         ! last hours endpoint of interest
  INTEGER,      INTENT(IN)    :: KTER         ! terrain height index  

!-------------------------------------------------------------------------------

  REAL, ALLOCATABLE :: tdiag(:)
  REAL              :: ylon,xlat,phour,hours,zhgt
  INTEGER           :: kk,kd,jf,jy,mg,kp,jm,jn,jh,jd,knum,kret 

!-------------------------------------------------------------------------------

! only first two diagnostic variables are currently used for plotting
  ALLOCATE (tdiag(ndiag), STAT=kret)

! dimensions set based upon first file opened
  KNUM=SIZE(NP,1)

! zero end point counter
  DO KP=1,KNUM  
     NP(KP,1)=0
     NP(KP,2)=0
     IF(KPTS.NE.0)OLVL(KP)=-1.0
  END DO

! first pass through data to determine dimensions
  IF(KPTS.EQ.0)THEN
     PHOUR=0.0
10   READ(10,'(I6,42X,F8.1)',END=20) KP,HOURS
     IF(HOURS.NE.PHOUR)THEN
        PHOUR=HOURS
!       array space should contain all time periods
!       11 May 2004 modification for mult traj
        NP(1,1)=NP(1,1)+1
     END IF
     IF(LHRS.GT.0.AND.ABS(HOURS).GT.LHRS) GOTO 10		! BS - 10-26-04
     NP(KP,2)=NP(KP,2)+1
     KPTS=MAX(NP(KP,2),NP(1,1),KPTS) 
     GOTO 10
20   REWIND(10)
     RETURN
  END IF

!-------------------------------------------------------------------------------

! initialize maximum forecast hour
  KFHR=0

  kret=0
  tloop : DO WHILE (kret.EQ.0)

  IF(tfmt.EQ.0)THEN
     READ(10,'(8I6,F8.1,2F8.3,9F8.1)',IOSTAT=kret)                          &
          KP,MG,JY,JM,JD,JH,JN,JF,HOURS,XLAT,YLON,ZHGT,(TDIAG(KD),KD=1,NDIAG)
  ELSE
     READ(10,'(8I6,F8.1,2F9.3,F8.1,9(1X,F8.1))',IOSTAT=kret)                &
          KP,MG,JY,JM,JD,JH,JN,JF,HOURS,XLAT,YLON,ZHGT,(TDIAG(KD),KD=1,NDIAG)
  END IF
  IF(kret.NE.0) EXIT tloop
  JY=MOD(JY,100)

! check numb of traj array exceedance
  IF(KP.GT.KNUM) CYCLE tloop

! check if HOURS greater than that wanted
  IF(LHRS.GT.0.AND.ABS(HOURS).GT.LHRS) CYCLE tloop           ! BS - 10-26-04

! check if this is the first endpoint
  IF(NP(KP,2).EQ.0)THEN
!    check for splitting multiple trajectories (unkown origin)
     IF(OLAT(KP).EQ.99.0.AND.OLON(KP).EQ.99.0)THEN
!       origin given by antecedent index
        KK=NP(MG,2)-1
        OLAT(KP)=TLAT(KK,MG)
        OLON(KP)=TLON(KK,MG)
     END IF 
  END IF

! starting position saved for plotting at initial time for all start points

  IF(OLVL(KP).LT.0.0)THEN
     IF(KAGL.EQ.0)THEN
        OLVL(KP)=TDIAG(1)
     ELSEIF(KAGL.EQ.1)THEN
        IF(KTER.NE.0)THEN
           OTER(KP)=TDIAG(KTER)                     ! rrd - 10/2/2007
           OLVL(KP)=ZHGT+TDIAG(KTER)
        ELSE
           OLVL(KP)=ZHGT
        END IF
     ELSEIF(KAGL.EQ.2)THEN
        OLVL(KP)=TDIAG(2)
     ELSE
        OLVL(KP)=TDIAG(NDIAG)
     END IF
  END IF

! skip starting location when in output listing
! increased dimension 10/30/2008
  IF(HOURS.EQ.0.0) THEN
     KM(SIZE(KM))=JN
     CYCLE tloop
  END IF

! set pointers for start and end for each trajectory
  NP(KP,2)=NP(KP,2)+1
  IF(NP(KP,1).EQ.0) THEN
     NP(KP,1)=MAXVAL(NP(:,2))
     NP(KP,2)=MAXVAL(NP(:,2))
  END IF

! test again because multiple file dimensions based upon the first file
  IF(NP(KP,2).GT.KPTS)THEN
     NP(KP,2)=KPTS           
     CYCLE tloop
  END IF

! save time and position
  KK=MAXVAL(NP(:,2))
  YY(KK)=JY
  MM(KK,KP)=JM
  DD(KK,KP)=JD
  HH(KK,KP)=JH
  KM(KK)=JN
  FF(KK)=JF
  RR(KK)=HOURS
  TLAT(KK,KP)=XLAT
  TLON(KK,KP)=YLON

! diagnostic variable for display purposes
  IF(KAGL.EQ.0)THEN
     THGT(KK,KP)=TDIAG(1)
  ELSEIF(KAGL.EQ.1)THEN
     IF(KTER.NE.0)THEN
!       terrain height in file and height is selected for display
!       then the AGL height is converted to MSL before plotting
!       WARNING: not valid if KAGL=0 selected in HYSPLIT setup.cfg
        TERR(KK,KP)=TDIAG(KTER)
        THGT(KK,KP)=ZHGT+TDIAG(KTER)
     ELSE
        THGT(KK,KP)=ZHGT
     END IF
  ELSEIF(KAGL.EQ.2)THEN
     THGT(KK,KP)=TDIAG(2)
  ELSE
     THGT(KK,KP)=TDIAG(NDIAG)
  END IF

! find the maximum forecast time (<12 hours = analysis)
  KFHR=MAX(KFHR,JF)

! find the forecast initialization time: last-forward first-backward
! IF(DIRCTN(1:4).EQ.'FORW'.OR.NP(1,1).EQ.1)THEN (07 Jun 2004 - BS)
  IF(DIRCTN(1:4).EQ.'FORW'.OR.DIRCTN(1:4).EQ.'BACK'.AND.NP(1,2).EQ.1)THEN
     JFYR=JY
     JFMO=JM
     JFDA=JD
     JFHR=JH
     CALL TMPLUS(JFYR,JFMO,JFDA,JFHR,-JF)
  END IF

  END DO tloop
  CLOSE (10)

END SUBROUTINE trjinp
