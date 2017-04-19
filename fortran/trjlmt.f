!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TRJLMT           DETERMINE TRAJCTORY LIMITS
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            ANALYSES TRAJECTORY ARRAY TO DETERMINE MAP LIMITS THAT
!            WILL PRODUCE THE OPTIMIUM MAP FOR ALL TRAJECTORIES
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 08 Dec 1998 (RRD)
!                 28 Mar 2000 (RRD) - added centroid calculation
!                 12 Dec 2000 (RRD) - fortran90 upgrade
!                 03 Jan 2001 (RRD) - vertical coordinate option
!                 11 Feb 2002 (RRD) - multiple trajectories in time
!                 27 Feb 2002 (RRD) - higher resolution selection
!                 19 Mar 2002 (RRD) - modified limits test
!                 17 Jul 2002 (RRD) - major revision 
!                 29 Aug 2002 (RRD) - test for no points in traj
!                 23 Oct 2002 (RRD) - exclude virtual split sources
!                 04 Dec 2002 (RRD) - convert include to module
!
! USAGE:  CALL TRJLMT(KAGL,NTRAJ,OLAT,OLON,PBOT,PTOP,TLAT,TLON,THGT,NP)
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

SUBROUTINE TRJLMT(KAGL,NTRAJ,OLAT,OLON,PBOT,PTOP,TLAT,TLON,THGT,NP)

  USE mapbox

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  INTEGER,  INTENT(IN)  :: KAGL            ! vertical coordinate flag
  INTEGER,  INTENT(IN)  :: NTRAJ           ! number of trajectories this file
  REAL,     INTENT(IN)  :: OLAT(:)         ! starting points
  REAL,     INTENT(IN)  :: OLON(:)         ! starting points
  REAL,     INTENT(IN)  :: TLAT(:,:)       ! latitude endpoints array
  REAL,     INTENT(IN)  :: TLON(:,:)       ! longitude endpoints array
  REAL,     INTENT(IN)  :: THGT(:,:)       ! height endpoints array
  INTEGER,  INTENT(IN)  :: NP(:,:)         ! number of endpoints per trajectory
  REAL,     INTENT(OUT) :: PBOT,PTOP       ! height limits

!-------------------------------------------------------------------------------

  INTEGER :: ic,jc,ig,jg,kt,kn,nhit,ipass,kret 
  REAL    :: plat,plon 

!-------------------------------------------------------------------------------
! two pass loop to refine trajectory domain

  maploop : DO IPASS=1,2

     IF(.NOT.ALLOCATED(LATLON))THEN
!       define grid for trajectory map optimization with default on first pass
        ALLOCATE(LATLON(KLAT,KLON),STAT=KRET)
        IF(KRET.NE.0)THEN
           WRITE(*,*)'Memory allocation error for internal grid:',KLAT,KLON
           WRITE(*,*)'Please report problem ... include endpoints file'
           STOP 900
        END IF
     END IF

     LATLON=0 ! zero trajectory position counter for mapping

!    add source points  
     sloop : DO KT=1,NTRAJ
        IF(OLAT(KT).EQ.99.0.AND.OLON(KT).EQ.99.0) CYCLE sloop
!       convert to global coordinates
        PLAT=OLAT(KT)
        PLON=OLON(KT)
        IF(PLON.LT.0.0)PLON=PLON+360.0
        IG=MIN(KLAT,INT((PLAT-GCLAT)/GINC)+1)
        JG=MIN(KLON,INT((PLON-GCLON)/GINC)+1)
!       count hits in cells
        LATLON(IG,JG)=LATLON(IG,JG)+1
     END DO sloop

!    find trajectory hits
     NHIT=0
     hloop : DO KT=1,NTRAJ
        IF(NP(KT,1).EQ.0.OR.NP(KT,2).EQ.0) CYCLE hloop

        DO KN=NP(KT,1),NP(KT,2)
!          convert to global coordinates
           PLAT=TLAT(KN,KT)
           PLON=TLON(KN,KT)
           IF(PLON.LT.0.0)PLON=PLON+360.0
           IG=MIN(KLAT,INT((PLAT-GCLAT)/GINC)+1)
           JG=MIN(KLON,INT((PLON-GCLON)/GINC)+1)
!          count hits in cells
           LATLON(IG,JG)=LATLON(IG,JG)+1
           NHIT=NHIT+1
        END DO
     END DO hloop

     IF(NHIT.EQ.0)THEN
        WRITE(*,*)'ERROR: no trajectories to plot'
        STOP 900
     END IF

!    first pass only refine grid for small plumes
     IF(IPASS.EQ.1)THEN
!       determine maximum plume extent in grid points
        IC=KLAT
        PLAT=0.0
        DO IG=1,KLAT
           IF(SUM(LATLON(IG,:)).NE.0)THEN
              PLAT=PLAT+GINC
              IC=MIN(IC,IG) 
           END IF
        END DO

        JC=KLON
        PLON=0.0
        DO JG=1,KLON
           IF(SUM(LATLON(:,JG)).NE.0)THEN
              PLON=PLON+GINC
              JC=MIN(JC,JG)  
           END IF
        END DO

!       for very small maps set up a finer grid
        IF(PLAT.LE.2.0.AND.PLON.LE.2.0)THEN
!          new corner point based upon minimum
           GCLAT=(IC-1)*GINC+GCLAT
           GCLON=(JC-1)*GINC+GCLON
           GINC=0.10
           KLAT=PLAT/GINC
           KLON=PLON/GINC
           DEALLOCATE (LATLON)
        ELSE
           EXIT maploop
        END IF
     END IF

  END DO maploop

!-------------------------------------------------------------------------------
! find vertical trajectory limits

  IF(KAGL.EQ.0)THEN
     PTOP=1000.0
     PBOT=0.0
  ELSE
     PTOP=0.0
     PBOT=100000.0
  END IF

  vloop : DO KT=1,NTRAJ

     IF(NP(KT,1).EQ.0.OR.NP(KT,2).EQ.0) CYCLE vloop

     DO KN=NP(KT,1),NP(KT,2)
!       save top height
        IF(KAGL.EQ.0)THEN
           PTOP=AMIN1(PTOP, THGT(KN,KT))
           PBOT=AMAX1(PBOT, THGT(KN,KT))
        ELSE
           PTOP=AMAX1(PTOP, THGT(KN,KT))
           PBOT=AMIN1(PBOT, THGT(KN,KT))
        END IF
     END DO

  END DO vloop

END SUBROUTINE trjlmt
