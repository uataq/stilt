!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONDSK           CONcentration output to DiSK
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONCENTRATION TO DISK WRITES OUT THE CONCENTRATION MATRIX TO
!   DISK FILES AT PRESELECTED SAMPLING INTERVALS.  AFTER OUTPUT
!   ELEMENTS ARE SET TO ZERO FOR THE NEXT CYLCE ACCUMULATION.
!   APPROPRIATE TIME INDEX RECORDS ARE WRITTEN TO EACH FILE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 14 Apr 2000 (RRD) - generalized start stop test
!                 21 Sep 2000 (RRD) - fortran90 update
!                 21 Dec 2000 (RRD) - concentration packing
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 05 Aug 2003 (RRD) - ichem=4 option for meteo grid
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 22 Dec 2004 (RRD) - non regular concentration grid
!                 01 Jun 2005 (RRD) - added minutes to the message
!                 17 Oct 2007 (RRD) - revised start/stop test  
!                 08 Jan 2008 (RRD) - single grid do not deallocate
!                 24 Jan 2008 (BS)  - revised start/stop test
!
! USAGE:  CALL CONDSK(CONC,DIRT,ICHEM,KGM,KTM,NUMGRD,NUMTYP,DT,JET,IFHR, 
!                     CPACK,CSUM)
!
!   INPUT ARGUMENT LIST:   see below
!   OUTPUT ARGUMENT LIST:  see below
!   INPUT FILES:           none
!   OUTPUT FILES:          units 21,22,etc 
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CONDSK(conc,dirt,ichem,kgm,ktm,numgrd,numtyp,dt,jet,ifhr,cpack,csum)

  USE funits
  use module_defgrid ! meteorological grid and file definitions

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC' ! pollutant and concentration grid

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  TYPE(cset),  INTENT(IN)     :: conc (:)            ! each conc grid 
  TYPE(pset),  INTENT(IN)     :: dirt (:)            ! each pollutant type 
  INTEGER,     INTENT(IN)     :: ichem               ! chemistry parameter   
  INTEGER,     INTENT(IN)     :: kgm                 ! current meteo grid number
  INTEGER,     INTENT(IN)     :: ktm                 ! current meteo time number
  INTEGER,     INTENT(IN)     :: numgrd              ! number of concen grids
  INTEGER,     INTENT(IN)     :: numtyp              ! number of pollutants
  REAL,        INTENT(IN)     :: dt                  ! time step
  INTEGER,     INTENT(IN)     :: jet                 ! current elapsed time
  INTEGER,     INTENT(IN)     :: ifhr                ! current forecast hour
  INTEGER,     INTENT(IN)     :: cpack               ! concen packing flag
  REAL,        INTENT(IN)     :: csum (:,:,:,:,:)    ! summation matrix

!-------------------------------------------------------------------------------
! internal variable definitions
!-------------------------------------------------------------------------------

  INTEGER                     :: ksb,ksd
  REAL                        :: plat,plon,xp,yp
  INTEGER                     :: ng,ki,kj,kg,kt,kl,kk
  INTEGER                     :: iyr,imo,ida,ihr,imn
  INTEGER                     :: mtime,kunit,kret,kntr
  INTEGER                     :: ii,jj,nxp,nyp,nzp

  TYPE conrec
     INTEGER(2)               :: ipnt
     INTEGER(2)               :: jpnt
     REAL                     :: conc
  END TYPE
  TYPE(conrec), ALLOCATABLE   :: condat (:)          ! nonzero concen values

!-------------------------------------------------------------------------------
! external variables definitions
!-------------------------------------------------------------------------------


  SAVE condat

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE TM2DAY(MACM,IY,IM,ID,IH,MN)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)    :: macm            ! accumulate minutes
  INTEGER,   INTENT(OUT)   :: iy,im,id,ih,mn  ! current date/time
  END SUBROUTINE tm2day
  END INTERFACE
!-------------------------------------------------------------------------------

! check for consistency
  IF(CPACK.EQ.0.AND.ICHEM.EQ.4)THEN
     WRITE(KF21,*)'*ERROR* condsk: ichem=4 requires cpack=1'
     STOP 900
  END IF

! convert current time as end sampling period
  CALL TM2DAY(JET,IYR,IMO,IDA,IHR,IMN)

!-------------------------------------------------------------------------------
! go through each grid
!-------------------------------------------------------------------------------

  NG=NUMGRD
  IF(CPACK.EQ.2)NG=1

  gloop : DO KG=1,NG

!=========================================================================
!    schematic timeline and description of variables
!    JET starts at "run start" and moves right(forward) or left(backward)
!
!                                   run                           
!                                  start                         
!                     (backward)     |      (forward)           
!    |----------<------o----->-------|----<-------o------->---------
!  (zero)      STOP   JET  START     |  START    JET     STOP     (time in min)
!
!  START = start of first sample (grid sample start time, in CONTROL file)
!  STOP  = end   of last  sample (grid sample end   time, in CONTROL file)
!  JET   = elapsed (current) time
!  DELTA = sample duration (can be many samples between START and STOP)
!  KSB   = number of minutes between START and JET  
!  KSD   = number of minutes between START and STOP (positive, by def)
!  MTIME = number of minutes between START and JET  (positive, by def)
!
!  KSB>0 --> current time is 'after' START (sampling began at least 1 min ago)
!               KSB>0 for both fwrd and back cases shown
!  KSB=0 --> current time equals START       (sampling starts now) 
!  KSB<0 --> current time is 'before ' START (sampling has not begun)
!
!  Conc is output when KSB>0, KSB<=KSD, and KSB is a multiple of DELTA
!=========================================================================

!    sample start stop
     IF(INT(DT).GT.0)THEN
        KSB=JET-CONC(KG)%START%MACC
     ELSE
        KSB=CONC(KG)%START%MACC-JET
     END IF
     KSD=ABS(CONC(KG)%STOP%MACC-CONC(KG)%START%MACC)

!    test for time within sampling interval
     IF(KSB.LE.0.OR.KSB.GT.KSD) CYCLE gloop		! 1-24-08 (BS)

!    current time must be at the output interval
     MTIME=ABS(JET-CONC(KG)%START%MACC)
     IF(MOD(MTIME,CONC(KG)%DELTA%MACC).NE.0) CYCLE gloop

     KUNIT=CONC(KG)%UNIT
     
!    write index record for sampling start time
     WRITE(KUNIT)                                            &
        CONC(KG)%NOW%YR, CONC(KG)%NOW%MO, CONC(KG)%NOW%DA,   &
        CONC(KG)%NOW%HR, CONC(KG)%NOW%MN, CONC(KG)%NOW%IC

!    write index record for sampling stop time
     WRITE(KUNIT)IYR,IMO,IDA,IHR,IMN,IFHR

!    determine loop indicies
     NXP=CONC(KG)%NUMB_LON
     NYP=CONC(KG)%NUMB_LAT
     NZP=CONC(KG)%LEVELS

!    when concentration packing
     IF(CPACK.EQ.1)THEN
        IF(.NOT.ALLOCATED(condat))THEN 
           ALLOCATE (condat(nxp*nyp),STAT=kret) 
           IF(kret.NE.0)THEN      
              WRITE(KF21,*)'*ERROR* condsk: temporary memory allocation'
              STOP 900
           END IF
        END IF
     END IF

!    write data record by pollutant type and level
     DO KT=1,NUMTYP
     DO KL=1,NZP

!       only write non-zero concentrations
        IF(CPACK.EQ.1)THEN
           KNTR=0
           DO JJ=1,NYP
    loop : DO II=1,NXP

              IF(ICHEM.NE.4)THEN
!                internal concentration grid is a lat/lon grid
                 KI=II
                 KJ=JJ
              ELSE
!                internal concentration grid matches meteorology grid
!                convert position to meteorological grid units
                 PLON=FLOAT(II-1)*CONC(KG)%DELT_LON+CONC(KG)%X1Y1_LON
                 PLAT=FLOAT(JJ-1)*CONC(KG)%DELT_LAT+CONC(KG)%X1Y1_LAT
                 IF(GRID(KGM,KTM)%LATLON)THEN
                    CALL GBL2XY(KGM,KTM,PLAT,PLON,XP,YP)
                 ELSE
                    CALL CLL2XY_wps(GRID(KGM,KTM)%GBASE,PLAT,PLON,XP,YP,GRID(KGM,KTM)%proj)
                 END IF
                 KI=NINT(XP)
                 KJ=NINT(YP)
                 IF(KI.LT.1.OR.KI.GT.GRID(KGM,KTM)%NX.OR.      &
                    KJ.LT.1.OR.KJ.GT.GRID(KGM,KTM)%NY) CYCLE loop
              END IF

              IF(CSUM(KI,KJ,KL,KT,KG).NE.0.0)THEN
                 KNTR=KNTR+1
                 CONDAT(KNTR)%IPNT=INT(II,KIND=2)
                 CONDAT(KNTR)%JPNT=INT(JJ,KIND=2)
                 CONDAT(KNTR)%CONC=CSUM(KI,KJ,KL,KT,KG)
              END IF

           END DO loop
           END DO

           IF(KNTR.GT.0)THEN
              WRITE(KUNIT)DIRT(KT)%IDENT, CONC(KG)%HEIGHT(KL), KNTR,   &
                   (CONDAT(KK)%IPNT,CONDAT(KK)%JPNT,CONDAT(KK)%CONC,KK=1,KNTR)
           ELSE
              WRITE(KUNIT)DIRT(KT)%IDENT, CONC(KG)%HEIGHT(KL), KNTR
           END IF

        ELSEIF(CPACK.EQ.2)THEN
!          station format file
           WRITE(KUNIT)DIRT(KT)%IDENT, CONC(1)%HEIGHT(KL), NUMGRD,           &
              (CONC(KNTR)%X1Y1_LON,CONC(KNTR)%X1Y1_LAT,CSUM(1,1,KL,KT,KNTR), &
               KNTR=1,NUMGRD)

        ELSE
!          write all grid points according to the old format
           WRITE(KUNIT)DIRT(KT)%IDENT, CONC(KG)%HEIGHT(KL),     &
               ((CSUM(II,JJ,KL,KT,KG),II=1,NXP),JJ=1,NYP)
        END IF

     END DO
     END DO

!    when concentration packing
     IF(CPACK.EQ.1.AND.NUMGRD.GT.1) DEALLOCATE (condat) 

!### diagnostic output for testing
!### WRITE(KF21,*)' NOTICE condsk: ',KG,' START:',CONC(KG)%START%MACC,  &
!###  ' JET:',JET,' STOP:',CONC(KG)%STOP%MACC
     WRITE(KF21,*)' NOTICE condsk: ',KG,' DELTA:',CONC(KG)%DELTA%MACC,  &
      ' KSB:',KSB,' KSD:',KSD,' MTIME:',MTIME

!    standard output
     WRITE(KF21,*)' NOTICE condsk: ',KG,' output at ',IMO,IDA,IHR,IMN

  END DO gloop

END SUBROUTINE condsk
