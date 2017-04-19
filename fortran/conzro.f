!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONZRO           CONcentration array set to ZeRO
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   MATRIX ELEMENTS ARE SET TO ZERO FOR THE NEXT CYLCE ACCUMULATION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 14 Apr 2000 (RRD) - generalized start stop tests
!                 21 Sep 2000 (RRD) - fortran90 update
!                 16 Dec 2001 (RRD) - argument list change
!                 02 Nov 2001 (RRD) - shortened averaging time options
!                 06 Mar 2002 (RRD) - reversed sample start time test
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 17 Jan 2007 (RRD) - support SNAP=2 option
!                 27 Feb 2007 (RRD) - drop MTIME=0 test
!                 17 Oct 2007 (RRD) - revised start stop test 
!
! USAGE:  CALL CONZRO(CONC,NUMGRD,DT,JET,IFHR,CSUM)
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

SUBROUTINE CONZRO(conc,numgrd,dt,jet,ifhr,csum)
  
  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'           ! pollutant and concentration grid

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  TYPE(cset), INTENT(INOUT)   :: conc(:)             ! each concentration grid 
  INTEGER,     INTENT(IN)     :: numgrd              ! number of concen grids
  REAL,        INTENT(IN)     :: dt                  ! time step (min)
  INTEGER,     INTENT(IN)     :: jet                 ! current elapsed time
  INTEGER,     INTENT(IN)     :: ifhr                ! current forecast hour
  REAL,        INTENT(INOUT)  :: csum (:,:,:,:,:)    ! summation matrix

!-------------------------------------------------------------------------------
! internal variables definitions
!-------------------------------------------------------------------------------

  INTEGER                     :: ksb,ksd
  INTEGER                     :: kg,mtime
  INTEGER                     :: iyr,imo,ida,ihr,imn

  INTERFACE
!-------------------------------------------------------------------------------
  SUBROUTINE TM2DAY(MACM,IY,IM,ID,IH,MN)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)    :: macm            ! accumulate minutes
  INTEGER,   INTENT(OUT)   :: iy,im,id,ih,mn  ! current date/time
  END SUBROUTINE tm2day
!-------------------------------------------------------------------------------
  END INTERFACE

!-------------------------------------------------------------------------------
! go through each grid
!-------------------------------------------------------------------------------

  gloop : DO KG=1,NUMGRD

!   convert current time as end sampling period (3/1/2002)
    IF(CONC(KG)%SNAP.EQ.1)THEN
       CALL TM2DAY(JET+INT(DT),IYR,IMO,IDA,IHR,IMN)
    ELSE
       CALL TM2DAY(JET,IYR,IMO,IDA,IHR,IMN)
    END IF

!   sample start stop
    IF(INT(DT).GT.0)THEN
       KSB=JET-CONC(KG)%START%MACC
    ELSE
       KSB=CONC(KG)%START%MACC-JET
    END IF
    KSD=ABS(CONC(KG)%STOP%MACC-CONC(KG)%START%MACC)

!   test for time within sampling interval
    IF(KSB.LT.0.OR.KSB.GE.KSD) CYCLE gloop

!   zero data array according to the type of averaging for each grid 

    IF(CONC(KG)%SNAP.EQ.0)THEN
!      with averaging (SNAP=0) the current time must be at the output interval
       MTIME=ABS(JET-CONC(KG)%START%MACC)
       IF(MOD(MTIME,CONC(KG)%DELTA%MACC).NE.0) CYCLE gloop

    ELSEIF(CONC(KG)%SNAP.EQ.2)THEN
!      for maximum values current time must be at the output interval
       MTIME=ABS(JET-CONC(KG)%START%MACC)
       IF(MOD(MTIME,CONC(KG)%DELTA%MACC).NE.0) CYCLE gloop

    ELSEIF(CONC(KG)%SNAP.LT.0)THEN
!      with shortened averaging time option zero out summation array
!      at ABS(SNAP) minutes prior to the end of the averaging period
       MTIME=ABS(JET+ABS(CONC(KG)%SNAP)-CONC(KG)%START%MACC)
       IF(MOD(MTIME,CONC(KG)%DELTA%MACC).NE.0) CYCLE gloop

    ELSE
!      with snapshot maps (SNAP=1) always zero out array
       CONTINUE
    END IF

!   zero out that level pollutant (3/6/2002)
    csum(:,:,:,:,KG) = 0.0

!   save end time as start of next sampling period
    CONC(KG)%NOW%YR=IYR
    CONC(KG)%NOW%MO=IMO
    CONC(KG)%NOW%DA=IDA
    CONC(KG)%NOW%HR=IHR
    CONC(KG)%NOW%MN=IMN
    CONC(KG)%NOW%IC=IFHR

  END DO gloop

END SUBROUTINE conzro 
