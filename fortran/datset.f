!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DATSET           DATe SETup for basic model parameters
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:00-10-11
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DATE SETUP SETS THE MOST BASIC MODEL SIMULATION PARAMETERS, BEFORE
!   ANY OTHER INFORMATION IS KNOWN. IF THERE EXISTS A FILE NAMED CONTROL
!   THEN ALL INPUT PARAMETERS WILL BE READ FROM THAT FILE.  OTHERWISE
!   INPUT IS EXPECTED ON STANDARD INPUT (UNIT 5).  IN THAT CASE AN OUTPUT
!   FILE CALLED STARTUP WILL BE CREATED THAT WILL CONTAIN ALL DATA ENTRIES
!   IT MAY BE USED IN SUBSEQUENT SIMULATIONS AS THE CONTROL FILE.
!   ENTRIES FROM STANDARD INPUT HAVE DEFAULT VALUES WHICH MAY BE SELECTED
!   BY JUST ENTERING "/".  THIS IS NOT VALID FOR CHARACTER STRINGS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 Oct 2000 (RRD) - initial version
!                 16 Mar 2001 (RRD) - interface intent change
!                 11 Feb 2002 (RRD) - formatted STARTUP file
!                 01 Apr 2004 (RRD) - generic file unit numbers
!                 15 Jun 2005 (RRD) - switched KF21 to KF25
!                 19 Aug 2008 (RRD) - optional minutes starting time
!
! USAGE:  CALL DATSET(NLOC,IBYR,IBMO,IBDA,IBHR,IBMN,IUNIT,KLEN,FNAME)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            units 5 or KF25 depending if CONTROL file found
!   OUTPUT FILES:           unit KF22 to STARTUP when input defined on unit 5
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE DATSET(NLOC,IBYR,IBMO,IBDA,IBHR,IBMN,IUNIT,KLEN,FNAME)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,      INTENT(OUT)   :: nloc   ! number of starting locations
  INTEGER,      INTENT(OUT)   :: ibyr,ibmo,ibda,ibhr,ibmn
  INTEGER,      INTENT(OUT)   :: iunit  ! unit for file of input parameters
  INTEGER,      INTENT(IN)    :: klen   ! length of unique file name string
  CHARACTER(*) ,INTENT(IN)    :: fname  ! unique file name to input data

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  LOGICAL           :: cntl
  CHARACTER(80)     :: label 

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE DECODI(IUNIT,VAR1,VAR2,VAR3,VAR4,VAR5)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IUNIT   ! unit number
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR1
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR2
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR3
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR4
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR5
  END SUBROUTINE decodi
  END INTERFACE

!-------------------------------------------------------------------------------
! create special file name for input
!-------------------------------------------------------------------------------

  IF(KLEN.LE.0)THEN
     LABEL='CONTROL'
  ELSE
     LABEL='CONTROL.'//FNAME(1:KLEN)
  END IF

! check for input data control file in local directory

  INQUIRE(FILE=LABEL,EXIST=CNTL)
  IF(CNTL)THEN
     IUNIT=KF25
     OPEN(IUNIT,FILE=LABEL)
  ELSE
     IUNIT=5
     OPEN(KF22,FILE='STARTUP')
  END IF

!-------------------------------------------------------------------------------
! starting time and location information
!-------------------------------------------------------------------------------

  IBMN=0
  IF(IUNIT.EQ.5)THEN
      IBYR=0
      IBMO=0
      IBDA=0
      IBHR=0
      WRITE(*,*)'Enter starting time (year, month, day, hour)'
      WRITE(*,*)IBYR,IBMO,IBDA,IBHR
   END IF

!# READ(IUNIT,*)IBYR,IBMO,IBDA,IBHR
   CALL DECODI(IUNIT,IBYR,IBMO,IBDA,IBHR,IBMN)

   IF(IUNIT.EQ.5)THEN
      WRITE(KF22,'(4I3)')IBYR,IBMO,IBDA,IBHR
      NLOC=1
      WRITE(*,*)'Enter number of starting locations'
      WRITE(*,*)NLOC
   END IF
   READ(IUNIT,*)NLOC
   IF(IUNIT.EQ.5)WRITE(KF22,'(I2.2)')NLOC

END SUBROUTINE datset
