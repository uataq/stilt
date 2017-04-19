!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TRJSET           TRaJectory output SETup
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TRAJECTORY OUTPUT FILE SETUP READS FILE AND DIRECTORY INFORMATION
!   OPENS END-POINTS FILE AND WRITES INDEX INFORMATION TO OUTPUT FILE
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 11 Jun 1997 (RRD)
!                  14 Dec 1998 (RRD) - multiple diagnostic output
!                  05 Jul 2000 (RRD) - generalized diagnostic output
!                  04 Sep 2000 (RRD) - fortran90 upgrade
!                  26 Sep 2001 (RRD) - simultaneous multiple meteorol
!                  21 May 2002 (RRD) - divergence option
!                  23 Jul 2002 (RRD) - terrain option
!                  02 Apr 2004 (RRD) - generic file unit numbers
!                  31 May 2005 (RRD) - updated trajectory file format
!                  29 Jan 2008 (RRD) - multi-processor option
!                  03 Jun 2008 (RRD) - embedded blanks dir/file
!
! USAGE:  CALL TRJSET(SPOT,NGRD,NTIM,NLOC,BACK,NDIAG,KVEL,IUNIT,LABEL,
!                     job_id,num_job)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             units 5 or KF21 if input comes from CONTROL file
!   OUTPUT FILES:            unit KF22 if input from unit 5
!                            unit KF11 for trajectory header information
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE TRJSET(SPOT,NGRD,NTIM,NLOC,BACK,NDIAG,KVEL,IUNIT,LABEL,   &
                  job_id,num_job)

  USE funits
  use module_defgrid ! meteorology file and grid

  IMPLICIT NONE

  INCLUDE 'DEFSPOT.INC' ! source information

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(rset),   INTENT(IN) :: spot (:)     ! source location characteristics
  INTEGER,      INTENT(IN) :: ngrd         ! number of grids defined
  INTEGER,      INTENT(IN) :: ntim         ! number of times defined
  INTEGER,      INTENT(IN) :: nloc         ! number of trajectories in file
  LOGICAL,      INTENT(IN) :: back         ! direction indicator
  INTEGER,      INTENT(IN) :: ndiag        ! number of diagnostic variables
  INTEGER,      INTENT(IN) :: kvel         ! vertical velocity method
  INTEGER,      INTENT(IN) :: iunit        ! model input
  CHARACTER(8), INTENT(IN) :: label(ndiag) ! diagnostic variable names
  INTEGER,      INTENT(IN) :: job_id       ! MPI job id (0 to proc-1)
  INTEGER,      INTENT(IN) :: num_job      ! MPI number of jobs

!-------------------------------------------------------------------------------
! internal variables      
!-------------------------------------------------------------------------------

  INTEGER, PARAMETER :: tfmt = 1           ! trajectory file formats:
                                           ! none = old format
                                           ! 0    = old format
                                           ! 1    = defined 5/31/2005

  CHARACTER(6)  :: myjob 
  CHARACTER(8)  :: dirctn(2), motion(0:6)
  CHARACTER(80) :: dir    =  '/main/trajectory/output/'
  CHARACTER(80) :: name   =  'file_name'

  INTEGER       :: kd = 1 
  INTEGER       :: n,nl,k,kk,kg,kt

!-------------------------------------------------------------------------------
! external variables      
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  DATA dirctn /'FORWARD','BACKWARD'/ 
  DATA motion /'OMEGA','ISOBARIC','THETA','DENSITY','ISOSIGMA','DIVERGE', &
               'TERRAIN'/
!-------------------------------------------------------------------------------

! set the proper file names and directory
  IF(IUNIT.EQ.5)THEN
     WRITE(*,'(A,I2,A)')' Enter trajectory output directory'
     WRITE(*,'(A)')DIR
  END IF
  READ(IUNIT,'(A)')DIR
  DIR=ADJUSTL(DIR)
  IF(IUNIT.EQ.5)WRITE(KF22,'(A)')DIR

  IF(IUNIT.EQ.5)THEN
     WRITE(*,'(A,I2,A)')' Enter trajectory output filename'
     WRITE(*,'(A)')NAME
  END IF
  READ(IUNIT,'(A)')NAME
  NAME=ADJUSTL(NAME)
  IF(IUNIT.EQ.5)WRITE(KF22,'(A)')NAME

! output to designated device
  KG=LEN_TRIM(DIR)
  KT=LEN_TRIM(NAME)
  IF(num_job.GT.1)THEN
     WRITE(MYJOB,'(A,I3.3)')'.',(job_id+1)
     OPEN(KF11,FILE=DIR(1:KG)//NAME(1:KT)//MYJOB)
  ELSE
     OPEN(KF11,FILE=DIR(1:KG)//NAME(1:KT))
  END IF

! determine number of meteorological files
  KK=0
  DO KG=1,NGRD
  DO KT=1,NTIM
     IF(GRID(KG,KT)%NUMBER.GE.0)KK=KK+1
  END DO
  END DO

! meteo file information index records
  WRITE(KF11,'(2I6)')KK,TFMT
  DO KG=1,NGRD
  DO KT=1,NTIM
     IF(GRID(KG,KT)%NUMBER.GE.0)                                           &
        WRITE(KF11,'(A8,5I6)')                                             &
        GRID(KG,KT)%MODEL_ID, FILE(KG,KT)%FIRST%YR, FILE(KG,KT)%FIRST%MO,  &
        FILE(KG,KT)%FIRST%DA, FILE(KG,KT)%FIRST%HR, FILE(KG,KT)%FIRST%IC
  END DO
  END DO

! trajectory information index records
  IF(BACK)KD=2

  IF(num_job.GT.1)THEN
!    starting points in CONTROL split amoung number of jobs
     NL=0
     DO N=1,NLOC
        IF(MOD(N,num_job).EQ.job_id)NL=NL+1
     END DO
     WRITE(KF11,'(I6,2(1X,A8))')NL,DIRCTN(KD),MOTION(KVEL)

     DO N=1,NLOC
        IF(MOD(N,num_job).EQ.job_id)THEN 
        WRITE(KF11,'(4I6,2F9.3,F8.1)') SPOT(N)%IBYR,SPOT(N)%IBMO,      &
        SPOT(N)%IBDA,SPOT(N)%IBHR,SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL
        END IF
     END DO

  ELSE
!    standard single processor trajectory simulation
     WRITE(KF11,'(I6,2(1X,A8))')NLOC,DIRCTN(KD),MOTION(KVEL)
     DO N=1,NLOC
        WRITE(KF11,'(4I6,2F9.3,F8.1)') SPOT(N)%IBYR,SPOT(N)%IBMO,      &
        SPOT(N)%IBDA,SPOT(N)%IBHR,SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL
     END DO
  END IF


! additional diagnostic meteorological information
  WRITE(KF11,'(I6,10(1X,A8))')NDIAG,(LABEL(K),K=1,NDIAG)

END SUBROUTINE trjset
