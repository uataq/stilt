!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSTMP           EMiSsion TeMPoral from a point
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION INITIALIZATION STARTS A NEW PUFF OR PARTICLE IF THE
!   CURRENT MODEL TIME IS WITHIN THE TIME LIMITS SPECIFIED FOR START
!   A POLLUTANT RELEASE. EMISSION SPECIFICATIONS ARE DEFINED BY
!   AN INPUT FILE, EACH SOURCE CAN EMIT DIFFERENT POLLUTANTS EACH
!   AT DIFFERNT RATES STARTING AND ENDING AT UNIQUE TIMES.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 07 Mar 2006 (RRD) - initial version from emspnt
!                 14 Apr 2006 (RRD) - associate with initial grid
!                 03 May 2007 (RRD) - added forward/backward option
!                 17 Jan 2008 (RRD) - multi-species on one particle 
!                 19 Mar 2008 (AS)  - emission duration for rise
!                 07 Aug 2008 (RRD) - exceeding array return code
!                 07 Aug 2008 (RRD) - refined INITD test
!
! USAGE:  CALL EMSTMP(SPRT,KG,NLOC,NUMTYP,KPM,INITD,DT,JET,NUMPAR,MAXPAR,   
!                 NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,       
!                 PAGE,PTYP,PGRD,job_id,num_job,kret)
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

SUBROUTINE EMSTMP(SPRT,KG,NLOC,NUMTYP,KPM,INITD,DT,JET,NUMPAR,MAXPAR,      &
                  NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,      &
                  PAGE,PTYP,PGRD,job_id,num_job,kret)

  USE funits

  IMPLICIT NONE

  INCLUDE 'DEFSPRT.INC' ! source emissions matrix 

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(qset),INTENT(IN)    :: sprt(:,:) ! source location characteristics
  INTEGER,   INTENT(IN)    :: kg        ! active grid number
  INTEGER,   INTENT(IN)    :: nloc      ! total number of source locations
  INTEGER,   INTENT(IN)    :: numtyp    ! number of pollutant types
  INTEGER,   INTENT(INOUT) :: kpm       ! number of puffs or particles
  INTEGER,   INTENT(IN)    :: initd     ! initial distribution type
  REAL,      INTENT(IN)    :: dt        ! time step (min)
  INTEGER,   INTENT(IN)    :: jet       ! current elapsed time (min)
  INTEGER,   INTENT(IN)    :: numpar    ! maximum number of particles permitted
  INTEGER,   INTENT(IN)    :: maxpar    ! maximum particle number

  INTEGER,   INTENT(INOUT) :: nsort (:)   ! index of sorted elements
  REAL,      INTENT(INOUT) :: mass  (:,:) ! mass of pollutant (arbitrary units)
  REAL,      INTENT(INOUT) :: xpos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: ypos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: zpos  (:)   ! puff center height (sigma)
  REAL,      INTENT(INOUT) :: sigh  (:)   ! horizontal puff sigma 
  REAL,      INTENT(INOUT) :: sigu  (:)   ! turbulence u'2    
  REAL,      INTENT(INOUT) :: sigv  (:)   ! turbulence v'2
  REAL,      INTENT(INOUT) :: sigw  (:)   ! turbulence w'2 or vertical puff sigma
  INTEGER,   INTENT(INOUT) :: hdwp  (:)   ! Horizontal distribution pollutant
  INTEGER,   INTENT(INOUT) :: page  (:)   ! pollutant age since release (min)
  INTEGER,   INTENT(INOUT) :: ptyp  (:)   ! pollutant type index number
  INTEGER,   INTENT(INOUT) :: pgrd  (:)   ! meteorological grid of puff position
  INTEGER,   INTENT(IN)    :: job_id      ! mpi implementation 
  INTEGER,   INTENT(IN)    :: num_job     ! mpi implementation
  INTEGER,   INTENT(OUT)   :: kret        ! exceeding array return code  

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL, PARAMETER  :: sigr   = 1.54      ! top-hat radius
  REAL, PARAMETER  :: PI     = 3.14159265358979

  LOGICAL          :: ETEST
  REAL             :: qsum,qval,qstep
  INTEGER          :: n,m,np,npar,ktp,kp,nphr,maxdim,numpol,initk

!-------------------------------------------------------------------------------
! Pollutants are emitted on different particles or all on the same particle.
! Currently this is the only non-chemistry associated emission routine that
! permits multiple species to be emitted on the same particle by increasing the
! size of the namelist variable MAXDIM=NUMTYP. Other combinations are not 
! permitted. The EMITIMES emission file requires NUMTYP entries for each
! emission location. Other emission characteristics are assigned from the first
! species in the list as defined in the CONTROL file.
!-------------------------------------------------------------------------------

  MAXDIM = SIZE(mass,1)  ! number of pollutants on single particle

  IF(MAXDIM.NE.1.AND.MAXDIM.NE.NUMTYP)THEN
     WRITE(KF21,*)'*ERROR*: emstmp - max mass dimension other than 1 or numtyp'
     WRITE(KF21,*)' maxdim = ',maxdim,'   numtyp = ',numtyp
     WRITE(*,*)'Abnormal termination (emstmp) - see MESSAGE file'
     STOP
  END IF
  KRET=0

  IF(MAXDIM.EQ.1)THEN
!    one pollutant per particle
     NUMPOL=NUMTYP
  ELSE
!    all pollutants on one particle, then position and other characteristics
!    of species #1 is used for all the other species on the same particle
     NUMPOL=1
  END IF

!-------------------------------------------------------------------------------
! set the number of units to emit
!-------------------------------------------------------------------------------

  INITK=INITD
  IF(INITK.GE.100) INITK=MOD(INITD/10,10)

  IF(INITK.EQ.1.OR.INITK.EQ.2)THEN
!    gaussian or top-hat emissions
     NPAR=1

  ELSE
     QSUM=0.0
     DO N=1,NLOC
     DO M=1,NUMTYP
!       maximum number of emission hours over all sources and pollutants  
        QSUM=MAX(QSUM,ABS((SPRT(N,M)%STOP-SPRT(N,M)%START)/60.0),ABS(DT/60.0))
     END DO
     END DO

!    particle emission rate: number particles per hour per source
     IF(MAXDIM.EQ.1)THEN
        NPHR=NINT(FLOAT(NUMPAR/NLOC/NUMTYP)/QSUM)
     ELSE
        NPHR=NINT(FLOAT(NUMPAR/NLOC)/QSUM)
     END IF

!    particle emissions per time step (need at least one per job)
     NPAR=MAX(num_job,INT(NPHR*ABS(DT/60.0)))
  END IF

!-------------------------------------------------------------------------------
! loop through the number of independent source locations
!-------------------------------------------------------------------------------

  NP=0 ! internal particle counter

nloop : DO N=1,NLOC
tloop : DO M=1,NUMPOL

     IF(MAXDIM.EQ.1)THEN
!       check if current pollutant requires emission 
        IF(DT.GT.0)THEN
!          forward/backward integration option (03 May 2007)
           ETEST=(JET.GE.SPRT(N,M)%START.AND.JET.LT.SPRT(N,M)%STOP)
        ELSE
           ETEST=(JET.LE.SPRT(N,M)%START.AND.JET.GT.SPRT(N,M)%STOP)
        END IF
     ELSE
!       check all pollutants, if any require emissions then continue    
        ETEST=.FALSE.
        DO KTP=1,NUMTYP
        IF(DT.GT.0)THEN
           IF(JET.GE.SPRT(N,KTP)%START.AND.JET.LT.SPRT(N,KTP)%STOP)ETEST=.TRUE.
        ELSE
           IF(JET.LE.SPRT(N,KTP)%START.AND.JET.GT.SPRT(N,KTP)%STOP)ETEST=.TRUE.
        END IF
        END DO
     END IF

!    check if this pollutant requires a start
     IF(ETEST)THEN

!       multiple emissions only for particles
        ploop : DO KP=1,NPAR

!          Use for simulations with multiple processors to only emit particles 
!          when process id even multiple of particle number. In a single process
!          environment num_job=1 and job_id=0.

           NP=NP+1
           IF(MOD(NP,num_job).NE.job_id)CYCLE ploop

           KPM=KPM+1
           IF(KPM.GT.MAXPAR)THEN
              KPM=MAXPAR
              WRITE(KF21,*)'Warning: emstmp - exceeding puff limit'
              KRET=1
              RETURN    
           END IF
           NSORT(KPM)=KPM

!          initial position from main program
           XPOS(KPM)=SPRT(N,M)%XP
           YPOS(KPM)=SPRT(N,M)%YP
           ZPOS(KPM)=SPRT(N,M)%QLVL

!          horizontal variances all start at zero
           SIGU(KPM)=0.0
           SIGV(KPM)=0.0
           SIGW(KPM)=0.0 

           IF(SPRT(N,M)%AREA.LE.0.0)THEN
!             points source defined for all species
              SIGH(KPM)=0.0
           ELSE
!             defined by source compute sigma for uniform radius
              SIGH(KPM)=SQRT(SPRT(N,M)%AREA/PI)/SIGR
           END IF

!          initial distribution (see main for definitions)
           HDWP(KPM)=INITD
!          initial age at zero
           PAGE(KPM)=0
!          pollutant type=1 for multiple species on one particle  
           PTYP(KPM)=M 
!          initial grid is the default startup grid from main
           PGRD(KPM)=SPRT(N,M)%KG
!          number of time steps per hour
           QSTEP=MAX(1.0,ABS(60.0/DT))

           DO KTP=1,MAXDIM
              IF(MAXDIM.EQ.1)THEN
!                one species per particle emit per time step
                 QVAL=SPRT(N,M)%RATE/QSTEP   
              ELSE
!                multiple species on one particle emit per time step
                 QVAL=SPRT(N,KTP)%RATE/QSTEP  
              END IF
              IF(DT.GT.0)THEN
                 ETEST=(JET.GE.SPRT(N,KTP)%START.AND.JET.LT.SPRT(N,KTP)%STOP)
              ELSE
                 ETEST=(JET.LE.SPRT(N,KTP)%START.AND.JET.GT.SPRT(N,KTP)%STOP)
              END IF

              IF(ETEST)THEN
!                emit per particle
                 MASS(KTP,KPM)=QVAL/NPAR 
              ELSE
                 MASS(KTP,KPM)=0.0
              END IF
           END DO

!          special case where non-zero value indicates a plume rise calculation
!          sigw and page used as temporary variables which are then set to zero
!          in the main program after the particles are emitted before advection
           IF(SPRT(N,M)%HEAT.GT.0.0)THEN
              SIGW(KPM)=SPRT(N,M)%HEAT
              PAGE(KPM)=ABS(SPRT(N,M)%STOP-SPRT(N,M)%START)
           END IF

!      particle loop
       END DO ploop

!  start time test
   END IF

!  pollutant type loop
   END DO tloop

!  number of sources loop
   END DO nloop

END SUBROUTINE emstmp
