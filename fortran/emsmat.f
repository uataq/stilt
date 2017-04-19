!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSMAT           EMiSsion MATrix configuration
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SETS THE EMISSION MATRIX, A MORE COMPLEX VERSION OF THE POINT SOURCE
!   ROUTINES THAT SETS RATES DIFFERENTLY BY LOCATION AND POLLUTANT
!   SPECIES AND MAY ALSO VARY IN TIME ACCORDING TO THE VALUES FOUND
!   IN THE EMISSIONS INPUT FILE READ BY THIS PROGRAM
!
! PROGRAM HISTORY LOG:
!   LAST REVISED:  07 Mar 2006 (RRD) - initial version from emsset  
!                  13 Apr 2006 (RRD) - associate initial meteo griD
!                  24 May 2006 (RRD) - format change
!                  03 May 2007 (RRD) - forward/backward option
!                  27 Jun 2008 (RRD) - continue with warning if N<NLOC
!                  15 Jul 2008 (RRD) - correction to index test
!
! USAGE:  CALL EMSMAT(SPRT,NGRD,JET,KG,KT,NLOC,NUMTYP,EFILE,QTEMP,BACK)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             unit 5 or unit KF21 if input from file CONTROL
!   OUTPUT FILES:            unit KF22 if input from unit 5
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE EMSMAT(SPRT,NGRD,JET,KG,KT,NLOC,NUMTYP,EFILE,QTEMP,BACK)

  USE funits
  use module_defgrid         ! meteorology grid and file

  IMPLICIT NONE

  INCLUDE 'DEFSPRT.INC'         ! special source-emission matrix

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(qset),    INTENT(INOUT) :: sprt(:,:) ! for each source and pollutant
  INTEGER,       INTENT(IN)    :: ngrd      ! number of meteo grids 
  INTEGER,       INTENT(IN)    :: jet       ! current time      
  INTEGER,       INTENT(IN)    :: kg,kt     ! active grid number
  INTEGER,       INTENT(IN)    :: nloc      ! number of source locations
  INTEGER,       INTENT(IN)    :: numtyp    ! number of pollutant types
  CHARACTER(80), INTENT(INOUT) :: efile     ! temporal emission file name
  LOGICAL,       INTENT(INOUT) :: qtemp     ! temporal file flag
  LOGICAL,       INTENT(IN)    :: back      ! integration direction

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER       :: j,k,n,m,kk,macc,kqinc,kret,nrec
  INTEGER       :: KHR,KMN,KDH,KDM,KYR,KMO,KDA,KDUR
  REAL          :: QLAT,QLON 
  LOGICAL       :: OFFG,POFFG

  SAVE MACC,KQINC

!-------------------------------------------------------------------------------

  INTERFACE
  SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)   :: iy,im,id,ih,mn       ! date and time
  INTEGER,  INTENT(OUT)  :: macc                 ! accumulated minutes     
  END SUBROUTINE tm2min
!-------------------------------------------------------------------------------
  SUBROUTINE GBL2XY(KG,KT,CLAT,CLON,X,Y)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: kg          ! active grid number
  INTEGER, INTENT(IN)  :: kt          ! active time number
  REAL,    INTENT(IN)  :: clat,clon   ! latlon location
  REAL,    INTENT(OUT) :: x,y         ! grid position
  END SUBROUTINE GBL2XY
  END INTERFACE

!-------------------------------------------------------------------------------
! open temporal emission file if it exists to replace values from control file 
! each group of records is written for a valid time period
!-------------------------------------------------------------------------------

  IF(EFILE.EQ.'') RETURN

  IF(.NOT.QTEMP)THEN     
     INQUIRE(FILE=EFILE,EXIST=QTEMP)
     IF(.NOT.QTEMP)THEN
        EFILE=''
        WRITE(KF21,*)'WARNING emsmat: emission file specified but not found!'
        WRITE(KF21,*)'     File name: ',EFILE(1:40)
        RETURN
     END IF
     OPEN(KF32,FILE=EFILE,ACTION='READ')
     READ(KF32,*)
     READ(KF32,*)
     WRITE(KF21,*)' NOTICE emsmat: opened emission file - ',EFILE(1:40)
     MACC=0
  END IF


  eloop : DO WHILE (MACC.LE.JET)

!    define the base starting time for each emission group and where 
!    kqinc represents the emission group time increment in hours 
     READ(KF32,*,IOSTAT=KRET)KYR,KMO,KDA,KHR,KQINC,NREC

     IF(KRET.NE.0)THEN
        IF(BACK) EXIT eloop
        WRITE(KF21,*)'WARNING emsmat: emissions termination'
        EFILE=''
        QTEMP=.FALSE.
        RETURN
     END IF

!    set time for next emission cycle start time (to read more data)
     CALL TM2MIN(KYR,KMO,KDA,KHR,0,MACC)
     MACC=MACC+KQINC*60

     K=0
     N=1
!    read the records for one emission time group, individual sources
!    can start at anytime within the period of the time increment

     DO J=1,NREC/NUMTYP
     OFFG=.FALSE.

     DO M=1,NUMTYP

!       non-zero values trigger replacement of control file input values
        READ(KF32,*,IOSTAT=KRET) KYR,KMO,KDA,KHR,KMN,KDUR,QLAT,QLON,   &
             SPRT(N,M)%QLVL,SPRT(N,M)%RATE,SPRT(N,M)%AREA,SPRT(N,M)%HEAT

        IF(KRET.NE.0)THEN
           WRITE(KF21,*)'*ERROR* emsmat: emissions file premature data termination'
           WRITE(KF21,*)'Required number of data records per group: ',NREC
           WRITE(KF21,*)'Data records processed in this time group: ',K
           STOP 900
        END IF

!       emission duration hours and minutes
        KDH=KDUR/100
        KDM=KDUR-KDH*100
        KDUR=KDH*60+KDM
        IF(BACK)THEN
           CALL TM2MIN(KYR,KMO,KDA,KHR,KMN, SPRT(N,M)%STOP)
           SPRT(N,M)%START=SPRT(N,M)%STOP+KDUR
        ELSE
           CALL TM2MIN(KYR,KMO,KDA,KHR,KMN, SPRT(N,M)%START)
           SPRT(N,M)%STOP=SPRT(N,M)%START+KDUR
        END IF

        POFFG=.FALSE.
!       keep emission point if it falls on any meteorological grid
        gloop : DO KK=1,NGRD
!          emission location in grid units
           IF(GRID(KK,KT)%LATLON)THEN
              CALL GBL2XY (KK,KT,QLAT,QLON,SPRT(N,M)%XP,SPRT(N,M)%YP)
           ELSE
              CALL CLL2XY_wps (GRID(KK,KT)%GBASE,QLAT,QLON,SPRT(N,M)%XP,SPRT(N,M)%YP,GRID(KK,KT)%proj)
           END IF

!          off-grid locations should have been removed in the main program
!          they also need to be removed in this section to match CONTROL file
           IF(SPRT(N,M)%XP.LT.2.0.OR.SPRT(N,M)%XP.GT.FLOAT(GRID(KK,KT)%NX)-1.0.OR. &
              SPRT(N,M)%YP.LT.2.0.OR.SPRT(N,M)%YP.GT.FLOAT(GRID(KK,KT)%NY)-1.0) THEN
              POFFG=.TRUE.
              SPRT(N,M)%KG=0
           ELSE
              POFFG=.FALSE.
              SPRT(N,M)%KG=KK
              EXIT gloop
           END IF
        END DO gloop

!       if any pollutant type offgrid then reject all pollutants for the N location
        IF(POFFG) THEN
          OFFG=.TRUE.
          SPRT(N,M)%KG=0
        ELSE
           K=K+1
        END IF

     END DO

!    increment source counter ... any sources off grid reuse array space
     IF(.NOT.OFFG) N=N+1

     END DO

     IF(NUMTYP.EQ.1.AND.K.LT.NLOC)THEN
        WRITE(KF21,*)'WARNING emsmat: emission file has fewer locations than CONTROL'
        WRITE(KF21,*)'Emission file: ',N,'   CONTROL file: ',NLOC      
     ELSEIF(K.NE.NLOC*NUMTYP)THEN
        WRITE(KF21,*)'*ERROR* emsmat: emission file mismatch with CONTROL'
        WRITE(KF21,*)'Emission file data records per group: ',NREC
        WRITE(KF21,*)'Actual number of data records  input: ',K 
        WRITE(KF21,*)'CONTROL file  locations * pollutants: ',NLOC,' * ',NUMTYP
        STOP 900
     ELSE
        CONTINUE
     END IF
     WRITE(KF21,*)' NOTICE emsmat: emissions updated - ',KYR,KMO,KDA,KHR

  END DO eloop

  IF(BACK)THEN
     REWIND(KF32)
     READ(KF32,*)
     READ(KF32,*)
     MACC=0
  END IF

END SUBROUTINE emsmat
