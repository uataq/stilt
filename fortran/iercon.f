!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  IERCON           Computes diagnostic smog species
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!     Integrated Empirical Rate model for CONcentrations uses the
!     instantaneous concentrations of smog, nox, and vocs,
!     summed over all puffs/particles, to compute the current
!     concentration (ppm) of O3, NO, and NO2.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 30 Nov 1998 (RRD)
!                 24 Jul 2003 (RRD) - fortran90 upgrade
!                 07 Aug 2003 (RRD) - conc grid = meteo grid
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 17 Oct 2007 (RRD) - forced positive time step
!
! USAGE:  CALL IERCON(CONC,DIRT,JET,NUMGRD,NUMTYP,DT,RK1,CSUM)
!   INPUT ARGUMENT LIST: 	see below
!   OUTPUT ARGUMENT LIST:	see below
!   INPUT FILES:		none
!   OUTPUT FILES:		none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE IERCON(CONC,DIRT,JET,NUMGRD,NUMTYP,DT,CSUM)

  USE funits
  USE metval
  USE iercfg
  use module_defgrid

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'

  TYPE(cset), INTENT(IN)    :: conc(:)           ! for each concentration grid
  TYPE(pset), INTENT(IN)    :: dirt(:)           ! for each pollutant type
  INTEGER,    INTENT(IN)    :: jet               ! elapsed time (min)
  INTEGER,    INTENT(IN)    :: numgrd            ! number of concentration grids
  INTEGER,    INTENT(IN)    :: numtyp            ! number of pollutants 
  REAL,       INTENT(IN)    :: dt                ! time step (min)
  REAL,       INTENT(INOUT) :: csum (:,:,:,:,:)  ! concentration sum matrix
!                                                ! (x,y,z,grids,species)

  INTEGER :: ii,jj,kk,kp,ks,ka,kg,kt
  INTEGER :: kvoc,knob,kspt,knox,kno2,knot,ko3t = 0

  REAL    :: ea,sea,dswf,rk1,temp,tfact
  REAL    :: tlat,tlon,csp,cnx,bno,cno,cn2,co3,rk3,const,term,quad


!-------------------------------------------------------------------------------
! check definitions for consistency

  KS=1  ! snapshot grid with last time period concentration
  KA=2  ! average concentration for results
  IF(CONC(KA)%SNAP.EQ.1.OR.CONC(KS)%SNAP.EQ.0)THEN
     WRITE(KF21,*)'*ERROR iercon: redefine snapshot and averaging grids'
     WRITE(KF21,*)' Grid #1 must be snapshot and Grid #2 is output average'
     STOP
  END IF

! averaging to grid valid this time
  IF((JET.LT.CONC(KA)%START%MACC).OR.(JET.GE.CONC(KA)%STOP%MACC))RETURN

! define species index numbers
  DO KK=1,NUMTYP
     IF(DIRT(KK)%IDENT.EQ.'VOC ')KVOC=KK
     IF(DIRT(KK)%IDENT.EQ.'NO# ')KNOB=KK
     IF(DIRT(KK)%IDENT.EQ.'SP  ')KSPT=KK
     IF(DIRT(KK)%IDENT.EQ.'NOx ')KNOX=KK
     IF(DIRT(KK)%IDENT.EQ.'NO2 ')KNO2=KK
     IF(DIRT(KK)%IDENT.EQ.'NO  ')KNOT=KK
     IF(DIRT(KK)%IDENT.EQ.'O3  ')KO3T=KK
  END DO

  IF(KVOC.EQ.0.OR.KNOB.EQ.0.OR.KSPT.EQ.0.OR.KNOX.EQ.0.OR.  &
     KNO2.EQ.0.OR.KNOT.EQ.0.OR.KO3T.EQ.0)THEN
     WRITE(KF21,*)'*ERROR* iercon: inappropriate application'
     WRITE(KF21,*)' Species not properly defined: ',DIRT(:)%IDENT 
     STOP
  END IF

!-------------------------------------------------------------------------------
! compute dependent chemistry as grid cell averages

! grid time summation averaging factor
  TFACT=ABS(DT)/CONC(KA)%DELTA%MACC
  KP=POINT(2)   ! set temporal pointer to the next time period
  KG=1          ! initial grid index set in main program
  KT=1          ! initial temporal index

  DO JJ=1,GRID(KG,KT)%NY
  DO II=1,GRID(KG,KT)%NX

!    temperature and downward flux at that grid cell
     TEMP=A(II,JJ,1,KP,KG)
     DSWF=DS(II,JJ,KP,KG)

!    convert grid units
     IF(GRID(KG,KT)%LATLON)THEN
!       CALL GBL2XY(KG,KT,TLAT,TLON,XP,YP)
        CALL GBL2LL(KG,KT,FLOAT(II),FLOAT(JJ),TLAT,TLON)
     ELSE
!       CALL CLL2XY(GRID(KG,KT)%GBASE,TLAT,TLON,XP,YP)
        CALL CXY2LL(GRID(KG,KT)%GBASE,FLOAT(II),FLOAT(JJ),TLAT,TLON)
     END IF

!    NO2 photolysis rate coefficient
     CALL SUNANG(JET,TLAT,TLON,EA,SEA)
     CALL IERNO2(EA,DSWF,RK1)

     level : DO KK=1,CONC(KA)%LEVELS

!       smog produced in grid cell
        CSP=CSUM(II,JJ,KK,KSPT,KS)
!       current value for NOx
        CNX=CSUM(II,JJ,KK,KNOX,KS)
!       background value for NO
        BNO=CSUM(II,JJ,KK,KNOB,KS)
!       skip computations when all species are zero
        IF(CSP+CNX+BNO.EQ.0.0) CYCLE level

!       photostationary constant NO+O3 -> NO2  (ppm-1 min-1)
        RK3=(9.4E+05/TEMP)*EXP(-1450.0/TEMP)

!       NOx losses from component converted to stable non-gaseous
!       nitrogen (nitrate particles) as linear proportion to SP
        CNX=CNX-MIN(PSTB*CSP,CNX)

!       known constant term
        CONST=CSP+BO3-BNO

!       quadaratic solution for NO
        TERM=RK1/RK3+CONST
        QUAD=0.5*SQRT(TERM*TERM+4.0*CNX*RK1/RK3)
        CNO=-0.5*TERM+QUAD

!       NO2 is taken from definition of NOx
        CN2=CNX-CNO

!       Ozone comes from photostationary state
        IF(CNO.LE.0.0.OR.RK1.LE.0.0)THEN
           CO3=AO3
        ELSE
           CO3=RK1*CN2/RK3/CNO+AO3
        END IF

!       sum dependent concentrations to averaging grid
        CSUM(II,JJ,KK,KNO2,KA)=CSUM(II,JJ,KK,KNO2,KA)+CN2*TFACT
        CSUM(II,JJ,KK,KNOT,KA)=CSUM(II,JJ,KK,KNOT,KA)+CNO*TFACT
        CSUM(II,JJ,KK,KO3T,KA)=CSUM(II,JJ,KK,KO3T,KA)+CO3*TFACT

     END DO level

  END DO
  END DO

END SUBROUTINE iercon
