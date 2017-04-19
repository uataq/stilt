!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GRSCON           Computes diagnostic smog species
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!     Calls Generic Reaction Set equations for ozone solution    
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 26 Aug 2003 (RRD) - initial version from iercon
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 17 Oct 2007 (RRD) - forced positive time step
!
! USAGE:  CALL GRSCON(CONC,DIRT,NUMGRD,NUMTYP,JET,DT,CSUM)
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

SUBROUTINE GRSCON(CONC,DIRT,NUMGRD,NUMTYP,JET,DT,CSUM)

  USE funits
  USE metval
  use module_defgrid

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'

  TYPE(cset), INTENT(IN)    :: conc(:)           ! for each concentration grid
  TYPE(pset), INTENT(IN)    :: dirt(:)           ! for each pollutant type
  INTEGER,    INTENT(IN)    :: numgrd            ! number of concentration grids
  INTEGER,    INTENT(IN)    :: numtyp            ! number of pollutants 
  INTEGER,    INTENT(IN)    :: jet               ! current time
  REAL,       INTENT(IN)    :: dt                ! time step (min)
  REAL,       INTENT(INOUT) :: csum (:,:,:,:,:)  ! concentration sum matrix
!                                                ! (x,y,z,grids,species)

  INTEGER, PARAMETER :: NEQN = 7  ! number of reactions
  REAL               :: CON0(NEQN), CON1(NEQN)
  REAL               :: KSR(NEQN)

  INTEGER :: ii,jj,kk,kp,kg,kt,kz,kintr,kchem,nlvl

  REAL    :: aa,bb,cc,ea,sea,dswf,rk1,temp,diftm,zenith
  REAL    :: zx,xp,yp,tlat,tlon,zbot,ztop,zagl,dist,ciso

  COMMON /ZZTOKK/ AA,BB,CC

!-------------------------------------------------------------------------------
! check definitions for consistency

  IF(CONC(1)%SNAP.EQ.0.OR.CONC(2)%SNAP.EQ.1.OR.CONC(3)%SNAP.EQ.0)THEN
     WRITE(KF21,*)'*ERROR gsrcon: redefine snapshot and averaging grids'
     WRITE(KF21,*)' Grid #1,#3 must be snapshot and Grid #2 is output average'
     STOP
  END IF

  KK=0
! define species index numbers
  IF(DIRT(1)%IDENT.NE.'NO2 ')KK=-1
  IF(DIRT(2)%IDENT.NE.'NO# ')KK=-1
  IF(DIRT(3)%IDENT.NE.'O3  ')KK=-1
  IF(DIRT(4)%IDENT.NE.'SGN ')KK=-1
  IF(DIRT(5)%IDENT.NE.'SNGN')KK=-1
  IF(DIRT(6)%IDENT.NE.'VOC ')KK=-1
  IF(DIRT(7)%IDENT.NE.'ISOP')KK=-1

  IF(KK.LT.0)THEN 
     WRITE(KF21,*)'*ERROR* grscon: inappropriate application'
     WRITE(KF21,*)' Species not defined in the proper order' 
     WRITE(KF21,*)' 1=no2 2=no 3=o3 4=sgn 5=sngm 6=voc 7=isop'
     STOP
  END IF

!-------------------------------------------------------------------------------
! compute dependent chemistry as grid cell averages

! meteorology file pointers             
  NLVL=SIZE(A,3)  ! number of meteo levels
  KP=POINT(2)     ! set temporal pointer to the next time period
  KG=1            ! initial grid index set in main program
  KT=1            ! initial temporal index

! number of chemistry time steps per meteo time step
  DIFTM=MIN(ABS(DT),6.0)  
  DO WHILE (MOD(INT(ABS(DT)),INT(DIFTM)).NE.0.AND.INT(DIFTM).GT.1)
     DIFTM=DIFTM-1.0
  END DO
  KCHEM=ABS(DT)/DIFTM

  jgrid : DO JJ=1,GRID(KG,KT)%NY
  igrid : DO II=1,GRID(KG,KT)%NX

!    skip if all species and levels are zero
     IF(SUM(CSUM(II,JJ,:,:,1)).EQ.0.0) CYCLE igrid

!    downward solar short wave flux at that grid cell
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

!    convert elevation angle to zenith angle
     ZENITH=90.0-EA

     ZBOT=0.0
     DO KK=1,CONC(1)%LEVELS

!       determine vertical cell midpoint height
        ZTOP=FLOAT(CONC(1)%HEIGHT(KK))
        ZAGL=0.5*(ZBOT+ZTOP)
        ZBOT=ZTOP

!       copy all species concentrations to processing array
        CON0=CSUM(II,JJ,KK,:,1)

!       run chemistry if any non-zero species exist
        IF(SUM(CON0).GT.0.0)THEN

!          compute the vertical meteorological index
           DIST=(BB*BB-4.0*AA*(CC-ZAGL))
           IF(DIST.GE.0.0)THEN
              ZX=(-BB+SQRT(DIST))/(2.0*AA)
           ELSE
              ZX=1.0
           END IF
           KZ=NINT(MIN(MAX(1.0,ZX),FLOAT(NLVL)))

!          temperature at cell level
           TEMP=A(II,JJ,KZ,KP,KG)

!          isoprene concentration
           CISO=CSUM(II,JJ,KK,7,1)

!          determine remaining reaction rates
           CALL GRATES(ZENITH,TEMP,KSR,RK1,1.0)

           DO KINTR=1,KCHEM
              CALL GCHEMS(DIFTM,CON0,CON1,KSR,1.0,24.97,CISO)
           END DO

!          new predicted values to array=3
           CSUM(II,JJ,KK,:,3)=CON1
           CSUM(II,JJ,KK,7,3)=CISO

        END IF

     END DO

  END DO igrid
  END DO jgrid

END SUBROUTINE grscon
