!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  IERSUM           Integration of VOC to smog products
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!     Integrated Empirical Rate model SUMmation - integration of the
!     conversion of HydroCarbons (HC) to Smog Products (SP).  Calculation
!     is performed independently on each particle or puff. Conversion to
!     air concentration is performed when mass is summed to the output grid.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 30 Nov 1998 (RRD)
!                 20 Apr 1999 (RRD) - terrain compression factor
!                 24 Jul 2003 (RRD) - fortran90 upgrade
!                 12 Aug 2003 (RRD) - conc grid = meteo grid
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 17 Oct 2007 (RRD) - forced positive time step
!
! USAGE:  CALL IERSUM(CONC,DIRT,JET,NUMTYP,XP,YP,ZMDL,ZSFC,ZPOS,DT, 
!                     MASS,TEMP,RK1,CNOW)
!
!   INPUT ARGUMENT LIST:	see below
!   OUTPUT ARGUMENT LIST:	see below
!   INPUT FILES:		none
!   OUTPUT FILES:		none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE IERSUM(CONC,DIRT,JET,NUMTYP,XP,YP,ZMDL,ZSFC,ZPOS,DT,  &
                  MASS,TEMP,RK1,CNOW)

  USE funits
  USE iercfg

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'                  ! concentration structure

  TYPE(cset), INTENT(IN)    :: conc(:)   ! for each concentration grid
  TYPE(pset), INTENT(IN)    :: dirt(:)   ! for each pollutant type
  INTEGER,    INTENT(IN)    :: JET       ! integration time (min)
  INTEGER,    INTENT(IN)    :: NUMTYP    ! number of pollutants  
  REAL,       INTENT(IN)    :: XP,YP     ! position of element
  REAL,       INTENT(IN)    :: ZMDL      ! model top scaling height (m)
  REAL,       INTENT(IN)    :: ZSFC      ! surface terrain height (m)
  REAL,       INTENT(IN)    :: ZPOS      ! particle sigma position
  REAL,       INTENT(IN)    :: DT        ! time step (min)
  REAL,       INTENT(IN)    :: TEMP      ! temperature at position (K)
  REAL,       INTENT(IN)    :: RK1       ! rate coefficient for NO2 photolysis

  REAL,       INTENT(IN)    :: CNOW (:,:,:,:)  ! snapshot concentration array
  REAL,       INTENT(INOUT) :: MASS (:)        ! mass of pollutant 

  INTEGER :: ii,jj,ks,kl,nxp,nyp
  INTEGER :: kvoc,kspt,knox,kiso = 0
  REAL    :: cnx,csp,rcn,zbot,ztop,tmass,zpar

!-------------------------------------------------------------------------------
! check definitions for consistency

  NXP = SIZE (cnow,1)   ! 1st index
  NYP = SIZE (cnow,2)   ! 2nd index
  KS  = 1               ! snapshot grid with last time period concentration

  IF(CONC(KS)%SNAP.NE.1)THEN
     WRITE(KF21,*)'*ERROR iersum: define snapshot grid'
     WRITE(KF21,*)' Grid #1 must be snapshot and Grid #2 is the temporary array'
     STOP
  END IF

! averaging to grid valid this time
  IF((JET.LT.CONC(KS)%START%MACC).OR.(JET.GE.CONC(KS)%STOP%MACC))RETURN

! define species index numbers
  DO KL=1,NUMTYP
     IF(DIRT(KL)%IDENT.EQ.'SP  ')KSPT=KL
     IF(DIRT(KL)%IDENT.EQ.'VOC ')KVOC=KL
     IF(DIRT(KL)%IDENT.EQ.'NOx ')KNOX=KL
     IF(DIRT(KL)%IDENT.EQ.'ISOP')KISO=KL
  END DO

  IF(KVOC.EQ.0.OR.KSPT.EQ.0.OR.KNOX.EQ.0.OR.KISO.EQ.0)THEN
     WRITE(KF21,*)'*ERROR* iersum: inappropriate application'
     WRITE(KF21,*)' Species not properly defined: ',DIRT(:)%IDENT
     STOP
  END IF

!-------------------------------------------------------------------------------
! determine local concentrations to test for NOx limitation

! find the grid position on current concentration grid
  II=NINT(XP)
  JJ=NINT(YP)
  IF(II.GT.NXP.OR.II.LT.1.OR.   &
     JJ.GT.NYP.OR.JJ.LT.1)THEN
     MASS=0.0
     RETURN
  END IF

! convert particle position to meters
  ZPAR=(ZMDL-ZSFC)*(1.0-ZPOS)
  ZPAR=MAX(ZPAR,1.0)

! all variations use particles in vertical
  ZBOT=0
  RCN=1.0

! sum through all levels within plume
  DO KL=1,CONC(KS)%LEVELS

!    determine vertical cell sizes (input defines top)
     ZTOP=CONC(KS)%HEIGHT(KL)

!    particle falls within vertical extent of cell
     IF(ZPAR.GE.ZBOT.AND.ZPAR.LE.ZTOP)THEN

!       smog produced in grid cell
        CSP=CNOW(II,JJ,KL,KSPT)
!       current value for NOx
        CNX=CNOW(II,JJ,KL,KNOX)

!       smog is minimum between light and NOX limited regimes
!       ratio of CSP/CNX must be less than BETA for integration
!       of hydrocarbons to smog products; ie CSP=MIN(CSP,BETA*CNX)
        IF(CNX.GT.0.0)THEN
           RCN=CSP/CNX
        ELSE
           RCN=BETA
        END IF

     END IF
     ZBOT=ZTOP
  END DO

!-------------------------------------------------------------------------------
! integrate light dependent term for smog production

  IF(RCN.LT.BETA)THEN
     TMASS=(RVOC*MASS(KVOC)+RISO*MASS(KISO))*                 &
            RK1*ABS(DT)*EXP(-1000.0*GAMMA*(1.0/TEMP-1.0/316.0))
     MASS(KSPT)=MASS(KSPT)+TMASS
  END IF

END SUBROUTINE iersum
