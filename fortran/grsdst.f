SUBROUTINE GRSDST(CONC,MASS,CSUM,KPM,ZPOS,XPOS,YPOS,PGRD,ZMDL,NUMTYP)

!------------------------------------------------------------------------------
! Mass redistribution Last Revised: 18 May 1999 
!                                   25 Aug 2003 (RRD) - fit into version 4.6
!------------------------------------------------------------------------------

  USE metval
  use module_defgrid

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  TYPE(cset),INTENT(IN)    :: conc(:)           ! for each concentration grid 
  REAL,      INTENT(INOUT) :: mass (:,:)        ! mass of pollutant 
  REAL,      INTENT(IN)    :: csum (:,:,:,:,:)  ! concentration summation matrix
  INTEGER,   INTENT(IN)    :: kpm               ! total number of puffs or particles
  REAL,      INTENT(IN)    :: zpos (:)          ! puff center height (sigma)
  REAL,      INTENT(IN)    :: xpos (:)          ! puff center positions (grid units)
  REAL,      INTENT(IN)    :: ypos (:)          ! puff center positions (grid units)
  INTEGER,   INTENT(IN)    :: pgrd (:)          ! meteorological grid of puff position
  REAL,      INTENT(IN)    :: zmdl              ! model domain top
  INTEGER,   INTENT(IN)    :: numtyp            ! number of pollutant types

!-------------------------------------------------------------------------------
! internal variable definitions
!-------------------------------------------------------------------------------

  INTEGER, ALLOCATABLE :: tnum (:,:,:)     ! particle total by cell
  INTEGER, ALLOCATABLE :: znum (:,:,:,:)   ! particles number by species

! define zero mass and concentration limits
! assume volume 50x50x1 km = 1.0E+12 m^3               
! and ppm = moles x 1.0E+03 / 1.0E+12 = moles x 1.0E-09
! If Zero ppm = 1.0E-15 THEN Zero moles = 4.0E-08
! If Zero ppm = 1.0E-12 THEN Zero moles = 4.0E-05

  REAL,      PARAMETER     ::   ZEROM   =   0.0
  REAL,      PARAMETER     ::   ZEROC   =   1.0E-12

  REAL       :: aa,bb,cc,temp,pres,dist,zagl
  REAL       :: zx,zpar,zbot,ztop,xmpd,carea,zver,diff
  INTEGER    :: l,kk,kp,kg,kt,ixp,jyp,kzp,lsp,nlvl


! common block to pass parameters for vertical grid
  COMMON /ZZTOKK/ AA,BB,CC

!-------------------------------------------------------------------------------
! compute particle summation

  IF(.NOT.ALLOCATED(tnum))THEN
     IXP=SIZE(CSUM,1)
     JYP=SIZE(CSUM,2)
     KZP=SIZE(CSUM,3)
     LSP=SIZE(CSUM,4)
     ALLOCATE (TNUM(IXP,JYP,KZP))
     ALLOCATE (ZNUM(IXP,JYP,KZP,LSP))
  END IF

! meteorology file pointers
  NLVL=SIZE(A,3)  ! number of meteo levels
  KP=POINT(2)     ! set temporal pointer to the next time period
  KG=1            ! initial grid index set in main program
  KT=1            ! initial temporal index

  TNUM=0
  ZNUM=0

  tloop : DO KP=1,KPM

!    skip terminated particles
     IF(PGRD(KP).EQ.0) CYCLE tloop 

!    find the grid position on current concentration grid
     IXP=NINT(XPOS(KP))
     JYP=NINT(YPOS(KP))

!    convert particle position to meters
     ZPAR=(ZMDL-ZT(IXP,JYP,KG))*(1.0-ZPOS(KP))
     ZPAR=MAX(ZPAR,1.0)

!    find the vertical index on concentration grid
     KZP=CONC(1)%LEVELS
     ZBOT=0.0
     DO KK=1,CONC(1)%LEVELS
        ZTOP=CONC(1)%HEIGHT(KK)
        IF(ZPAR.GT.ZBOT.AND.ZPAR.LE.ZTOP)KZP=KK
        ZBOT=ZTOP
     END DO

!    sum zero mass
     DO L=1,NUMTYP
        IF(MASS(L,KP).LE.ZEROM)ZNUM(IXP,JYP,KZP,L)=ZNUM(IXP,JYP,KZP,L)+1 
     END DO
     TNUM(IXP,JYP,KZP)=TNUM(IXP,JYP,KZP)+1

  END DO tloop 

!-------------------------------------------------------------------------------
! mass redistribution ... divide chemistry into particle mass           

  ploop : DO KP=1,KPM

!    skip terminated particles
     IF(PGRD(KP).EQ.0) CYCLE ploop 
          
!    find the grid position on current concentration grid
     IXP=NINT(XPOS(KP))
     JYP=NINT(YPOS(KP))

!    convert particle position to meters
     ZPAR=(ZMDL-ZT(IXP,JYP,KG))*(1.0-ZPOS(KP))
     ZPAR=MAX(ZPAR,1.0)

!    compute the vertical meteorological index
     ZAGL=ZMDL*(1.0-MIN(1.0,ZPAR))
     DIST=(BB*BB-4.0*AA*(CC-ZAGL))
     IF(DIST.GE.0.0)THEN
        ZX=(-BB+SQRT(DIST))/(2.0*AA)
     ELSE
        ZX=1.0
     END IF
     KZP=NINT(MIN(MAX(1.0,ZX),FLOAT(NLVL)))
     PRES=P(IXP,JYP,KZP,KG,KT)
     TEMP=A(IXP,JYP,KZP,KG,KT)

!    cell area
     CAREA=GX(IXP,JYP,KG)*GY(IXP,JYP,KG)

!    find the vertical index on concentration grid
     KZP=CONC(1)%LEVELS
     ZBOT=0.0
     DO KK=1,CONC(1)%LEVELS
        ZTOP=CONC(1)%HEIGHT(KK)
        IF(ZPAR.GT.ZBOT.AND.ZPAR.LE.ZTOP)KZP=KK
        ZBOT=ZTOP
     END DO

!    cell depth
     IF(KZP.EQ.1)THEN
        ZVER=CONC(1)%HEIGHT(KZP)
     ELSE
        ZVER=CONC(1)%HEIGHT(KZP)-CONC(1)%HEIGHT(KZP-1)
     ENDIF

!    mass adjustment ratios
     DO L=1,NUMTYP
	 
!       snapshot difference before and after chemistry
        DIFF=CSUM(IXP,JYP,KZP,L,3)-CSUM(IXP,JYP,KZP,L,1)
	 
        IF(DIFF.GT.0.0)THEN     
!          increasing concentration

           IF(MASS(L,KP).LE.ZEROM.OR.CSUM(IXP,JYP,KZP,L,1).LE.ZEROC)THEN
!             convert concentration to cell mass and divide mass equally
!             over all particles that have zero mass in the ratio of particles
!             with zero mass to the total number of particles

              MASS(L,KP)=CSUM(IXP,JYP,KZP,L,3)*ZVER*CAREA*(PRES/TEMP)/         &
                        (TNUM(IXP,JYP,KZP)*8.314E-02*1.0E+06)
           ELSE
              MASS(L,KP)=(1.0-(ZNUM(IXP,JYP,KZP,L)/TNUM(IXP,JYP,KZP)))*        &
                          MASS(L,KP)*CSUM(IXP,JYP,KZP,L,3)/CSUM(IXP,JYP,KZP,L,1)
           END IF
         
        ELSE
!          decreasing concentration
	 
           IF(CSUM(IXP,JYP,KZP,L,1).GT.0.0.AND.CSUM(IXP,JYP,KZP,L,3).GT.ZEROC)THEN 
              MASS(L,KP)=MASS(L,KP)*(CSUM(IXP,JYP,KZP,L,3)/CSUM(IXP,JYP,KZP,L,1))
           ELSE
              MASS(L,KP)=0.0
           END IF    
        END IF
	   
      END DO
           
  END DO ploop 

END SUBROUTINE grsdst
