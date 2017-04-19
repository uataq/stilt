!#########################################################################
! GEMPAR - Load particles or puffs into the global model as a source term
! Last revised: 29 May 2008 (RRD) - initial version
!               01 Jul 2008 (RRD) - changed conversion test to age  
!               09 Sep 2008 (RRD) - half grid point correction
!               01 Oct 2008 (RRD) - added poslon test
!-------------------------------------------------------------------------
!   HDWP     Initial Distribution        Conversion Property
!   0        3D particle                 None
!   1        Gaussian Puff               None
!   2        Top-Hat Puff                None
!   3        Gaussian Puff/Particle      None
!   4        Top-Hat Puff/Particle       None
!   103      3D particle                 Gaussian Puff/Particle
!   104      3D particle                 Top-Hat Puff/Particle
!   130      Gaussian Puff/Particle      3D particle
!   140      Top-Hat Puff/Particle       3D particle
!--------------------------------------------------------------------------

SUBROUTINE GEMPAR(KPM,HDWP,PTYP,PAGE,XPOS,YPOS,ZPOS,PGRD,MASS,CONAGE,ZMDL)

  USE gemkon  
  USE gemvar  
  USE gemcfg
  use module_defgrid

  IMPLICIT NONE

  INTEGER,  INTENT(IN)    :: kpm        ! total number of puffs or particles
  INTEGER,  INTENT(IN)    :: hdwp (:)   ! horizontal distribution of pollutant
  INTEGER,  INTENT(IN)    :: ptyp (:)   ! pollutant type
  INTEGER,  INTENT(IN)    :: page (:)   ! particle age (minutes) 
  REAL,     INTENT(IN)    :: xpos (:)   
  REAL,     INTENT(IN)    :: ypos (:)   
  REAL,     INTENT(IN)    :: zpos (:)  
  INTEGER,  INTENT(INOUT) :: pgrd (:)   ! meteorological grid index 
  REAL,     INTENT(INOUT) :: mass (:,:)
  INTEGER,  INTENT(IN)    :: conage     ! conversion age (minutes)  
  REAL,     INTENT(IN)    :: zmdl       ! vertical index scaling height

  LOGICAL   :: poslon = .true.          ! all positive longitudes

  REAL      :: zlevel, depth, volume
  INTEGER   :: i,j,k,m,mm,mp,nx,ny,nz,np,kp,kgrd,maxdim
  REAL      :: clat,clon,clat1,clon1,dlat,dlon

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  IF(KINIT.LT.0)RETURN

  MAXDIM = SIZE (mass,1)  ! number of species per single particle

! set the longitude system flag
  IF(clon1.LT.0.0)poslon=.false.

!----------------------------------------------------------------------------

  ploop : DO KP=1,KPM

!    check for on-grid and puff/particle conversion property 
     IF(PGRD(KP).LE.0.OR.HDWP(KP).LT.100) CYCLE ploop

!    determine if large enough for conversion (never true when conage<0) 
     IF(PAGE(KP).LT.CONAGE) CYCLE ploop

!    if particle on global grid and multiple meteo grids have been defined
!    and conage<0 ... move particle to the global model
     IF(PAGE(KP).LT.ABS(CONAGE).AND.PGRD(KP).NE.KGRD) CYCLE ploop

!    convert positions to lat/lon before dump
     IF(GRID(PGRD(KP),1)%LATLON)THEN
        CALL GBL2LL(PGRD(KP),1,XPOS(KP),YPOS(KP),CLAT,CLON)
     ELSE
        CALL CXY2LL(GRID(PGRD(KP),1)%GBASE,XPOS(KP),YPOS(KP),CLAT,CLON)
     END IF

!    adjust the longitude system to reflect meteo data
     IF(poslon)THEN
        IF(clon.LT.0.0)  clon=clon+360.0
     ELSE
        IF(clon.GT.180.0)clon=360.0-clon
     END IF

!    convert sigma to agl assuming terrain height = 0
     ZLEVEL=(1.0-ZPOS(KP))*ZMDL

!    compute horizontal particle index on GEM grid 
     j=1+(clat-clat1+dlat/2.0)/dlat
     i=1+(clon-clon1+dlon/2.0)/dlon
     IF(i.LT. 1) i=nx+i
     IF(i.GT.nx) i=i-nx

!    find vertical index
     k=1
     DO WHILE (HHH(i,j,k).LE.zlevel.AND.k.LT.nz)
        k=k+1
     END DO

     IF(k.LT.nz)THEN
        depth=HHH(i,j,k+1)-HHH(i,j,k)
     ELSE
        depth=HHH(i,j,k)-HHH(i,j,k-1)
     END IF

!    grid cell volume
     volume=GXY(i,j)*depth

     DO M=1,MAXDIM
        IF(maxdim.EQ.1)THEN
!          only one species per particle
           MM=1
           MP=PTYP(KP) 
        ELSE
!          all species on one particle
           MP=M 
           MM=M
        END IF

!       add mass to the particle position grid cell  
        IF(mp.LE.np)THEN
           XXX(i,j,k,mp)=XXX(i,j,k,mp)+MASS(MM,KP)/(volume*RRR(i,j,k))
           MASS(MM,KP)=0.0
        ELSE
!          undefined mass (should never occur)
           CYCLE ploop
        END IF
     END DO

!    turn off particle
     PGRD(KP)=0

  END DO ploop

END SUBROUTINE gempar
