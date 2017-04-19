!###############################################################################
! GEMDEP - Global Eulerian Model DEPosition routines                       
!------------------------------------------------------------------------------
! From sub depset the following logical definitions apply to this routine 
! dogas = true : particle parameters are NOT defined 
! dodry = true : deposition velocity is defined as none zero
! dogrv = true : particle parameters are defined
! note that when dodry = true, sets dogrv = false, hence no calculation of
! settling velocity, but if dogas = false, particles still settle with VG=VD 
! 
! LUS - LANDUSE CATEGORY        index
! RGH - ROUGHNESS LENGTH        m
! FRV - FRICTION VELOCITY       m/s
! PSI - INTEGRATED STABILITY
! SWF - SHORT WAVE FLIX         w/m2
! U10 - 10 meter wind           m/s
! V10 - 10 meter wind           m/s
! T02 - 2  meter temperature    deg-K
! P6H - precipitation rate      m/min
!-------------------------------------------------------------------------------
! LAST REVISED: 29 May 2008 (RRD) - initial version
!               09 Sep 2008 (RRD) - half grid point correction
!               02 Feb 2009 (RRD) - restructure loop for use in gemeqn
!-------------------------------------------------------------------------------

  SUBROUTINE gemdep (DIRT,N,IMO)

  USE gemkon   
  USE gemvar
  USE gemcfg

  IMPLICIT NONE

  INCLUDE     'DEFCONC.INC'           ! concentration and pollutant structure
  TYPE(pset),  INTENT(IN) :: DIRT (:) ! characteristics for each pollutant type
  INTEGER*4,   INTENT(IN) :: N        ! pollutant index number 
  INTEGER*4,   INTENT(IN) :: IMO      ! current month                      

  REAL*8    :: beta,betad,betaw
  REAL*8    :: flux,fluxd,fluxw
  REAL*8    :: dzi,decay
  REAL*4    :: vc,vd,vg,dens,pdia,frea

  REAL*4    :: clat,clat1,clon1,dlat,dlon

  INTEGER*4 :: i,j,k                  ! basic loop indicies for 3D grid
  INTEGER*4 :: kbot,ktop              ! cloud layer indicies 
  INTEGER*4 :: nx,ny,nz,np,kgrd       ! 3D grid dimensions

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE DEPDRY(DIRT,OLAT,IBMO,KPOL,LAND,ROUG,SFCL,USTR,PSI,SFLX,AIRD,   &
                  TEMP,PDIA,VG,VD)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'                 ! concentration and pollutant structure
  TYPE(pset), INTENT(IN)  :: dirt(:)    ! for each pollutant type
  REAL,     INTENT(IN)    :: olat       ! origin latitude
  INTEGER,  INTENT(IN)    :: ibmo       ! computational month
  INTEGER,  INTENT(IN)    :: kpol       ! polluant index number
  INTEGER,  INTENT(IN)    :: land       ! land-use category
  REAL,     INTENT(IN)    :: roug       ! aerodynamic roughness length (m)
  REAL,     INTENT(IN)    :: sfcl       ! height of constant flux layer (m)
  REAL,     INTENT(IN)    :: ustr       ! friction velocity (m/s)
  REAL,     INTENT(IN)    :: psi        ! integrated stability function heat
  REAL,     INTENT(IN)    :: sflx       ! solar irradiation at sfc (watts/m2)
  REAL,     INTENT(IN)    :: aird       ! ambient air density (g/m3)
  REAL,     INTENT(IN)    :: temp       ! canopy air temperature (deg K)
  REAL,     INTENT(IN)    :: pdia       ! particle diameter (meters)
  REAL,     INTENT(IN)    :: vg         ! gravitational settling velocity (m/s)
  REAL,     INTENT(OUT)   :: vd         ! total deposition velocity (m/s)
  END SUBROUTINE depdry
  END INTERFACE
!-------------------------------------------------------------------------------

  IF(KINIT.LT.0)RETURN

! radioactive decay
  IF(DIRT(n)%DORAD)THEN
!    convert half-life (days) to time constant (1/min)
     DECAY=LOG(0.5)/(DIRT(n)%RHALF*1440.0)
     DECAY=EXP(DELTA*DECAY)
  END IF

! default dry deposition (vd) / gravitational settling (vg)
  VD=0.0
  VG=0.0
  IF(DIRT(n)%DODRY)THEN
!    explicit definition of the dry deposition velocity
     VD=DIRT(n)%DRYVL
     VG=VD 
!    turn off gravitational settling if gas 
     IF(DIRT(n)%DOGAS) VG=0.0
  END IF

! horizonal grid loop
  DO J=1,ny
  DO I=1,nx

!    analyze for cloud layers to compute wet removal
     IF(DIRT(n)%DOWET.AND.P6H(i,j).GT.0.0)THEN
!       determine bottom and top of the precip layer (80% to 60%)
        kbot=0
        ktop=nz   
        DO k=1,nz    
           IF(kbot.EQ.0.AND.MMM(i,j,k).GE.80.0)kbot=k 
           IF(kbot.NE.0.AND.ktop.EQ.nz.AND.MMM(i,j,k).LE.60.0)ktop=k 
        END DO

!       default layer (to be consistent with rain even if no cloud)
        IF(kbot.EQ.0)THEN
           DO k=1,nz   
              IF(HHH(i,j,k).LT. 500.0)kbot=k
              IF(HHH(i,j,k).LT.3000.0)ktop=k
           END DO
        END IF
    END IF

! zero out the deposition flux profile
  ZZZ = 0.0

! profile loop
  DO K=nz,1,-1

!    radioactive decay
     IF(DIRT(n)%DORAD) QQQ(i,j,k)=QQQ(i,j,k)*DECAY

!----------------------------------------------------
!    compute gravitational settling for particles
!    when not already explicitly defined 

     IF(DIRT(n)%DOGRV)THEN
!       particle density from g/cc g/m3
        DENS=DIRT(n)%PDENS*1.0E+06
!       particle diameter (um) to (m)
        PDIA=DIRT(n)%PDIAM*1.0E-06
!       settling velocity where air density (kg/m3) --> (g/m3)
        VG=PDIA*PDIA*GRAV*(DENS-RRR(i,j,k)*1000.0)/(18.0*DMVC)
!       slip correction (mean free path = particle diameter)
        FREA=FREP*(ROW/RRR(i,j,k))
        VC=1.0+(2.0*FREA/PDIA)*(1.26+0.4*EXP(-0.55*PDIA/FREA))
!       final value apply shape correction factor
        VG=VG*VC/DIRT(n)%SHAPE
     END IF

!----------------------------------------------------
!    resistance deposition method

     IF(DIRT(n)%DORES.AND.k.EQ.1)THEN
!       non-settling deposition only occurs in the lowest cell
        clat=(j-1)*dlat+clat1-dlat/2.0
        CALL DEPDRY(DIRT,CLAT,IMO,N,LUS(i,j),RGH(i,j),REAL(HHH(i,j,k)),   &
                    FRV(i,j),PSI(i,j),SWF(i,j),REAL(RRR(i,j,k)*1000.0),   &
                    REAL(TTT(i,j,k)),PDIA,VG,VD)
     END IF

!    grid cell height
     IF (k.GT.1.AND.k.LT.nz) THEN
        DZI=MAX(DZMIN,0.5*(HHH(i,j,k+1)-HHH(i,j,k-1)))
     ELSEIF (k.EQ.1) THEN
        DZI=MAX(DZMIN,HHH(i,j,2)-HHH(i,j,1))
     ELSE
        DZI=MAX(DZMIN,HHH(i,j,nz)-HHH(i,j,nz-1))
     END IF

!    use layer depth to compute removal fraction
!    compute rate constant 1/min
     BETAD=0.0
!    settling applies over all levels
     IF(VG.GT.0.0)            BETAD=60.0*VG/DZI                            
!    dry deposition only from the lowest layer 
     IF(VD.GT.0.0.AND.K.EQ.1) BETAD=60.0*VD/DZI                            

!-------------------------------------------------------
!    simplified wet removal (layer either below or in-cloud)

     BETAW=0.0
     IF(DIRT(n)%DOWET.AND.P6H(i,j).GT.0.0)THEN
        IF(DIRT(n)%DOGAS)THEN
!          equilibrium concentration (mass units) for gases
!          applies as long as material is below cloud top
           IF(k.LT.ktop) BETAW=DIRT(n)%WETGAS*RGAS   &
                              *TTT(i,j,k)*P6H(i,j)/DZI

        ELSE
!          particle wet removal
           IF(k.lt.kbot)THEN
              BETAW=+DIRT(n)%WETLO*60.0        !  below cloud 
           ELSEIF(k.LE.ktop)THEN
              BETAW=DIRT(n)%WETIN*P6H(i,j)/DZI !  within cloud 
           ELSE
              CONTINUE                         !  above cloud  
           END IF
        END IF
     END IF

!----------------------------------------------------
!    mass removal fluxes

     FLUX= 0.0  
     FLUXD=0.0
     FLUXW=0.0
     BETA=BETAD+BETAW
     IF(BETA.EQ.0.0)THEN
        CONTINUE
     ELSEIF(DELTA*BETA.LT.0.01)THEN
!       small removal values assume linear approximation
        FLUX =QQQ(i,j,k)*DELTA*BETA
        FLUXW=QQQ(i,j,k)*DELTA*BETAW
        FLUXD=FLUX-FLUXW
     ELSE
!       apply exponential for removal > 1%
        FLUX =QQQ(i,j,k)*(1.0-EXP(-DELTA*BETA))
        FLUXW=QQQ(i,j,k)*(1.0-EXP(-DELTA*BETAW))
        FLUXD=FLUX-FLUXW
     END IF

     IF(FLUX.GT.0.0)THEN
        FLUX=MIN(QQQ(i,j,k),FLUX)  ! can't remove more than exists
!       mass decrease due to both dry and wet removal
        QQQ(i,j,k)=QQQ(i,j,k)-FLUX

        IF(k.GT.1)THEN
!          save dry flux (gravitational settling adds to lower cell)  
           ZZZ(k)=FLUXD
!          only wet removal accumulates when above the lowest cell
           DDD(i,j,n)=DDD(i,j,n)+FLUXW*DZI
        ELSE
!          all removal from lowest box accumulates
           DDD(i,j,n)=DDD(i,j,n)+FLUX*DZI
        END IF
     END IF
  END DO

  DO k=1,(nz-1)
!    what is removed at cell(k+1) is added to cell(k)
     QQQ(i,j,k)=QQQ(i,j,k)+ZZZ(k+1)
  END DO  

  END DO
  END DO

END SUBROUTINE gemdep
