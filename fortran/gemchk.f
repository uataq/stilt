!###############################################################################
! GEMCHK - Analysis of the meteorological data file for the global model    
!-------------------------------------------------------------------------------
! LAST REVISED: 28 May 2008 (RRD) - initial version 
!               25 Jul 2008 (AS)  - vertical integral output
!               09 Sep 2008 (RRD) - half grid point correction
!               13 Nov 2008 (RRD) - improved namelist backup format
!-------------------------------------------------------------------------------
 
SUBROUTINE gemchk (ngrd,ntyp)

  USE gemkon
  USE gemvar     
  USE gemcfg
  USE funits
  use module_defgrid

  IMPLICIT NONE
 

  INTEGER*4, INTENT(IN)  :: NGRD  ! number of meteo grids 
  INTEGER*4, INTENT(IN)  :: NTYP  ! number of pollutants  

  LOGICAL   :: ftest
  INTEGER*4 :: i,j,k,n,kgrd
  INTEGER*4 :: nx,ny,nz,np
  REAL*4    :: clat,clon,clat1,clon1,dlat,dlon,hgts

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

!--------------------------------------
! namelist initialization 

  NAMELIST/GEMPARM/HMIX,CFACT,CKON,CMIN,CMAX,VMIX,WFACT,KDIVG,KMASS,KINIT, &
                   GEMDUMP,GEMCONC,GEMZINT,KZBOT,KZTOP,KZAVG,WFREQ,IHOUR 

  INQUIRE(FILE='GEMPARM.CFG',EXIST=FTEST)
  IF(FTEST)THEN
     OPEN(KF26,FILE='GEMPARM.CFG')
     READ(KF26,GEMPARM)
     CLOSE(KF26)
     WRITE(*,*)' NOTICE gemchk: set GEM namelist parameters from input file'  
  ELSE
     OPEN(KF26,FILE='GEMPARM.CFG')
     WRITE(KF26,'(A)')' &GEMPARM'
     WRITE(KF26,'(A,F8.1,A)')' hmix = ',hmix,','
     WRITE(KF26,'(A,F3.1,A)')' cfact = ',cfact,','
     WRITE(KF26,'(A,F3.1,A)')' ckon = ',ckon,',' 
     WRITE(KF26,'(A,F3.1,A)')' cmin = ',cmin,','
     WRITE(KF26,'(A,E7.1,A)')' cmax = ',cmax,','
     WRITE(KF26,'(A,F4.1,A)')' vmix = ',vmix,','
     WRITE(KF26,'(A,F3.1,A)')' wfact = ',wfact,','
     WRITE(KF26,'(A,I1  ,A)')' kdivg = ',kdivg,','
     WRITE(KF26,'(A,I1  ,A)')' kmass = ',kmass,','
     WRITE(KF26,'(A,I1  ,A)')' kinit = ',kinit,','
     WRITE(KF26,'(A,I1  ,A)')' kzbot = ',kzbot,','
     WRITE(KF26,'(A,I1  ,A)')' kztop = ',kztop,','
     WRITE(KF26,'(A,I1  ,A)')' kzavg = ',kzavg,','
     WRITE(KF26,'(A,I1  ,A)')' wfreq = ',wfreq,','
     WRITE(KF26,'(A,I1  ,A)')' ihour = ',ihour,','
     K=INDEX(GEMDUMP,' ')-1
     WRITE(KF26,'(3A)')' gemdump = ',''''//gemdump(:k),''','
     K=INDEX(GEMCONC,' ')-1
     WRITE(KF26,'(3A)')' gemconc = ',''''//gemconc(:k),''','
     K=INDEX(GEMZINT,' ')-1
     WRITE(KF26,'(3A)')' gemzint = ',''''//gemzint(:k),''','
     WRITE(KF26,'(A)')' /'
     CLOSE(KF26)
     WRITE(*,*)' *ERROR* gemchk: GEMPARM.CFG namelist file not found!' 
     STOP 800 
  END IF
  IF(KINIT.LT.0)RETURN

!-------------------------------------
! pollutant initialization 
! currently each pollutant is defined in the array space XXX(:,:,:,np) 

  NP=NTYP

!------------------------------------
! meteorology initialization


! find the first global meteorology grid
  KGRD=0
  loop1 : DO N=1,NGRD
!     IF(GRID(N,1)%GBLDAT)THEN
     IF(GRID(N,1)%GBLDAT.EQ."gl")THEN
        KGRD=N
        EXIT loop1
     END IF
  END DO loop1

  IF(KGRD.EQ.0)THEN
     WRITE(KF21,*)'*ERROR* gemchk: global meteorology file not defined'
     STOP 800 
  END IF

  IF(DREC(KGRD,1)%Z_FLAG.NE.2)THEN
     WRITE(KF21,*)'*ERROR* gemchk: meteorological coordinate not pressure'
     STOP 800
  END IF

! define various grid parameters 
  NX   =GRID(KGRD,1)%NX
  NY   =GRID(KGRD,1)%NY
  NZ   =GRID(KGRD,1)%NZ-1
  CLAT1=GRID(KGRD,1)%SYNC_LAT
  CLON1=GRID(KGRD,1)%SYNC_LON
  DLAT =GRID(KGRD,1)%REF_LAT
  DLON =GRID(KGRD,1)%REF_LON

! meteorology data file unit number
  KUMET=FILE(KGRD,1)%KUNIT 

  IF(DLON.NE.1.0.AND.DLON.NE.2.5)THEN
     WRITE(KF21,*)'*ERROR* gemchk: use only 1.0 or 2.5 deg lat-lon grids'
     STOP 800 
  END IF

! allocate all the meteorological data array space
  CALL gemset 
      
! note that level=0 represents surface fields not used in the calculation
  DO N=1,(NZ+1)
     PPP(N-1) =DREC(KGRD,1)%HEIGHT(N)
     NVAR(N-1)=DREC(KGRD,1)%NUM_VARB(N)
  END DO

  DO K=1,NZ
!    define STP sigma coordinate for each GEM level
     SIG(K)=(PPP(K)-PTOP)/(PSFC-PTOP)

!    set the standard atm index value for each level
     KNDX(K)=0
     loop2 : DO N=1,38
        IF(PPP(K).EQ.STDPPP(N))THEN
           KNDX(K)=N
           EXIT loop2 
        END IF
     END DO loop2  

     IF(KNDX(K).EQ.0)THEN
        WRITE(KF21,*)'*ERROR* gemchk: input level not in std atm - ',K,PPP(K)
        STOP 800 
     END IF
  END DO

!------------------------------------
! surface parameters initialization

  DO j=1,ny
  DO i=1,nx
     clat=dlat*(j-1)+clat1-dlat/2.0
     clon=dlon*(i-1)+clon1-dlon/2.0
     IF(clon.GT.180.0)clon=clon-360.0
     CALL SFCINP(CLAT,CLON,RGH(i,j),LUS(i,j),.FALSE.,HGTS)
  END DO
  END DO

END SUBROUTINE gemchk
