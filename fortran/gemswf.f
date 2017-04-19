!###############################################################################
! GEMSWF - Computes the incident solar short-wave flux 
!-------------------------------------------------------------------------------
! LAST REVISED: 28 May 2008 (RRD) - initial version
!               09 Sep 2008 (RRD) - half grid point correction
!-------------------------------------------------------------------------------

SUBROUTINE gemswf (jet)               

  USE gemcfg
  USE gemvar  
  USE funits 

  IMPLICIT NONE

  INTEGER*4, INTENT(IN) :: jet  

  REAL*4,   ALLOCATABLE :: qq(:) ! fractional RH profile

  REAL*4                :: tr,ea,sea
  REAL*4                :: clat,clon,clat1,clon1,dlat,dlon
  INTEGER*4             :: i,j,nx,ny,nz,np,kgrd,kret

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  SAVE qq

!-------------------------------------------------------------------------------
  INTERFACE  
  SUBROUTINE SUNFLX(NLVL,SEA,QQ,SWF,TR)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)   :: nlvl      ! number of vertical levels in profile
  REAL,     INTENT(IN)   :: sea       ! sine of the solar elevation angle
  REAL,     INTENT(IN)   :: qq   (:)  ! RH fraction profile
  REAL,     INTENT(OUT)  :: swf       ! incident short wave flux (w/m2)
  REAL,     INTENT(OUT)  :: tr        ! transmissivity
  END SUBROUTINE sunflx
!-------------------------------------------------------------------------------
  SUBROUTINE SUNANG(JET,OLAT,OLON,EA,SEA)
  INTEGER,  INTENT(IN)   :: jet       ! elapsed minutes since January 1st 1970
  REAL,     INTENT(IN)   :: olat      ! latitude (+ = north)
  REAL,     INTENT(IN)   :: olon      ! longitude (- = west)
  REAL,     INTENT(OUT)  :: ea        ! solar elevation angle in deg at time
  REAL,     INTENT(OUT)  :: sea       ! sine of EA
  END SUBROUTINE sunang
  END INTERFACE
!-------------------------------------------------------------------------------

  IF(KINIT.LT.0)RETURN

  IF(.NOT.ALLOCATED(qq))THEN
     ALLOCATE (qq(nz), STAT=kret) 
     IF(kret.NE.0)THEN
        WRITE(KF21,'(A,I3)')'*ERROR gemswf: memory allocation - ',kret
        STOP 800       
     END IF
  END IF

  DO j=1,ny
  DO i=1,nx

     clat=dlat*(j-1)+clat1-dlat/2.0
     clon=dlon*(i-1)+clon1-dlon/2.0
     qq=mmm(i,j,:)/100.0
  
     CALL SUNANG(jet,clat,clon,ea,sea)  
     CALL SUNFLX(nz,sea,qq,swf(i,j),tr)

  END DO
  END DO

END SUBROUTINE gemswf 
