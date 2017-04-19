!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  STBHOR           STaBility HORizonal
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:98-12-17
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   COMPUTES THE HORIZONTAL MIXING COEFFICIENT BASED UPON THE
!   SMAGORINSKY FORMULATION OF THE HORIZONTAL WIND FIELD DEFORMATION
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 12 Dec 1998 (RRD)
!                  29 Sep 2000 (RRD) - fortran90 upgrade
!                  09 Mar 2001 (RRD) - global lat lon grid option
!                  06 Nov 2003 (RRD) - convert to turbulent velocity
!                  11 Dec 2003 (RRD) - zero velocity test
!
! USAGE:  CALL STBHOR(NXS,NYS,NLVL,GX,GY,U,V,H,E)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

!dwen(20090806) SUBROUTINE STBHOR(NXS,NYS,NLVL,GX,GY,U,V,H,E)
SUBROUTINE STBHOR(NXS,NYS,NLVL,GX,GY,U,V,H,hm,E)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,     INTENT(IN)    :: nxs,nys          ! horizontal subgrid dimensions
  INTEGER,     INTENT(IN)    :: nlvl             ! number of output levels
  REAL,        INTENT(IN)    :: gx(:,:)          ! horizontal grid spacing
  REAL,        INTENT(IN)    :: gy(:,:)          ! horizontal grid spacing
  REAL,        INTENT(IN)    :: u (:,:,:)        ! horizontal wind component
  REAL,        INTENT(IN)    :: v (:,:,:)        ! horizontal wind component
  REAL,        INTENT(OUT)   :: h (:,:,:)        ! u component turbulence
  REAL,        INTENT(OUT)   :: e (:,:,:)        ! v component turbulence

!dwen(20090806) ***********************
!              add hm as horizontal mixing coefficient
  real,        intent(out)   :: hm(:,:,:)
!dwen           ***********************

  REAL,        PARAMETER     :: hmin = 432.0     ! Mixing minimum
  REAL,        PARAMETER     :: hmax = 43200.0   ! Mixing maximum
  REAL,        PARAMETER     :: hscale = 10800.0 ! horz time scale (sec)

!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------

  INTEGER                    :: i,j,k,ip,jp,im,jm
  REAL                       :: uvm2,uvp2,hmix,dudy,dvdx,dudx,dvdy,deft

!-------------------------------------------------------------------------------

! process each node on subgrid
  DO K=1,NLVL
  DO J=1,NYS
  DO I=1,NXS

!    to obtain the horizontal diffusivity (from Smagorinsky, 1963)
!    from velocity deformation (shearing and tension stress)
!    centered differences and all horizontal velocities in grid/min
!    (u2-u1)/(2 delta-x), where x = 1 grid length

!    check edge limits (P-plus and M-minus)
     IP=MIN(I+1,NXS)
     JP=MIN(J+1,NYS)
     IM=MAX(I-1,1)
     JM=MAX(J-1,1)

!    sample code for divergence (min^-1) correction
!    DELU=0.5*(U(IP,J,K)-U(IM,J,K))
!    DELV=0.5*(V(I,JP,K)-V(I,JM,K))
!    DIVG=DELU+DELV

!    deformation (1/min)
     DUDY=U(I,JP,K)-U(I,JM,K)
     DVDX=V(IP,J,K)-V(IM,J,K)
     DUDX=U(IP,J,K)-U(IM,J,K)
     DVDY=V(I,JP,K)-V(I,JM,K)
     DEFT=0.5*SQRT((DUDY+DVDX)*(DUDY+DVDX)+(DUDX-DVDY)*(DUDX-DVDY))

!    deformation in 1/min => mixing m2/sec
     HMIX=0.014*GX(I,J)*GY(I,J)*DEFT/60.0

!    check limits
     HMIX=MAX(HMIN, MIN(HMAX, HMIX))

!dwen(20090806) *******************
     hm(i,j,k)=hmix
!dwen  ***************************
!    Convert mixing coefficient to turbulent velocity variance
     UVP2=HMIX/HSCALE

!    then proportion according to the mean velocity components
     UVM2=U(I,J,K)*U(I,J,K)+V(I,J,K)*V(I,J,K)
     IF(UVM2.GT.0.0)THEN
!dwen(20090806) **************
        hm(i,j,k)=u(i,j,k)*u(i,j,k)*hmix/uvm2
!dwen    *********************
        H(I,J,K)=U(I,J,K)*U(I,J,K)*UVP2/UVM2
        E(I,J,K)=V(I,J,K)*V(I,J,K)*UVP2/UVM2
     ELSE
!dwen(20090806) *************
        hm(i,j,k)=0.0
!dwen ***********************
        H(I,J,K)=0.0
        E(I,J,K)=0.0
     END IF

  END DO
  END DO
  END DO

END SUBROUTINE stbhor
