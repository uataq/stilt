!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METSLP           METeorological SLoPe of a surface
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL SLOPE COMPUTES THE SLOPE OF METEOROLOGICAL VARIABLE
!   AND THEIR LOCAL TIME DERIVATIVE AND DETERMINES THE VERTICAL VELOCITY
!   REQUIRED TO STAY ON THE SURFACE FOR THAT VARIABLE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 29 Sep 2000 (RRD) - fortran90 upgrade
!                 14 Aug 2001 (RRD) - boundary values set
!                 09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL METSLP(K1,K2,NXS,NYS,NLVL,DM,ZSG,U,V,X,W)
!
!   INPUT ARGUMENT LIST:      see below
!   OUTPUT ARGUMENT LIST:     see below
!   INPUT FILES:              none
!   OUTPUT FILES:             none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE METSLP(K1,K2,NXS,NYS,NLVL,DM,ZSG,U,V,X,W)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,    INTENT(IN)    :: k1,k2            ! last and next time indicies
  INTEGER,    INTENT(IN)    :: nxs,nys          ! dimensions of sub-grid
  INTEGER,    INTENT(IN)    :: nlvl             ! number of data levels
  REAL,       INTENT(IN)    :: dm               ! time between obs
  REAL,       INTENT(IN)    :: zsg (:)          ! sigma levels
  REAL,       INTENT(IN)    :: u (:,:,:)        ! wind component
  REAL,       INTENT(IN)    :: v (:,:,:)        ! wind component
  REAL,       INTENT(IN)    :: x (:,:,:,:)      ! data mapping variable  
  REAL,       INTENT(OUT)   :: w (:,:,:)        ! vertical motion (ds/dt)

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER                   :: i,j,k
  REAL                      :: dmt,delx,delz,vert,dxd,dyd,time,wvel

!-------------------------------------------------------------------------------
  DMT = 1.0 / DM
!-------------------------------------------------------------------------------

  DO K=1,NLVL
  DO J=2,NYS-1
  DO I=2,NXS-1

!    compute horizontal gradients according to level
     IF(K.EQ.1)THEN                
        DELX=X(I,J,K+1,K2)-X(I,J,K,K2)
        DELZ=ZSG(K+1)-ZSG(K)
     ELSEIF(K.EQ.NLVL)THEN
        DELX=X(I,J,K,K2)-X(I,J,K-1,K2)
        DELZ=ZSG(K)-ZSG(K-1)
     ELSE
        DELX=X(I,J,K+1,K2)-X(I,J,K-1,K2)
        DELZ=ZSG(K+1)-ZSG(K-1)
     END IF 

!    vertical slope of surface (per sigma)
     VERT=DELX/DELZ

     IF(VERT.NE.0.0)THEN

!       horizontal slope of surface (grid/min)
        DXD=0.5*U(I,J,K)*(X(I+1,J,K,K2)-X(I-1,J,K,K2))
        DYD=0.5*V(I,J,K)*(X(I,J+1,K,K2)-X(I,J-1,K,K2))

!       local derivative (per min)
        TIME=(X(I,J,K,K2)-X(I,J,K,K1))*DMT

!       w velocity in sigma/min to maintain surface
        WVEL=-(TIME+DXD+DYD)/VERT
!       restrict maximum to 0.005 sig/min (about 1 m/s)
        W(I,J,K) = SIGN( MIN(0.005,ABS(WVEL)) , WVEL)

     END IF

  END DO
  END DO

! set boundary values    

  DO I=2,NXS-1
     W(I,1,K)  =W(I,2,K)
     W(I,NYS,K)=W(I,NYS-1,K)
  END DO

  DO J=1,NYS
     W(1,J,K)  =W(2,J,K)
     W(NXS,J,K)=W(NXS-1,J,K)
  END DO

  END DO

END SUBROUTINE metslp
