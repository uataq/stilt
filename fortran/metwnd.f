!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METWND           METeorologically derived WiND
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICALLY DERIVED WIND DETERMINES WHICH VARIABLE IS TO BE
!   ANALZED FOR ITS SLOPE DEPENDING UPON THE VERTICAL VELOCITY OPTION
!   SELECTED ON INPUT - THESE OPTIONS ASSUME VARIOUS PROPERTIES ARE
!   CONSTANT ALONG THE TRAJECTORY.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 29 Sep 2000 (RRD) - fortran90 upgrade
!                 21 May 2002 (RRD) - added divergence option
!                 23 Jul 2002 (RRD) - added eta terrain correction
!                 14 Oct 2003 (RRD) - removed density from argument list
!                 02 Apr 2004 (RRD) - generic file unit numbers
!
! USAGE: CALL METWND(K1,K2,KVEL,NXS,NYS,NLVL,DM,ZMDL,ZSG,ZT,U,V,W,P,T,A)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            none
!   OUTPUT FILES:           none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE METWND(K1,K2,KVEL,NXS,NYS,NLVL,DM,ZMDL,ZSG,ZT,U,V,W,P,T,A)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  INTEGER,     INTENT(IN)    :: k1,k2         ! time index pointers
  INTEGER,     INTENT(IN)    :: kvel          ! vertical motion flag
  INTEGER,     INTENT(IN)    :: nxs,nys       ! dimensions of sub-grid
  INTEGER,     INTENT(IN)    :: nlvl          ! number of data levels
  REAL,        INTENT(IN)    :: dm            ! time between obs
  REAL,        INTENT(IN)    :: zmdl          ! vertical model domain top (m)
  REAL,        INTENT(IN)    :: zsg (:)       ! sigma levels
  REAL,        INTENT(IN)    :: zt(:,:)       ! terrain heights 
  REAL,        INTENT(IN)    :: u (:,:,:)     ! velocity components
  REAL,        INTENT(IN)    :: v (:,:,:)     ! velocity components
  REAL,        INTENT(INOUT) :: w (:,:,:)     ! velocity components
  REAL,        INTENT(IN)    :: p (:,:,:,:)   ! pressure             
  REAL,        INTENT(IN)    :: t (:,:,:,:)   ! potential temperature          
  REAL,        INTENT(IN)    :: a (:,:,:,:)   ! ambient temperature          

!-------------------------------------------------------------------------------

  REAL,  ALLOCATABLE :: d (:,:,:,:)     ! density              
  REAL,    PARAMETER :: rdry  = 287.04  ! dry air (J/Kg-K)
  REAL,    PARAMETER :: p2jm  = 100.0   ! mb to j/m3
  INTEGER            :: i,j,k,l         ! array indicies
  INTEGER            :: kret            ! error flag

!-------------------------------------------------------------------------------

  INTERFACE
  SUBROUTINE METSLP(K1,K2,NXS,NYS,NLVL,DM,ZSG,U,V,X,W)
  IMPLICIT NONE
  INTEGER,    INTENT(IN)    :: k1,k2            ! last and next time indicies
  INTEGER,    INTENT(IN)    :: nxs,nys          ! dimensions of sub-grid
  INTEGER,    INTENT(IN)    :: nlvl             ! number of data levels
  REAL,       INTENT(IN)    :: dm               ! time between obs
  REAL,       INTENT(IN)    :: zsg (:)          ! sigma levels
  REAL,       INTENT(IN)    :: u (:,:,:)        ! wind component
  REAL,       INTENT(IN)    :: v (:,:,:)        ! wind component
  REAL,       INTENT(IN)    :: x (:,:,:,:)      ! data mapping variable  
  REAL,       INTENT(OUT)   :: w (:,:,:)        ! vertical motion (ds/dt)
  END SUBROUTINE metslp

  SUBROUTINE METDIV(NXS,NYS,NLVL,ZMDL,ZSG,U,V,W)
  IMPLICIT NONE
  INTEGER,    INTENT(IN)    :: nxs,nys          ! dimensions of sub-grid
  INTEGER,    INTENT(IN)    :: nlvl             ! number of data levels
  REAL,       INTENT(IN)    :: zmdl             ! vertical model domain top (m)
  REAL,       INTENT(IN)    :: zsg (:)          ! sigma levels
  REAL,       INTENT(IN)    :: u (:,:,:)        ! wind component
  REAL,       INTENT(IN)    :: v (:,:,:)        ! wind component
  REAL,       INTENT(OUT)   :: w (:,:,:)        ! vertical motion (ds/dt)
  END SUBROUTINE metdiv

  SUBROUTINE METTER(NXS,NYS,NLVL,ZMDL,ZSG,ZT,U,V,W)
  IMPLICIT NONE
  INTEGER,    INTENT(IN)    :: nxs,nys          ! dimensions of sub-grid
  INTEGER,    INTENT(IN)    :: nlvl             ! number of data levels
  REAL,       INTENT(IN)    :: zmdl             ! vertical model domain top (m)
  REAL,       INTENT(IN)    :: zsg (:)          ! sigma levels
  REAL,       INTENT(IN)    :: zt(:,:)          ! terrain height
  REAL,       INTENT(IN)    :: u (:,:,:)        ! wind component
  REAL,       INTENT(IN)    :: v (:,:,:)        ! wind component
  REAL,       INTENT(INOUT) :: w (:,:,:)        ! vertical motion (ds/dt)
  END SUBROUTINE metter

  END INTERFACE

!-------------------------------------------------------------------------------

  SELECT CASE (KVEL)
  CASE (0)
     RETURN

!--------------------------------------------------
! constant pressure
  CASE (1)
     CALL METSLP(K1,K2,NXS,NYS,NLVL,DM,ZSG,U,V,P,W)

!--------------------------------------------------
! constant temperature (isentropic)
  CASE (2)
     CALL METSLP(K1,K2,NXS,NYS,NLVL,DM,ZSG,U,V,T,W)

!--------------------------------------------------
! constant density surfaces 
  CASE (3)
!    Previous versions density computed in prf??? routines and saved in array
!    However the constant density option is not frequently used and therefore
!    the density value is now computed from pressure and temperature when needed

!    allocate array
     I=SIZE(P,1)
     J=SIZE(P,2)
     K=SIZE(P,3)
     L=SIZE(P,4)
     ALLOCATE (D(I,J,K,L), STAT=KRET)
     IF(KRET.NE.0)THEN
        WRITE(KF21,*)'*ERROR* metwnd: allocation of density array - ',KRET
        STOP 900
     END IF

!    use pressure and temperature to find density 
     D=(P2JM*P)/(A*RDRY)

!    constant density
     CALL METSLP(K1,K2,NXS,NYS,NLVL,DM,ZSG,U,V,D,W)

!    free up space
     DEALLOCATE (D)

!--------------------------------------------------
! constant sigma level
  CASE (4)
     W = 0.0

!--------------------------------------------------
! divergence
  CASE (5)
     CALL METDIV(NXS,NYS,NLVL,ZMDL,ZSG,U,V,W)

!--------------------------------------------------
! terrain correction should only be called for ETA 
  CASE (6)
     CALL METTER(NXS,NYS,NLVL,ZMDL,ZSG,ZT,U,V,W)
  END SELECT 

END SUBROUTINE metwnd
