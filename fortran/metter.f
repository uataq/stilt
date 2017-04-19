!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METTER           METeorological TERrain factor
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:2002-07-23
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL TERRAIN CORRECTS THE VERTICAL VELCOITY FIELD
!   TO ACCOUNT FOR THE SLOPE OF THE TERRAIN.  THE CORRECTION IS  
!   REQUIRED FOR INPUT DATA WHERE THE VERTICAL VELOCITIES WERE
!   COMPUTED ON HORIZONTAL SURFACES SO THAT THEY REFLECT THE
!   VERTICAL MOTION OF AIR MOVING OVER CHANGING TERRAIN ELEVATION.
!   THE HYSPLIT INTERNAL COORIDINATE IS TERRAIN FOLLOWING, HENCE
!   THIS COMPONENT MUST BE REMOVED FROM THE FIELD, LEAVING ONLY
!   THE SYNOPTIC CONTRIBUTION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 23 Jul 2002 (RRD) - initial version
!
! USAGE:  CALL METTER(NXS,NYS,NLVL,ZMDL,ZSG,ZT,U,V,W)
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

SUBROUTINE METTER(NXS,NYS,NLVL,ZMDL,ZSG,ZT,U,V,W)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,    INTENT(IN)    :: nxs,nys          ! dimensions of sub-grid
  INTEGER,    INTENT(IN)    :: nlvl             ! number of data levels
  REAL,       INTENT(IN)    :: zmdl             ! vertical model domain top (m)
  REAL,       INTENT(IN)    :: zsg (:)          ! sigma levels
  REAL,       INTENT(IN)    :: zt(:,:)          ! terrain height 
  REAL,       INTENT(IN)    :: u (:,:,:)        ! wind component
  REAL,       INTENT(IN)    :: v (:,:,:)        ! wind component
  REAL,       INTENT(INOUT) :: w (:,:,:)        ! vertical motion (ds/dt)

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER                   :: i,j,k
  REAL                      :: delu,delv,wvel

!-------------------------------------------------------------------------------

  DO J=2,NYS-1
  DO I=2,NXS-1

!    velocity contributon (min^-1) due to slope of terrain surface
     DELU=0.5*(ZT(I+1,J)-ZT(I-1,J))*U(I,J,1)
     DELV=0.5*(ZT(I,J+1)-ZT(I,J-1))*V(I,J,1)

!    [ZMDL-ZT] converts from 1/min to sigma/min with implicit sign reversal.
!    Positive wvel is upward motion while positive sigma is downward motion.
     WVEL=(DELU+DELV)/(ZMDL-ZT(I,J))

     DO K=1,NLVL

!       This downward motion is added to the w field to remove the unwanted  
!       upward component.  The correction is scaled by sigma such that a full
!       correction is applied at the surface and no correction at the model top.
        W(I,J,K)=W(I,J,K)+WVEL*ZSG(K)

     END DO

  END DO
  END DO

! set boundary values at half adjacent    

  DO I=2,NXS-1
     W(I,1,:)  =0.5*W(I,2,:)
     W(I,NYS,:)=0.5*W(I,NYS-1,:)
  END DO

  DO J=1,NYS
     W(1,J,:)  =0.5*W(2,J,:)
     W(NXS,J,:)=0.5*W(NXS-1,J,:)
  END DO

END SUBROUTINE metter
