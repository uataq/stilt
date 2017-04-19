!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METDIV           METeorological DIVergence
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:02-05-21
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL DIVERGENCE COMPUTES THE VERTICAL VELOCITY FROM THE 
!   VERTICAL INTEGRATION OF THE HORIZONTAL VELOCOTY DIVERGENCE. 
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 21 May 2002 (RRD) - initial version
!                 24 Jun 2002 (RRD) - modified dS/dZ
!                 06 Feb 2003 (RRD) - correction to integration
!
! USAGE:  CALL METDIV(NXS,NYS,NLVL,ZMDL,ZSG,U,V,W)
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

SUBROUTINE METDIV(NXS,NYS,NLVL,ZMDL,ZSG,U,V,W)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,    INTENT(IN)    :: nxs,nys          ! dimensions of sub-grid
  INTEGER,    INTENT(IN)    :: nlvl             ! number of data levels
  REAL,       INTENT(IN)    :: zmdl             ! vertical model domain top (m)
  REAL,       INTENT(IN)    :: zsg (:)          ! sigma levels
  REAL,       INTENT(IN)    :: u (:,:,:)        ! wind component
  REAL,       INTENT(IN)    :: v (:,:,:)        ! wind component
  REAL,       INTENT(OUT)   :: w (:,:,:)        ! vertical motion (ds/dt)

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER                   :: i,j,k
  REAL                      :: divg,dels,delu,delv,wvel

!-------------------------------------------------------------------------------

  W=0.0 

  DO J=2,NYS-1
  DO I=2,NXS-1

!    initialize vertical velocity at ground level
     WVEL=0.0

     DO K=1,NLVL

!       vertical integration distance in sigma units
        IF(K.EQ.1)THEN
           DELS=ZSG(K)-1.0 
        ELSEIF(K.EQ.NLVL)THEN
           DELS=ZSG(K)-ZSG(K-1)
        ELSE
           DELS=(ZSG(K+1)-ZSG(K-1))/2.0
        END IF

!       divergence in min^-1
        DELU=0.5*(U(I+1,J,K)-U(I-1,J,K))
        DELV=0.5*(V(I,J+1,K)-V(I,J-1,K))
        DIVG=DELU+DELV

!       integrate column for w velocity
        WVEL=WVEL+DIVG

!       conversion from 1/min to sigma/min
!       positive divergence implies downward motion
        W(I,J,K)=-WVEL*DELS

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

END SUBROUTINE metdiv
