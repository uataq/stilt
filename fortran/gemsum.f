!###############################################################################
! GEMSUM - Select pollutant target and call equation integration subroutine 
!-------------------------------------------------------------------------------
! LAST REVISED: 28 May 2008 (RRD) - initial version
!               02 Feb 2009 (RRD) - resturcture loops
!-------------------------------------------------------------------------------

SUBROUTINE gemsum (DIRT,DT,BACK,CDEP,IMO)

  USE gemcfg
  USE gemvar

  IMPLICIT NONE

  INCLUDE     'DEFCONC.INC'           ! concentration and pollutant structure
  TYPE(pset),  INTENT(IN) :: DIRT (:) ! for each pollutant type

  INTEGER*4,   INTENT(IN) :: DT       ! external time step
  LOGICAL,     INTENT(IN) :: BACK     ! integration direction
  LOGICAL,     INTENT(IN) :: CDEP     ! deposition flag      
  INTEGER*4,   INTENT(IN) :: IMO      ! current month

  REAL*8    :: qsum
  INTEGER*4 :: k,n,nx,ny,nz,np,kgrd,loop

  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

!-------------------------------------
  INTERFACE
    SUBROUTINE gemdep (DIRT,N,IMO)
    IMPLICIT NONE
    INCLUDE     'DEFCONC.INC'           
    TYPE(pset),  INTENT(IN) :: DIRT (:) 
    INTEGER*4,   INTENT(IN) :: N   
    INTEGER*4,   INTENT(IN) :: IMO 
    END SUBROUTINE gemdep
  END INTERFACE
!------------------------------------

! align internal equation timestep (delta) with external
! time step (dt) for proper integration time step  

  IF(dt.LT.delta)THEN
     delta=dt
     loop=1
  ELSEIF(dt.EQ.delta)THEN
     loop=1
  ELSE
     DO WHILE (MOD(dt,delta).NE.0.AND.delta.GT.1)
        delta=delta-1
     END DO
     loop=dt/delta
  END IF
  IF(BACK)delta=-delta

! perform integration for each pollutant species
! and keep track of total system mass for diagnostics

  qsum = 0.0
  DO n=1,np
     QQQ => XXX(:,:,:,n)

     DO k=1,loop
        CALL gemeqn 

!       call any deposition or chemistry routines after the
!       advection and diffusion because they do not interact with
!       adjacent grid cells hence can be applied to master array
        IF(CDEP) CALL gemdep(dirt,n,imo)

     END DO
     qsum = qsum + qtot
  END DO
  NULLIFY (QQQ)
  qtot = qsum

END SUBROUTINE gemsum
