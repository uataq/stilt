!-------------------------------------------------------------------------------
! Configure subroutine interface argument lists for global routines       
! Last Revised: 30 Aug 2001 (RRD) - simultaneous multiple meteorology
!-------------------------------------------------------------------------------

INTERFACE

!-------------------------------------------------------------------------------
  SUBROUTINE GBL2XY(KG,KT,CLAT,CLON,X,Y)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: clat,clon   ! latlon location       
  REAL,    INTENT(OUT) :: x,y         ! grid position         
  END SUBROUTINE GBL2XY
!-------------------------------------------------------------------------------
  SUBROUTINE GBL2LL(KG,KT,X,Y,CLAT,CLON)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: x,y         ! grid position         
  REAL,    INTENT(OUT) :: clat,clon   ! latlon location       
  END SUBROUTINE GBL2LL
!-------------------------------------------------------------------------------

END INTERFACE
