!###############################################################################
! GEMTST - Routine to dump GEM variables at specific lat-lon locations
!-------------------------------------------------------------------------------
! LAST REVISED: 29 May 2008 (RRD) - initial version
!               09 Sep 2008 (RRD) - half grid point correction
!-------------------------------------------------------------------------------

SUBROUTINE gemtst (meto)

  USE gemcfg
  USE gemvar  

  IMPLICIT NONE

  INCLUDE 'DEFMETO.INC'
  TYPE(aset),INTENT(IN) :: meto  ! surface and advection variables

  REAL*4                :: clat1,clon1,dlat,dlon
  INTEGER*4             :: i,j,nx,ny,nz,np,kgrd,kret

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  IF(KINIT.LT.0)RETURN

  j=1+(meto%plat-clat1+dlat/2.0)/dlat
  i=1+(meto%plon-clon1+dlon/2.0)/dlon
  IF(i.LT. 1) i=nx+i
  IF(i.GT.nx) i=i-nx

  WRITE(*,100)'Variable',   'HYSPLIT',        'GEM-VALUE'
  WRITE(*,200)'Latitude',    METO%PLAT,        j
  WRITE(*,200)'Longitude',   METO%PLON,        i
  WRITE(*,400)'Land Use',    METO%LAND,        LUS(i,j)
  WRITE(*,300)'Solar Flux',  METO%DSWF,        SWF(i,j)
  WRITE(*,300)'Precip Rate', METO%RAIN*1000.0, P6H(i,j)*1000.0
  WRITE(*,300)'Roughness',   METO%AERO,        RGH(i,j)
  WRITE(*,300)'Friction V',  METO%USTR,        FRV(i,j)
  WRITE(*,300)'Stability',   METO%PSI,         PSI(i,j)
  WRITE(*,300)'Wind Speed',  METO%UBAR,        SQRT(U10(i,j)**2.0+V10(i,j)**2.0)
  WRITE(*,300)'Temperature', -1.0,             T02(i,j)
  READ(*,*)
  
  100 FORMAT(A20,2A10)
  200 FORMAT(A20,F10.4,I10)
  300 FORMAT(A20,2F10.4)
  400 FORMAT(A20,2I10)

END SUBROUTINE gemtst 
