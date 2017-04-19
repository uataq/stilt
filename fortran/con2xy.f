!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CON2XY           REMAP LAT LON VALUES TO POSITION ON MAP
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            MAPS THE INPUT CONCENTRATION ARRAY FROM A LAT/LON GRID TO
!            AND X/Y CONFORMAL MAP PROJECTION USING NEAREST NEIGHBOR
!            CONCEPT.  LINEAR INTERPOLATION IS PERFORMED AT THE PLUME
!            EDGES BETWEEN NON-ZERO VALUES. ADDITIONAL SMOOTHING OPTIONAL.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 (RRD)
!                 15 Mar 2000 (RRD) - improved high latitude interpolation
!                 03 Nov 2000 (RRD) - longitude correction
!                 21 Nov 2000 (RRD) - fortran90 upgrade
!                 12 Dec 2000 (RRD) - optional smoothing factor
!                 21 Feb 2001 (RRD) - grid span correction
!                 05 May 2004 (RRD) - correction for global grid cut
!                 12 Jan 2007 (RRD) - cylindrical equidistant option
!                 24 Jun 2008 (RRD) - fixed error in distance weighting
!                 05 Aug 2008 (RRD) - concentration smoothing argument
!
! USAGE: CALL CON2XY(KPROJ,NXP,NYP,XCON,TFACT,SCAN,CONC,PARMAP,NLAT,NLON, 
!                    DLAT,DLON,CLAT,CLON)
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

SUBROUTINE CON2XY(KPROJ,NXP,NYP,XCON,TFACT,SCAN,CONC,PARMAP,NLAT,NLON,         &
                  DLAT,DLON,CLAT,CLON)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER, INTENT(IN)   :: kproj       ! map projection         
  INTEGER, INTENT(IN)   :: nxp,nyp     ! output grid dimensions
  REAL,    INTENT(IN)   :: tfact       ! concentration units multiplier
  INTEGER, INTENT(IN)   :: scan        ! map scan index for smoothing       
  REAL,    INTENT(IN)   :: conc (:,:)  ! concentration grid on lat/lon system
  REAL,    INTENT(IN)   :: parmap (9)  ! concentration grid on lat/lon system
  INTEGER, INTENT(IN)   :: nlat,nlon   ! size of input concentration grid
  REAL,    INTENT(IN)   :: dlat,dlon   ! grid resolution in lat/lon 
  REAL,    INTENT(IN)   :: clat,clon   ! center of input concentration grid
  REAL,    INTENT(OUT)  :: xcon (:,:)  ! concentration array on conformal grid 

!-------------------------------------------------------------------------------
! internally defined variables
!-------------------------------------------------------------------------------

  REAL    :: yj,xi,rd,pcon,plat,plon,wght,csum,cmax,xdis,ydis 
  INTEGER :: i,j,ii,jj,ic,jc 
  INTEGER :: maxdim (2)

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE CONSMT( RVAL,SCAN )
  IMPLICIT NONE
  REAL,    INTENT(INOUT) :: rval (:,:)      ! array to be smoothed
  INTEGER, INTENT(IN)    :: scan            ! smoothing factor
  END SUBROUTINE consmt
  END INTERFACE
!-------------------------------------------------------------------------------

  CMAX=0.0

! fill output array
  DO II=1,NXP
  DO JJ=1,NYP

!    find lat/lon of conformal grid point
     IF(KPROJ.EQ.4)THEN
        CALL CYL2LL(FLOAT(II),FLOAT(JJ),PLAT,PLON)
     ELSE
        CALL CXY2LL(PARMAP,FLOAT(II),FLOAT(JJ),PLAT,PLON)
     END IF

!    longitude correction
     IF((CLON.GE.0.0.AND.PLON.LT.0.0).OR.(PLON-CLON).LE.0.0)PLON=PLON+360.0

!    convert to grid point number
     XI=1.0+(PLON-CLON)/DLON
     YJ=1.0+(PLAT-CLAT)/DLAT

!    compute index of lower left point
     IC=INT(XI)
     JC=INT(YJ)

!    point is within the lat-lon concentration grid
     IF(IC.GE.1.AND.IC.LT.NLON.AND.JC.GE.1.AND.JC.LT.NLAT)THEN
        CSUM=0.0
        WGHT=0.0
        DO I=0,1
        DO J=0,1
!          compute contribution from surrounding concentrations
           PCON=CONC(IC+I,JC+J)
           CMAX=MAX(CMAX,PCON)  
           XDIS=(XI-(IC+I))*(XI-(IC+I))
           YDIS=(YJ-(JC+J))*(YJ-(JC+J))
           IF(XDIS+YDIS.GT.0.0)THEN
              RD=1.0/(XDIS+YDIS) 
           ELSE
              RD=1.0E+03
           END IF
           CSUM=CSUM+PCON*RD
           WGHT=WGHT+RD
        END DO
        END DO

!    point is along the eastern boundary
     ELSEIF(IC.EQ.NLON.AND.JC.GE.1.AND.JC.LT.NLAT)THEN
        CSUM=0.0
        WGHT=0.0
        DO J=0,1
!          compute value of surrounding concentrations
           PCON=CONC(IC,JC+J)
           YDIS=(YJ-(JC+J))*(YJ-(JC+J))
           IF(YDIS.GT.0.0)THEN
              RD=1.0/YDIS
           ELSE
              RD=1.0E+03
           END IF
           CSUM=CSUM+PCON*RD
           WGHT=WGHT+RD
        END DO

!    point is along the northern boundary
     ELSEIF(JC.EQ.NLAT.AND.IC.GE.1.AND.IC.LT.NLON)THEN
        CSUM=0.0
        WGHT=0.0
        DO I=0,1
!          compute value of surrounding concentrations
           PCON=CONC(IC+I,JC)
           XDIS=(XI-(IC+I))*(XI-(IC+I))
           IF(XDIS.GT.0.0)THEN
              RD=1.0/XDIS
           ELSE
              RD=1.0E+03
           END IF
           CSUM=CSUM+PCON*RD
           WGHT=WGHT+RD
        END DO

     ELSEIF(IC.EQ.NLON.AND.JC.EQ.NLAT)THEN
        WGHT=1.0
        CSUM=CONC(IC,JC)

     ELSE
!       point outside of the concentration grid
        WGHT=0.0

     END IF

     IF(WGHT.GT.0.0)THEN
        XCON(II,JJ)=TFACT*CSUM/WGHT
     ELSE
        XCON(II,JJ)=0.0
     END IF

  END DO
  END DO

! optional smoothing (scan=1 gives 1,2,1 weights)   
  IF( scan.GT.0 ) CALL CONSMT(xcon,scan)

! replace maximum value with orignal data maximum
  MAXDIM=MAXLOC(XCON)
  IF(CMAX.GT.0.0) XCON(MAXDIM(1),MAXDIM(2))=CMAX

END SUBROUTINE con2xy
