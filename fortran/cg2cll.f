	SUBROUTINE CG2CLL (STCPRM, XLAT,XLONG, UG,VG, UE,VN)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      PARAMETER (PI=3.14159265358979,RADPDG=PI/180,DGPRAD=180/PI)
      DOUBLE PRECISION XPOLG,YPOLG,ALONG,SLONG,CLONG,ROT
      REAL STCPRM(9)
        ALONG = CSPANF( XLONG - STCPRM(2), -180., 180.)
        ROT = - STCPRM(1) * ALONG
!* Revised 2/12/02 to allow cartographic wind vector transformations everywhere
!* with rotation to nominal longitudes at the poles, to match U,V values on a
!* Lat-Lon grid.
        SLONG = SIN( RADPDG * ROT )
        CLONG = COS( RADPDG * ROT )
        XPOLG = SLONG * STCPRM(5) + CLONG * STCPRM(6)
        YPOLG = CLONG * STCPRM(5) - SLONG * STCPRM(6)
!*  CHANGE MADE 3/9/99 TO ALLOW UG,VG TO HAVE SAME STORAGE AS UE,VN
        TEMP = YPOLG * UG - XPOLG * VG
        VN = YPOLG * VG + XPOLG * UG
        UE = TEMP
!* PERMITTING ROTATE OF WINDS IN PLACE
        RETURN
      END
