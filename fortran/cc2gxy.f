!	SUBROUTINE CC2GXY (STCPRM, X,Y, UE,VN, UG,VG,proj)
	SUBROUTINE CC2GXY (STCPRM, X,Y, UE,VN, UG,VG)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL

!      use map_utils
      use module_defgrid, only : vmiss, vmissle
      
!      type (proj_info), intent(in) :: proj

      PARAMETER (REARTH=6371.2)
      PARAMETER (PI=3.14159265358979,RADPDG=PI/180,DGPRAD=180/PI)
      REAL STCPRM(9)
      DOUBLE PRECISION XPOLG,YPOLG,TEMP,XI0,ETA0

      if (all(stcprm .le. vmissle)) then
!     TBD: add calls to map_utils routines for wind rotation
         write (*,*) 'unsupported cc2gxy for wps projection'
         stop 'unsupported cc2gxy for wps projection'
      end if

      XI0 = ( X - STCPRM(3) ) * STCPRM(7) / REARTH
      ETA0 = ( Y - STCPRM(4) ) * STCPRM(7) /REARTH
      XPOLG = STCPRM(6) - STCPRM(1) * XI0
      YPOLG = STCPRM(5) - STCPRM(1) * ETA0
      TEMP = SQRT ( XPOLG ** 2 + YPOLG ** 2 )
!* Revised 2/12/02 to allow cartographic wind vector transformations everywhere
!* except at the poles, with WMO conventions only at the poles.
      IF (TEMP .LE. 1.0e-3) THEN
!* Close to either pole, convert to XLAT,XLON for compatibility
!* with CXY2LL usasge
          CALL CXY2LL(STCPRM, X,Y, XLAT,XLONG)
          CALL CC2GLL(STCPRM, XLAT,XLONG, UE,VN, UG,VG)
        ELSE
!* Elsewhere use vector algebra and avoid time consuming trigonometry
          XPOLG = XPOLG / TEMP
          YPOLG = YPOLG / TEMP
!*  CHANGE MADE 3/9/99 TO ALLOW UG,VG TO HAVE SAME STORAGE AS UE,VN
          TEMP = YPOLG * UE + XPOLG * VN
          VG = YPOLG * VN - XPOLG * UE
          UG = TEMP
!* PERMITTING ROTATE OF WINDS IN PLACE
        ENDIF
      RETURN
      END
