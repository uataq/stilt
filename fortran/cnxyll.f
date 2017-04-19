      SUBROUTINE CNXYLL (STCPRM, XI,ETA, XLAT,XLONG)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
!  MAIN TRANSFORMATION ROUTINE FROM CANONICAL (EQUATOR-CENTERED,
!  RADIAN UNIT) COORDINATES
      PARAMETER (PI=3.14159265358979,RADPDG=PI/180,DGPRAD=180/PI)
      REAL STCPRM(9)
      DOUBLE PRECISION GAMMA,TEMP,ARG1,ARG2,YMERC,ALONG,GXI,CGETA
      GAMMA = STCPRM(1)
      CGETA = 1.D0 - GAMMA * ETA
      GXI = GAMMA * XI
!  CALCULATE EQUIVALENT MERCATOR COORDINATE
      ARG2 = ETA + (ETA * CGETA - GXI * XI)
      ARG1 = GAMMA * ARG2
      IF (ARG1 .GE. 1.0) THEN
!  DISTANCE TO NORTH (OR SOUTH) POLE IS ZERO (OR IMAGINARY ;) )
        XLAT = SIGN(90.,STCPRM(1))
!        XLONG = STCPRM(2)
        XLONG = 90. + XLAT
! Change made 02/12/02 to acommodate WMO reporting conventions.  North
! pole is longitude 180., so "North" points to the Greenwich Meridian,
! South Pole is longitude 0. so "North" again points to the Greenwich
! Meridian.
        RETURN
      ENDIF
      IF (ABS(ARG1) .LT. .01) THEN
!  CODE FOR GAMMA SMALL OR ZERO.  THIS AVOIDS ROUND-OFF ERROR OR DIVIDE-
!  BY ZERO IN THE CASE OF MERCATOR OR NEAR-MERCATOR PROJECTIONS.
        TEMP = (ARG1 / (2. - ARG1) )**2
        YMERC = ARG2 / (2. - ARG1) * (1.    + TEMP *                           &
                                     (1./3. + TEMP *                           &
                                     (1./5. + TEMP *                           &
                                     (1./7. ))))
      ELSE
! CODE FOR MODERATE VALUES OF GAMMA
        YMERC = - LOG ( 1. - ARG1 ) /2. / GAMMA
      ENDIF
!  CONVERT YMERC TO LATITUDE
      TEMP = EXP( - ABS(YMERC) )
      XLAT = SIGN(ATAN2((1. - TEMP) * (1. + TEMP), 2. * TEMP), YMERC)
!  FIND LONGITUDES
      IF ( ABS(GXI) .LT. .01*CGETA ) THEN
!  CODE FOR GAMMA SMALL OR ZERO.  THIS AVOIDS ROUND-OFF ERROR OR DIVIDE-
!  BY ZERO IN THE CASE OF MERCATOR OR NEAR-MERCATOR PROJECTIONS.
        TEMP = ( GXI /CGETA )**2
        ALONG = XI / CGETA * (1.    - TEMP *                                   &
                               (1./3. - TEMP *                                 &
                               (1./5. - TEMP *                                 &
                               (1./7.   ))))
      ELSE
! CODE FOR MODERATE VALUES OF GAMMA
        ALONG = ATAN2( GXI, CGETA) / GAMMA
      ENDIF
      XLONG = SNGL(STCPRM(2) + DGPRAD * ALONG)
      XLAT = XLAT * DGPRAD
      RETURN
      END
