	SUBROUTINE CPOLXY (STCPRM, X,Y, ENX,ENY,ENZ)
!*  WRITTEN ON 11/26/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      PARAMETER (REARTH=6371.2)
      PARAMETER (PI=3.14159265358979,RADPDG=PI/180,DGPRAD=180/PI)
      REAL STCPRM(9)
      DOUBLE PRECISION XPOL,YPOL,TEMP,XI0,ETA0,XI,ETA,RADIAL
      DOUBLE PRECISION TEMP2,YMERC,ARG,OARG,CLAT
      XI0 = ( X - STCPRM(3) ) * STCPRM(7) / REARTH
      ETA0 = ( Y - STCPRM(4) ) * STCPRM(7) /REARTH
      XI = XI0 * STCPRM(5) - ETA0 * STCPRM(6)
      ETA = ETA0 * STCPRM(5) + XI0 * STCPRM(6)
      RADIAL = 2. * ETA -  STCPRM(1) * (XI*XI + ETA*ETA)
      TEMP = STCPRM(1) * RADIAL
      IF (TEMP .GE. 1.) THEN
        ENX = 0.
        ENY = 0.
        ENZ = SIGN(1.,STCPRM(1))
        RETURN
      ENDIF
      IF (ABS(TEMP).LT.1.E-2) THEN
        TEMP2 = (TEMP / (2. - TEMP))**2
        YMERC = RADIAL / (2. - TEMP) * (1. + TEMP2 *                           &
                                     (1./3. + TEMP2 *                          &
                                     (1./5. + TEMP2 *                          &
                                     (1./7.))))
      ELSE
        YMERC = -.5 * LOG(1. - TEMP) / STCPRM(1)
      ENDIF
      ARG = EXP( YMERC )
      OARG = 1./ARG
      CLAT = 2./(ARG + OARG)
      ENZ = (ARG - OARG) * CLAT /2.
      TEMP = CLAT / SQRT(1. - TEMP)
      XPOL = - XI * STCPRM(1) * TEMP
      YPOL = (1. - ETA * STCPRM(1) ) * TEMP
      ENX = XPOL * STCPRM(5) + YPOL * STCPRM(6)
      ENY = YPOL * STCPRM(5) - XPOL * STCPRM(6)
      RETURN
      END
