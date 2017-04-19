      REAL FUNCTION CGSZXY (STRCMP, X,Y)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      PARAMETER (REARTH=6371.2,ALMST1=.99999)
      REAL STRCMP(9)
      DOUBLE PRECISION YMERC,EFACT
      DOUBLE PRECISION XI0,ETA0,XI,ETA
      XI0 = ( X - STRCMP(3) ) * STRCMP(7) / REARTH
      ETA0 = ( Y - STRCMP(4) ) * STRCMP(7) /REARTH
      XI = XI0 * STRCMP(5) - ETA0 * STRCMP(6)
      ETA = ETA0 * STRCMP(5) + XI0 * STRCMP(6)
      RADIAL = 2. * ETA - STRCMP(1) * (XI*XI + ETA*ETA)
      EFACT = STRCMP(1) * RADIAL
      IF (EFACT .GT. ALMST1) THEN
        IF (STRCMP(1).GT.ALMST1) THEN
          CGSZXY = 2. * STRCMP(7)
        ELSE
          CGSZXY = 0.
        ENDIF
        RETURN
      ENDIF
      IF (ABS(EFACT) .LT. 1.E-2) THEN
        TEMP = (EFACT / (2. - EFACT) )**2
        YMERC = RADIAL / (2. - EFACT) * (1.    + TEMP *                        &
                                        (1./3. + TEMP *                        &
                                        (1./5. + TEMP *                        &
                                        (1./7. ))))
      ELSE
        YMERC = - LOG( 1. - EFACT ) /2. /STRCMP(1)
      ENDIF
      IF (YMERC .GT. 6.) THEN
        IF (STRCMP(1) .GT. ALMST1) THEN
          CGSZXY = 2. * STRCMP(7)
        ELSE
          CGSZXY = 0.
        ENDIF
      ELSE IF (YMERC .LT. -6.) THEN
        IF (STRCMP(1) .LT. -ALMST1) THEN
          CGSZXY = 2. * STRCMP(7)
        ELSE
          CGSZXY = 0.
        ENDIF
      ELSE
        EFACT = EXP(YMERC)
        CGSZXY = 2. * STRCMP(7) * EXP (STRCMP(1) * YMERC)                      &
                               / (EFACT + 1./EFACT)
      ENDIF
      RETURN
      END
