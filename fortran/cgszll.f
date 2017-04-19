
      REAL FUNCTION CGSZLL (STRCMP, XLAT,XLONG)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      PARAMETER (PI=3.14159265358979,RADPDG=PI/180,DGPRAD=180D0/PI)
      REAL STRCMP(9)
      DOUBLE PRECISION SLAT,YMERC,EFACT
      IF (XLAT .GT. 89.995) THEN
!* CLOSE TO NORTH POLE
        IF (STRCMP(1) .GT. 0.9999) THEN
!* AND TO GAMMA == 1.
          CGSZLL = 2. * STRCMP(7)
          RETURN
        ENDIF
        EFACT = COS(RADPDG * XLAT)
        IF (EFACT .LE. 0.) THEN
          CGSZLL = 0.
          RETURN
        ELSE
          YMERC = - LOG( EFACT /(1. + SIN(RADPDG * XLAT)))
        ENDIF
      ELSE IF (XLAT .LT. -89.995) THEN
!* CLOSE TO SOUTH POLE
        IF (STRCMP(1) .LT. -0.9999) THEN
!* AND TO GAMMA == -1.0
          CGSZLL = 2. * STRCMP(7)
          RETURN
        ENDIF
        EFACT = COS(RADPDG * XLAT)
        IF (EFACT .LE. 0.) THEN
          CGSZLL = 0.
          RETURN
        ELSE
          YMERC = LOG( EFACT /(1. - SIN(RADPDG * XLAT)))
        ENDIF
      ELSE
      SLAT = SIN(RADPDG * XLAT)
      YMERC = LOG((1. + SLAT) / (1. - SLAT))/2.
!      EFACT = EXP(YMERC)
!      CGSZLL = 2. * STRCMP(7) * EXP (STRCMP(1) * YMERC)
!     C                   / (EFACT + 1./EFACT)
      ENDIF
      CGSZLL = STRCMP(7) * COS(RADPDG * XLAT) * EXP(STRCMP(1) *YMERC)
      RETURN
      END
