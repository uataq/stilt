      SUBROUTINE CCRVXY (STRCMP, X,Y, GX,GY)
!*  WRITTEN ON 9/20/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      PARAMETER (REARTH=6371.2)
      PARAMETER (PI=3.14159265358979,RADPDG=PI/180,DGPRAD=180/PI)
      REAL STRCMP(9)
      DOUBLE PRECISION XPOLG,YPOLG,TEMP,YMERC,EFACT,CURV
      TEMP = STRCMP(1) * STRCMP(7) /REARTH
      XPOLG = STRCMP(6) + TEMP * (STRCMP(3) - X)
      YPOLG = STRCMP(5) + TEMP * (STRCMP(4) - Y)
      TEMP = SQRT ( XPOLG ** 2 + YPOLG ** 2 )
      IF (TEMP.GT.0.) THEN
        YMERC = - LOG( TEMP) /STRCMP(1)
        EFACT = EXP(YMERC)
         CURV = ( (STRCMP(1) - 1.D0) * EFACT +                                 &
                  (STRCMP(1) + 1.D0) / EFACT )                                 &
                 * .5D0 / REARTH
        GX = XPOLG * CURV / TEMP
        GY = YPOLG * CURV / TEMP
      ELSE
        IF (ABS(STRCMP(1)) .EQ. 1.) THEN
          GX = 0.
          GY = 0.
        ELSE
          GX = 1./REARTH
          GY = 1./REARTH
        ENDIF
      ENDIF
      RETURN
      END
