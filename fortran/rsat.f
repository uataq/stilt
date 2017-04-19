!#######################################################################
      REAL FUNCTION RSAT_T(PP,TP)
!
!  RETURNS SATURATION SPECIFIC HUMIDITY (g/g) FROM
!  PRESSURE (mbar) AND TEMPERATURE (K)
!  AFTER STULL
!
!  $Id: rsat.f,v 1.1 2009/10/26 15:36:55 jel Exp $
!
      IMPLICIT REAL (A-H,O-Z)
!
! CHG:(12/3/01)
      REAL,intent(in) :: PP,TP

!      DATA IFF /0/

!     get saturation vap. pres. (mbar) from T (K)
      IF(TP.GE.273.15)THEN
!     liquid
         ES=6.1078*exp(17.29694*(TP-273.16)/(TP-35.86))
      ELSE
!     ice
         ES=exp(23.33086-6111.72784/TP+0.15215*log(TP))
      END IF
!     get sat. spec. hum. rsat (g/g) from es (mbar) and p (mbar)
      RSAT_T=0.622*ES/(PP-ES)
      RETURN
      END
