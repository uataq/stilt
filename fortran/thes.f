!#######################################################################
! SUBPROGRAM:  DTHES_T
!

!
!  RETURNS DIFFERENCE OF PARCELS SATURATION EQUIVALENT
!  POTENTIAL TEMPERATURE (K) AT TEMPERATURE TP (K)
!  FROM TRUE SATURATION EQUIVALENT POTENTIAL TEMPERATURE (K)
!  ALSO THE DERIVATIVE OF THIS DIFFERENCE WITH RESPECT TO TEMPERATURE (K)
!  input variables:
!    PP:    PRESSURE (mbar),
!    TP:    TEMPERATURE (K),
!    TTEST: TRUE SATURATION EQUIVALENT POTENTIAL TEMPERATURE (K)
!
!  output variables:
!    TTES:  DIFFERENCE IN THETAES
!    DTTES: DERIVATIVE OF THETAES DIFFERENCE
!
!  FORMULAS AFTER STULL
!  by CHG (12/3/01)
!
! $Id: thes.f,v 1.1 2009/10/26 15:36:55 jel Exp $
!
      SUBROUTINE THES_T(PP,TP,TTEST,TTES,DTTES)
      IMPLICIT REAL (A-H,O-Z)

!     get saturation vap. pres. (mbar) from T (K)
      IF(TP.GE.273.16)THEN
!     liquid
         ES=6.1078*exp(17.29694*(TP-273.16)/(TP-35.86))
         DES=ES*17.29694*(1/(TP-35.86)-(TP-273.16)/((TP-35.86)**2))
      ELSE
!     ice
         ES=exp(23.33086-6111.72784/TP+0.15215*log(TP))
         DES=ES*(6111.72784/(TP*TP) + 0.15215/TP)
      END IF
!     get sat. spec. hum. rsat (g/g) from es (mbar) and p (mbar)
         RS=0.622*ES/(PP-ES)
         DRS=0.622*DES/(PP-ES)+DES*RS/(PP-RS)

      TTES=TP*(1000.0/PP)**0.286+RS*2500.0*(1000.0/PP)**0.286
      TTES=TTES-TTEST
      DTTES=(1000.0/PP)**0.286+DRS*2500.0*(1000.0/PP)**0.286
      RETURN
      END
