!$$$  SUBROUTINE DOCUMENTATION BLOCK
!
! Output subroutine that reads a 1-D array of 4-character strings identifying output variables
!        that will be written out to PARTICLE.DAT
!      02/28/2004 by Marcos Longo and John Lin
!
! $Id: output.f,v 1.1 2009/10/26 15:36:54 jel Exp $
!
!$$$
SUBROUTINE OUTPUT (IFILN,IVMAX,VARSIWANT,NTURB,IETIM,KP,ICNDX,XOUT, &
                   YOUT,ZOUT,ZSFC,SIGMAW,TEMP,SAMPTTCUM,FOOTCUM,    &
                   SWF,WWOUT,ZMLAVG,RAIN,CRAI,ZFXCUM,SHTF,WHTF,     &
                   TCLD,DMASSWT,TLGR,DENSLOCAL,RHFR,SPHU,SOLW,LCLD, &
                   ZLOCNEXT,PRESLOCAL,TEMPLOCAL)

   IMPLICIT NONE

! INTENT(IN) VARIABLES
   INTEGER        , INTENT(IN)  :: IFILN,IVMAX,NTURB,IETIM,KP,ICNDX
   REAL           , INTENT(IN)  :: XOUT,YOUT,ZOUT,ZSFC,SIGMAW,TEMP,       &
                                   SAMPTTCUM,FOOTCUM,SWF,WWOUT,ZMLAVG,    &
                                   RAIN,CRAI,ZFXCUM,SHTF,WHTF,TCLD,       &
                                   DMASSWT,TLGR,DENSLOCAL,RHFR,SPHU,SOLW, &
                                   LCLD,ZLOCNEXT,PRESLOCAL,TEMPLOCAL
   CHARACTER(LEN=*), INTENT(IN) :: VARSIWANT(IVMAX)

! INTERNAL VARIABLES
   REAL                         :: WRITING(IVMAX)
   CHARACTER(LEN=270)           :: FORMATO
   CHARACTER(LEN=12)            :: FORMATO2
   INTEGER                      :: IV,E2,EF


!---------------------------------------------------------------------------------------------------
   FORMATO=''

! JCL:(3/2/2004) this code doesn't work. . .
! CHECK IF THERE IS ANY NOT-ALLOWED VARIABLE WHEN NTURB = 1
!      IF (NTURB.EQ.1) THEN
!        WRITE(*,*) VARSIWANT
!        WHERE(VARSIWANT == 'sigw' .OR. VARSIWANT == 'SIGW' .OR.
!     &        VARSIWANT == 'tlgr' .OR. VARSIWANT == 'TLGR')
!          STOP 'You cannot output SIGW or TLGR when NTURB = 1!!!'
!        END WHERE
!      END IF


! FILLING THE OUTPUT VECTOR
   DO IV= 1, IVMAX
     SELECT CASE(VARSIWANT(IV))
     CASE ('time','TIME')  !time since start of simulation
       WRITE (FORMATO2,'(a) ') '1x,f10.0,'
       WRITING(IV) = IETIM
     CASE ('sigw','SIGW')  !standard deviation of vertical velocity [m/s]
       WRITE (FORMATO2,'(a) ') '1x,f10.3,'
       WRITING(IV) = SIGMAW
     CASE ('tlgr','TLGR')  !Lagrangian decorrelation timescale [s]
       WRITE (FORMATO2,'(a) ') '1x,f13.4,'
       WRITING(IV) = TLGR
     CASE ('long','LONG')  !particle longitude position [degrees]
       WRITE (FORMATO2,'(a) ') '1x,f10.4,'
       WRITING(IV) = XOUT
     CASE ('lati','LATI')  !particle latitude position [degrees]
       WRITE (FORMATO2,'(a) ') '1x,f10.4,'
       WRITING(IV) = YOUT
     CASE ('zagl','ZAGL')  !particle vertical position [m above-ground-level]
       WRITE (FORMATO2,'(a) ') '1x,f10.4,'
       WRITING(IV) = ZOUT
     CASE ('zsfc','ZSFC')  !terrain height [m]
       WRITE (FORMATO2,'(a) ') '1x,f10.4,'
       WRITING(IV) = ZSFC
     CASE ('indx','INDX')  !particle index
       WRITE (FORMATO2,'(a) ') '1x,f8.0,'
       WRITING(IV) = KP
     CASE ('icdx','ICDX')  !cloud index (1 = updraft, 2=environment, 3=downdraft)
       WRITE (FORMATO2,'(a) ') '3x,f2.0,'
       WRITING(IV) = ICNDX
     CASE ('temp','TEMP')  !air temperature at lowest model layer [K]
       WRITE (FORMATO2,'(a) ') '1x,f9.2,'
       WRITING(IV) = TEMP
     CASE ('temz','TEMZ')  !air temperature at particle location [K]
       WRITE (FORMATO2,'(a) ') '1x,f9.2,'
       WRITING(IV) = TEMPLOCAL
     CASE ('pres','PRES')  !air pressure at particle location [mbar]
       WRITE (FORMATO2,'(a) ') '1x,f9.2,'
       WRITING(IV) = PRESLOCAL
     CASE ('samt','SAMT')  !amount of time that particle spends below 'VEGHT' [min]
       WRITE (FORMATO2,'(a) ') '1x,f9.4,'
       WRITING(IV) = SAMPTTCUM
     CASE ('foot','FOOT')  !sensitivity of mixing ratio to surface fluxes [ppm/(micro-moles/m2/s)]
       WRITE (FORMATO2,'(a) ') '1x,f12.7,'
       WRITING(IV) = FOOTCUM
     CASE ('shtf','SHTF')  !sensible heat flux [W/m2]
       WRITE (FORMATO2,'(a) ') '1x,f9.2,'
       WRITING(IV) = SHTF
     CASE ('whtf','WHTF')  !latent heat flux [W/m2]
       WRITE (FORMATO2,'(a) ') '1x,f9.2,'
       WRITING(IV) = WHTF
     CASE ('tcld','TCLD')  !total cloud cover [%]
       WRITE (FORMATO2,'(a) ') '1x,f12.7,'
       WRITING(IV) = TCLD
     CASE ('dmas','DMAS')  !particle weight--changes due to mass violation [init value = 1.0]
       WRITE (FORMATO2,'(a) ') '1x,f16.3,'
       WRITING(IV) = DMASSWT
     CASE ('dens','DENS')  !air density [kg/m3]
       WRITE (FORMATO2,'(a) ') '1x,f8.5,'
       WRITING(IV) = DENSLOCAL
     CASE ('rhfr','RHFR')  !relative humidity fraction [0~1.0]
       WRITE (FORMATO2,'(a) ') '1x,f8.5,'
       WRITING(IV) = RHFR
     CASE ('sphu','SPHU')  !specific humidity [g/g]
       WRITE (FORMATO2,'(a) ') '1x,f13.10,'
       WRITING(IV) = SPHU
     CASE ('solw','SOLW')  !soil moisture [units ???]
       WRITE (FORMATO2,'(a) ') '1x,f12.5,'
       WRITING(IV) = SOLW
     CASE ('lcld','LCLD')  !low cloud cover [%]
       WRITE (FORMATO2,'(a) ') '1x,f12.7,'
       WRITING(IV) = LCLD
     CASE ('zloc','ZLOC')  !limit of convection heights [m]
       WRITE (FORMATO2,'(a) ') '1x,f9.2,'
       WRITING(IV) = ZLOCNEXT
     CASE ('dswf','DSWF')  !downward shortwave radiation [W/m2]
       WRITE (FORMATO2,'(a) ') '1x,f10.2,'
       WRITING(IV) = SWF
     CASE ('wout','WOUT')  !vertical mean wind [m/s]
       WRITE (FORMATO2,'(a) ') '1x,f12.6,'
       WRITING(IV) = WWOUT
     CASE ('mlht','MLHT')  !mixed-layer height [m]
       WRITE (FORMATO2,'(a) ') '1x,f12.2,'
       WRITING(IV) = ZMLAVG
     CASE ('rain','RAIN')  !total rain fall rate [m/min]
       WRITE (FORMATO2,'(a) ') '1x,g11.4,'
       WRITING(IV) = RAIN
     CASE ('crai','CRAI')  !convective rain fall rate [m/min]
       WRITE (FORMATO2,'(a) ') '1x,g11.4,'
       WRITING(IV) = CRAI
     CASE ('zcfx','ZFX1')  !vertical displacement due to
                           !convective flux (deep or shallow,
                           !up or downdraft) along trajectory [m]
       WRITE (FORMATO2,'(a) ') '1x,f9.2,'
       WRITING(IV) = ZFXCUM
     CASE DEFAULT
       FORMATO = 'Weird variable in VARSIWANT(SETUP.CFG):  '//VARSIWANT(IV)
       PRINT '(a) ', TRIM(FORMATO)
       STOP
     END SELECT
     EF = INDEX(FORMATO,' ',.FALSE.)-1
     E2 = INDEX(FORMATO2,' ',.FALSE.)-1
     FORMATO = FORMATO(1:EF)//FORMATO2(1:E2)
   END DO
   EF = INDEX(FORMATO,',',.TRUE.)-1
   WRITE (UNIT = IFILN,FMT='('//FORMATO(1:EF)//') ') WRITING(1:IVMAX)

END SUBROUTINE OUTPUT
