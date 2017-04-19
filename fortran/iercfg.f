MODULE iercfg

  REAL :: BO3=0.000    ! plume background values for O3 (ppm)
  REAL :: AO3=0.040    ! composite ambient background for O3 (ppm)
  REAL :: PSTB= 0.125  ! proportion of smog produced that is converted to stable nitrate
  REAL :: RVOC= 0.0067 ! weighted hydrocarbon reactivity (ppm/ppm-c)
  REAL :: RISO= 0.0117 ! spiked mixure isoprene reactivity (ppm/ppm-c)
  REAL :: GAMMA= 4.71  ! smog chamber empirical constant
  REAL :: BETA= 4.09   ! NOx stochiometric coefficient for NOx limited regime

  SAVE BO3,AO3,PSTB,RVOC,RISO,GAMMA,BETA

END MODULE iercfg
