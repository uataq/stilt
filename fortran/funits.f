!-----------------------------------------------------------
! Module to define FORTRAN I/O Unit Numbers 
!-----------------------------------------------------------
! Last Revised: 01 Apr 2004 (RRD) - initial version 
!               15 Jun 2005 (RRD) - split control & message
!               25 Jul 2008 (RRD) - gem routines
!-----------------------------------------------------------

MODULE funits
   
  INTEGER, PARAMETER :: KF01 = 10 ! Meteorological base 
  INTEGER, PARAMETER :: KF02 = 65 ! Meteorological maximum

  INTEGER, PARAMETER :: KF11 = 68 ! Conc or Traj base    
  INTEGER, PARAMETER :: KF12 = 69 ! Concentration maximum

  INTEGER, PARAMETER :: KF21 = 70 ! MESSAGE 
  INTEGER, PARAMETER :: KF22 = 71 ! STARTUP  
  INTEGER, PARAMETER :: KF23 = 72 ! Particle Dump Input
  INTEGER, PARAMETER :: KF24 = 73 ! Particle Dump Output
  INTEGER, PARAMETER :: KF25 = 74 ! CONTROL 
  INTEGER, PARAMETER :: KF26 = 75 ! namelist: SETUP.CFG (in), CONC.CFG (out); ZICONTROL (in), WINDERR (in), ZIERR (in)
  INTEGER, PARAMETER :: KFPARDAT = 76 ! PARTICLE.DAT
  INTEGER, PARAMETER :: KFJCLMSG = 80 ! JCLMESSAGE (also, before JCLMESSAGE is opened: ZSG_LEVS.IN)

  INTEGER, PARAMETER :: KF27 = 87 ! GEM concentration output
  INTEGER, PARAMETER :: KF28 = 88 ! GEM concentration dump
  INTEGER, PARAMETER :: KF29 = 89 ! GEM integrated profile 

  INTEGER, PARAMETER :: KF41 = 77 ! Landuse & ASCDATA.CFG
  INTEGER, PARAMETER :: KF42 = 78 ! Roughness length
  INTEGER, PARAMETER :: KF43 = 79 ! Terrain height

  INTEGER, PARAMETER :: KF31 = 80 ! Gridded emission array
  INTEGER, PARAMETER :: KF32 = 81 ! Temporal emission values
  INTEGER, PARAMETER :: KF33 = 82 ! IER | CHEM | PRCHEM

  INTEGER, PARAMETER :: KF50 = 90 ! lagrangian: LAGSET.CFG
  INTEGER, PARAMETER :: KF51 = 91 ! lagrangian: min out file
  INTEGER, PARAMETER :: KF52 = 99 ! lagrangian: max out file

END MODULE funits 
