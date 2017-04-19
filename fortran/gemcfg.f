!###############################################################################
! GEMCFG - global model configuration variables                             
!-------------------------------------------------------------------------------
! LAST REVISED: 19 May 2008 (RRD) - initial version 
!               25 Jul 2008 (AS)  - vertical integrated output
!               13 Nov 2008 (RRD) - changed default kinit to =0
!-------------------------------------------------------------------------------
 
MODULE gemcfg 

!---------------------------------------------------------
! model integration values  

  REAL*8    :: HMIX   = 2.0E+05  ! Horizontal mixing coefficient at equator
  REAL*4    :: VMIX   = 50.0     ! Maximum vertical mixing (50 m2/s)
  REAL*4    :: WFACT  = 1.0      ! Vertical velocity scale factor (0.0 to 1.0)
  INTEGER*4 :: KDIVG  = 0        ! Vertical velocity method (0=data 1=diverg)
  INTEGER*4 :: KMASS  = 2        ! Global Mass (0=skip 1=show 2=conserve)
  INTEGER*4 :: DELTA  = 60       ! Integration time step


!---------------------------------------------------------
! model initialization options

!                                 <0 gem routines not called
!                                  0:zero 1:bands 2:dewpt 3:gemdump 4:gbl37
  INTEGER*4 :: KINIT  = 0        ! Initialization 
  REAL*4    :: CKON   = 0.0      ! Initialization value when kinit<>above
  REAL*4    :: CMIN   = 0.0      ! Concentration range for initialization
  REAL*4    :: CMAX   = 1.0E+25
  REAL*8    :: QTOT   = 0.0      ! Total system mass


!---------------------------------------------------------
! concentration output 


  CHARACTER(80) :: GEMDUMP = 'gemdump.bin'   ! daily 3D dump file
  CHARACTER(80) :: GEMCONC = 'gemconc.bin'   ! concentration output file
  CHARACTER(80) :: GEMZINT = 'gemzint.bin'   ! vertically integrated output file

  INTEGER*4 :: KZBOT = 1         ! bottom output level index number
  INTEGER*4 :: KZTOP = 1         ! top output level index number
  INTEGER*4 :: KZAVG = 0         ! output all levels (0) or average (1)     
  INTEGER*4 :: WFREQ = 3         ! write output frequency (hours <=24)
  INTEGER*4 :: IHOUR = 0         ! initial hour for output (0Z-23Z)

  REAL*4    :: CFACT  = 1.0      ! Concentration output units conversion factor

  SAVE

END MODULE gemcfg 
