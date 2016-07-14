!-------------------------------------------------------------------------------
!
!  Description
!  -----------
!
!  Interface for TM3 background CO2 fields as particle boundary conditions. The
!  plain TM3 binary file format "*.b" is assumed.
!  For more description see the calling R code "getTM3.bin.r".
!
!  The algorithm does not work when the receptor time minus the backtime crosses
!  to another year. For that time instants a default value of -999 is returned.
!  All internal times count in "julian hours".
!
!  Contains the Fortran-95 extension CONVERT='BIG_ENDIAN' in two OPEN statements.
!
!  Compilation e.g. by
!     pgf95 -shared -Mbounds date_sub.f90 TM3IO.f90 getTM3bin.f90 -o getTM3bin.so
!     rm *.mod
!     chmod g+rwx *.so
!  or
!     ifort -fPIC -shared -assume byterecl date_sub.f90 TM3IO.f90 getTM3bin.f90 -o getTM3bin.so
!     rm *.mod
!     chmod g+rwx *.so
!
!-------------------------------------------------------------------------------
!  History
!  -------
!
!  09.08.2007  Initial version. Stefan Koerner
!  23.10.2007  Modified to handle annual files. Stefan Koerner
!  02.11.2007  Format issue in error message resolved. Stefan Koerner
!  18.01.2008  Outsourcing some code pieces to new module TM3IO. Stefan Koerner
!
!  $Id: getTM3bin.f90,v 1.4 2008-01-21 16:08:11 skoerner Exp $
!-------------------------------------------------------------------------------

INCLUDE 'date_sub.f90'
INCLUDE 'TM3IO.f90'


MODULE moddata
   USE TM3IO
   IMPLICIT NONE
   SAVE


   !---  metadata of the raw binary TM3 file

   INTEGER, PARAMETER :: nx = 72, ny = 48, nz = 19          ! TM3 fg spatial field dimensions
   REAL   , PARAMETER :: sigma(nz+1) = (/1.00000, 0.99000, 0.97418, 0.95459, 0.93054, 0.86642, &
                                         0.77773, 0.66482, 0.53529, 0.40334, 0.28421, 0.23286, &
                                         0.18783, 0.14916, 0.11654, 0.08945, 0.06723, 0.04920, &
                                         0.02309, 0.00000/)

   !--- I/O

   TYPE(TM3fparamtype) :: fparams_mix, fparams_ps


   !--- CO2 inversion offset

   REAL :: offset


END MODULE moddata



SUBROUTINE getTM3bin (yr4, mon, day, hr, co2inifile, n, btime, lat, lon, p, &
                      xco2)
   USE moddata
   IMPLICIT NONE

   INTEGER       , INTENT(IN)  :: yr4, mon, day, hr, n
   CHARACTER(255), INTENT(IN)  :: co2inifile
   REAL          , INTENT(IN)  :: btime(n), lat(n), lon(n), p(n)
   REAL          , INTENT(OUT) :: xco2(n)


   !-- local variables

   CHARACTER(255)     :: fname
   INTEGER            :: k, i, j, l, &
                         rec_x(n), rec_ps(n), oldrec_x
   REAL               :: tr(n)                  ! requested times, unit: fractional "julian hours"


   ! file record data

   INTEGER            :: rl, tau, date(6)
   REAL               :: x(nx,ny,nz), ps(nx,ny)

!---------------------------------------------------------------------------------------------------
   !  default behaviour for uncaught errors
   xco2 = -999.


   IF (fparams_mix%lun == -1) THEN

      ! open the mixing ratio file, with possible replacement of YYYY by the receptor year
      fname = co2inifile
      k = INDEX(fname,'YYYY')
      IF (k > 0) WRITE (fname(k:k+3),'(i4)') yr4            ! year substitution
      CALL TM3IO_Open (fname, fparams_mix)

      ! construct surface pressure file name and open
      l = INDEX(co2inifile, '/', BACK=.TRUE.)
      fname = co2inifile(1:l)//'stagc_ps_YYYY_fg.b'
      k = INDEX(fname,'YYYY')
      WRITE (fname(k:k+3),'(i4)') yr4                       ! year substitution
      CALL TM3IO_Open (fname, fparams_ps)

      ! get constant CO2 offset
      fname = co2inifile(1:l)//'c0.d'
      k = getFreeLun()
      OPEN (k, FILE=fname, STATUS='OLD', ACTION='READ', ERR=90)
      READ (k,*) offset
      CLOSE (k)

      ! file time parameter initialization
      CALL TM3IO_getFileTimes (fparams_mix)
      CALL TM3IO_getFileTimes (fparams_ps)

   END IF


   !--- compute the requested file records, determined by the requested time instants

   tr = julian(yr4, mon, day)*24+hr - btime                 ! requested reference time instants
!   WHERE (tr < t0_x) tr = tr + CEILING((t0_x-tr)/(365*24))*365*24 ! shift by years if outside range
!   WHERE (tr > t1_x) tr = tr - CEILING((tr-t1_x)/(365*24))*365*24
   rec_x = NINT((tr-fparams_mix%t0)/fparams_mix%dt) + 1

!   WHERE (tr < t0_ps) tr = tr + CEILING((t0_ps-tr)/(365*24))*365*24! shift by years if outside range
!   WHERE (tr > t1_ps) tr = tr - CEILING((tr-t1_ps)/(365*24))*365*24
   rec_ps = NINT((tr-fparams_ps%t0)/fparams_ps%dt) + 1


   !--- main loop over the records of the mixing ratio file
   ! reference times are mostly backwards ordered, thus the -1 increment

   oldrec_x = 0
   DO k=n,1,-1

      IF (rec_x(k)  < 1 .OR. rec_x(k)  > fparams_mix%lastrec .OR. &
          rec_ps(k) < 1 .OR. rec_ps(k) > fparams_ps%lastrec) THEN
         PRINT '(a)', 'getTM3bin error: reference time not in file'
         PRINT '(4(a,i0),a,f8.4)', 'yr4=', yr4, ', mon=', mon, ', day=', day, ', hr=', hr, &
                                   ', btime=', btime(k)
         CYCLE
      END IF

      IF (rec_x(k) /= oldrec_x) THEN
         READ (fparams_mix%lun,REC=rec_x(k) ) rl, tau, date, x
         READ (fparams_ps%lun ,REC=rec_ps(k)) rl, tau, date, ps
         oldrec_x = rec_x(k)
      END IF

      ! spatial indices
      i = MIN(MAX( NINT( (lon(k)+180.)*    nx/360.+1 ), 1), nx)
      j = MIN(MAX( NINT( (lat(k)+ 90.)*(ny-1)/180.+1 ), 1), ny)
      IF (j == 1 .OR. j == ny) i=1                          ! pole box

      DO l=1,nz-1
         IF (sigma(l+1)*ps(i,j) <= p(k)*100.) EXIT
      END DO

      xco2(k) = x(i,j,l) + offset

   END DO


   RETURN

!---------------------------------------------------------------------------------------------------
90 CALL errhand ('file '//TRIM(fname)//' cannot be opened.')

CONTAINS

   SUBROUTINE errhand (msg)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: msg

   !------------------------------------------------------------------------------------------------
      PRINT '(2a)', 'getTM3bin error: ', msg
      PRINT '(a)' , 'Returning with dummy data.'
      RETURN

   END SUBROUTINE errhand

END SUBROUTINE getTM3bin



!PROGRAM test
!   IMPLICIT NONE
!
!   INTEGER, PARAMETER :: nmax=100
!   CHARACTER(255)     :: co2inifile
!   INTEGER            :: yr4, mon, day, hr, n, k
!   REAL               :: btime(nmax), lat(nmax), lon(nmax), p(nmax), xco2(nmax)
!
!!---------------------------------------------------------------------------------------------------
!   OPEN (10, FILE='input.txt', ACTION='READ')
!   READ (10,*) yr4, mon, day, hr
!   READ (10,*) co2inifile
!   READ (10,*) n
!   IF (n > nmax) STOP 'n > nmax. Stop.'
!   READ (10,*) (btime(k), lat(k), lon(k), p(k), k=1,n)
!   CLOSE (10)
!
!   CALL getTM3bin (yr4, mon, day, hr, co2inifile, n, btime, lat, lon, p, &
!                   xco2)
!
!END PROGRAM test
!
!-- next line begins 'input.txt'
!2005 5 1 14
!'/Net/Groups/BSY/STILT/fluxes_input/TM3/mu1.0_070_mix_2005_fg.b'
!100
!  43.00000 35.5461  -7.4472 610.99
!  70.50000 35.6342  -8.9635 587.89
!  72.50000 35.6237  -7.3059 656.12
!  72.50000 35.6183  -9.2256 585.82
!  73.00000 35.6083  -8.4006 636.70
!  73.50000 36.4041 -10.5025 529.58
!  75.50000 35.6085  -6.2940 666.05
!  75.50000 35.6285  -7.9757 645.61
!  80.50000 35.6224  -6.5629 657.03
!  82.50000 35.6242  -6.8885 665.01
!  87.33333 35.6350  -3.5168 760.42
!  91.66667 38.8603 -10.4604 473.75
!  93.66667 38.6618 -10.5370 501.21
!  94.00000 35.9956 -10.4512 549.12
!  95.33333 41.3264 -10.5612 623.07
!  95.33333 42.1623 -10.5739 666.91
!  96.00000 41.7450 -10.5880 601.43
!  97.00000 37.7990 -10.5207 525.64
!  98.00000 38.1561 -10.4663 535.44
!  99.00000 41.3030 -10.5599 617.68
!  99.33333 41.7066 -10.5313 684.21
!  99.66667 35.8167 -10.4476 582.01
! 100.00000 38.4309 -10.4453 472.81
! 100.00000 40.8518 -10.5044 622.49
! 101.00000 37.6652 -10.4940 568.78
! 101.00000 41.0603 -10.5312 590.63
! 101.00000 39.2297 -10.5104 549.35
! 101.66667 41.2104 -10.5722 637.79
! 101.66667 39.8479 -10.5098 577.47
! 102.33333 35.5783   2.0739 768.67
! 102.33333 41.2691 -10.4825 659.65
! 103.00000 37.9277 -10.4458 546.93
! 103.33333 40.9134 -10.4949 618.64
! 103.33333 41.4514 -10.5831 664.60
! 103.66667 37.3891 -10.4962 546.38
! 104.33333 41.1100 -10.4828 604.76
! 105.33333 40.9691 -10.6034 632.60
! 105.66667 39.7322 -10.5739 618.69
! 105.66667 38.2666 -10.5052 581.71
! 106.00000 35.5939   2.3787 695.24
! 106.33333 40.1950 -10.5509 627.25
! 107.50000 35.6194   5.4086 725.81
! 112.00000 39.0602 -10.5740 628.87
! 114.50000 35.5302  -5.6852 632.18
! 115.00000 35.6354   1.2462 619.41
! 115.50000 35.6316  -5.4384 619.38
! 119.00000 35.6161  -2.3457 616.20
! 121.00000 35.6305  -5.5561 662.34
! 121.66667 39.1387 -10.4695 669.39
! 124.66667 38.8036 -10.4905 660.91
! 127.00000 35.6127  -3.8208 650.08
! 128.33333 35.6235  -4.7474 748.43
! 132.00000 38.5037 -10.4739 684.58
! 135.00000 35.6008  -3.7553 680.98
! 135.33333 36.3061 -10.4776 569.81
! 137.00000 39.5976 -10.5526 629.69
! 137.33333 35.6383  -4.6102 765.41
! 137.33333 38.6444 -10.4455 608.61
! 137.66667 39.2497 -10.4775 633.49
! 140.00000 43.3496 -10.4444 631.93
! 140.33333 39.4573 -10.5248 626.55
! 140.66667 39.5334 -10.4641 645.69
! 141.00000 35.6384  -7.1377 658.52
! 141.33333 44.4134 -10.4762 669.09
! 141.33333 44.1915 -10.4409 676.63
! 141.33333 43.5797 -10.4671 631.98
! 141.66667 39.1198 -10.5201 602.81
! 143.00000 39.1371 -10.5108 620.56
! 143.66667 43.1335 -10.5431 735.24
! 144.00000 45.5374 -10.5453 721.53
! 144.00000 37.0939 -10.5366 573.47
! 144.66667 36.6631 -10.5076 623.16
! 144.66667 45.5829 -10.4639 692.58
! 145.00000 40.4227 -10.5190 588.28
! 145.33333 38.3080 -10.4835 598.23
! 145.66667 43.9259 -10.6189 685.01
! 146.00000 45.4223 -10.5099 614.22
! 147.00000 45.4028 -10.4710 641.64
! 147.00000 39.5165 -10.4678 623.14
! 147.66667 39.5850 -10.4509 626.21
! 148.00000 41.1687 -10.4895 675.42
! 148.33333 38.4709 -10.4455 644.27
! 148.33333 41.7955 -10.4887 701.17
! 148.66667 40.9388 -10.4895 691.26
! 148.66667 39.1222 -10.5180 627.47
! 148.66667 36.9991 -10.4825 614.87
! 149.33333 44.4657 -10.4647 840.94
! 149.66667 43.2367 -10.4755 804.67
! 150.33333 44.9523 -10.4659 765.09
! 150.33333 39.5546 -10.5278 616.84
! 150.33333 37.1189 -10.4528 626.04
! 150.66667 38.5497  -8.4985 661.49
! 150.66667 39.0444  -8.9349 619.73
! 150.66667 40.5079  -7.3471 847.21
! 150.66667 41.3194  -3.4585 771.75
! 150.66667 40.2970  -6.9485 812.63
! 150.66667 43.5896  -8.4579 800.00
! 150.66667 36.1144 -10.3086 619.14
! 150.66667 37.4456  -8.2302 687.92
! 150.66667 39.0805  -9.0701 662.91
