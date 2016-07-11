!-------------------------------------------------------------------------------
!
!  Description
!  -----------
!
!  Interface for TM3 background CH4 fields as particle boundary conditions. The
!  plain TM3 binary file format "*.b" is assumed.
!  For more description see the calling R code "get.TM3CH4.bin.r".
!
!  The algorithm does not work when the receptor time minus the backtime crosses
!  a year border. For that time instants a default value of -999 is returned.
!  All internal times count in "julian hours".
!
!  Contains the Fortran-95 extension CONVERT='BIG_ENDIAN' in two OPEN statements.
!
!  Compilation e.g. by
!     pgf95 -shared -Mbounds date_sub.f90 TM3IO.f90 getTM3CH4bin.f90 -o getTM3CH4bin.so
!     rm *.mod
!     chmod g+rwx *.so
!  or
!     ifort -fPIC -shared -assume byterecl date_sub.f90 TM3IO.f90 getTM3CH4bin.f90 -o getTM3CH4bin.so
!     rm *.mod
!     chmod g+rwx *.so
!
!-------------------------------------------------------------------------------
!  History
!  -------
!
!  09.08.2007  Initial version. Stefan Koerner
!  23.10.2007  Modified to handle annual files. Stefan Koerner
!  23.10.2007  Modified (and renamed) to handle CH4 files. Christoph Gerbig
!  18.01.2008  Outsourcing some code pieces to new module TM3IO. Stefan Koerner
!
!  $Id: getTM3CH4bin.f90,v 1.5 2008-01-21 16:08:12 skoerner Exp $
!-------------------------------------------------------------------------------

INCLUDE 'date_sub.f90'
INCLUDE 'TM3IO.f90'


MODULE moddata
   USE TM3IO
   IMPLICIT NONE
   SAVE


   !---  metadata of the raw binary TM3 file

   INTEGER, PARAMETER :: nx = 72, ny = 48, nz = 19, nt = 3          ! TM3 fg spatial field dimensions
   REAL   , PARAMETER :: sigma(nz+1) = (/1.00000, 0.99000, 0.97418, 0.95459, 0.93054, 0.86642, &
                                         0.77773, 0.66482, 0.53529, 0.40334, 0.28421, 0.23286, &
                                         0.18783, 0.14916, 0.11654, 0.08945, 0.06723, 0.04920, &
                                         0.02309, 0.00000/)

   !--- I/O

   TYPE(TM3fparamtype) :: fparams_mix, fparams_ps


END MODULE moddata



SUBROUTINE getTM3CH4bin (yr4, mon, day, hr, ch4inifile, n, btime, lat, lon, p, &
                         xch4)
   USE moddata
   USE date_sub, ONLY: jd, cdate
   IMPLICIT NONE

   INTEGER       , INTENT(IN)  :: yr4, mon, day, hr, n
   CHARACTER(255), INTENT(IN)  :: ch4inifile
   REAL          , INTENT(IN)  :: btime(n), lat(n), lon(n), p(n)
   REAL          , INTENT(OUT) :: xch4(n)


   !-- local variables

   CHARACTER(255)     :: fname
   INTEGER            :: k, i, j, l, &
                         rec_x(n), rec_ps(n), oldrec_x
   LOGICAL            :: warning_issued


   ! file record data

   INTEGER            :: rl, tau, date(6), yyyy, mm, dd
   REAL               :: x(nx,ny,nz), ps(nx,ny)

!---------------------------------------------------------------------------------------------------
   !  default behaviour for uncaught errors
   xch4 = -999.
   !--- for writing a file with test data
   ! OPEN (28,file='/Net/Groups/BSY/tpobw/input.txt')
   ! WRITE (28,*)  yr4, mon, day, hr
   ! WRITE (28,3a) ''', TRIM(ch4inifile), ''''
   ! WRITE (28,*)  n
   ! DO k=1,n
   !    WRITE (28,*) btime(k), lat(k), lon(k), p(k)
   ! ENDDO
   ! CLOSE (28)


   IF (fparams_mix%lun == -1) THEN

      ! open the files, with possible replacement of YYYY by the receptor year
      fname = ch4inifile
      k = INDEX(fname,'YYYY')
      IF (k > 0) WRITE (fname(k:k+3),'(i4)') yr4      ! year substitution
      CALL TM3IO_Open (fname, fparams_mix)

      ! construct surface pressure file name and open
      l = INDEX(ch4inifile, '/', BACK=.TRUE.)
      fname = ch4inifile(1:l)//'stagc_ps_mm_YYYY_fg.b'
      k = INDEX(fname,'YYYY')
      WRITE (fname(k:k+3),'(i4)') yr4      ! year substitution
      CALL TM3IO_Open (fname, fparams_ps)


      ! file time parameter initialization
      CALL TM3IO_getFileTimes (fparams_mix)
      CALL TM3IO_getFileTimes (fparams_ps)

   END IF


   !--- compute the requested file records, determined by the requested time instants

   warning_issued = .FALSE.
   DO k = 1, n
      CALL cdate(jd(yr4, mon, day) - NINT(btime(k)/24), yyyy, mm, dd)
      IF ((yyyy < fparams_mix%firstyear .OR. yyyy > fparams_mix%lastyear) &
           .AND. .NOT. warning_issued) THEN
         PRINT '(/a)'                , 'Notice (subroutine getTM3CH4bin): requested reference time not in file.'
         PRINT '(a,i0,a,i4,2(1x,i2))', 'Input item ', k, ', requested date is ', yyyy, mm, dd
         PRINT '(2(a,i0)/)'          , 'Resetting year to ', fparams_mix%firstyear, ' <= yyyy <= ', fparams_mix%lastyear
         warning_issued = .TRUE.
      END IF
      yyyy = MIN(MAX(yyyy,fparams_mix%firstyear),fparams_mix%lastyear)
      rec_x(k) = mm + (yyyy-2001) * 12
      rec_ps(k) = mm
   END DO



   !--- main loop over the records of the mixing ratio file
   ! reference times are mostly backwards ordered, thus the -1 increment

   oldrec_x = 0
   DO k=n,1,-1
      IF (rec_x(k)  < 1 .OR. rec_x(k)  > fparams_mix%lastrec .OR. &
          rec_ps(k) < 1 .OR. rec_ps(k) > fparams_ps%lastrec) THEN
         PRINT '(a)', 'getTM3CH4bin error: reference time not in file'
         PRINT '(3(a,i0),a,f9.4)', 'yr4=', yr4, ', mon=', mon, 'day=', day, 'btime=', btime(k)
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

      xch4(k) = x(i,j,l)

   END DO


   RETURN

!---------------------------------------------------------------------------------------------------
90 CALL errhand ('file '//TRIM(fname)//' cannot be opened.')


CONTAINS

   SUBROUTINE errhand (msg)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: msg

   !------------------------------------------------------------------------------------------------
      PRINT '(2a)', 'getTM3CH4bin error: ', msg
      PRINT '(a)' , 'Returning with dummy data.'
      RETURN

   END SUBROUTINE errhand

END SUBROUTINE getTM3CH4bin



! PROGRAM test
!    IMPLICIT NONE
!
!    INTEGER, PARAMETER :: nmax=100
!    CHARACTER(255)     :: ch4inifile
!    INTEGER            :: yr4, mon, day, hr, n, k
!    REAL               :: btime(nmax), lat(nmax), lon(nmax), p(nmax), xch4(nmax)
!
! !---------------------------------------------------------------------------------------------------
!    OPEN (10, FILE='input.txt', ACTION='READ')
!    READ (10,*) yr4, mon, day, hr
!    READ (10,*) ch4inifile
!    READ (10,*) n
!    IF (n > nmax) STOP 'n > nmax. Stop.'
!    READ (10,*) (btime(k), lat(k), lon(k), p(k), k=1,n)
!    CLOSE (10)
!
!    CALL getTM3CH4bin (yr4, mon, day, hr, ch4inifile, n, btime, lat, lon, p, &
!                       xch4)
!
! END PROGRAM test
