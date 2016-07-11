!*******************************************************************************
!  access to IER human emission data file with aggregation
!*******************************************************************************
!
!  On the MPI Jena cluster, the compile command reads (in one line)
!     ifort -I/Net/Groups/BSY/STILT/netCDF/netcdf-3.6.2/f90 -fPIC -shared -o getfossEUnetCDF.so getfossEUnetCDF.f90
!     /Net/Groups/BSY/STILT/netCDF/netcdf-3.6.2/fortran/.libs/libnetcdff.so /Net/Groups/BSY/STILT/netCDF/netcdf-3.6.2/libsrc/.libs/libnetcdf.so
!
!-------------------------------------------------------------------------------
!  History
!  -------
!
!  11.01.2007  Initial version. Stefan Koerner
!  25.01.2007  File access optimizations for increased speed. Stefan Koerner
!  30.01.2007  Allow access on multiple files. Stefan Koerner
!  11.12.2007  Handles time-independent emissions (dim(time)==1). Stefan Koerner
!
!  $Id: getfossEUnetCDF.f90,v 1.4 2009-01-27 15:45:12 gerbig Exp $
!-------------------------------------------------------------------------------

INCLUDE 'date_sub.f90'
INCLUDE 'mrgrnk.f90'

MODULE common_data

   USE netCDF
   IMPLICIT NONE

   SAVE

   INTEGER, PARAMETER :: maxfiles=10
   CHARACTER(255)     :: fnames(maxfiles)=' '
   INTEGER            :: ncstat, ncid(maxfiles)=-HUGE(1), varid(maxfiles)=-HUGE(1)
   LOGICAL            :: time_dependent(maxfiles)

CONTAINS


   !***********************************************************************
   !  error handler for netCDF
   !***********************************************************************

   SUBROUTINE nc_handle_err (stat)
      INTEGER, INTENT(IN) :: stat
      IF (stat /= nf90_noerr) THEN
         PRINT '(2a)', 'netCDF error: ', NF90_STRERROR(stat)
         STOP 'Stop.'
      ENDIF
   END SUBROUTINE nc_handle_err


END MODULE common_data



!***********************************************************************
!  assign a file index to the given file name
!***********************************************************************

SUBROUTINE fname2fi (fname, fi)
   USE common_data
   IMPLICIT NONE

   CHARACTER(*), INTENT(IN)  :: fname                       ! filename
   INTEGER     , INTENT(OUT) :: fi                          ! file index

   INTEGER :: k

!---------------------------------------------------------------------------------------------------
   DO k=1,maxfiles
      IF (fname == fnames(k)) EXIT
   END DO
   IF (k <= maxfiles) THEN
      fi = k                                                ! already in
      RETURN
   ELSE                                                     ! not in
      ! find in fnames the first index with empty string
      DO k=1,maxfiles
         IF (fnames(k) == ' ') EXIT
      END DO
      IF (k > maxfiles) THEN
         STOP 'Subroutine fname2fi: all file entries allocated. Stop.'
      ELSE
         fi = k
         fnames(fi) = fname
         CALL IER_Open (fi)
      END IF

   END IF


END SUBROUTINE fname2fi



SUBROUTINE IER_Open (fi)
   USE common_data
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: fi                                ! file index, file name is fnames(fi)

   INTEGER             :: dimid_time, k


!---------------------------------------------------------------------------------------------------
   IF (fi < 1 .OR. fi > maxfiles) THEN
      PRINT '(2(a,i0))', 'Subroutine IER_Open: argument fi must be >= 1 and <= ', maxfiles, &
         ' but is ', fi
      STOP 'Stop.'
   END IF
   WRITE (*,'(3a)', ADVANCE='NO') 'Opening file ', TRIM(fnames(fi)), '... '
   ncstat = NF90_OPEN(TRIM(fnames(fi)), nf90_nowrite, ncid(fi))
   CALL nc_handle_err (ncstat)
   WRITE (*,'(a)', ADVANCE='NO') 'done.'


   !--- get ID of variable

   ncstat = NF90_INQ_VARID(ncid(fi), 'emission', varid(fi))
   CALL nc_handle_err (ncstat)



   !--- check whether there is a time dependence in the emission

   ncstat = NF90_INQ_DIMID(ncid(fi), 'time', dimid_time)
   CALL nc_handle_err (ncstat)
   ncstat = NF90_INQUIRE_DIMENSION(ncid(fi), dimid_time, len=k)
   CALL nc_handle_err (ncstat)
   time_dependent(fi) = k > 1
   IF (time_dependent(fi)) THEN
      PRINT *
   ELSE
      PRINT '(a)', ' (time-independent emission field)'
   END IF


END SUBROUTINE IER_Open



SUBROUTINE getfossEUnetCDF (fname, n, date, i, j, ires, jres, f)

!-------------------------------------------------------------------------------
!
!  Algorithm
!  ---------
!
!  Consult the author :)
!
!-------------------------------------------------------------------------------
   USE common_data
   USE date_sub, ONLY: iday, izlr
   USE m_mrgrnk, ONLY: mrgrnk
   IMPLICIT NONE

   CHARACTER(255), INTENT(IN)  :: fname                     ! filename
   INTEGER       , INTENT(IN)  :: n                         ! number of requested values
   INTEGER       , INTENT(IN)  :: date(5,n),  &             ! YY MM DD hh mm - reference times
                                  i(n), j(n), &             ! indices of cells on aggregated grid
                                  ires(n), jres(n)          ! cell aggregation parameters
   REAL          , INTENT(OUT) :: f(n)                      ! averaged flux values


   !--- local variables

   INTEGER, PARAMETER :: imax=376, jmax=324                 ! horiz. field dimensions in file
   INTEGER            :: fi, k, l,                   &
                         i1(n), i2(n), j1(n), j2(n), &      ! field boundaries to read from file
                         day2000, dow2000, dow, delta_dow   ! day # in 2000, day of week
   INTEGER            :: ti_file(n), &                      ! time index on file (1..)
                         it_file(n)                         ! ranking array

   !  for file access
   INTEGER            :: ti_file_last                       ! time index of the last loaded field
   LOGICAL            :: mask_ti(n)
   INTEGER            :: i1_ti, i2_ti, j1_ti, j2_ti         ! envelope for all areas at current time
   REAL               :: a(imax, jmax)


!---------------------------------------------------------------------------------------------------
   CALL fname2fi (fname, fi)


   DO k=1,n

      !--- compute the spatial boundaries of the netCDF field

      i1(k) = (i(k)-1)*ires(k)+1
      i2(k) = MIN(i(k)*ires(k),imax)
      j1(k) = (j(k)-1)*jres(k)+1
      j2(k) = MIN(j(k)*jres(k),jmax)
      IF (i1(k) < 1 .OR. i1(k) > imax .OR. j1(k) < 1 .OR. j1(k) > jmax) THEN
         PRINT '(a,i0)'      , 'Subroutine getfossEUnetCDF: Field index exceedance at input item ', k
         PRINT '(a,5(1x,i0))', 'i, j, ires, jres : ', i(k), j(k), ires(k), jres(k)
         PRINT '(a,4(1x,i0))', 'i1, i2, j1, j2   : ', i1(k), i2(k), j1(k), j2(k)
         STOP 'Stop.'
      END IF


      IF (time_dependent(fi)) THEN

         !--- map the required time instants to those of the file
         !--- - the same weekday next to the month/day in 2000

         day2000   = iday(yyyy=2000, mm=date(2,k), dd=date(3,k))
         dow2000   = izlr(yyyy=2000, mm=date(2,k), dd=date(3,k))
         dow       = izlr(yyyy=date(1,k), mm=date(2,k), dd=date(3,k))
         delta_dow = dow2000-dow
         IF (delta_dow < -3) THEN
            day2000 = day2000 - (delta_dow+7)
         ELSE IF (delta_dow > +3) THEN
            day2000 = day2000 - (delta_dow-7)
         ELSE
            day2000 = day2000 - delta_dow
         END IF

         !  force into 2000
         IF (day2000 < 1  ) day2000 = day2000 + 7
         IF (day2000 > 366) day2000 = day2000 - 7

         !  time index on file (= hour number in 2000, starting with 1)
         ti_file(k) = (day2000-1)*24 + date(4,k) + 1

      ELSE

         ti_file(k) = 1

      END IF


   END DO


   !--- for file access optimization, do a ranking of the reference times

   CALL mrgrnk (XDONT=ti_file, IRNGT=it_file)               ! on output, ti_file(it_file(:)) is
                                                            !    ordered

   !--- read the grid cells, ordered in time

   ti_file_last = -HUGE(ti_file_last)
   DO k=1,n

      l = it_file(k)

      IF (ti_file(l) /= ti_file_last) THEN

         !--- determine the spatial envelope for this time instant

         mask_ti = ti_file==ti_file(l)
         i1_ti = MINVAL(i1, MASK=mask_ti)
         i2_ti = MAXVAL(i2, MASK=mask_ti)
         j1_ti = MINVAL(j1, MASK=mask_ti)
         j2_ti = MAXVAL(j2, MASK=mask_ti)

         ncstat = NF90_GET_VAR (ncid(fi), varid(fi), start=(/i1_ti, j1_ti, ti_file(l)/), &
                                values=a(i1_ti:i2_ti,j1_ti:j2_ti))
         CALL nc_handle_err (ncstat)
         ti_file_last = ti_file(l)

      END IF

      !  return the average
      f(l) = SUM(a(i1(l):i2(l),j1(l):j2(l)))/((i2(l)-i1(l)+1)*(j2(l)-j1(l)+1))

   END DO


END SUBROUTINE getfossEUnetCDF



SUBROUTINE IER_Close (fi)
   USE common_data
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: fi

!---------------------------------------------------------------------------------------------------
   ncstat = nf90_close(ncid(fi))
   CALL nc_handle_err (ncstat)
   fnames(fi) = ' '
   ncid(fi)   = -HUGE(1)
   varid(fi)  = -HUGE(1)

END SUBROUTINE IER_Close
