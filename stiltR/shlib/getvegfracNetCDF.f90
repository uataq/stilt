!*********************************************************************************
!  input subroutines for reading and aggregatin spatially vegetation fraction maps
!*********************************************************************************
!
!  This is a stripped-down version of VPRM_read.f90 to handle only a vegetation fraction map.
!
!  On the MPI Jena cluster, the compile command reads
!     ifort -I/Net/Groups/BSY/STILT/netCDF/netcdf-3.6.2/f90 -fPIC -shared -o getvegfracNetCDF.so getvegfracNetCDF.f90
!     /Net/Groups/BSY/STILT/netCDF/netcdf-3.6.2/fortran/.libs/libnetcdff.so /Net/Groups/BSY/STILT/netCDF/netcdf-3.6.2/libsrc/.libs/libnetcdf.so
!
!-------------------------------------------------------------------------------
!  History
!  -------
!
!  12.11.2007  Initial version. Stefan Koerner
!
!  $Id: getvegfracNetCDF.f90,v 1.2 2009-01-27 15:45:12 gerbig Exp $
!-------------------------------------------------------------------------------

MODULE common_data

   USE netCDF
   IMPLICIT NONE

   SAVE

   INTEGER     , PARAMETER :: nf=7                             ! number of files
   INTEGER     , PARAMETER :: imax=376, jmax=324               ! horiz. field dimensions in files
   CHARACTER(8), PARAMETER :: varname(nf)=(/'EVI     ', 'EVI_MAX ', 'EVI_MIN ', &   ! 1..3
                                            'LSWI    ', 'LSWI_MAX', 'LSWI_MIN', &   ! 4..6
                                            'VEG_FRA '/)                            ! 7
   INTEGER                 :: ncstat, ncid(nf)=-HUGE(1), varid(nf)=-HUGE(1)
   REAL                    :: VEG_FRA_file(8,imax,jmax)

   LOGICAL                 :: initialized=.FALSE.


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
!  assign file and variable IDs
!***********************************************************************

SUBROUTINE vegfrac_init (dpath, year)
   USE common_data
   IMPLICIT NONE

   CHARACTER(*), INTENT(IN) :: dpath
   INTEGER     , INTENT(IN) :: year

   INTEGER        :: k
   CHARACTER(256) :: fname
   CHARACTER(4)   :: YYYY

!---------------------------------------------------------------------------------------------------
   PRINT '(3a)', 'Looking for file in ''', TRIM(dpath), ''''
   DO k=7,nf

      WRITE (YYYY,'(i4)') year
      fname = 'VPRM_input_'//TRIM(varname(k))//'_d01_'//YYYY//'.nc'
      WRITE (*,'(3a)', ADVANCE='NO') 'Opening file ', TRIM(fname), '... '
      ncstat = NF90_OPEN(TRIM(dpath)//TRIM(fname), nf90_nowrite, ncid(k))
      CALL nc_handle_err (ncstat)

      ncstat = NF90_INQ_VARID(ncid(k), TRIM(varname(k)), varid(k))
      CALL nc_handle_err (ncstat)

      PRINT '(a)', 'done.'

   END DO


   !--- read the time-independent field

   ncstat = NF90_GET_VAR (ncid(7), varid(7), VEG_FRA_file)
   CALL nc_handle_err (ncstat)


   initialized = .TRUE.

END SUBROUTINE vegfrac_init



INCLUDE 'date_sub.f90'

SUBROUTINE getvegfracNetCDF (dpath, n, date, i, j, ires, jres, &
                             VEG_FRA)

!-------------------------------------------------------------------------------
!
!   dpath is the directory with trailing "/" where the following files are
!   expected:
!   VPRM_input_VEG_FRA_d01_YYYY.nc
!   Here, "YYYY" is expected to be the year in date(1)
!
!-------------------------------------------------------------------------------
   USE common_data
   USE date_sub, ONLY: iday
   IMPLICIT NONE

   CHARACTER(255), INTENT(IN)  :: dpath                     ! path to netCDF
   INTEGER       , INTENT(IN)  :: n                         ! number of requested values
   INTEGER       , INTENT(IN)  :: date(5,n),  &             ! YY MM DD hh mm - reference times
                                  i(n), j(n), &             ! indices of cells on aggregated grid
                                  ires(n), jres(n)          ! cell aggregation parameters
   REAL          , INTENT(OUT) :: VEG_FRA(8,n)

   !--- local variables

   INTEGER            :: k, p, l, m, &
                         i1(n), i2(n), j1(n), j2(n)         ! field boundaries to read from file

   EXTERNAL           :: vegfrac_init

!---------------------------------------------------------------------------------------------------
   ! error defaults
   VEG_FRA   = -999.

   IF (ANY(date(1,:) /= date(1,1))) STOP 'getvegfracNetCDF: subroutine does not handle multiple years. Stop.'

   IF (.NOT. initialized) CALL vegfrac_init (dpath, date(1,1))


   DO k=1,n

      !--- compute the spatial boundaries of the netCDF field

      i1(k) = (i(k)-1)*ires(k)+1
      i2(k) = MIN(i(k)*ires(k),imax)
      j1(k) = (j(k)-1)*jres(k)+1
      j2(k) = MIN(j(k)*jres(k),jmax)
      IF (i1(k) < 1 .OR. i1(k) > imax .OR. j1(k) < 1 .OR. j1(k) > jmax) THEN
         PRINT '(a,i0)'      , 'Subroutine getvegfracNetCDF: Field index exceedance at input item ', k
         PRINT '(a,5(1x,i0))', 'i, j, ires, jres : ', i(k), j(k), ires(k), jres(k)
         PRINT '(a,4(1x,i0))', 'i1, i2, j1, j2   : ', i1(k), i2(k), j1(k), j2(k)
         STOP 'Stop.'
      END IF

   END DO


   !--- compute the grid average

   DO k=n,1,-1                               ! backwards because this is the expected order on input

      !--- return the average

      VEG_FRA(:,k)   = 0.
      DO m=j1(k),j2(k)
         DO l=i1(k),i2(k)
            DO p=1,8
               VEG_FRA(p,k) = VEG_FRA(p,k) + VEG_FRA_file(p,l,m)
            ENDDO
         ENDDO
      ENDDO
      VEG_FRA(:,k) = VEG_FRA(:,k)/((i2(k)-i1(k)+1)*(j2(k)-j1(k)+1))

   END DO


END SUBROUTINE getvegfracNetCDF



SUBROUTINE vegfrac_Close
   USE common_data
   IMPLICIT NONE

   INTEGER :: fi

!---------------------------------------------------------------------------------------------------
   DO fi=7,nf
      ncstat = NF90_CLOSE(ncid(fi))
      CALL nc_handle_err (ncstat)
      ncid(fi)   = -HUGE(1)
      varid(fi)  = -HUGE(1)
   END DO

   initialized = .FALSE.

END SUBROUTINE vegfrac_Close
