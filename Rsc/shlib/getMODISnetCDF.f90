!*******************************************************************************
!  input subroutines for driving the vegetation model VPRM incl. aggregation
!*******************************************************************************
!
!  On the MPI Jena cluster, the compile command reads (in one line)
!     ifort -I/Net/Groups/BSY/STILT/netCDF/netcdf-3.6.2/f90 -fPIC -shared -o getMODISnetCDF.so getMODISnetCDF.f90
!     /Net/Groups/BSY/STILT/netCDF/netcdf-3.6.2/fortran/.libs/libnetcdff.so /Net/Groups/BSY/STILT/netCDF/netcdf-3.6.2/libsrc/.libs/libnetcdf.so
!
!-------------------------------------------------------------------------------
!  History
!  -------
!
!  08.11.2007  Initial version. Stefan Koerner
!  29.11.2007  Modified error handling when reference time is not in file. Stefan Koerner
!  14.01.2008  More descriptive message when reference time is not in file. Stefan Koerner
!  16.07.2008  Added nest so that outside of it vegetation fraction is zero. Christoph Gerbig
!  02.07.2008  Modification to deal with standard NetCDF convention (used after June 2009) Christoph Gerbig
!
!  $Id: getMODISnetCDF.f90,v 1.8 2009-09-10 17:59:19 gerbig Exp $
!-------------------------------------------------------------------------------

MODULE common_data

   USE netCDF
   IMPLICIT NONE

   SAVE

   INTEGER     , PARAMETER :: nf=7                             ! number of files
!   CHARACTER(8), PARAMETER :: varname(nf)=(/'EVI     ', 'EVI_MAX ', 'EVI_MIN ', &   ! 1..3
!                                            'LSWI    ', 'LSWI_MAX', 'LSWI_MIN', &   ! 4..6
!                                            'VEG_FRA '/)                            ! 7
   !new format for files processed after 2009
   CHARACTER(24), PARAMETER :: varname(nf)=(/'evi                     ', &
                                             'evi_max                 ', &
                                             'evi_min                 ', &
                                             'lswi                    ', &
                                             'lswi_max                ', &
                                             'lswi_min                ', &
                                             'vegetation_fraction_map '/)                            ! 7
   CHARACTER(24), PARAMETER :: varnamef(nf)=(/'EVI                     ', &
                                             'EVI_MAX                 ', &
                                             'EVI_MIN                 ', &
                                             'LSWI                    ', &
                                             'LSWI_MAX                ', &
                                             'LSWI_MIN                ', &
                                             'VEG_FRA                 '/)                            ! 7
   INTEGER                 :: ncstat, ncid(nf)=-HUGE(1), varid(nf)=-HUGE(1)
!   REAL                    :: EVI_MAX_file(8,imax,jmax), EVI_MIN_file(8,imax,jmax),   &
!                              LSWI_MAX_file(8,imax,jmax), LSWI_MIN_file(8,imax,jmax), &
!                              VEG_FRA_file(8,imax,jmax)
   !new format for files processed after 2009
   REAL , ALLOCATABLE      :: EVI_MAX_file(:,:,:), EVI_MIN_file(:,:,:),   &
                              LSWI_MAX_file(:,:,:), LSWI_MIN_file(:,:,:), &
                              VEG_FRA_file(:,:,:)
   INTEGER                 :: startday, timax                  ! timax = number of time instants
   INTEGER                 :: imax, jmax                       ! horiz. field dimensions in files
   INTEGER                 :: nveg                             ! dimension for number of vegetation types in files

   LOGICAL                 :: initialized=.FALSE.
   LOGICAL                 :: newf                             ! flag for new file format


CONTAINS


   !***********************************************************************
   !  error handler for netCDF
   !***********************************************************************

   SUBROUTINE nc_handle_err (stat)
      INTEGER, INTENT(IN) :: stat
      IF (stat /= nf90_noerr) THEN
         PRINT '(2a)', 'netCDF error: ', NF90_STRERROR(stat)
         PRINT '(a)' , 'in subroutine getMODISnetCDF'
         STOP 'Stop.'
      ENDIF
   END SUBROUTINE nc_handle_err


END MODULE common_data



!***********************************************************************
!  assign file and variable IDs
!***********************************************************************

SUBROUTINE VPRM_init (dpath, year)
   USE common_data
   IMPLICIT NONE

   CHARACTER(*), INTENT(IN) :: dpath
   INTEGER     , INTENT(IN) :: year

   INTEGER        :: k
   CHARACTER(256) :: fname
   CHARACTER(4)   :: YYYY

!---------------------------------------------------------------------------------------------------
   PRINT '(3a)', 'Looking for files in ''', TRIM(dpath), ''''
   DO k=1,nf

      WRITE (YYYY,'(i4)') year
!      fname = 'VPRM_input_'//TRIM(varname(k))//'_d01_'//YYYY//'.nc'
   !new format for files processed after 2009
      fname = 'VPRM_input_'//TRIM(varnamef(k))//'_'//YYYY//'.nc'
      INQUIRE(FILE=TRIM(dpath)//TRIM(fname),EXIST=newf)
      IF (.NOT. newf) THEN
         fname = 'VPRM_input_'//TRIM(varnamef(k))//'_d01_'//YYYY//'.nc'
      END IF

      WRITE (*,'(3a)', ADVANCE='NO') 'Opening file ', TRIM(fname), '... '
      ncstat = NF90_OPEN(TRIM(dpath)//TRIM(fname), nf90_nowrite, ncid(k))
      CALL nc_handle_err (ncstat)

      IF (newf) THEN
         ncstat = NF90_INQ_VARID(ncid(k), TRIM(varname(k)), varid(k))
      ELSE
         ncstat = NF90_INQ_VARID(ncid(k), TRIM(varnamef(k)), varid(k))
      END IF
      CALL nc_handle_err (ncstat)

      PRINT '(a)', 'done.'

   END DO

   !--- allocate time-independent using dimensions from NetCDF files

   IF (newf) THEN
   PRINT '(a)'      , 'NEW FILE FORMAT'
      ncstat = NF90_INQ_DIMID(ncid(1), 'lon', k)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQUIRE_DIMENSION(ncid(1), k, len=imax)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQ_DIMID(ncid(1), 'lat', k)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQUIRE_DIMENSION(ncid(1), k, len=jmax)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQ_DIMID(ncid(1), 'vprm_classes', k)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQUIRE_DIMENSION(ncid(1), k, len=nveg)
      CALL nc_handle_err (ncstat)
      allocate(EVI_MAX_file(imax,jmax,nveg))
      allocate(EVI_MIN_file(imax,jmax,nveg))
      allocate(LSWI_MAX_file(imax,jmax,nveg))
      allocate(LSWI_MIN_file(imax,jmax,nveg))
      allocate(VEG_FRA_file(imax,jmax,nveg))
   ELSE
   PRINT '(a)'      , 'OLD FILE FORMAT'
      ncstat = NF90_INQ_DIMID(ncid(1), 'west_east', k)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQUIRE_DIMENSION(ncid(1), k, len=imax)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQ_DIMID(ncid(1), 'south_north', k)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQUIRE_DIMENSION(ncid(1), k, len=jmax)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQ_DIMID(ncid(1), 'vprm_classes', k)
      CALL nc_handle_err (ncstat)
      ncstat = NF90_INQUIRE_DIMENSION(ncid(1), k, len=nveg)
      CALL nc_handle_err (ncstat)
      allocate(EVI_MAX_file(nveg,imax,jmax))
      allocate(EVI_MIN_file(nveg,imax,jmax))
      allocate(LSWI_MAX_file(nveg,imax,jmax))
      allocate(LSWI_MIN_file(nveg,imax,jmax))
      allocate(VEG_FRA_file(nveg,imax,jmax))
   END IF

   PRINT '(a)'      , 'Subroutine getMODISnetCDF, VPRM_init: dimensions for static fields:'
   PRINT '(a,3(1x,i0))', 'lon, lat, veg.classes   : ', imax, jmax, nveg

   !--- read the time-independent fields

   ncstat = NF90_GET_VAR (ncid(2), varid(2), EVI_MAX_file)
   CALL nc_handle_err (ncstat)
   ncstat = NF90_GET_VAR (ncid(3), varid(3), EVI_MIN_file)
   CALL nc_handle_err (ncstat)
   ncstat = NF90_GET_VAR (ncid(5), varid(5), LSWI_MAX_file)
   CALL nc_handle_err (ncstat)
   ncstat = NF90_GET_VAR (ncid(6), varid(6), LSWI_MIN_file)
   CALL nc_handle_err (ncstat)
   ncstat = NF90_GET_VAR (ncid(7), varid(7), VEG_FRA_file)
   CALL nc_handle_err (ncstat)


   !--- get the start day of the year for computations time -> file indices

!   ncstat = NF90_INQ_VARID(ncid(1), 'START_DAY_OF_YEAR', k)
   !new format for files processed after 2009
   IF (newf) THEN
      ncstat = NF90_INQ_VARID(ncid(1), 'start_day_of_year', k)
   ELSE
      ncstat = NF90_INQ_VARID(ncid(1), 'START_DAY_OF_YEAR', k)
   END IF
   CALL nc_handle_err (ncstat)
   ncstat = NF90_GET_VAR (ncid(1), k, startday)
   CALL nc_handle_err (ncstat)


   !--- get the number of 8-day intervals in the files for time checking

!   ncstat = NF90_INQ_DIMID(ncid(1), 'timesteps', k)
   !new format for files processed after 2009
   IF (newf) THEN
      ncstat = NF90_INQ_DIMID(ncid(1), 'time', k)
   ELSE
      ncstat = NF90_INQ_DIMID(ncid(1), 'timesteps', k)
   END IF
   CALL nc_handle_err (ncstat)
   ncstat = NF90_INQUIRE_DIMENSION(ncid(1), k, len=timax)
   CALL nc_handle_err (ncstat)

   initialized = .TRUE.

END SUBROUTINE VPRM_init



INCLUDE 'date_sub.f90'

SUBROUTINE getMODISnetCDF (dpath, n, date, i, j, ires, jres, subill, subjll, subiur, subjur, nv, &
                           EVI, EVI_amax, EVI_amin, LSWI, LSWI_amax, LSWI_amin, VEG_FRA)

!-------------------------------------------------------------------------------
!
!   dpath is the directory with trailing "/" where the following files are
!   expected:
!   VPRM_input_EVI_d01_YYYY.nc
!   VPRM_input_EVI_MAX_d01_YYYY.nc
!   VPRM_input_EVI_MIN_d01_YYYY.nc
!   VPRM_input_LSWI_d01_YYYY.nc
!   VPRM_input_LSWI_MAX_d01_YYYY.nc
!   VPRM_input_LSWI_MIN_d01_YYYY.nc
!   VPRM_input_VEG_FRA_d01_YYYY.nc
!   Here, "YYYY" is expected to be the year in date(1)
!
!
!-------------------------------------------------------------------------------
   USE common_data
   USE date_sub, ONLY: iday
   IMPLICIT NONE

   CHARACTER(255), INTENT(IN)  :: dpath                     ! path to netCDFs
   INTEGER       , INTENT(IN)  :: n                         ! number of requested values
   INTEGER       , INTENT(IN)  :: date(5,n),  &             ! YY MM DD hh mm - reference times
                                  i(n), j(n), &             ! indices of cells on aggregated grid
                                  ires(n), jres(n), &       ! cell aggregation parameters
                                  subill, subjll, &         ! subgrid specification - lower left indices
                                  subiur, subjur            ! subgrid specification - upper right indices
   INTEGER       , INTENT(IN)  :: nv                        ! number of expected number of vegetation classes
                                                            ! (needed to pass on from R so can define EVI array extent here)
   REAL          , INTENT(OUT) :: EVI(nv,n),                     &    ! averaged species values
                                  EVI_amin(nv,n), EVI_amax(nv,n), &    ! annual min and max
                                  LSWI(nv,n),                    &
                                  LSWI_amin(nv,n), LSWI_amax(nv,n), VEG_FRA(nv,n)
!   REAL , ALLOCATABLE , INTENT(OUT) :: EVI(:,:),                     &    ! averaged species values
!                                  EVI_amin(:,:), EVI_amax(:,:), &    ! annual min and max
!                                  LSWI(:,:),                    &
!                                  LSWI_amin(:,:), LSWI_amax(:,:), VEG_FRA(:,:)

   !--- local variables

   INTEGER            :: k, p, l, m, &
                         i1(n), i2(n), j1(n), j2(n)         ! field boundaries to read from file
   INTEGER            :: ti_file(n)                         ! time index on file (1..)
   REAL   , PARAMETER :: missing=-999.                      ! result at undefined space/time points or error
   LOGICAL            :: warning_issued=.FALSE.

   !  for file access
   INTEGER            :: ti_file_last                       ! time index of the last loaded field
   LOGICAL            :: mask_ti(n)
   INTEGER            :: i1_ti, i2_ti, j1_ti, j2_ti         ! envelope for all areas at current time
   REAL , ALLOCATABLE :: EVI_file(:,:,:), LSWI_file(:,:,:)


   EXTERNAL           :: VPRM_init

!---------------------------------------------------------------------------------------------------
   IF (ANY(date(1,:) /= date(1,1))) STOP 'getMODISnetCDF: subroutine does not handle multiple years. Stop.'

   IF (.NOT. initialized) CALL VPRM_init (dpath, date(1,1))

   ! error defaults
   EVI       = missing
   EVI_amin  = missing
   EVI_amax  = missing
   LSWI      = missing
   LSWI_amin = missing
   LSWI_amax = missing
   VEG_FRA   = missing



   DO k=1,n

      !--- compute the spatial boundaries of the netCDF fields

      i1(k) = (i(k)-1)*ires(k)+1
      i2(k) = MIN(i(k)*ires(k),imax)
      j1(k) = (j(k)-1)*jres(k)+1
      j2(k) = MIN(j(k)*jres(k),jmax)
      IF (i1(k) < 1 .OR. i1(k) > imax .OR. j1(k) < 1 .OR. j1(k) > jmax) THEN
         PRINT '(a,i0)'      , 'Subroutine getMODISnetCDF: Field index exceedance at input item ', k
         PRINT '(a,5(1x,i0))', 'i, j, ires, jres : ', i(k), j(k), ires(k), jres(k)
         PRINT '(a,4(1x,i0))', 'i1, i2, j1, j2   : ', i1(k), i2(k), j1(k), j2(k)
         STOP 'Stop.'
      END IF


      !--- map the required time instants to those of the file

      !  time index on file
      ti_file(k) = (iday(date(1,k),date(2,k),date(3,k)) - startday)/8 + 1
      IF (ti_file(k) < 1 .OR. ti_file(k) > timax) THEN
         ti_file(k) = -1          ! "not in file" value
         IF (.NOT. warning_issued) THEN
            PRINT '(/a)'                , 'Notice (subroutine getMODISnetCDF): requested reference time not in file.'
            PRINT '(a,i0,a,i4,4(1x,i2))', 'Input item ', k, ', requested date is ', date(:,k)
            PRINT '(a,f5.0,a/)'         , 'Missing values ', missing, ' assigned at these points.'
            ! PRINT *, iday(date(1,k),date(2,k),date(3,k)), startday
            ! PRINT *, ti_file(k), timax
            warning_issued = .TRUE.
         END IF
      END IF

   END DO




   !--- read at the given space/time coordinates
   IF (newf) THEN
      allocate(EVI_file(imax,jmax,nveg))
      allocate(LSWI_file(imax,jmax,nveg))
   ELSE
      allocate(EVI_file(nveg,imax,jmax))
      allocate(LSWI_file(nveg,imax,jmax))
   END IF

   ti_file_last = -HUGE(ti_file_last)
   DO k=n,1,-1                               ! backwards because this is the expected order on input

      IF (ti_file(k) == -1) THEN

         EVI_file  = -999.                   ! "not in file" values
         LSWI_file = -999.

      ELSE IF (ti_file(k) /= ti_file_last) THEN

         !--- determine the spatial envelope for this time instant

         mask_ti = ti_file==ti_file(k)
         i1_ti = MINVAL(i1, MASK=mask_ti)
         i2_ti = MAXVAL(i2, MASK=mask_ti)
         j1_ti = MINVAL(j1, MASK=mask_ti)
         j2_ti = MAXVAL(j2, MASK=mask_ti)

   !new format for files processed after 2009
         IF (newf) THEN
            ncstat = NF90_GET_VAR (ncid(1), varid(1), start=(/i1_ti, j1_ti, ti_file(k),1/), &
                                   count=(/i2_ti-i1_ti+1,j2_ti-j1_ti+1,1,nveg/), &
                                   values=EVI_file(i1_ti:i2_ti,j1_ti:j2_ti,:))
            CALL nc_handle_err (ncstat)
            ncstat = NF90_GET_VAR (ncid(4), varid(4), start=(/i1_ti, j1_ti, ti_file(k),1/), &
                                   count=(/i2_ti-i1_ti+1,j2_ti-j1_ti+1,1,nveg/), &
                                   values=LSWI_file(i1_ti:i2_ti,j1_ti:j2_ti,:))
            CALL nc_handle_err (ncstat)
         ELSE
            ncstat = NF90_GET_VAR (ncid(1), varid(1), start=(/1,i1_ti, j1_ti, ti_file(k)/), &
                                   values=EVI_file(:,i1_ti:i2_ti,j1_ti:j2_ti))
            CALL nc_handle_err (ncstat)
            ncstat = NF90_GET_VAR (ncid(4), varid(4), start=(/1,i1_ti, j1_ti, ti_file(k)/), &
                                   values=LSWI_file(:,i1_ti:i2_ti,j1_ti:j2_ti))
            CALL nc_handle_err (ncstat)
         END IF

         ti_file_last = ti_file(k)

      END IF


      !--- return the weighted average, except of VEG_FRA

      EVI(:,k)       = 0.
      EVI_amin(:,k)  = 0.
      EVI_amax(:,k)  = 0.
      LSWI(:,k)      = 0.
      LSWI_amin(:,k) = 0.
      LSWI_amax(:,k) = 0.
      VEG_FRA(:,k)   = 0.
!      set VEG_FRA to zero outside of nest
      if ( subill>i1(7) ) VEG_FRA_file(:,i1(7):(subill-1),:) = 0.
      if ( subiur<imax ) VEG_FRA_file(:,(subiur+1):imax,:) = 0.
      if ( subjll>j1(7) ) VEG_FRA_file(:,:,j1(7):(subjll-1)) = 0.
      if ( subjur<jmax ) VEG_FRA_file(:,:,(subjur+1):jmax) = 0.

      if (newf) THEN
         DO p=1,nveg
            DO m=j1(k),j2(k)
               DO l=i1(k),i2(k)
                  EVI(p,k)       = EVI(p,k)       + EVI_file(l,m,p)     *VEG_FRA_file(l,m,p)
                  EVI_amax(p,k)  = EVI_amax(p,k)  + EVI_MAX_file(l,m,p) *VEG_FRA_file(l,m,p)
                  EVI_amin(p,k)  = EVI_amin(p,k)  + EVI_MIN_file(l,m,p) *VEG_FRA_file(l,m,p)
                  LSWI(p,k)      = LSWI(p,k)      + LSWI_file(l,m,p)    *VEG_FRA_file(l,m,p)
                  LSWI_amax(p,k) = LSWI_amax(p,k) + LSWI_MAX_file(l,m,p)*VEG_FRA_file(l,m,p)
                  LSWI_amin(p,k) = LSWI_amin(p,k) + LSWI_MIN_file(l,m,p)*VEG_FRA_file(l,m,p)
                  VEG_FRA(p,k)   = VEG_FRA(p,k)   + VEG_FRA_file(l,m,p)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO m=j1(k),j2(k)
            DO l=i1(k),i2(k)
               DO p=1,8
                  EVI(p,k)       = EVI(p,k)       + EVI_file(p,l,m)     *VEG_FRA_file(p,l,m)
                  EVI_amax(p,k)  = EVI_amax(p,k)  + EVI_MAX_file(p,l,m) *VEG_FRA_file(p,l,m)
                  EVI_amin(p,k)  = EVI_amin(p,k)  + EVI_MIN_file(p,l,m) *VEG_FRA_file(p,l,m)
                  LSWI(p,k)      = LSWI(p,k)      + LSWI_file(p,l,m)    *VEG_FRA_file(p,l,m)
                  LSWI_amax(p,k) = LSWI_amax(p,k) + LSWI_MAX_file(p,l,m)*VEG_FRA_file(p,l,m)
                  LSWI_amin(p,k) = LSWI_amin(p,k) + LSWI_MIN_file(p,l,m)*VEG_FRA_file(p,l,m)
                  VEG_FRA(p,k)   = VEG_FRA(p,k)   + VEG_FRA_file(p,l,m)
               ENDDO
            ENDDO
         ENDDO
      END IF
      WHERE (ABS(VEG_FRA(:,k)) > TINY(1.))
         EVI(:,k)       = EVI(:,k)      /VEG_FRA(:,k)
         EVI_amax(:,k)  = EVI_amax(:,k) /VEG_FRA(:,k)
         EVI_amin(:,k)  = EVI_amin(:,k) /VEG_FRA(:,k)
         LSWI(:,k)      = LSWI(:,k)     /VEG_FRA(:,k)
         LSWI_amax(:,k) = LSWI_amax(:,k)/VEG_FRA(:,k)
         LSWI_amin(:,k) = LSWI_amin(:,k)/VEG_FRA(:,k)
      END WHERE                              ! where VEG_FRA==0, every weighted quantity is also 0
                                             ! because all contributions to VEG_FRA must then be 0
      VEG_FRA(:,k) = VEG_FRA(:,k)/((i2(k)-i1(k)+1)*(j2(k)-j1(k)+1))

   END DO


END SUBROUTINE getMODISnetCDF



SUBROUTINE VPRM_Close
   USE common_data
   IMPLICIT NONE

   INTEGER :: fi

!---------------------------------------------------------------------------------------------------
   DO fi=1,nf
      ncstat = NF90_CLOSE(ncid(fi))
      CALL nc_handle_err (ncstat)
      ncid(fi)   = -HUGE(1)
      varid(fi)  = -HUGE(1)
   END DO

   initialized = .FALSE.

END SUBROUTINE VPRM_Close
