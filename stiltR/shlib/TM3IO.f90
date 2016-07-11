!-------------------------------------------------------------------------------
!
!  Purpose
!
!  Collects several subroutines around random access of TM3 output files.
!
!-------------------------------------------------------------------------------
!  History
!  -------
!
!  01/18/2008  Initial version. Stefan Koerner
!  04/08/2008  Speed improvements and TM3IO_Close added. Stefan Koerner
!
!  $Id: TM3IO.f90,v 1.2 2008-04-08 14:08:35 skoerner Exp $
!-------------------------------------------------------------------------------

MODULE TM3IO

   IMPLICIT NONE

   TYPE TM3fparamtype
      CHARACTER(255) :: fname=''
      INTEGER        :: lun=-1                  ! use this to determine whether a file is opened
      INTEGER        :: t0, t1, dt, lastrec     ! start time, end time and time delta [julian hour]
      INTEGER        :: firstyear, lastyear
   END TYPE TM3fparamtype


   PRIVATE :: errhand                           ! subroutine at the end


CONTAINS


!***************************************************************************************************
!  Opens a TM3-structured binary file for random access and assigns partial information to fparam
!***************************************************************************************************

SUBROUTINE TM3IO_Open (fname, fparam)
   IMPLICIT NONE

   CHARACTER(*)       , INTENT(IN ) :: fname
   TYPE(TM3fparamtype), INTENT(OUT) :: fparam

   LOGICAL :: op
   INTEGER :: k, l
   REAL    :: z

!---------------------------------------------------------------------------------------------------
   !--- first check whether the file is already open
   INQUIRE (FILE=fname, OPENED=op)
   IF (op) CALL errhand ('TM3IO_Open error: file '''//TRIM(fname)//''' is already open.')


   !--- preparation of direct (non-sequential) file access
   ! check whether INTEGERs are 4 bytes long, the I/O record length is in bytes and REALs have
   ! the same I/O length as INTEGERs
   IF (BIT_SIZE(k) /= 32) CALL errhand ('TM3IO_Open error: INTEGER data type does not have 4 bytes (subr. relies on that).')
   INQUIRE (IOLENGTH=k) l
   INQUIRE (IOLENGTH=l) z
   IF (k /= 4 .OR. k /= l) CALL errhand ('TM3IO_Open error: wrong I/O length(s) of scalar data type(s).')

   ! open the TM3 file
   WRITE (*,'(3a)', ADVANCE='NO') 'Opening file ', TRIM(fname), '... '
   l = getFreeLun()
   OPEN (l, FILE=fname, ACCESS='DIRECT', RECL=4, &
         CONVERT='BIG_ENDIAN', STATUS='OLD', ACTION='READ', ERR=90)
   READ (l, REC=1) k                            !  size of the sequential record w/out length markers
   CLOSE (l)
   OPEN (l, FILE=fname, ACCESS='DIRECT', RECL=4+k+4, &
         CONVERT='BIG_ENDIAN', STATUS='OLD', ACTION='READ', ERR=90)
   PRINT '(a)', 'ok.'

   fparam%lun = l
   fparam%fname=fname
   RETURN

!---------------------------------------------------------------------------------------------------
90 CALL errhand ('TM3IO_Open error: file '//TRIM(fname)//' cannot be opened.')

END SUBROUTINE TM3IO_Open



!***************************************************************************************************
!  Fill the structure "fparam"
!***************************************************************************************************
SUBROUTINE TM3IO_getFileTimes (fparam, verbose)
   IMPLICIT NONE

   TYPE(TM3fparamtype), INTENT(INOUT) :: fparam
   LOGICAL            , OPTIONAL      :: verbose

   INTEGER :: l, tau, date(6), &
              k

!---------------------------------------------------------------------------------------------------
   !--- get first time
   READ (fparam%lun,REC=1) l, tau, date
   fparam%t0 = julian(date(1), date(2), date(3))*24+date(4)
   fparam%firstyear = date(1)

   !--- get delta, assuming constant time-step
   READ (fparam%lun,REC=2) l, tau, date
   fparam%dt = julian(date(1), date(2), date(3))*24+date(4) - fparam%t0

   !--- get last time, assuming constant time-step
   k = getFreeLun()
   OPEN (k, FILE=fparam%fname, FORM='UNFORMATTED', CONVERT='BIG_ENDIAN', STATUS='OLD', &
      ACTION='READ', POSITION='APPEND')
   BACKSPACE (k)
   READ (k) tau, date
   CLOSE (k)
   fparam%t1 = julian(date(1), date(2), date(3))*24+date(4)
   fparam%lastyear = date(1)
   fparam%lastrec = (fparam%t1-fparam%t0)/fparam%dt + 1


   !--- printout of structure fparam
   IF (PRESENT(verbose)) THEN
      IF (verbose) THEN
         PRINT '(/a)'  , 'getFileTimes: printout of structure fparam'
         PRINT '(3a)'  , '   fname = ''', TRIM(fparam%fname), ''''
         PRINT '(a,i0)', '   lun = ', fparam%lun
         PRINT '(a,i0)', '   t0 = ', fparam%t0
         PRINT '(a,i0)', '   t1 = ', fparam%t1
         PRINT '(a,i0)', '   dt = ', fparam%dt
         PRINT '(a,i0)', '   lastrec = ', fparam%lastrec
         PRINT '(a,i0)', '   firstyear = ', fparam%firstyear
         PRINT '(a,i0)', '   lastyear = ', fparam%lastyear
         PRINT *
      END IF
   END IF


END SUBROUTINE TM3IO_getFileTimes



FUNCTION julian(yy, mm, dd)
! function to get small julian day values because of usage of fractional julian days
! (less rounding errors)
   USE date_sub, ONLY: jd
   IMPLICIT NONE

   INTEGER             :: julian
   INTEGER, INTENT(IN) :: yy, mm, dd

   julian = jd(yy,mm,dd)-jd(1990,1,1)

END FUNCTION julian



FUNCTION getFreeLun()
!
!  get a none-occupied logical unit number
!
   IMPLICIT NONE

   INTEGER            :: getFreeLun

   INTEGER, PARAMETER :: maxlun=99
   LOGICAL            :: occ

   DO getFreeLun=10,maxlun                   ! look from 10
      INQUIRE(UNIT=getFreeLun, OPENED=occ)
      IF (.NOT. occ) RETURN
   END DO
   PRINT '(a,i0)', 'getFreeLun: no free unit number less or equal ', maxlun
   STOP 'getFreeLun: STOP.'

END FUNCTION getFreeLun



SUBROUTINE errhand (errmsg)
   IMPLICIT NONE

   CHARACTER(*), INTENT(IN) :: errmsg

!---------------------------------------------------------------------------------------------------
   PRINT '(/2a)', TRIM(errmsg), ' Stop.'
   STOP

END SUBROUTINE errhand



!***************************************************************************************************
!  Closes a TM3-structured binary file and sets the file handler to "not connected"
!***************************************************************************************************

SUBROUTINE TM3IO_Close (fparam)
   IMPLICIT NONE

   TYPE(TM3fparamtype), INTENT(INOUT) :: fparam

!---------------------------------------------------------------------------------------------------
   CLOSE (fparam%lun)
   fparam%lun   = -1
   fparam%fname = ''

END SUBROUTINE TM3IO_Close



END MODULE TM3IO
