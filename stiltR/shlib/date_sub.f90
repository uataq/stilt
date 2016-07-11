MODULE date_sub

!      COLLECTED AND PUT TOGETHER JANUARY 1972, H. D. KNOBLE .
!      ORIGINAL REFERENCES ARE CITED IN EACH ROUTINE.

! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-22  Time: 10:23:47

IMPLICIT NONE

CONTAINS

!        ARITHMETIC FUNCTIONS 'IZLR' AND 'IDAY' ARE TAKEN FROM REMARK ON
!        ALGORITHM 398, BY J. DOUGLAS ROBERTSON, CACM 15(10):918.

FUNCTION iday(yyyy, mm, dd) RESULT(ival)
!------IDAY IS A COMPANION TO CALEND; GIVEN A CALENDAR DATE, YYYY, MM,
!           DD, IDAY IS RETURNED AS THE DAY OF THE YEAR.
!           EXAMPLE: IDAY(1984, 4, 22) = 113

INTEGER, INTENT(IN) :: yyyy, mm, dd
INTEGER             :: ival

ival = 3055*(mm+2)/100 - (mm+10)/13*2 -91 + (1-(MOD(yyyy, 4)+3)/4 +  &
       (MOD(yyyy, 100) + 99)/100 - (MOD(yyyy, 400)+399)/400)*(mm+10)/13 + dd

RETURN
END FUNCTION iday


FUNCTION izlr(yyyy, mm, dd) RESULT(ival)
!------IZLR(YYYY, MM, DD) GIVES THE WEEKDAY NUMBER 0 = SUNDAY, 1 = MONDAY,
!      ... 6 = SATURDAY.  EXAMPLE: IZLR(1970, 1, 1) = 4 = THURSDAY

INTEGER, INTENT(IN) :: yyyy, mm, dd
INTEGER             :: ival

ival = MOD((13*(mm+10-(mm+10)/13*12)-1)/5 + dd + 77 + 5*(yyyy+(mm-14)/12 -  &
           (yyyy+(mm-14)/12)/100*100)/4 + (yyyy+(mm-14)/12)/400 -  &
           (yyyy+(mm-14)/12)/100*2, 7)

RETURN
END FUNCTION izlr


SUBROUTINE calend(yyyy, ddd, mm, dd)
!=============CALEND WHEN GIVEN A VALID YEAR, YYYY, AND DAY OF THE YEAR, DDD,
!             RETURNS THE MONTH, MM, AND DAY OF THE MONTH, DD.
!             SEE ACM ALGORITHM 398, TABLELESS DATE CONVERSION, BY
!             DICK STONE, CACM 13(10):621.

INTEGER, INTENT(IN)   :: yyyy
INTEGER, INTENT(IN)   :: ddd
INTEGER, INTENT(OUT)  :: mm
INTEGER, INTENT(OUT)  :: dd

INTEGER :: t

t = 0
IF(MOD(yyyy, 4) == 0) t = 1

!-----------THE FOLLOWING STATEMENT IS NECESSARY IF YYYY IS < 1900 OR > 2100.
IF(MOD(yyyy, 400) /= 0 .AND. MOD(yyyy, 100) == 0) t = 0

dd = ddd
IF(ddd > 59+t) dd = dd + 2 - t
mm = ((dd+91)*100)/3055
dd = (dd+91) - (mm*3055)/100
mm = mm - 2
!----------MM WILL BE CORRECT IFF DDD IS CORRECT FOR YYYY.
IF(mm >= 1 .AND. mm <= 12) RETURN
WRITE(*,1) ddd
1 FORMAT('0$$$CALEND: DAY OF THE YEAR INPUT =', i11,  ' IS OUT OF RANGE.')

STOP
END SUBROUTINE calend


SUBROUTINE cdate(jd, yyyy, mm, dd)
!=======GIVEN A JULIAN DAY NUMBER, NNNNNNNN, YYYY,MM,DD ARE RETURNED AS THE
!              CALENDAR DATE. JD = NNNNNNNN IS THE JULIAN DATE FROM AN EPOCH
!              IN THE VERY DISTANT PAST.  SEE CACM 1968 11(10):657,
!              LETTER TO THE EDITOR BY FLIEGEL AND VAN FLANDERN.
!    EXAMPLE CALL CDATE(2440588, YYYY, MM, DD) RETURNS 1970 1 1 .

INTEGER, INTENT(IN)   :: jd
INTEGER, INTENT(OUT)  :: yyyy
INTEGER, INTENT(OUT)  :: mm
INTEGER, INTENT(OUT)  :: dd

INTEGER :: l, n

l = jd + 68569
n = 4*l/146097
l = l - (146097*n + 3)/4
yyyy = 4000*(l+1)/1461001
l = l - 1461*yyyy/4 + 31
mm = 80*l/2447
dd = l - 2447*mm/80
l = mm/11
mm = mm + 2 - 12*l
yyyy = 100*(n-49) + yyyy + l
RETURN
END SUBROUTINE cdate


SUBROUTINE daysub(jd, yyyy, mm, dd, wd, ddd)
!========GIVEN JD, A JULIAN DAY # (SEE ASF JD), THIS ROUTINE CALCULATES DD,
!        THE DAY NUMBER OF THE MONTH; MM, THE MONTH NUMBER; YYYY THE YEAR;
!        WD THE WEEKDAY NUMBER, AND DDD THE DAY NUMBER OF THE YEAR.

!   EXAMPLE: CALL DAYSUB(2440588, YYYY, MM, DD, WD, DDD) YIELDS 1970 1 1 4 1.

INTEGER, INTENT(IN)   :: jd
INTEGER, INTENT(OUT)  :: yyyy
INTEGER, INTENT(OUT)  :: mm
INTEGER, INTENT(OUT)  :: dd
INTEGER, INTENT(OUT)  :: wd
INTEGER, INTENT(OUT)  :: ddd

CALL cdate(jd, yyyy, mm, dd)
wd = izlr(yyyy, mm, dd)
ddd = iday(yyyy, mm, dd)

RETURN
END SUBROUTINE daysub


FUNCTION jd(yyyy, mm, dd) RESULT(ival)

INTEGER, INTENT(IN)  :: yyyy
INTEGER, INTENT(IN)  :: mm
INTEGER, INTENT(IN)  :: dd
INTEGER              :: ival

!              DATE ROUTINE JD(YYYY, MM, DD) CONVERTS CALENDER DATE TO
!              JULIAN DATE.  SEE CACM 1968 11(10):657, LETTER TO THE
!              EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
!    EXAMPLE JD(1970, 1, 1) = 2440588

ival = dd - 32075 + 1461*(yyyy+4800+(mm-14)/12)/4 +  &
       367*(mm-2-((mm-14)/12)*12)/12 - 3*((yyyy+4900+(mm-14)/12)/100)/4

RETURN
END FUNCTION jd


FUNCTION ndays(mm1, dd1, yyyy1, mm2, dd2, yyyy2) RESULT(ival)

INTEGER, INTENT(IN)  :: mm1
INTEGER, INTENT(IN)  :: dd1
INTEGER, INTENT(IN)  :: yyyy1
INTEGER, INTENT(IN)  :: mm2
INTEGER, INTENT(IN)  :: dd2
INTEGER, INTENT(IN)  :: yyyy2
INTEGER              :: ival

!==============NDAYS IS RETURNED AS THE NUMBER OF DAYS BETWEEN TWO
!              DATES; THAT IS  MM1/DD1/YYYY1 MINUS MM2/DD2/YYYY2,
!              WHERE DATEI AND DATEJ HAVE ELEMENTS MM, DD, YYYY.
!-------NDAYS WILL BE POSITIVE IFF DATE1 IS MORE RECENT THAN DATE2.

ival = jd(yyyy1, mm1, dd1) - jd(yyyy2, mm2, dd2)

RETURN
END FUNCTION ndays


SUBROUTINE date_stamp( string, want_day, short )
CHARACTER (LEN=*), INTENT(OUT)  :: string
LOGICAL, INTENT(IN), OPTIONAL   :: want_day, short

! Returns the current date as a character string
! e.g.
! want_day     short      string
! .TRUE.      .TRUE.      Thursday, 23 Dec 1999
! .TRUE.      .FALSE.     Thursday, 23 December 1999  <- default
! .FALSE.     .TRUE.      23 Dec 1999
! .FALSE.     .FALSE.     23 December 1999

INTEGER :: val(8), pos
LOGICAL :: want_d, sh
CHARACTER (LEN=9) :: day(0:6) = (/ 'Sunday   ', 'Monday   ', 'Tuesday  ',  &
                      'Wednesday', 'Thursday ', 'Friday   ', 'Saturday ' /)
CHARACTER (LEN=9) :: month(1:12) = (/ 'January  ', 'February ', 'March    ',  &
                         'April    ', 'May      ', 'June     ', 'July     ',  &
                         'August   ', 'September', 'October  ', 'November ',  &
                         'December ' /)

want_d = .TRUE.
IF (PRESENT(want_day)) want_d = want_day
sh = .FALSE.
IF (PRESENT(short)) sh = short

CALL DATE_AND_TIME(VALUES=val)

IF (want_d) THEN
  pos = izlr(val(1), val(2), val(3))
  string = TRIM( day(pos) ) // ','
  pos = LEN_TRIM( string ) + 2
ELSE
  pos = 1
  string = ' '
END IF

WRITE(string(pos:pos+1), '(i2)') val(3)
IF (sh) THEN
  string(pos+3:pos+5) = month(val(2))(1:3)
  pos = pos + 7
ELSE
  string(pos+3:) = month(val(2))
  pos = LEN_TRIM( string ) + 2
END IF

WRITE( string(pos:pos+3), '(i4)') val(1)

RETURN
END SUBROUTINE date_stamp

END MODULE date_sub
!
!
!
!PROGRAM test_datesub
!
!!======DATESUB.FOR with Sample Drivers.
!
!USE date_sub
!IMPLICIT NONE
!INTEGER            :: yyyy, mm, dd, wd, ddd, mma, dda, ndiff, val(8)
!CHARACTER (LEN=50) :: message
!
!! Test date_stamp
!message = ' date_stamp = '
!CALL date_stamp( message(15:) )
!WRITE(*, '(a)') message
!message = ' date_stamp = '
!CALL date_stamp( message(15:), want_day=.FALSE.)
!WRITE(*, '(a)') message
!message = ' date_stamp = '
!CALL date_stamp( message(15:), short=.TRUE.)
!WRITE(*, '(a)') message
!message = ' date_stamp = '
!CALL date_stamp( message(15:), want_day=.FALSE., short=.TRUE.)
!WRITE(*, '(a)') message
!
!!  Is this a leap year? I.e. is 12/31/yyyy the 366th day of the year?
!CALL DATE_AND_TIME(VALUES=val)
!yyyy = val(1)
!
!IF(iday(yyyy, 12, 31) == 366) THEN
!  WRITE(*,*) yyyy, ' is a Leap Year'
!ELSE
!  WRITE(*,*) yyyy, ' is not a Leap Year'
!END IF
!
!!  DAYSUB SHOULD RETURN: 1970, 1, 1, 4, 1
!CALL daysub(jd(1970, 1, 1), yyyy, mm, dd, wd, ddd)
!IF(yyyy /= 1970 .OR. mm /= 1 .OR. dd /= 1 .OR. wd /= 4 .OR. ddd /= 1) THEN
!  WRITE(*,*) 'DAYSUB Failed; YYYY, MM, DD, WD, DDD = ', yyyy, mm, dd, wd, ddd
!  STOP
!END IF
!
!!  DIFFERENCE BETWEEN TO SAME MONTHS AND DAYS OVER 1 LEAP YEAR IS 366.
!ndiff = ndays(5, 22, 1984, 5, 22, 1983)
!IF(ndiff /= 366) THEN
!  WRITE(*,*) 'NDAYS FAILED; NDIFF = ', ndiff
!ELSE
!!  RECOVER MONTH AND DAY FROM YEAR AND DAY NUMBER.
!  CALL calend(yyyy, ddd, mma, dda)
!  IF(mma /= 1 .AND. dda /= 1) THEN
!    WRITE(*,*) 'CALEND FAILED; MMA, DDA = ', mma, dda
!  ELSE
!    WRITE(*,*) '** DATE MANIPULATION SUBROUTINES SIMPLE TEST OK.'
!  END IF
!END IF
!
!STOP
!END PROGRAM test_datesub
