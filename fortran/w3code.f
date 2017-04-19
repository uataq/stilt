!$$$  SUBPROGRAM DOCUMENTATION  BLOCK
! SUBPROGRAM:  W3FI63        UNPACK GRIB FIELD TO GRIB GRID
!   PRGMMR: DRAXLER          ORG: R/ARL       DATE:99-02-19
!
! ABSTRACT: UNPACK A GRIB (EDITION 1) FIELD TO THE EXACT GRID
!   SPECIFIED IN THE GRIB MESSAGE, ISOLATE THE BIT MAP, AND MAKE
!   THE VALUES OF THE PRODUCT DESCRIPTON SECTION (PDS) AND THE
!   GRID DESCRIPTION SECTION (GDS) AVAILABLE IN RETURN ARRAYS.
!
!   WHEN DECODING IS COMPLETED, DATA AT EACH GRID POINT HAS BEEN
!          RETURNED IN THE UNITS SPECIFIED IN THE GRIB MANUAL.
!
! PROGRAM HISTORY LOG:
!   19 Feb 1999 (RRD) - initial version
!   07 Mar 2001 (RRD) - added additional error messages
!   12 Feb 2002 (RRD) - ecmwf compatibility 
!   18 Jun 2002 (RRD) - time range indicator correction
!   01 Oct 2002 (RRD) - support for nonhydrostatic meso model
!   25 Feb 2003 (RRD) - replaced ichar function
!   04 Apr 2003 (RRD) - reversed lat-lon direction increment
!
! USAGE:    CALL W3FI63(BUFF,KPDS,KGDS,KBMS,RVAL,KPTR,KRET)
!   INPUT ARGUMENT LIST:
!     BUFF     - GRIB FIELD - "GRIB" THRU "7777"   CHAR*1
!
!   OUTPUT ARGUMENT LIST:
!     RVAL     - ARRAY CONTAINING DATA ELEMENTS
!     KPDS     - ARRAY CONTAINING PDS ELEMENTS.  (EDITION 1)
!          (1)   - ID OF CENTER
!          (2)   - GENERATING PROCESS ID NUMBER
!          (3)   - GRID DEFINITION
!          (4)   - GDS/BMS FLAG (RIGHT ADJ COPY OF OCTET 8)
!          (5)   - INDICATOR OF PARAMETER
!          (6)   - TYPE OF LEVEL
!          (7)   - HEIGHT/PRESSURE , ETC OF LEVEL
!          (8)   - YEAR INCLUDING (CENTURY-1)
!          (9)   - MONTH OF YEAR
!          (10)  - DAY OF MONTH
!          (11)  - HOUR OF DAY
!          (12)  - MINUTE OF HOUR
!          (13)  - INDICATOR OF FORECAST TIME UNIT
!          (14)  - TIME RANGE 1
!          (15)  - TIME RANGE 2
!          (16)  - TIME RANGE FLAG
!          (17)  - NUMBER INCLUDED IN AVERAGE
!          (18)  - VERSION NR OF GRIB SPECIFICATION
!          (19)  - VERSION NR OF PARAMETER TABLE
!          (20)  - NR MISSING FROM AVERAGE/ACCUMULATION
!          (21)  - CENTURY OF REFERENCE TIME OF DATA
!          (22)  - UNITS DECIMAL SCALE FACTOR
!          (23)  - SUBCENTER NUMBER
!          (24)  - PDS BYTE 29, FOR NMC ENSEMBLE PRODUCTS
!                  128 IF FORECAST FIELD ERROR
!                   64 IF BIAS CORRECTED FCST FIELD
!                   32 IF SMOOTHED FIELD
!                  WARNING: CAN BE COMBINATION OF MORE THAN 1
!          (25)  - PDS BYTE 30, NOT USED
!       (26-35)  - RESERVED
!       (36-N)   - CONSECUTIVE BYTES EXTRACTED FROM PROGRAM
!                  DEFINITION SECTION (PDS) OF GRIB MESSAGE
!     KGDS     - ARRAY CONTAINING GDS ELEMENTS.
!          (1)   - DATA REPRESENTATION TYPE
!          (19)  - NUMBER OF VERTICAL COORDINATE PARAMETERS
!          (20)  - OCTET NUMBER OF THE LIST OF VERTICAL COORDINATE
!                  PARAMETERS
!                  OR
!                  OCTET NUMBER OF THE LIST OF NUMBERS OF POINTS
!                  IN EACH ROW
!                  OR
!                  255 IF NEITHER ARE PRESENT
!          (21)  - FOR GRIDS WITH PL, NUMBER OF POINTS IN GRID
!          (22)  - NUMBER OF WORDS IN EACH ROW
!       LATITUDE/LONGITUDE GRIDS
!          (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
!          (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
!          (4)   - LA(1) LATITUDE OF ORIGIN
!          (5)   - LO(1) LONGITUDE OF ORIGIN
!          (6)   - RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)
!          (7)   - LA(2) LATITUDE OF EXTREME POINT
!          (8)   - LO(2) LONGITUDE OF EXTREME POINT
!          (9)   - DI LONGITUDINAL DIRECTION INCREMENT
!          (10)  - DJ LATITUDINAL DIRECTION OF INCREMENT
!          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
!       POLAR STEREOGRAPHIC GRIDS
!          (2)   - N(I) NR POINTS ALONG LAT CIRCLE
!          (3)   - N(J) NR POINTS ALONG LON CIRCLE
!          (4)   - LA(1) LATITUDE OF ORIGIN
!          (5)   - LO(1) LONGITUDE OF ORIGIN
!          (6)   - RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
!          (7)   - LOV GRID ORIENTATION
!          (8)   - DX - X DIRECTION INCREMENT
!          (9)   - DY - Y DIRECTION INCREMENT
!          (10)  - PROJECTION CENTER FLAG
!          (11)  - SCANNING MODE (RIGHT ADJ COPY OF OCTET 28)
!       MERCATOR GRIDS
!          (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
!          (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
!          (4)   - LA(1) LATITUDE OF ORIGIN
!          (5)   - LO(1) LONGITUDE OF ORIGIN
!          (6)   - RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)
!          (7)   - LA(2) LATITUDE OF LAST GRID POINT
!          (8)   - LO(2) LONGITUDE OF LAST GRID POINT
!          (9)   - LATIT - LATITUDE OF PROJECTION INTERSECTION
!          (10)  - RESERVED
!          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
!          (12)  - LONGITUDINAL DIR GRID LENGTH
!          (13)  - LATITUDINAL DIR GRID LENGTH
!       LAMBERT CONFORMAL GRIDS
!          (2)   - NX NR POINTS ALONG X-AXIS
!          (3)   - NY NR POINTS ALONG Y-AXIS
!          (4)   - LA1 LAT OF ORIGIN (LOWER LEFT)
!          (5)   - LO1 LON OF ORIGIN (LOWER LEFT)
!          (6)   - RESOLUTION (RIGHT ADJ COPY OF OCTET 17)
!          (7)   - LOV - ORIENTATION OF GRID
!          (8)   - DX - X-DIR INCREMENT
!          (9)   - DY - Y-DIR INCREMENT
!          (10)  - PROJECTION CENTER FLAG
!          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
!          (12)  - LATIN 1 - FIRST LAT FROM POLE OF SECANT CONE INTER
!          (13)  - LATIN 2 - SECOND LAT FROM POLE OF SECANT CONE INTER
!     KBMS       - BITMAP DESCRIBING LOCATION OF OUTPUT ELEMENTS.
!                            (ALWAYS CONSTRUCTED)
!     KPTR       - 20 WORD ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!          (1)   - TOTAL LENGTH OF GRIB MESSAGE
!          (2)   - LENGTH OF INDICATOR (SECTION  0)
!          (3)   - LENGTH OF PDS       (SECTION  1)
!          (4)   - LENGTH OF GDS       (SECTION  2)
!          (5)   - LENGTH OF BMS       (SECTION  3)
!          (6)   - LENGTH OF BDS       (SECTION  4)
!          (7)   - VALUE OF CURRENT BYTE
!          (8)   - BIT POINTER
!          (9)   - GRIB START BIT NR
!         (10)   - GRIB/GRID ELEMENT COUNT
!         (11)   - NR UNUSED BITS AT END OF SECTION 3
!         (12)   - BIT MAP FLAG (COPY OF BMS OCTETS 5,6)
!         (13)   - NR UNUSED BITS AT END OF SECTION 2
!         (14)   - BDS FLAGS (RIGHT ADJ COPY OF OCTET 4)
!         (15)   - NR UNUSED BITS AT END OF SECTION 4
!         (16)   - RESERVED
!         (17)   - RESERVED
!         (18)   - RESERVED
!         (19)   - RESERVED
!         (20)   - RESERVED
!     KRET       - FLAG INDICATING QUALITY OF COMPLETION
!
! REMARKS: WHEN DECODING IS COMPLETED, DATA AT EACH GRID POINT HAS BEEN
!          RETURNED IN THE UNITS SPECIFIED IN THE GRIB MANUAL.
!
!          VALUES FOR RETURN FLAG (KRET)
!     KRET = 0 - NORMAL RETURN, NO ERRORS
!          = 1 - 'GRIB' NOT FOUND IN FIRST 100 CHARS
!          = 2 - '7777' NOT IN CORRECT LOCATION
!          = 3 - UNPACKED FIELD IS LARGER THAN 260000
!          = 4 - GDS/ GRID NOT ONE OF CURRENTLY ACCEPTED VALUES
!          = 5 - GRID NOT CURRENTLY AVAIL FOR CENTER INDICATED
!          = 8 - TEMP GDS INDICATED, BUT GDS FLAG IS OFF
!          = 9 - GDS INDICATES SIZE MISMATCH WITH STD GRID
!          =10 - INCORRECT CENTER INDICATOR
!          =11 - BINARY DATA SECTION (BDS) NOT COMPLETELY PROCESSED.
!                PROGRAM IS NOT SET TO PROCESS FLAG COMBINATIONS
!                SHOWN IN OCTETS 4 AND 14.
!          =12 - BINARY DATA SECTION (BDS) NOT COMPLETELY PROCESSED.
!                PROGRAM IS NOT SET TO PROCESS FLAG COMBINATIONS
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  GENERIC
!
!$$$

      SUBROUTINE W3FI63(BUFF,KPDS,KGDS,KBMS,RVAR,KPTR,KRET)

      CHARACTER(1), INTENT(IN)  :: BUFF(:)
      LOGICAL(1),   INTENT(OUT) :: KBMS(:)
      INTEGER,      INTENT(OUT) :: KPDS(:)
      INTEGER,      INTENT(OUT) :: KGDS(:)
      REAL,         INTENT(OUT) :: RVAR(:)
      INTEGER,      INTENT(OUT) :: KPTR(:)
      INTEGER,      INTENT(OUT) :: KRET

      INTERFACE
         SUBROUTINE SECT0(buff,klen,kptr,kret)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         INTEGER,      INTENT(OUT) :: KLEN
         INTEGER,      INTENT(OUT) :: KPTR(:)
         INTEGER,      INTENT(OUT) :: KRET
         END SUBROUTINE sect0

         SUBROUTINE SECT1(buff,klen,kpds,kptr)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         INTEGER,      INTENT(OUT) :: KLEN
         INTEGER,      INTENT(OUT) :: KPDS(:)
         INTEGER,      INTENT(OUT) :: KPTR(:)
         END SUBROUTINE sect1

         SUBROUTINE SECT2(buff,klen,kgds,kptr,kret)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         INTEGER,      INTENT(OUT) :: KLEN
         INTEGER,      INTENT(OUT) :: KGDS(:)
         INTEGER,      INTENT(OUT) :: KPTR(:)
         INTEGER,      INTENT(OUT) :: KRET    
         END SUBROUTINE sect2

         SUBROUTINE SECT3(buff,klen,kptr)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         INTEGER,      INTENT(OUT) :: KLEN
         INTEGER,      INTENT(OUT) :: KPTR(:)
         END SUBROUTINE sect3

         SUBROUTINE SECT4(buff,klen,kptr,mfact,refval,nbpv,kret)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         INTEGER,      INTENT(OUT) :: KLEN
         INTEGER,      INTENT(OUT) :: KPTR(:)
         INTEGER,      INTENT(OUT) :: MFACT
         REAL,         INTENT(OUT) :: REFVAL
         INTEGER,      INTENT(OUT) :: NBPV
         INTEGER,      INTENT(OUT) :: KRET
         END SUBROUTINE sect4

         SUBROUTINE SECT5(buff,klen,kret)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         INTEGER,      INTENT(OUT) :: KLEN
         INTEGER,      INTENT(OUT) :: KRET
         END SUBROUTINE sect5

         SUBROUTINE DECODE(buff,rvar,nxyp,kfact,mfact,refval,nbpv)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         REAL,         INTENT(OUT) :: RVAR(:)
         INTEGER,      INTENT(IN)  :: NXYP
         INTEGER,      INTENT(IN)  :: KFACT
         INTEGER,      INTENT(IN)  :: MFACT
         REAL,         INTENT(IN)  :: REFVAL
         INTEGER,      INTENT(IN)  :: NBPV
         END SUBROUTINE decode
      END INTERFACE

      kbyte=1

!---->indicator "section 0"
      CALL SECT0(buff(kbyte:),klen,kptr,kret)
      kbyte=kbyte+klen
      if(kret.ne.0)return

!---->product definition "section 1"
      CALL SECT1(buff(kbyte:),klen,kpds,kptr)
      kbyte=kbyte+klen

!---->grid description "section 2"
      if(kpds(4).eq.128.or.kpds(4).eq.192)then
         CALL SECT2(buff(kbyte:),klen,kgds,kptr,kret)
         if(kret.ne.0)return
         kbyte=kbyte+klen
      end if

!---->bit map "section 3"
      if(kpds(4).eq.64.or.kpds(4).eq.192)then
         CALL SECT3(buff(kbyte:),klen,kptr)
         kbyte=kbyte+klen
         write(*,*)'Bit map decoding not supported: ', klen
         kret=99
         kbms(1)=.true.  
         return
      else
         kbms(1)=.false.
      end if

!---->binary data "section 4"
      CALL SECT4(buff(kbyte:),klen,kptr,mfact,refval,nbpv,kret)
      if(kret.ne.0)return
      kfact=kpds(22)
      nxyp=kgds(2)*kgds(3)
      CALL DECODE(buff(kbyte+11:),rvar,nxyp,kfact,mfact,refval,nbpv)
      kbyte=kbyte+klen

!---->grib end "section 5"
      CALL SECT5(buff(kbyte:),klen,kret)
      kbyte=kbyte+klen
      if(kret.ne.0)return

      END SUBROUTINE w3fi63

!-----------------------------------------------------------------
! light version of w3fi63 only decodes pds and gds

      SUBROUTINE W3FI00(BUFF,KPDS,KGDS,KBMS,RVAR,KPTR,KRET)

      CHARACTER(1), INTENT(IN)  :: BUFF(:)
      LOGICAL(1),   INTENT(OUT) :: KBMS(:)
      INTEGER,      INTENT(OUT) :: KPDS(:)
      INTEGER,      INTENT(OUT) :: KGDS(:)
      REAL,         INTENT(OUT) :: RVAR(:)
      INTEGER,      INTENT(OUT) :: KPTR(:)
      INTEGER,      INTENT(OUT) :: KRET

      INTERFACE
         SUBROUTINE SECT0(buff,klen,kptr,kret)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         INTEGER,      INTENT(OUT) :: KLEN
         INTEGER,      INTENT(OUT) :: KPTR(:)
         INTEGER,      INTENT(OUT) :: KRET
         END SUBROUTINE sect0

         SUBROUTINE SECT1(buff,klen,kpds,kptr)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         INTEGER,      INTENT(OUT) :: KLEN
         INTEGER,      INTENT(OUT) :: KPDS(:)
         INTEGER,      INTENT(OUT) :: KPTR(:)
         END SUBROUTINE sect1

         SUBROUTINE SECT2(buff,klen,kgds,kptr,kret)
         CHARACTER(1), INTENT(IN)  :: BUFF(:)
         INTEGER,      INTENT(OUT) :: KLEN
         INTEGER,      INTENT(OUT) :: KGDS(:)
         INTEGER,      INTENT(OUT) :: KPTR(:)
         INTEGER,      INTENT(OUT) :: KRET    
         END SUBROUTINE sect2
      END INTERFACE

      kbyte=1

!---->indicator "section 0"
      CALL SECT0(buff(kbyte:),klen,kptr,kret)
      kbyte=kbyte+klen
      if(kret.ne.0)return

!---->product definition "section 1"
      CALL SECT1(buff(kbyte:),klen,kpds,kptr)
      kbyte=kbyte+klen

!---->grid description "section 2"
      if(kpds(4).eq.128.or.kpds(4).eq.192)then
         CALL SECT2(buff(kbyte:),klen,kgds,kptr,kret)
         if(kret.ne.0)return
         kbyte=kbyte+klen
      end if

!---->bit map "section 3"
      kbms(1)=.false.

!---->binary data "section 4"
      rvar(1)=0.0

      END SUBROUTINE w3fi00

!-----------------------------------------------------------------

      SUBROUTINE SECT0(buff,klen,kptr,kret)
      CHARACTER(1), INTENT(IN)  :: BUFF(:)
      INTEGER,      INTENT(OUT) :: KLEN
      INTEGER,      INTENT(OUT) :: KPTR(:)
      INTEGER,      INTENT(OUT) :: KRET

!     ------------------------------------------------------------
!     Use the jchar function in all subs to replace of ichar
!     due to undefined nature of ichar for values greater than 127
!     in some Fortran compilers.  Other solutions available ...
!     see pakinp.f in ../source directory

      character*1 mychr
      jchar(mychr)=iand(ichar(mychr),255)
!     ------------------------------------------------------------

      IF(buff(1)//buff(2)//buff(3)//buff(4).NE.'GRIB')THEN
         WRITE(*,*)'ERROR w3lib: record start not GRIB'
         KRET=1
         RETURN
      END IF
      KPTR(1)=jchar(buff(7))+jchar(buff(6))*256+jchar(buff(5))*65536
      klen=8
      KPTR(2)=klen
      KRET=0
      END SUBROUTINE sect0

!-----------------------------------------------------------------

      SUBROUTINE SECT1(buff,klen,kpds,kptr)
      CHARACTER(1), INTENT(IN)  :: BUFF(:)
      INTEGER,      INTENT(OUT) :: KLEN
      INTEGER,      INTENT(OUT) :: KPDS(:)
      INTEGER,      INTENT(OUT) :: KPTR(:)

!     ------------------------------------------------------------
      character*1 mychr
      jchar(mychr)=iand(ichar(mychr),255)
!     ------------------------------------------------------------

      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      kptr(3)=klen

      kpds(1)=jchar(buff(5))
      kpds(2)=jchar(buff(6))
      kpds(3)=jchar(buff(7))
      kpds(4)=jchar(buff(8))
      kpds(5)=jchar(buff(9))
      kpds(6)=jchar(buff(10))
      kpds(7)=jchar(buff(12))+jchar(buff(11))*256
      kpds(8)=jchar(buff(13))
      kpds(9)=jchar(buff(14))
      kpds(10)=jchar(buff(15))
      kpds(11)=jchar(buff(16))
      kpds(12)=jchar(buff(17))
      kpds(13)=jchar(buff(18))

      IF(jchar(buff(21)).EQ.10)THEN
         kpds(15)=jchar(buff(20))+jchar(buff(19))*256
         kpds(14)=kpds(15)
      ELSEIF(jchar(buff(18)).GE.8)THEN
         kpds(14)=jchar(buff(19))*jchar(buff(18))
         kpds(15)=jchar(buff(20))*jchar(buff(18))
      ELSE
         kpds(14)=jchar(buff(19))
         kpds(15)=jchar(buff(20))
      END IF

      kpds(16)=jchar(buff(21))
      kpds(17)=jchar(buff(23))+jchar(buff(22))*256
      kpds(18)= 1
      kpds(19)=jchar(buff(4))
      kpds(20)=jchar(buff(24))
      kpds(21)=jchar(buff(25))


!     decimal scale factor
      kfact=jchar(buff(28))+jchar(buff(27))*256
      if(iand(kfact,32768).ne.0)kfact=-iand(kfact,32767)
      kpds(22)=kfact

      kpds(23)=jchar(buff(26))

      END SUBROUTINE sect1

!-----------------------------------------------------------------

      SUBROUTINE SECT2(buff,klen,kgds,kptr,kret)
      CHARACTER(1), INTENT(IN)  :: BUFF(:)
      INTEGER,      INTENT(OUT) :: KLEN
      INTEGER,      INTENT(OUT) :: KGDS(:)
      INTEGER,      INTENT(OUT) :: KPTR(:)
      INTEGER,      INTENT(OUT) :: KRET    

      COMMON /VCOORD/ LEVELS, OFFSET(100), SIGMA(100)

!     ------------------------------------------------------------
      character*1 mychr
      jchar(mychr)=iand(ichar(mychr),255)
!     ------------------------------------------------------------

      kret=0  ! assume normal return unless otherwise 

      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      kptr(4)=klen

      kgds(1)=jchar(buff(6))
      kgds(2)=jchar(buff(8))+jchar(buff(7))*256
      kgds(3)=jchar(buff(10))+jchar(buff(9))*256

      lat1=jchar(buff(13))+jchar(buff(12))*256+jchar(buff(11))*65536
      if(iand(lat1,8388608).ne.0)lat1=-iand(lat1,8388607)
      lon1=jchar(buff(16))+jchar(buff(15))*256+jchar(buff(14))*65536
      if(iand(lon1,8388608).ne.0)lon1=-iand(lon1,8388607)
      kgds(4)=lat1
      kgds(5)=lon1

      kgds(6)=jchar(buff(17))

      lov=jchar(buff(20))+jchar(buff(19))*256+jchar(buff(18))*65536
      if(iand(lov,8388608).ne.0)lov=-iand(lov,8388607)
      kgds(7)=lov

      mxds=jchar(buff(23))+jchar(buff(22))*256+jchar(buff(21))*65536
      if(iand(mxds,8388608).ne.0)mxds=-iand(mxds,8388607)
      kgds(8)=mxds

      myds=jchar(buff(26))+jchar(buff(25))*256+jchar(buff(24))*65536

      if(kgds(1).eq.0)then
!        lat/lon grid
         kgds(9) =jchar(buff(25))+jchar(buff(24))*256
         kgds(10)=jchar(buff(27))+jchar(buff(26))*256

      elseif(kgds(1).eq.1)then
!        mercator projection
         kgds(9)=myds
         kgds(10)=jchar(buff(27))
         lon=jchar(buff(31))+jchar(buff(30))*256+jchar(buff(29))*65536
         if(iand(lon,8388608).ne.0)lon=-iand(lon,8388607)
         lat=jchar(buff(34))+jchar(buff(33))*256+jchar(buff(32))*65536
         if(iand(lat,8388608).ne.0)lat=-iand(lat,8388607)
         kgds(12)=lon
         kgds(13)=lat

      elseif(kgds(1).eq.3)then
!        lambert projection
         kgds(9)=myds
         kgds(10)=jchar(buff(27))
         lt1=jchar(buff(31))+jchar(buff(30))*256+jchar(buff(29))*65536
         if(iand(lt1,8388608).ne.0)lt1=-iand(lt1,8388607)
         lt2=jchar(buff(34))+jchar(buff(33))*256+jchar(buff(32))*65536
         if(iand(lt2,8388608).ne.0)lt2=-iand(lt2,8388607)
         kgds(12)=lt1
         kgds(13)=lt2

      elseif(kgds(1).eq.5)then
!        polar sterographic
         kgds(9)=myds
         kgds(10)=jchar(buff(27))

      elseif(kgds(1).eq.203)then
!        non-hydrostatic mesoscale model
         kgds(9) =jchar(buff(27))+jchar(buff(26))*256
         kgds(10)=jchar(buff(25))+jchar(buff(24))*256

      else
!        grid not currently accepted
         kret=4
         return
      end if

      kgds(11)=jchar(buff(28))
      kgds(19)=jchar(buff(4))
      kgds(20)=jchar(buff(5))

!     ecmwf vertical coordinate parameters
      
      if(kgds(19).eq.0)return
      levels=kgds(19)/2

      do kv=1,kgds(19)
         k=kgds(20)+(kv-1)*4
         mexp=iand(jchar(buff(k)),127)
         mant=jchar(buff(k+3))+jchar(buff(k+2))*256+jchar(buff(k+1))*65536
         rvar=16.0**(mexp-64)*(mant/16777216.0)
         if(iand(jchar(buff(k)),128).ne.0)rvar=-rvar
         kk=mod(kv-1,levels)+1
         if(kv.le.levels)then
            offset(kk)=rvar/100.0
         else
            sigma(kk)=rvar 
         end if
      end do

      END SUBROUTINE sect2

!-----------------------------------------------------------------

      SUBROUTINE SECT3(buff,klen,kptr)
      CHARACTER(1), INTENT(IN)  :: BUFF(:)
      INTEGER,      INTENT(OUT) :: KLEN
      INTEGER,      INTENT(OUT) :: KPTR(:)

!     ------------------------------------------------------------
      character*1 mychr
      jchar(mychr)=iand(ichar(mychr),255)
!     ------------------------------------------------------------

      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      kptr(5)=klen
      END SUBROUTINE sect3

!-----------------------------------------------------------------

      SUBROUTINE SECT4(buff,klen,kptr,mfact,refval,nbpv,kret)
      CHARACTER(1), INTENT(IN)  :: BUFF(:)
      INTEGER,      INTENT(OUT) :: KLEN
      INTEGER,      INTENT(OUT) :: KPTR(:)
      INTEGER,      INTENT(OUT) :: MFACT
      REAL,         INTENT(OUT) :: REFVAL
      INTEGER,      INTENT(OUT) :: NBPV
      INTEGER,      INTENT(OUT) :: KRET

!     ------------------------------------------------------------
      character*1 mychr
      jchar(mychr)=iand(ichar(mychr),255)
!     ------------------------------------------------------------

      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      kptr(6)=klen

      kpack=jchar(buff(4))
      if(iand(kpack,128).ne.0.or.                                              &
         iand(kpack, 64).ne.0.or.                                              &
         iand(kpack, 16).ne.0)then
         WRITE(*,*)'ERROR w3lib: unsupported packing format' 
         kret=11
         return
      end if

!     binary scale factor
      mfact=jchar(buff(6))+jchar(buff(5))*256
      if(iand(mfact,32768).ne.0)mfact=-iand(mfact,32767)

!     floating point value consists of 7 bit exponent and 24 bit mantissa
      mexp=iand(jchar(buff(7)),127)
      mant=jchar(buff(10))+jchar(buff(9))*256+jchar(buff(8))*65536

!     value equals 2^(-24)*MANT*16^(MEXP-64)
      refval=16.0**(mexp-64)*(mant/16777216.0)
      if(iand(jchar(buff(7)),128).ne.0)refval=-refval

!     number of bits per variable
      nbpv=jchar(buff(11))

      kret=0

      END SUBROUTINE sect4

!-----------------------------------------------------------------

      SUBROUTINE SECT5(buff,klen,kret)
      CHARACTER(1), INTENT(IN)  :: BUFF(:)
      INTEGER,      INTENT(OUT) :: KLEN
      INTEGER,      INTENT(OUT) :: KRET

      klen=4
      IF(buff(1)//buff(2)//buff(3)//buff(4).NE.'7777')THEN
         WRITE(*,*)'-->',buff(1)//buff(2)//buff(3)//buff(4),'<--'
         WRITE(*,*)'ERROR w3lib: record termination not 7777'
         KRET=2
         RETURN
      END IF
      kret=0

      END SUBROUTINE sect5

!-----------------------------------------------------------------

      SUBROUTINE DECODE(buff,rvar,nxyp,kfact,mfact,refval,nbpv)
      CHARACTER(1), INTENT(IN)  :: BUFF(:)
      REAL,         INTENT(OUT) :: RVAR(:)
      INTEGER,      INTENT(IN)  :: NXYP
      INTEGER,      INTENT(IN)  :: KFACT
      INTEGER,      INTENT(IN)  :: MFACT
      REAL,         INTENT(IN)  :: REFVAL
      INTEGER,      INTENT(IN)  :: NBPV

!     ------------------------------------------------------------
      character*1 mychr
      jchar(mychr)=iand(ichar(mychr),255)
!     ------------------------------------------------------------

!     go through each element of the real array
      do kv=1,nxyp

!        start with bit sum as zero
         rsum=0.0

!        constant field packed number of bits per variable nbpv=0
         if(nbpv.gt.0)then

!           integer sum is converted to real sum after exponent scaling
            ksum=0

!           compute starting (1) and ending (2) bit position
            kbit1=(kv-1)*nbpv+1
            kbit2=kv*nbpv

!           loop through each bit left (high) to right (low)
            do kb=kbit1,kbit2


!              compute which byte contains the bit "kb"
               kbyte=(kb-1)/8+1

!              determine the relative bit position for "kb" in kbyte
               mbit=kb-(kbyte-1)*8

!              leftmost bit starts with value of zero and then
!              increments by power of two for each shift to the right
               ksum=2*ksum

!              add that bit value to sum if mbit in kbyte is on
               kval=jchar(buff(kbyte))
               if(iand(kval,2**(8-mbit)).ne.0)ksum=ksum+1
!              if(iand(jchar(buff(kbyte)),2**(8-mbit)).ne.0)ksum=ksum+1

            end do

!           compute real data point adjusted for exponential scaling
!           using the binary scale factor from BDS section octets 5-6
            if(ksum.gt.0) rsum = float(ksum) * 2.0**mfact

         end if

!        real value
         rvar(kv)= (refval+rsum) / 10.0**kfact

      end do

      END SUBROUTINE decode
