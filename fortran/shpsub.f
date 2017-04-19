!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SHPSUB           SHaPe file SUBroutines
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:08-12-08
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONSISTS OF A SERIES OF SUBROUTINES DESIGNED TO DECODE POLY-
!   LINE AND POLYGON SHAPEFILES.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 12 Aug 2008 (RRD) - initial version 
!
! USAGE:  CALL SHPSUB
!   INPUT ARGUMENT LIST:
!     See individual routines      
!     CALL shphead(buff,file_length)
!     CALL shpdat0(buff,rec_numb,rec_length)
!     CALL shprec5(buff,num_parts,num_points)
!     CALL shploc5(buff,jbyte,num_parts,parts)
!     CALL shpvec5(buff,jbyte,Xpnt,Ypnt)
!     CALL asc2int (x,y)
!     CALL asc2dbl (x,y)
!     CALL flip4 (inp,out)
!     CALL flip8 (inp,out)
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     NONE                                       
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE shphead(buff,file_length)

  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)  :: BUFF(*)
  INTEGER(4),   INTENT(OUT) :: file_length
  
  CHARACTER(1)  :: buff4(4),buff8(8)
  REAL(8)       :: Xmin,Ymin,Xmax,Ymax
  INTEGER(4)    :: file_code,version,shape_type

  LOGICAL(4)    :: big_endian, diag
  COMMON /shpdesc/ big_endian, diag

  CALL asc2int(buff(1:4),file_code)      
  IF(file_code.NE.9994)THEN
     big_endian=.FALSE.
  ELSE
     big_endian = .TRUE.
  END IF

  IF(.NOT.big_endian)THEN
     CALL flip4(buff(25:28),buff4)
  ELSE
     buff4=buff(25:28)
  END IF
  CALL asc2int(buff4,file_length)      
  file_length=file_length*2 ! convert words to bytes

  IF(big_endian)THEN
     CALL flip4(buff(29:32),buff4)
  ELSE
     buff4=buff(29:32)
  END IF
  CALL asc2int(buff4,version) 

  IF(big_endian)THEN
     CALL flip4(buff(33:36),buff4)
  ELSE
     buff4=buff(33:36)
  END IF
  CALL asc2int(buff4,shape_type) 

  IF(diag)THEN
     write(*,*)'Big_endian:  ',big_endian
     write(*,*)'File length: ',file_length
     write(*,*)'Version num: ',version
     write(*,*)'Shape type:  ',shape_type
  END IF

!-----------------------------------

  IF(big_endian)THEN
     CALL flip8(buff(37:44),buff8)
  ELSE
     buff8=buff(37:44)
  END IF
  CALL asc2dbl(buff8,Xmin) 

  IF(big_endian)THEN
     CALL flip8(buff(45:52),buff8)
  ELSE
     buff8=buff(45:52)
  END IF
  CALL asc2dbl(buff8,Ymin) 

  IF(big_endian)THEN
     CALL flip8(buff(53:60),buff8)
  ELSE
     buff8=buff(53:60)
  END IF
  CALL asc2dbl(buff8,Xmax) 

  IF(big_endian)THEN
     CALL flip8(buff(61:68),buff8)
  ELSE
     buff8=buff(61:68)
  END IF
  CALL asc2dbl(buff8,Ymax) 

  IF(diag) write(*,*)'Bounding box: ', Xmin,Xmax,Ymin,Ymax 

END SUBROUTINE shphead 

!###################################################################

SUBROUTINE shpdat0(buff,rec_numb,rec_length)

  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)  :: buff(*)
  INTEGER(4),   INTENT(OUT) :: rec_numb
  INTEGER(4),   INTENT(OUT) :: rec_length
  
  CHARACTER(1)  :: buff4(4)

  LOGICAL(4)    :: big_endian, diag
  COMMON /shpdesc/ big_endian, diag

  IF(.NOT.big_endian)THEN
     CALL flip4(buff(1:4),buff4)
  ELSE
     buff4=buff(1:4)
  END IF
  CALL asc2int(buff4,rec_numb)      

  IF(.NOT.big_endian)THEN
     CALL flip4(buff(5:8),buff4)
  ELSE
     buff4=buff(5:8)
  END IF
  CALL asc2int(buff4,rec_length) 
  rec_length=rec_length*2 ! convert words to bytes

  IF(diag)THEN
     write(*,*)'Record number: ',rec_numb
     write(*,*)'Record length: ',rec_length
  END IF

END SUBROUTINE shpdat0 

!###################################################################

SUBROUTINE shprec5(buff,num_parts,num_points)

  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)  :: buff(*)
  INTEGER,      INTENT(OUT) :: num_parts
  INTEGER,      INTENT(OUT) :: num_points
  
  CHARACTER(1)  :: buff4(4)
  INTEGER(4)    :: rec_type

  LOGICAL(4)    :: big_endian, diag
  COMMON /shpdesc/ big_endian, diag

  IF(big_endian)THEN
     CALL flip4(buff(1:4),buff4)
  ELSE
     buff4=buff(1:4)
  END IF
  CALL asc2int(buff4,rec_type)      

  IF(big_endian)THEN
     CALL flip4(buff(37:40),buff4)
  ELSE
     buff4=buff(37:40)
  END IF
  CALL asc2int(buff4,num_parts) 

  IF(big_endian)THEN
     CALL flip4(buff(41:44),buff4)
  ELSE
     buff4=buff(41:44)
  END IF
  CALL asc2int(buff4,num_points) 

  IF(diag)THEN
     write(*,*)'Shape type:    ',rec_type
     write(*,*)'Number parts:  ',num_parts
     write(*,*)'Number points: ',num_points
  END IF

END SUBROUTINE shprec5 

!###################################################################

SUBROUTINE shploc5(buff,jbyte,num_parts,parts)

  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)    :: buff(*)
  INTEGER(4),   INTENT(INOUT) :: jbyte     
  INTEGER(4),   INTENT(IN)    :: num_parts
  INTEGER(4),   INTENT(OUT)   :: parts(*)
  
  CHARACTER(1)  :: buff4(4)
  INTEGER(4)    :: k

  LOGICAL(4)    :: big_endian, diag
  COMMON /shpdesc/ big_endian, diag

  DO k=1,num_parts
     IF(big_endian)THEN
        CALL flip4(buff(jbyte:jbyte+3),buff4)
     ELSE
        buff4=buff(jbyte:jbyte+3)
     END IF
     CALL asc2int(buff4,parts(k))      
     IF(diag) write(*,*)'Part (num,loc): ',k,parts(k)
     jbyte=jbyte+4
  END DO

END SUBROUTINE shploc5 

!###################################################################

SUBROUTINE shpvec5(buff,jbyte,Xpnt,Ypnt)

  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)    :: buff(*)
  INTEGER,      INTENT(INOUT) :: jbyte     
  REAL(8),      INTENT(OUT)   :: Xpnt,Ypnt
  
  CHARACTER(1)  :: buff8(8)
  INTEGER(4)    :: k

  LOGICAL(4)    :: big_endian, diag
  COMMON /shpdesc/ big_endian, diag

  IF(big_endian)THEN
     CALL flip8(buff(jbyte:jbyte+7),buff8)
  ELSE
     buff8=buff(jbyte:jbyte+7)
  END IF
  CALL asc2dbl(buff8,Xpnt)      
  jbyte=jbyte+8

  IF(big_endian)THEN
     CALL flip8(buff(jbyte:jbyte+7),buff8)
  ELSE
     buff8=buff(jbyte:jbyte+7)
  END IF
  CALL asc2dbl(buff8,Ypnt)      
  jbyte=jbyte+8

END SUBROUTINE shpvec5 

!###################################################################

subroutine asc2int (x,y)
   integer(4), intent(in)  :: x
   integer(4), intent(out) :: y
   y=x
end subroutine asc2int

subroutine asc2dbl (x,y)
   real(8), intent(in)  :: x
   real(8), intent(out) :: y
   y=x
end subroutine asc2dbl

subroutine flip4 (inp,out)
   integer      :: k
   character(4) :: inp,out
   DO k=1,4  
      out(5-k:5-k)=inp(k:k)
   END DO
end subroutine flip4

subroutine flip8 (inp,out)
   integer      :: k
   character(8) :: inp,out
   DO k=1,8  
      out(9-k:9-k)=inp(k:k)
   END DO
end subroutine flip8
