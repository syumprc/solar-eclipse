C
C     File: writestr.f
C     Purpose: Subroutine to write a string on a FORTRAN unit
C              (For writing from C or C++ code.)
C
      SUBROUTINE WRITESTR (UNITNO, STRING)
      INTEGER UNITNO
      CHARACTER*(*) STRING

 1    FORMAT (A)
      WRITE (UNITNO, 1) STRING
      END
