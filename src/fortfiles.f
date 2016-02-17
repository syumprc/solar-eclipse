C
C     File: fortfiles.f
C     Purpose: Subroutines relating to opening or closing fortran units
C              from C
C
      SUBROUTINE OPENOUT (STATUS)
C
C     Set up units for PREPED
C     (Code partly from the now obsolete INPUT.F, before call to PREPED)
C
      INTEGER STATUS
C
 10   FORMAT (A)

C      OPEN (1, STATUS='SCRATCH', FORM='UNFORMATTED', ERR=910)
C      OPEN (2, STATUS='SCRATCH', FORM='UNFORMATTED', ERR=920)
      CALL RESETALLSC
      OPEN (3, FILE='solar.hst', STATUS='UNKNOWN', ERR=930)
      OPEN (4, FILE='solar.smp', STATUS='UNKNOWN', ERR=940)
C
C Successful completion
C
      STATUS = 0
      RETURN
C
C Errors.  Return status code and let C handle them
C
 930  STATUS = 3
      RETURN
 940  STATUS = 4
      CLOSE (3)
      RETURN
C
      END

C OPENPHEN opens phenotype file in unit 7 and skips to old part
C This used to be part of preprep, but now has to be done earlier
C because of directory change
C
      SUBROUTINE OPENPHEN (PHENFILE, STATUS)
      CHARACTER*(*) PHENFILE
      INTEGER STATUS
      CHARACTER*256 LINE

      OPEN (7, FILE=PHENFILE, STATUS='OLD', ERR=910)
C
C Now we must skip over new header sections.  Skip until we find the
C !Data mark
C
 10   FORMAT (A)
 100  READ (7, 10, ERR=920) LINE
      IF (LINE .NE. '!Data') GOTO 100
C
C Successful completion
C
      STATUS = 0
      RETURN
C
C Errors.      
C
 910  STATUS = 1
      RETURN
 920  STATUS = 2
      CLOSE (7)
      RETURN
      END

C CLOSEPHEN closes phenotype file, or whatever happens to be open on unit 2

      SUBROUTINE CLOSEPHEN
      CLOSE (7)
      RETURN
      END
C
C     Close output units 1, 2, 3, 4
C
      SUBROUTINE CLOSEOUT

      CALL RESETALLSC
      CLOSE (3)
      CLOSE (4)
      RETURN
      END
