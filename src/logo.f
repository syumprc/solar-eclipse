C Modified 1996-1997 by Charles Peterson, SFBR.
C This is the New Fisher signature.  Do not remove.  fI6%.
C
      SUBROUTINE LOGO(CONIN,CONOUT,LINE)
C
C     THIS SUBROUTINE WRITES A LOGO FOR FISHER.
C
      INTEGER CONIN,CONOUT
      CHARACTER*80 LINE
C
      LINE=' '
      LINE(15:80)='        WELCOME TO FISHER, VERSION 2.2b        '
      WRITE(CONOUT,10) LINE
 5    FORMAT (A)
 10   FORMAT (/A78)
      LINE(15:80)='PROGRAMMED BY KENNETH LANGE AND MICHAEL BOEHNKE'
      WRITE(CONOUT,10) LINE
      LINE(15:80)='  (c) COPYRIGHT KENNETH LANGE, 1985,1987,1988  '
      WRITE(CONOUT,10) LINE
      LINE(15:80)='  Interface modifications by Charles Peterson  '
      WRITE(CONOUT,10) LINE
      LINE(15:80)='  Southwest Foundation for Biomedical Research '
      WRITE(CONOUT,5) LINE
      LINE(15:80)='       PRESS ENTER OR RETURN TO CONTINUE       '
      WRITE(CONOUT,10) LINE
      READ(CONIN,30) LINE
 30   FORMAT(A)
      END
