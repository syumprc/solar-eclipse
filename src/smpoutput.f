C Modified 1996-1997 by Charles Peterson, SFBR.
C This is the New Fisher signature.  Do not remove.  fI6%.
C
      SUBROUTINE SMPOUTPUT(MODFIL,NCNSTR,NPAR,CNSTR,PNAME,PAR,IBIG,
     &                      ITER,LOGLIK,QDFORM,NSAMP,NTRAIT,LAST,UNIT4,
     &                      WORK2,MAXTAB,MAXPAR,INASYC,TARRAY,NPED,
     &                      PARMAX,PARMIN,
     &                      PAREST,UNIT3,CONOUT,IROBST,IDISCRETE,BIG)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DOUBLE PRECISION PAR(MAXPAR),LOGLIK,
     &                 CNSTR(NCNSTR,NPAR),
     &                 WORK2(MAXTAB),
     &                 PARMIN(*),PARMAX(*),PAREST(*),BIG,MF,SDVAL
      INTEGER          I,J,UNIT4,NCNSTR,IBIG,ITER,UNIT3
      INTEGER          CONOUT,CNCOUNT,HOWCLOSE,DIGITS
      INTEGER          GETVIRT,NVTRAIT,SDPAR
      CHARACTER        PNAME(MAXPAR)*(*)
      CHARACTER*40     MODFIL
      LOGICAL          INASYC,IROBST,IDISCRETE
      LOGICAL          VERBOSE,VRESULT,VLOGLIKE,SDLOW
      REAL             TARRAY(2)
      CHARACTER*320    CNSTRING
      CHARACTER*300    CNSTRING2
 
      VRESULT = VERBOSE ("RESULTS")
      VLOGLIKE = VERBOSE ("LOGLIKE")

      IF (ITER .GT. LAST) ITER = ITER - 1

      IF (VRESULT) WRITE(CONOUT,40)
      WRITE(UNIT4,40) 
 40   FORMAT (//T25,'Model Parameter Optimizations'/)
      CALL FLUSH (UNIT4)
      CALL FLUSH (CONOUT)

       CALL WRITEPHEAD (UNIT4, 2)
       IF (VRESULT) CALL WRITEPHEAD (CONOUT, 2)

       II = 0
       DO I = 1, NPAR

C Set up constraint description string
C   Note to old-timers: this replaces the old M array
         
          CNCOUNT = 0
          CNSTRING = ' '
          DO 5 J = 1, NCNSTR
             IF (CNSTR(J,I) .NE. 0) THEN
                IF (CNCOUNT .EQ. 0) THEN
                   WRITE (CNSTRING, 1) J
 1                 FORMAT ('#',I2)
                ELSE
                   CNSTRING2 = CNSTRING(1:LNBLNK(CNSTRING))
                   IF (J.LE.9) THEN
                      WRITE(CNSTRING,2) CNSTRING2(1:LNBLNK(CNSTRING2)),J
 2                    FORMAT (A,',',I1)
                   ELSE
                      WRITE(CNSTRING,3) CNSTRING2(1:LNBLNK(CNSTRING2)),J
 3                    FORMAT (A,',',I2)
                   ENDIF
                ENDIF
                CNCOUNT = CNCOUNT + 1
             ENDIF
 5        CONTINUE

          IF (INASYC.AND.WORK2(I).NE.0) THEN

C We got standard error (WORK2(I))

             CALL WRITERESULTS (UNIT4,I,PAR(I),WORK2(I),
     &            CNSTRING(1:LNBLNK(CNSTRING)),
     &            PAREST(I), PARMIN(I), PARMAX(I))
             IF (VRESULT) THEN
                CALL WRITERESULTS (CONOUT,I,PAR(I),WORK2(I),
     &               CNSTRING(1:LNBLNK(CNSTRING)),
     &               PAREST(I), PARMIN(I), PARMAX(I))
             ENDIF
          ELSE

C We didn't get standard error, use -1.0 to signal that

             CALL WRITERESULTS (UNIT4,I,PAR(I),-1.0D0,
     &            CNSTRING(1:LNBLNK(CNSTRING)),
     &            PAREST(I), PARMIN(I), PARMAX(I))
             IF (VRESULT) THEN
                CALL WRITERESULTS (CONOUT,I,PAR(I),-1.0D0,
     &               CNSTRING(1:LNBLNK(CNSTRING)),
     &               PAREST(I), PARMIN(I), PARMAX(I))
             ENDIF
          ENDIF
      ENDDO

      IF (ITER .NE. IBIG) THEN
         DIGITS = HOWCLOSE(BIG,LOGLIK)
         IF (DIGITS .LT. 9) THEN
            IF (VLOGLIKE) THEN
               WRITE(CONOUT,15) IBIG
               WRITE(UNIT4,15) IBIG
 15            FORMAT(/,'** Max Likelihood at Iteration ',I5,
     &              ', Not the Last Iteration **')
            ENDIF
         ENDIF
      ENDIF

      WRITE(UNIT4,35) LOGLIK
      IF (VLOGLIKE) WRITE(CONOUT,35) LOGLIK
 35   FORMAT(/T28,'Loglikelihood =',F15.6)

      IF (.NOT.IDISCRETE) THEN
         NVTRAIT = GETVIRT()
         IF (NVTRAIT.EQ.0) THEN
            NVTRAIT = NTRAIT
         ENDIF
         IF (VRESULT) WRITE(CONOUT,45) QDFORM/(NSAMP*NVTRAIT) 
         WRITE(UNIT4,45) QDFORM/(NSAMP*NTRAIT) 
 45      FORMAT(/T21,'Normalized Quadratic =',F15.6)
      ENDIF

      IF (VRESULT) WRITE(CONOUT,25) ITER
      WRITE(UNIT4,25) ITER
   25 FORMAT(/T31,'Iterations = ',I7)

C Check standard deviations

      SDLOW = .FALSE.
      SDVAL = 0.5
      IF (SDPAR(1).GT.0) THEN
         IF (PAR(SDPAR(1)) .LT. 0.5) THEN
            SDLOW = .TRUE.
            SDVAL = PAR(SDPAR(1))
         ENDIF
      ENDIF
      IF (SDPAR(2).GT.0) THEN
         IF (PAR(SDPAR(2)) .LT. 0.5) THEN
            SDLOW = .TRUE.
            IF (PAR(SDPAR(2)).LT.SDVAL) THEN
               SDVAL = PAR(SDPAR(2))
            ENDIF
         ENDIF
      ENDIF
      IF (SDLOW) THEN
         IF (SDVAL .NE. 0) THEN
            MF = 1.0 / SDVAL
            WRITE (UNIT4,55) MF
            IF (VRESULT) WRITE (CONOUT,55) MF
 55   FORMAT (/,T8,'Warning!  Trait Standard Deviation is below 0.5'
     */T8,'Multiplying trait values by ',F6.1,' is recommended to'
     *,' ensure accuracy'/)
         ELSE
            WRITE (UNIT4,56)
            IF (VRESULT) WRITE (CONOUT,56)
 56   FORMAT (/,T8,'Warning!  Trait Standard Deviation is below 0.5'
     */T8,'Multiplying trait by a factor is recommended to'
     *,' ensure accuracy'/)
         ENDIF
      ENDIF

      END
