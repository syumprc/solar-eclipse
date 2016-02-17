C Modified 1996-1997 by Charles Peterson, SFBR.
C This is the New Fisher signature.  Do not remove.  fI6%.
C
      SUBROUTINE OUTPUT(PAR,PNAME,LOGLIK,QUADNORM,ITER
     1,LAST,NPAR,NSTEP,UNIT3,STAND
     2,BIG_ITER,CONOUT,DISCRETE,BIG)
C
C     THIS SUBROUTINE OUTPUTS THE LOGLIKELIHOOD AND PARAMETERS
C     AT EACH ITERATION.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION PAR(NPAR),LOGLIK,QUADNORM,BIG
      INTEGER I,UNIT3,CONOUT,BIG_ITER
      CHARACTER PNAME(NPAR)*(*)
      LOGICAL STAND, VERBOSE, VITERATE, DISCRETE
      SAVE FIRST,IBIG
      character*32 modeltype
      logical notevd
      integer evdphase
C
      notevd = .true.
      viterate = verbose ("ITERATE")

      call soption ("ModelType",modeltype)
      call ioption ("EvdOptimizations",evdphase)
      if (modeltype.eq."Evd") then
         notevd = .false.
         viterate = .false.
      end if

      IF (ITER.EQ.1) THEN
      BIG=-1.0D20
      FIRST=LOGLIK
      IF (NOTEVD) WRITE(UNIT3,5)
 5    FORMAT(/T31,'Iteration History')

      IF (NOTEVD) WRITE(UNIT3,10) (PNAME(I),I=1,NPAR) 
      IF (VITERATE) WRITE(CONOUT,10) (PNAME(I),I=1,NPAR)
 10   FORMAT(/,'Iter Step Loglikelihood',(T24,4(1X,A13),:))
      END IF

      IF (STAND) LOGLIK=LOGLIK-FIRST
      IF (LOGLIK.GE.BIG) THEN
      IBIG=ITER
      BIG=LOGLIK
      END IF

 20   FORMAT(/,I3,1X,I3,2X,D14.7,(T24,4(1X,D13.6),:))
 30   FORMAT(//' WARNING: The maximum loglikelihood does not occur at'
     1,' the last iteration.'/
     2' It occurs at iteration ',I4,'.'/)
 40   FORMAT(/,'    The normalized quadratic form at the preceding'
     1,' iteration is',D11.4,'.')

      IF (NOTEVD) THEN
         WRITE(UNIT3,20) ITER,NSTEP,LOGLIK,(PAR(I),I=1,NPAR)
         IF (VITERATE) 
     1        WRITE(CONOUT,20) ITER,NSTEP,LOGLIK,(PAR(I),I=1,NPAR)
         IF (.NOT. DISCRETE) THEN
            WRITE(UNIT3,40) QUADNORM
            IF (VITERATE)
     1           WRITE(CONOUT,40) QUADNORM
         ENDIF
      END IF

C  This now done in smp_output.f
C      IF (ITER.EQ.LAST .AND. ITER .NE. IBIG) THEN
C         WRITE(UNIT3,30) IBIG
C         WRITE(CONOUT,30) IBIG
C      END IF

      BIG_ITER=IBIG
      END
