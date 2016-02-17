      SUBROUTINE MVNCDF(N, MAX, AFFECT, DISC, MU, RHO, LOGLIKE)
C-----------------------------------------------------------------
C
C       8/21/96
C       John Blangero SFBR
C
C       9/13/06
C       approximation improved
C       Thomas Dyer SFBR
C
C       This subroutine is designed to calculate the MVN CDF
C       using the Mendell-Elston procedure as described
C       by Hasstedt.
C
C-----------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                 
      INTEGER N, MAX
      DOUBLE PRECISION AFFECT(N), MU(N), LOGLIKE, PROB
      DOUBLE PRECISION S(N), T(N), A
      DOUBLE PRECISION BIGV(N), SMV
      DOUBLE PRECISION RHO(MAX,MAX)
      DOUBLE PRECISION ALNORM, PHIDENS, ZERO
      DOUBLE PRECISION PHIT, PHIS, ZMIN
      DOUBLE PRECISION CDFNT, CDFNS
      CHARACTER*1 DISC(N)
      EXTERNAL ALNORM, PHIDENS

      ZMIN = -1.D37
      ZERO = 0

      DO 2 I = 1, N
        IF (DISC(I) .EQ. 'd') THEN
          IF (AFFECT(I) .EQ. 1) THEN
            S(I) = MU(I)             
            T(I) = 100
          ELSE
            S(I) = -100
            T(I) = MU(I)
          ENDIF  
        ELSEIF (DISC(I) .EQ. 'c') THEN
          S(I) = MU(I)
          SMV = 1
        ENDIF
        BIGV(I) = 1
2     CONTINUE

      LOGLIKE = 0
      DO 3 I = 1, N
        IF (DISC(I) .EQ. 'd') THEN
          CDFNS = ALNORM(S(I), .FALSE.)
          CDFNT = ALNORM(T(I), .FALSE.)
          PHIS = PHIDENS(S(I), ZERO)
          PHIT = PHIDENS(T(I), ZERO)
          A = (PHIS - PHIT)/(CDFNT - CDFNS)
          PROB = CDFNT - CDFNS
          IF (PROB .GT. 0) THEN
            LOGLIKE = LOGLIKE + DLOG(CDFNT - CDFNS)
          ELSE
            LOGLIKE = ZMIN
            RETURN
          ENDIF
          SMV = A**2 - (S(I)*PHIS - T(I)*PHIT)/(CDFNT - CDFNS)
        ELSEIF (DISC(I) .EQ. 'c') THEN
          A = S(I)
          PROB = PHIDENS(S(I), ZERO)/DSQRT(BIGV(I))
          IF (PROB .GT. 0) THEN
            LOGLIKE = LOGLIKE +
     &                DLOG(PHIDENS(S(I), ZERO)) - 0.5*DLOG(BIGV(I))
          ELSE
            LOGLIKE = ZMIN
            RETURN
          ENDIF
        ENDIF

        DO 4 J = I+1, N
          S(J) = (S(J) - RHO(I,J)*A)/DSQRT(1 - RHO(I,J)**2*SMV)
          T(J) = (T(J) - RHO(I,J)*A)/DSQRT(1 - RHO(I,J)**2*SMV)
          BIGV(J) = BIGV(J)*(1 - RHO(I,J)**2*SMV) 

          DO 5 K = J+1, N
            RHO(J,K) = RHO(J,K) - RHO(I,J)*RHO(I,K)*SMV
            RHO(J,K) = RHO(J,K)/(DSQRT(1 - RHO(I,J)**2*SMV)
     &                           *DSQRT(1 - RHO(I,K)**2*SMV))
            RHO(K,J) = RHO(J,K) 
5         CONTINUE
4       CONTINUE
3     CONTINUE

      RETURN
      END
