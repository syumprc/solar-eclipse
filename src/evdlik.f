      SUBROUTINE EVDLIK(N, MAX, MU, COV, LOGLIKE, H2, EVAL, EVEC)
C-----------------------------------------------------------------
C
C       10/12/09
C       John Blangero SFBR
C       Thomas Dyer SFBR
C
C       This subroutine implements ...
C
C-----------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                 
      INTEGER N, MAX
      DOUBLE PRECISION MU(N), LOGLIKE
      DOUBLE PRECISION COV(MAX,MAX)
      DOUBLE PRECISION EVAL(N), EVEC(N,N)
      DOUBLE PRECISION Z(N), TAU(N)
      DOUBLE PRECISION SD(N), DP(N)

C      print *,"entered evdlik"
C      print *,"n=",n
C      print *,"max=",max
C      do I = 1,N
C         print *,"mu(",i,")=",mu(I)
C      enddo
C      do I=1,N
C         do j=1,N
C            print *,"cov(",i,",",j,")=",cov(i,j)
C         enddo
C      enddo
C      do I = 1,N
C         print *,"eval(",i,")=",eval(i)
C      enddo
C      do I=1,N
C         do j=1,N
C            print *,"evec(",i,",",j,")=",evec(i,j)
C         enddo
C      enddo
     
      DO I = 1, N
         SD(I) = DSQRT(COV(I,I))
         Z(I) = MU(I)/SD(I)
         DP(I) = EVAL(I)*H2 + 1 - H2
      ENDDO

      LOGLIKE = 0
      DO I = 1, N
         TAU(I) = 0
         DO J = 1, N
            TAU(I) = TAU(I) + EVEC(J,I)*Z(J)
         ENDDO
         TAU(I) = TAU(I)/DSQRT(DP(I))
      ENDDO

      DO I = 1, N
         LOGLIKE = LOGLIKE - .5*TAU(I)**2 - .5*DLOG(DP(I)) - DLOG(SD(I))
      ENDDO

      RETURN
      END
