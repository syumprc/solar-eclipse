      SUBROUTINE HESSIN(DF,HESS,OBS,PAR,F,ITER,KASE,NCASE,
     &NOBS,NPAR,DIFFER)                                             
C                                                                    
C     THIS SUBROUTINE PERMITS RECOMPUTATION OF THE OBJECTIVE FUNCTION
C     F AND ITS DIFFERENTIAL DF.  IF DIFFER(2) IS TRUE, THEN PROVIDE 
C     AN APPROXIMATION TO THE SECOND DIFFERENTIAL, I.E, HESSIAN OF F.
C                                                                    
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                             
      DOUBLE PRECISION DF(NPAR),HESS(NPAR,NPAR),OBS(NOBS),
     &PAR(NPAR)                                                     
      LOGICAL DIFFER(2)                                              
C                                                                    
      END                                                            
