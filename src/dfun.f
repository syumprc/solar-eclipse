      SUBROUTINE DFUN(DF,DG,GHESS,HESS,OBS,PAR,PARMAX,WORK1,DP,F
     1,ITER,MAXPAR,NCASE,NDERIV,NOBS,NPASS,DIFFER,FORWRD,MALE,TRTTYPE
     2,VTRAITS,VARDATA,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT
     3,MAXPEO,PARMIN)
C                                                                        
C     THIS SUBROUTINE CONTROLS THE COMPUTATION OF F AND ITS FIRST        
C     NDERIV PARTIAL DERIVATIVES.  WHEN THE DERIVATIVES ARE              
C     COMPUTED NUMERICALLY, FORWARD OR CENTRAL DIFFERENCES ARE           
C     USED DEPENDING ON THE LOGICAL VARIABLE FORWRD.  DIFFER(1)          
C     SHOULD BE INPUT AS TRUE WHEN ANALYTIC DERIVATIVES ARE              
C     AVAILABLE.  THE NUMERICAL DIFFERENTIATION INTERVAL IS              
C     ADJUSTED TO TAKE INTO ACCOUNT THE MAGNITUDE OF A PARAMATER         
C     AND WHETHER IT LIES ON ITS UPPER BOUND.  DIFFER(2) IS TRUE         
C     WHEN AN APPROXIMATE HESSIAN IS AVAILABLE.                          
C                                                                        
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                 
      DOUBLE PRECISION DF(MAXPAR),DG(MAXPAR),GHESS(MAXPAR,MAXPAR)
     1,HESS(MAXPAR,MAXPAR),OBS(NOBS,NCASE),PAR(MAXPAR)
     2,PARMAX(MAXPAR),WORK1(MAXPAR),COV(MAXPEO,MAXPEO),MU(MAXPEO)
     3,AFFECT(MAXPEO),VARDATA(NVAR,NTOT),PARMIN(MAXPAR),GPLUS
     4,GMINUS,GMINUS2,BOUNDIFF
      INTEGER PASS,NVAR,NTOT,NIND(NPED),NASCER(NPED),NPED,MAXPEO,VTRAITS
      INTEGER NCUMIND(0:NPED)
      LOGICAL DIFFER(2),FORWRD,MALE(NTOT)
      CHARACTER*1 TRTTYPE(NTOT)
C                                                                        
      CALL DOPTION ('BounDiff',BOUNDIFF)

      F=0.                                                               
      DO 10 I=1,NDERIV                                                   
 10   DF(I)=0.                                                           
      IF (DIFFER(2)) THEN                                                
      DO 20 J=1,NDERIV                                                   
      DO 30 I=1,NDERIV                                                   
 30   HESS(I,J)=0.                                                       
 20   CONTINUE                                                           
      END IF                                                             
C                                                                        
C     DO NOTHING ON THE FIRST FEW PASSES.                                
C                                                                        
      DO 40 PASS=1,NPASS-1                                               
      DO 50 KASE=1,NCASE                                                 
 50   CALL DDFUN(DG,OBS(1,KASE),PAR,G,ITER,KASE,NCASE,NOBS
     1,MAXPAR,NPASS,PASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
 40   CONTINUE                                                           
C                                                                        
C     COMPUTE THE FUNCTION AND ITS DERIVATIVES ON THE LAST PASS.         
C                                                                        
      DO 60 KASE=1,NCASE                                                 
      IF (DIFFER(1)) THEN                                                
      DO 70 I=1,NDERIV                                                   
 70   DG(I)=0.                                                           
      END IF                                                             
      CALL DDFUN(DG,OBS(1,KASE),PAR,G,ITER,KASE,NCASE,NOBS
     1,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)

      IF (.NOT.DIFFER(1)) THEN                                           
      DO 80 I=1,NDERIV                                                   
      D=DP*MAX(ABS(PAR(I)),1.0D0)                                        
      PTEMP=PAR(I)                                                       
      IF (FORWRD) THEN                                                   
C                                                                        
C     COMPUTE THE PARTIAL DERIVATIVE BY A FORWARD DIFFERENCE.            
C                                                                        
         IF (PTEMP+D.LE.PARMAX(I)
     1.OR. BOUNDIFF.EQ.-1 .OR. BOUNDIFF.EQ.2
     2.OR.(BOUNDIFF.EQ.0.AND.(PARMAX(I).NE.0.AND.PARMAX(I).NE.1))) THEN

            IF (PTEMP+D.GT.PARMAX(I).AND.BOUNDIFF.EQ.2) D = -D

            PAR(I)=PTEMP+D
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GPLUS,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
            DG(I)=(GPLUS-G)/D

         ELSE

C If we would otherwise pass over an upper boundary, then
C Compute two lower differences to estimate what the forward
C difference would be expected to be

            PAR(I) = PTEMP - 2*D
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GMINUS2,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
            PAR(I) = PTEMP - D
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GMINUS,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
            PAR(I) = PTEMP
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GPLUS,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)

C Mestimate = M01 + (M01 - M12) = 2*M01 - M12

            DG(I) = 2.0*(GPLUS-GMINUS)/D - ((GMINUS-GMINUS2)/D)

         ENDIF

      ELSE                                                               
C                                                                        
C     COMPUTE THE PARTIAL DERIVATIVE BY A CENTRAL DIFFERENCE.            
C
C        but check upper and lower bounds, computing an estimate if
C        either bound would be exceeded.  The estimate is based on
C        a smaller delta than used for forward differences for
C        empirical reasons.
C
         IF (PTEMP-D.GE.PARMIN(I).AND.(PTEMP+D.GT.PARMAX(I)
     1.AND.BOUNDIFF.NE.-1.AND.BOUNDIFF.NE.2.AND.(.NOT.(BOUNDIFF.EQ.0
     2.AND.(.NOT.(PARMAX(I).EQ.0.OR.PARMAX(I).EQ.1)))))) THEN

C Upper bound would be crossed

            PAR(I) = PTEMP - (D*2)
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GMINUS2,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
            PAR(I) = PTEMP - D
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GMINUS,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
            PAR(I) = PTEMP
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GPLUS,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)

C Mcen = Mneg + Mpos
C Mpos estimate is Mneg + (Mneg - Mneg2)
C therefore, Mcen = 3*Mneg - Mneg2

            DG(I) = ((3.0*(GPLUS-GMINUS)/D)-((GMINUS-GMINUS2)/D))/2.0

         ELSE IF (PTEMP-D.LT.PARMIN(I).AND.BOUNDIFF.NE.-1
     1.AND.BOUNDIFF.NE.2.AND.(.NOT.(BOUNDIFF.EQ.0
     2.AND.(.NOT.(PARMIN(I).EQ.-1.OR.PARMIN(I).EQ.0))))) THEN

C Lower bound would be crossed

            PAR(I) = PTEMP
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GMINUS2,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
            PAR(I) = PTEMP + D
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GMINUS,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
            PAR(I) = PTEMP + (D*2)
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GPLUS,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)

            DG(I) = ((3.0*(GMINUS-GMINUS2)/D)-((GPLUS-GMINUS)/D))/2.0

C Earlier incorrect formula that also
C surprisingly it worked quite well
C            DG(I) = (3.0*(GPLUS-GMINUS)/D) - ((GMINUS-GMINUS2)/D)


         ELSE

C No boundary crossed.  This is the default case

            PAR(I)=PTEMP+D
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GPLUS,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
            PAR(I)=PTEMP-D
            CALL DDFUN(WORK1,OBS(1,KASE),PAR,GMINUS,ITER,KASE,NCASE
     1,NOBS,MAXPAR,NPASS,NPASS,DIFFER,VARDATA,MALE,TRTTYPE,VTRAITS
     2,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO)
            DG(I)=(GPLUS-GMINUS)/(2*D)


         END IF
      END IF   

      PAR(I)=PTEMP                                                       


 80   CONTINUE                                                           
      END IF                                                             
C                                                                        
C     UPDATE THE FUNCTION AND ITS DIFFERENTIAL.  APPROXIMATE THE         
C     HESSIAN IF DESIRED.                                                
C                                                                        
      IF (DIFFER(2)) THEN                                                
      DO 90 J=1,NDERIV                                                   
      DO 100 I=1,NDERIV                                                  
 100  GHESS(I,J)=0.                                                      
 90   CONTINUE                                                           
      END IF                                                             
      CALL HESSIN(DG,GHESS,OBS(1,KASE),PAR,G,ITER,KASE,NCASE       
     :,NOBS,MAXPAR,DIFFER)                                               
      F=F+G                                                              
      DO 110 I=1,NDERIV                                                  
 110  DF(I)=DF(I)+DG(I)                                                  
      IF (DIFFER(2)) THEN                                                
      DO 120 J=1,NDERIV                                                  
      DO 130 I=1,NDERIV                                                  
 130  HESS(I,J)=HESS(I,J)+GHESS(I,J)                                     
 120  CONTINUE                                                           
      END IF                                                             
 60   CONTINUE                                                           
      END                                                                
