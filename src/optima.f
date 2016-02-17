      SUBROUTINE OPTIMA(CNSTR,CVALUE,DELTA,DF,DFOLD,GRID,HESS,OBS        
     1,PAR,PARMAX,PARMIN,PAROLD,TABLE,WORK1,WORK2,NSWEEP,PNAME,CONV      
     2,DP,F,SMALL,TOL,IOUNIT,MAXPAR,MAXTAB,MCNSTR,MXITER,MXSTEP    
     3,NCASE,NCNSTR,NCONV,NOBS,NPAR,NPASS,NPOINT,STDERR,TRAVEL,DIFFER   
     4,VMEAN,VVAR,VMIN,VMAX,NVAR,NTRAIT,UNIT3,CONOUT,PAREST
     5,VARDATA,NTOT,NIND,NASCER,NPED,STATUS,NCUMIND,COV,MU,AFFECT
     6,MAXPEO,MALE,IBIG,LOGLIK,INASYCV,ITER,BIG,TRTTYPE,VTRAITS)
C                                                                        
C     THIS SUBROUTINE COMPUTES THE CONSTRAINED MINIMUM OF A FUNCTION F   
C     BY THE VARIABLE METRIC METHOD OF BIGGS, HAN, AND POWELL.  SEE:     
C     M.J.D. POWELL(1978) "A FAST ALGORITHM FOR NONLINEARLY CONSTRAINED  
C     OPTIMIZATION CALCULATIONS."  PROCEEDINGS OF THE 1977 DUNDEE        
C     CONFERENCE ON NUMERICAL ANALYSIS. G.A. WATSON EDITOR. SPRINGER     
C     VERLAG.  FOR LEAST SQUARES PROBLEMS IT USES THE CLASSICAL GAUSS    
C     NEWTON METHOD.  IN THE LIST BELOW ITEMS MARKED BY * OR ** MUST     
C     BE PROVIDED BY THE USER.  THOSE ITEMS MARKED BY * SHOULD BE        
C     PASSED TO SUBROUTINE SEARCH.  THOSE ITEMS MARKED BY ** SHOULD BE   
C     INITIALIZED IN SUBROUTINE INITAL.                                  
C                                                                        
C     CNSTR**          MATRIX OF LINEAR EQUALITY CONSTRAINTS             
C     CVALUE**         CONSTANTS FOR PARAMETER/CONSTRAINT INNER PRODUCTS 
C     DELTA            UPDATE DIRECTION FOR PARAMETERS                   
C     DF               CURRENT DIFFERENTIAL OF THE FUNCTION F            
C     DFOLD            PREVIOUS DIFFERENTIAL OF THE FUNCTION F           
C     F                VALUE OF FUNCTION TO BE MINIMIZED                 
C     GRID**           GRID OF POINTS TO EVALUATE F ON                   
C     HESS**           CURRENT APPROXIMATE HESSIAN OF F                  
C     OBS**            MATRIX OF OBSERVATIONS                            
C     PAR**            CURRENT PARAMETERS                                
C     PARMAX,PARMIN**  PARAMETER MAXIMA AND MINIMA                       
C     PAROLD           PREVIOUS PARAMETERS                               
C     TABLE            TABLEAU FOR QUADRATIC PROGRAMMING PROBLEM         
C     WORK1,WORK2      WORK VECTORS                                      
C     NSWEEP           INDICATOR FOR WHICH PARAMETERS HAVE BEEN SWEPT    
C                        IN TABLE                                        
C     PNAME**          PARAMETER NAMES                                   
C     CONV*            CONVERGENCE CRITERION FOR CHANGE IN F             
C     DP*              NUMERICAL DIFFERENTIATION INTERVAL                
C     SMALL*           SMALL POSITIVE NUMBER FOR CHECKING BOUNDS         
C     TOL*             TOLERANCE FOR MATRIX SWEEPING                     
C     IOUNIT*          OUTPUT UNIT NUMBER                                
C     MAXPAR           MAXIMUM(NPAR,1)                                   
C     MAXTAB           NCNSTR+NPAR+1                                     
C     MCNSTR           MAXIMUM(NCNSTR,1)                                 
C     MODEL*           USER MODEL NUMBER                                 
C     MXITER*          MAXIMUM NUMBER OF ITERATIONS                      
C     MXSTEP*          MAXIMUM NUMBER OF STEPS PER ITERATION             
C     NCASE*           NUMBER OF CASES IN A PROBLEM                      
C     NOBS*            NUMBER OF OBSERVATIONS PER CASE                   
C     NPASS*           NUMBER OF PASSES PER ITERATION                    
C     NCONV*           NUMBER OF TIMES CONVERGENCE CRITERION MUST BE MET 
C     NCNSTR*          NUMBER OF LINEAR EQUALITY CONSTRAINTS             
C     NPAR*            NUMBER OF PARAMETERS                              
C     NPOINT*          NUMBER OF POINTS FOR 'GRID' OPTION                
C     TRAVEL*          'SEARCH' OR 'GRID' OPTION                         
C     STDERR*          0 FOR NO STANDARD ERRORS, 1 FOR STANDARD ERRORS   
C                      BASED ON THE OBSERVED INFORMATION, 2 FOR STANDARD 
C                      ERRORS BASED ON THE FINAL APPROXIMATE HESSIAN, 3  
C                      FOR STANDARD ERRORS IN A LEAST SQUARES PROBLEM    
C     DIFFER*          DIFFER(1) TRUE FOR EXACT FIRST DIFFERENTIAL       
C                      DIFFER(2) TRUE FOR APPROXIMATE SECOND DIFFERENTIAL
C     MODFIL*          MODEL FILE NAME
C     OUTFIL*          OUTPUT FILE NAME
C                                                                        
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                 
      DOUBLE PRECISION CNSTR(MCNSTR,MAXPAR),CVALUE(MCNSTR)               
     1,DELTA(MAXPAR),DF(MAXPAR),DFOLD(MAXPAR),GRID(NPOINT,MAXPAR)        
     2,HESS(MAXPAR,MAXPAR),OBS(NOBS,NCASE),PAR(MAXPAR),PARMAX(MAXPAR)    
     3,PARMIN(MAXPAR),PAROLD(MAXPAR),TABLE(MAXTAB,MAXTAB),WORK1(MAXPAR)  
     4,WORK2(MAXTAB),VMEAN(NVAR),VVAR(NVAR),VMIN(NVAR),VMAX(NVAR)
     5,PAREST(MAXPAR),VARDATA(NVAR,NTOT),COV(MAXPEO,MAXPEO),MU(MAXPEO)
     6,AFFECT(MAXPEO),CONV,DP,F,SMALL,TOL,BIG
      INTEGER NSWEEP(MAXPAR),PROBLM,SEED,STDERR,NVAR,NTRAIT,UNIT3
     1,CONOUT,IERROR,NTOT,NPED,NIND(NPED),NASCER(NPED),STATUS
     2,NCUMIND(0:NPED),MAXPEO,IBIG,IOUNIT,VTRAITS
      CHARACTER*8 PNAME(MAXPAR),TRAVEL                                   
      CHARACTER*32 MODELTYPE
      CHARACTER*1 TRTTYPE(NTOT)
      LOGICAL DIFF(2),DIFFER(2),FORWRD,INASYCV,STAND,MALE(NTOT)
      LOGICAL DISCRETE,VERBOSE,CHANGED,LCHANGED,VITERATE
      INTEGER ISNAN_F,MOMEGA,OMEGATYP,IFDISC,EVDPHASE
      DOUBLE PRECISION QDFORM,LOGLIK,BCLIFF,MXCLIFFS
      SAVE PROBLM,SEED,DISCRETE
      DATA PROBLM,SEED,DISCRETE/0,25431,.TRUE./

      VITERATE = VERBOSE('ITERATE')
      WRITE (UNIT3,1)
      IF (VITERATE) WRITE (CONOUT,1)
C
C WARNING!  DO NOT CHANGE THE FOLLOWING MESSAGE WITHOUT EDITING
C THE LINE IN SOLAR.TCL THAT TESTS THIS MESSAGE FOR RESIDUAL
C
 1    FORMAT(/
     1,'                   '
     2,'Using SOLAR Discrete and Mixed Trait Modeling')

C
C If high verbosity, output model type and evd info to terminal
C If not high verbosity, only output to solar.out file
C
      CALL SOPTION ('MODELTYPE',MODELTYPE)
      WRITE (UNIT3,2) MODELTYPE
 2    FORMAT(/
     1,'                   Model Type: ',A)
      IF (VITERATE) WRITE (CONOUT,2) MODELTYPE

      call ioption ('EVDPHASE',evdphase)
      call ifanydisc (ifdisc)
      momega= omegatyp()
      if (evdphase.eq.0) then
         if (modeltype.eq."Evd".and.((momega.lt.1.or.momega.gt.2).or.
     &        ifdisc.eq.1)) then
C
C evd1 modeltype set but cannot process as evd
C give warning as to why evd cannot be used
C
            write (unit3, 4)
            if (viterate) write (conout,4)
 4    format('                   * Evd Processing is Not Applicable')
            if (ifdisc.eq.0) then
               write (unit3, 5)
               if (viterate) write (conout,5)
 5    format('                   * Because of Model Omega')
            else
               write (unit3, 6)
               if (viterate) write (conout,6)
 6    format('                   * Because of discrete trait(s)')
            end if
         end if 
      end if

C tobit now obsolescent
C      IF (MODELTYPE .EQ. 'Tobit') THEN
C         WRITE (CONOUT,3)
C      ENDIF
C 3    FORMAT (/
C     1,'                   Using tobit feature')


      STATUS = 0
      QDFORM=0.0
      STAND=.FALSE.
      CALL DOPTION ('BCliff',BCLIFF)
      CALL DOPTION ('MaxCliffs',MXCLIFFS)
C                                                                        
C     INITIALIZE THE PARAMETER VALUES, THEIR BOUNDS, THEIR CONSTRAINTS,  
C     AND THEIR NAMES.                                                   
C                                                                        
      IF (NPAR.GT.0) THEN                                                
      DO 10 I=1,NPAR                                                     
C     PAR(I)=1.0D-6                                                      
C     PARMAX(I)=1.0D20                                                   
C     PARMIN(I)=-1.0D20                                                  
C     WRITE(PNAME(I),'(I6)') I                                           
C     PNAME(I)(2:4)='PAR'                                                
      DO 10 J=1,NCNSTR                                                   
      CNSTR(J,I)=0.0D0                                                   
 10   CVALUE(J)=0.0D0                                                    
C                                                                        
C     SET THE INITIAL HESSIAN TO THE IDENTITY.                           
C                                                                        
      IF (TRAVEL.EQ.'SEARCH') THEN                                       
      DO 20 I=1,NPAR                                                     
      DO 30 J=1,NPAR                                                     
 30   HESS(J,I)=0.0D0                                                    
 20   HESS(I,I)=1.0D0                                                    
C                                                                        
C     FILL IN THE GRID WITH A SMALL POSITIVE VALUE.                      
C                                                                        
      ELSE IF (TRAVEL.EQ.'GRID') THEN                                    
      DO 40 I=1,NPAR                                                     
      DO 40 J=1,NPOINT                                                   
 40   GRID(J,I)=1.0D-6                                                   
      END IF                                                             
C                                                                        
C     LET C++ CHANGE THE INITIAL SETTINGS.                          
C                                                                        
      IERROR=0
      CALL INITAL(CNSTR,CVALUE,PAR,PARMAX,PARMIN,VMEAN,VVAR,PNAME
     1,MCNSTR,NPAR,NTRAIT,NVAR,UNIT3,VMIN,VMAX,CONOUT,PAREST,IERROR)
      IF (IERROR.NE.0) THEN
         STATUS = 1
         RETURN
      END IF
      END IF                                                             
C                                                                        
C     COMPUTE FUNCTION VALUES OVER A USER DEFINED GRID OF POINTS.        
C                                                                        
      IF (TRAVEL.EQ.'GRID') THEN                                         
      LAST=NPOINT                                                        
      DO 50 ITER=1,NPOINT                                                
      DO 60 J=1,NPAR                                                     
 60   PAROLD(J)=GRID(ITER,J)                                             
      CALL DFUN(DF,WORK2,TABLE,HESS,OBS,PAROLD,PARMAX,WORK1,DP,FOLD
     1,ITER,MAXPAR,NCASE,0,NOBS,NPASS,DIFFER,FORWRD,MALE,TRTTYPE
     2,VTRAITS
     3,VARDATA,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO
     4,PARMIN)
      IF (ITER.EQ.1) F=FOLD                                              
      IF (FOLD.LE.F) THEN                                                
      F=FOLD                                                             
      DO 70 J=1,NPAR                                                     
 70   PAR(J)=PAROLD(J)                                                   
      END IF                                                             
      LOGLIK = -FOLD
      CALL OUTPUT(PAROLD,PNAME,LOGLIK,QDFORM,ITER,LAST
     1,NPAR,0,IOUNIT,STAND,IBIG,CONOUT,DISCRETE,BIG)
 50   CONTINUE
C                                                                        
C     OTHERWISE SEARCH THE FUNCTION SURFACE.  FIRST CALL PREOPT          
C     TO CHECK THAT THE PARAMETERS SATISFY THEIR BOUNDS AND THEIR        
C     LINEAR EQUALITY CONSTRAINTS.  PREOPT ALSO CHECKS THAT THE          
C     CONSTRAINTS ARE NOT REDUNDANT.  IF ANY OF THESE CHECKS FAIL,       
C     THEN STOP.                                                         
C                                                                        
      ELSE                                                               
      CALL PREOPT(CNSTR,CVALUE,PAR,PARMAX,PARMIN,TABLE,WORK2             
     1,PNAME,CNORM,TOL,IERROR,MAXPAR,MAXTAB,MCNSTR,NCNSTR,NPAR
     2,UNIT3,TRAVEL,CONOUT)
      IF (IERROR.GE.1) THEN
         STATUS = 1
         RETURN
      END IF
C                                                                        
C     INITIALIZE SOME VARIABLES.                                         
C                                                                        
      ITER=1                                                             
      LAST=MXITER                                                        
      NCRIT=0                                                            
      FORWRD=.TRUE.                                                      
C                                                                        
C     COMPUTE THE INITIAL FUNCTION VALUE AND DIFFERENTIAL, AND OUTPUT    
C     THE FIRST ITERATION.                                               
C                                                                        
      CALL DFUN(DF,WORK2,TABLE,HESS,OBS,PAR,PARMAX,WORK1,DP,F,ITER
     1,MAXPAR,NCASE,NPAR,NOBS,NPASS,DIFFER,FORWRD,MALE,TRTTYPE,VTRAITS
     2,VARDATA,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO
     3,PARMIN)
      LOGLIK = -F
      CALL OUTPUT(PAR,PNAME,LOGLIK,QDFORM,ITER,LAST
     1,NPAR,0,IOUNIT,STAND,IBIG,CONOUT,DISCRETE,BIG)
C                                                                        
C     ENTER THE MAIN ITERATION LOOP.                                     
C                                                                        
      IF (MXITER.GT.1) THEN                                              
      DO 80 ITER=2,MXITER                                                
C                                                                        
C     CREATE THE TABLEAU FOR THE QUADRATIC PROGRAMMING PROBLEM.          
C                                                                        
 110  CHANGED = .TRUE.

 111  LCHANGED = CHANGED
      CHANGED = .FALSE.

      CALL SETTAB(CNSTR,CVALUE,DF,HESS,PAR,TABLE,WORK2,CNORM             
     1,MAXPAR,MAXTAB,MCNSTR,NCNSTR,NPAR,NTAB)                            
C                                                                        
C     SOLVE THE QUADRATIC PROGRAMMING PROBLEM FOR THE NEXT               
C     STEP DIRECTION DELTA.                                              
C                                                                        
      CALL QDPROG(DELTA,PAR,PARMAX,PARMIN,TABLE,WORK1,WORK2,NSWEEP       
     1,SMALL,TOL,MAXPAR,MAXTAB,NCNSTR,NCYCLE,NPAR,NTAB)                  
C                                                                        
C     IF NCYCLE IS NEGATIVE, THE QUADRATIC PROGRAMMING PROBLEM IS        
C     SOLVED.  IF NCYCLE IS ZERO, IT IS IMPOSSIBLE TO ADEQUATELY         
C     SWEEP THE TABLEAU.  IF NCYCLE IS POSITIVE, THERE IS A              
C     POSSIBLE INFINITE LOOP IN THE QUADRATIC PROGRAMMING ALGORITHM.     
C     IN EITHER OF THE LAST TWO CASES, RESET THE HESSIAN AND TRY AGAIN.  
C                                                                        
      IF (NCYCLE.GE.0) THEN                                              
         HMIN=1.0D20                                                        
         DO 90 J=1,NPAR                                                     
 90         IF (HESS(J,J).GT.0.0D0) HMIN=MIN(HMIN,HESS(J,J))                   
         IF (HMIN.EQ.1.0D20) HMIN=1.0D0                                     
         DO 100 J=1,NPAR                                                    
            IF (HESS(J,J).LE.HMIN) THEN
               HESS(J,J)=HMIN*(1.0D0+RANDOM(SEED))
               CHANGED = .TRUE.
            END IF
 100     CONTINUE
C   If two passes and still no changes, call it an error
C      (CONVERR throws error) to prevent infinite loop
         IF ((.NOT. CHANGED) .AND. (.NOT. LCHANGED)) THEN
            CALL CONVERR (F,PAR)
         END IF
         GO TO 111
      END IF                                                             
C                                                                        
C     COMPUTE THE INNER PRODUCT D OF DELTA AND THE DIFFERENTIAL DF.      
C     IF D IS POSITIVE, THEN DELTA IS NOT A DESCENT DIRECTION.           
C     CONVERGENCE HAS OCCURRED, OR THE SEARCH IS IN DEEP TROUBLE.        
C                                                                        
      D=0.0D0                                                            
      DO 120 J=1,NPAR                                                    
 120  D=D+DF(J)*DELTA(J)                                                 
C                                                                        
C     WHEN THE DIFFERENTIAL IS COMPUTED NUMERICALLY, CHANGE FROM         
C     FORWARD DIFFERENCES TO CENTRAL DIFFERENCES AND VICE VERSA          
C     BASED ON THE QUANTITY D.  REDO THE QUADRATIC PROGRAMMING           
C     PROBLEM IF NECESSARY.                                              
C                                                                        
      IF (.NOT.DIFFER(1)) THEN                                           
      IF (FORWRD.AND.D.GE.0.0D0) THEN                                    
      FORWRD=.FALSE.                                                     
      CALL DFUN(DF,WORK2,TABLE,HESS,OBS,PAR,PARMAX,WORK1,DP,F,ITER
     1,MAXPAR,NCASE,NPAR,NOBS,NPASS,DIFFER,FORWRD,MALE,TRTTYPE,VTRAITS
     2,VARDATA,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO
     3,PARMIN)
      GO TO 110                                                          
      END IF                                                             
      FORWRD=D.LE.-CONV                                                  
      END IF                                                             
C                                                                        
C     ENTER THE STEP DECREMENTING LOOP.  T IS THE FRACTION OF DELTA      
C     TAKEN.                                                             
C                                                                        
      T=1.0D0                                                            
      NSTEP=0                                                            
      D=MIN(D,0.0D0)                                                     
C                                                                        
C     RECORD THE OLD DATA IN PREPARATION FOR THE NEXT STEP.              
C                                                                        
      FOLD=F                                                             
      DO 130 J=1,NPAR                                                    
      PAROLD(J)=PAR(J)                                                   
 130  DFOLD(J)=DF(J)                                                     
C                                                                        
C     COMPUTE A NEW POINT AND A NEW FUNCTION VALUE.                      
C                                                                        
 150  DO 140 J=1,NPAR                                                    
 140  PAR(J)=PAROLD(J)+T*DELTA(J)                                        
      CALL DFUN(DF,WORK2,TABLE,HESS,OBS,PAR,PARMAX,WORK1,DP,F,ITER
     1,MAXPAR,NCASE,NPAR,NOBS,NPASS,DIFFER,FORWRD,MALE,TRTTYPE,VTRAITS
     2,VARDATA,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO
     3,PARMIN)
C                                                                        
C     IF THERE IS NOT A SUFFICIENT DECREASE IN F, THEN TRY TO            
C     FIND A BETTER POINT ALONG THE DIRECTION DELTA.  COMPUTE            
C     THE MINIMUM POINT FOR THE QUADRATIC IN T WHICH PASSES              
C     THROUGH FOLD AND F AND HAS SLOPE D AT T=0.  IF THIS                
C     MINIMUM IS TOO CLOSE TO 0, DECREMENT T BY ONLY 90 PER CENT.        
C                                                                        
C CPP adds: If we've fallen off the cliff, decrement T by BCLIFF
C

      IF ((F.EQ.0D0.OR.ABS(F).GT.1.0E+35.OR.0.NE.ISNAN_F(ABS(F)))
     1     .AND.NSTEP.LT.MXCLIFFS) THEN
         T=BCLIFF*T
         NSTEP=NSTEP+1
         GO TO 150
      ELSE IF (F.GT.FOLD+0.1D0*T*D.AND.NSTEP.LT.MXSTEP) THEN
         T1=-0.5D0*D*T*T/(F-FOLD-T*D)
         T2=0.1D0*T
         T=MAX(T1,T2)
         NSTEP=NSTEP+1
         GO TO 150
      END IF
C                                                                        
C     QUIT WHEN THERE IS A SUFFICIENT DECREASE IN F OR TOO               
C     MANY STEPS.  CHECK THE CONVERGENCE CRITERION.  IF IT HAS           
C     BEEN SATISFIED NCONV TIMES, THEN EXIT THE MAIN LOOP.               
C     OTHERWISE, OUTPUT THE CURRENT ITERATION.                           
C                                                                        
      IF (ABS(FOLD-F).GT.CONV) NCRIT=-1                                  
      NCRIT=NCRIT+1                                                      
      IF (NCRIT.GE.NCONV) GO TO 160                                      
      LOGLIK = -F
      CALL OUTPUT(PAR,PNAME,LOGLIK,QDFORM,ITER,LAST,NPAR
     1,NSTEP,IOUNIT,STAND,IBIG,CONOUT,DISCRETE,BIG)
C                                                                        
C     IF THIS IS ORDINARY MINIMIZATION, THEN UPDATE THE CURRENT          
C     APPROXIMATION TO THE HESSIAN.                                      
C                                                                        
      IF (.NOT.DIFFER(2)) THEN                                           
C                                                                        
C     RESET DELTA SO THAT IT IS THE ACTUAL STEP TAKEN.                   
C                                                                        
      DO 170 J=1,NPAR                                                    
 170  DELTA(J)=T*DELTA(J)                                                
C                                                                        
C     PREPARE TO UPDATE THE HESSIAN.  STORE IN WORK1 THE                 
C     PRODUCT HESS*DELTA.  WORK1 APPROXIMATES THE DIFFERENCE             
C     IN DIFFERENTIALS DF-DFOLD.  STORE IN C1 AND C2 THE                 
C     INNER PRODUCT OF THESE TWO VECTORS WITH DELTA.                     
C                                                                        
      C1=0.0D0                                                           
      DO 180 J=1,NPAR                                                    
      S=0.0D0                                                            
      DO 190 K=1,NPAR                                                    
 190  S=S+HESS(J,K)*DELTA(K)                                             
      WORK1(J)=S                                                         
 180  C1=C1+DELTA(J)*S                                                   
      C2=0.0D0                                                           
      DO 200 J=1,NPAR                                                    
 200  C2=C2+(DF(J)-DFOLD(J))*DELTA(J)                                    
C                                                                        
C     IF C2 IS TOO SMALL, BIAS DF-DFOLD BY TAKING A CONVEX               
C     COMBINATION WITH WORK1.  STORE THE RESUTING VECTOR                 
C     IN WORK2.  ITS INNER PRODUCT WITH DELTA WILL BE C4.                
C                                                                        
      IF (C1.GT.0.0D0) THEN                                              
      IF (C2.GT.0.2D0*C1) THEN                                           
      C3=1.0D0                                                           
      ELSE                                                               
      C3=0.8D0*C1/(C1-C2)                                                
      END IF                                                             
      DO 210 J=1,NPAR                                                    
 210  WORK2(J)=C3*(DF(J)-DFOLD(J))+(1.0D0-C3)*WORK1(J)                   
      C4=C3*C2+(1.0D0-C3)*C1                                             
C                                                                        
C     NOW RESET THE HESSIAN USING THE RANK TWO BFGS UPDATE.              
C                                                                        
      DO 220 J=1,NPAR                                                    
      DO 220 K=1,NPAR                                                    
 220  HESS(K,J)=HESS(K,J)-WORK1(J)*WORK1(K)/C1+WORK2(J)*WORK2(K)/C4      
      END IF                                                             
      END IF                                                             

C     CALL UPD(CNSTR,CVALUE,PAR,PARMAX,PARMIN
C    1,PNAME,NCNSTR,NPAR,NPOINT,NSAMP,MODFIL)

 80   CONTINUE                                                           
      IF (STDERR.GT.0) THEN                                              
      DO 230 J=1,NPAR                                                    
      DO 230 I=1,NPAR                                                    
 230  HESS(I,J)=0.0D0                                                    
      END IF                                                             
      RETURN                                                             
      END IF                                                             
C                                                                        
C     CONVERGENCE HAS OCCURRED.  OUTPUT THE LAST ITERATION               
C                                                                        
 160  IF (ITER.NE.MXITER) THEN                                           
      LOGLIK = -F
      CALL OUTPUT(PAR,PNAME,LOGLIK,QDFORM,ITER,ITER
     1,NPAR,NSTEP,IOUNIT,STAND,IBIG,CONOUT,DISCRETE,BIG)
      END IF                                                             
C                                                                        
C     IF THE ASYMPTOTIC COVARIANCE MATRIX IS DESIRED, THEN RECOMPUTE     
C     THE HESSIAN AND CALL ASYCOV.  USE CENTRAL DIFFERENCES FOR THE      
C     FIRST PARTIALS AND FORWARD DIFFERENCES FOR THE SECOND PARTIALS.    
C     ADJUST THE DIFFERENTIATION INTERVAL FOR THE SECOND PARTIALS.       
C                                                                        
      IF (STDERR.EQ.1) THEN                                              
      IF (.NOT.DIFFER(1).AND.FORWRD) CALL DFUN(DF,WORK2,TABLE,HESS,OBS
     1,PAR,PARMAX,WORK1,DP,F,ITER,MAXPAR,NCASE,NPAR,NOBS,NPASS
     2,DIFFER,.FALSE.,MALE,TRTTYPE,VTRAITS
     3,VARDATA,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO
     4,PARMIN)
      IF (.NOT.DIFFER(1)) THEN                                           
      DP23=DP**0.66667D0                                                 
      ELSE                                                               
      DP23=DP                                                            
      END IF                                                             
C                                                                        
C     COMPUTE THE SECOND PARTIALS AND CALL ASYCOV.  NBOUND IS THE        
C     NUMBER OF PARAMETERS ON A BOUNDARY.  NOTE THAT IT IS               
C     UNNECESSARY TO COMPUTE THE HESSIAN AT THIS STAGE IF IT IS          
C     ALREADY GIVEN.                                                     
C                                                                        
      DIFF(1)=DIFFER(1)                                                  
      DIFF(2)=.FALSE.                                                    
      NBOUND=0                                                           
      DO 240 J=1,NPAR                                                    
      IF (PAR(J).LE.PARMIN(J)+SMALL.OR.PAR(J).GE.PARMAX(J)-SMALL) THEN   
      NBOUND=NBOUND+1                                                    
      NSWEEP(J)=0                                                        
      ELSE                                                               
      NSWEEP(J)=1                                                        
      DPJ=DP23*MAX(ABS(PAR(J)),1.0D0)                                    
      PAR(J)=PAR(J)+DPJ                                                  
      CALL DFUN(DFOLD,WORK2,TABLE,HESS,OBS,PAR,PARMAX,WORK1,DP,FOLD
     1,ITER,MAXPAR,NCASE,J,NOBS,NPASS,DIFF,FORWRD,MALE,TRTTYPE,VTRAITS
     2,VARDATA,NVAR,NTOT,NIND,NASCER,NPED,NCUMIND,COV,MU,AFFECT,MAXPEO
     3,PARMIN)
      PAR(J)=PAR(J)-DPJ                                                  
      DO 250 I=1,J                                                       
      IF (NSWEEP(I).EQ.0.OR.NSWEEP(J).EQ.0) THEN                         
      HESS(I,J)=0.0D0                                                    
      ELSE                                                               
      HESS(I,J)=(DFOLD(I)-DF(I))/DPJ                                     
      END IF                                                             
 250  HESS(J,I)=HESS(I,J)                                                
      END IF                                                             
 240  CONTINUE                                                           
      END IF                                                             
C                                                                        
C     ADJUST THE HESSIAN FOR A LEAST SQUARES PROBLEM BY DIVIDING         
C     BY THE RESIDUAL MEAN SQUARE.                                       
C                                                                        
      IF (STDERR.EQ.3) THEN                                              
      SIGMA=MAX(2.0D0*F/MAX(NCASE-NPAR+NBOUND+NCNSTR,1),1.0D-10)         
      DO 260 I=1,NPAR                                                    
      DO 260 J=1,NPAR                                                    
 260  HESS(J,I)=HESS(J,I)/SIGMA                                          
      END IF                                                             
      if(STDERR.GT.0)INASYCV=.FALSE.
      IF (STDERR.GT.0) CALL DASYCOV(CNSTR,CVALUE,DF,HESS,PAR,TABLE        
     :,WORK1,WORK2,NSWEEP,PNAME,CNORM,SMALL,TOL,IOUNIT,MAXPAR,MAXTAB     
     :,MCNSTR,NCNSTR,NPAR,INASYCV)
      END IF                                                             
      END                                                                

