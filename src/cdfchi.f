      SUBROUTINE cdfchi(which,p,q,x,df,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFCHI( WHICH, P, Q, X, DF, STATUS, BOUND )
C               Cumulative Distribution Function
C               CHI-Square distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the chi-square
C     distribution given values for the others.
C
C
C                              Arguments
C
C
C     WHICH --> Integer indicating which of the next three argument
C               values is to be calculated from the others.
C               Legal range: 1..3
C               iwhich = 1 : Calculate P and Q from X and DF
C               iwhich = 2 : Calculate X from P,Q and DF
C               iwhich = 3 : Calculate DF from P,Q and X
C                    INTEGER WHICH
C
C     P <--> The integral from 0 to X of the chi-square
C            distribution.
C            Input range: [0, 1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C     X <--> Upper limit of integration of the non-central
C            chi-square distribution.
C            Input range: [0, +infinity).
C            Search range: [0,1E100]
C                    DOUBLE PRECISION X
C
C     DF <--> Degrees of freedom of the
C             chi-square distribution.
C             Input range: (0, +infinity).
C             Search range: [ 1E-100, 1E100]
C                    DOUBLE PRECISION DF
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
C                3 if P + Q .ne. 1
C               10 indicates error returned from cumgam.  See
C                  references in cdfgam
C                    INTEGER STATUS
C
C     BOUND <-- Undefined if STATUS is 0
C
C               Bound exceeded by parameter number I if STATUS
C               is negative.
C
C               Lower search bound if STATUS is 1.
C
C               Upper search bound if STATUS is 2.
C
C
C                              Method
C
C
C     Formula    26.4.19   of Abramowitz  and     Stegun, Handbook  of
C     Mathematical Functions   (1966) is used   to reduce the chisqure
C     distribution to the incomplete distribution.
C
C     Computation of other parameters involve a seach for a value that
C     produces  the desired  value  of P.   The search relies  on  the
C     monotinicity of P with the other parameter.
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION tol
      PARAMETER (tol=1.0D-8)
      DOUBLE PRECISION atol
      PARAMETER (atol=1.0D-50)
      DOUBLE PRECISION zero,inf
      PARAMETER (zero=1.0D-100,inf=1.0D100)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,df,p,q,x
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ccum,cum,fx,porq,pq
      LOGICAL qhi,qleft,qporq
C     ..
C     .. External Functions ..
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL cumchi,dinvr,dstinv
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      IF (.NOT. ((which.LT.1).OR. (which.GT.3))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 3.0D0
   20 status = -1
C      RETURN
      GO TO 666

   30 IF (which.EQ.1) GO TO 70
      IF (.NOT. ((p.LT.0.0D0).OR. (p.GT.1.0D0))) GO TO 60
      IF (.NOT. (p.LT.0.0D0)) GO TO 40
      bound = 0.0D0
      GO TO 50

   40 bound = 1.0D0
   50 status = -2
C      RETURN
      GO TO 666

   60 CONTINUE
   70 IF (which.EQ.1) GO TO 110
      IF (.NOT. ((q.LE.0.0D0).OR. (q.GT.1.0D0))) GO TO 100
      IF (.NOT. (q.LE.0.0D0)) GO TO 80
      bound = 0.0D0
      GO TO 90

   80 bound = 1.0D0
   90 status = -3
C      RETURN
      GO TO 666

  100 CONTINUE
  110 IF (which.EQ.2) GO TO 130
      IF (.NOT. (x.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      status = -4
C      RETURN
      GO TO 666

  120 CONTINUE
  130 IF (which.EQ.3) GO TO 150
      IF (.NOT. (df.LE.0.0D0)) GO TO 140
      bound = 0.0D0
      status = -5
C      RETURN
      GO TO 666

  140 CONTINUE
  150 IF (which.EQ.1) GO TO 190
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 180
      IF (.NOT. (pq.LT.0.0D0)) GO TO 160
      bound = 0.0D0
      GO TO 170

  160 bound = 1.0D0
  170 status = 3
C      RETURN
      GO TO 666

  180 CONTINUE
  190 IF (which.EQ.1) GO TO 220
      qporq = p .LE. q
      IF (.NOT. (qporq)) GO TO 200
      porq = p
      GO TO 210

  200 porq = q
  210 CONTINUE
  220 IF ((1).EQ. (which)) THEN
          status = 0
          CALL cumchi(x,df,p,q)
          IF (porq.GT.1.5D0) THEN
              status = 10
C              RETURN
              GO TO 666

          END IF

      ELSE IF ((2).EQ. (which)) THEN
          x = 5.0D0
          CALL dstinv(0.0D0,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,x,fx,qleft,qhi)
  230     IF (.NOT. (status.EQ.1)) GO TO 270
          CALL cumchi(x,df,cum,ccum)
          IF (.NOT. (qporq)) GO TO 240
          fx = cum - p
          GO TO 250

  240     fx = ccum - q
  250     IF (.NOT. ((fx+porq).GT.1.5D0)) GO TO 260
          status = 10
C          RETURN
          GO TO 666

  260     CALL dinvr(status,x,fx,qleft,qhi)
          GO TO 230

  270     IF (.NOT. (status.EQ.-1)) GO TO 300
          IF (.NOT. (qleft)) GO TO 280
          status = 1
          bound = 0.0D0
C          GO TO 290
          GO TO 666

  280     status = 2
          bound = inf
          GO TO 666
C  290     CONTINUE
  300     CONTINUE

      ELSE IF ((3).EQ. (which)) THEN
          df = 5.0D0
          CALL dstinv(zero,inf,0.5D0,0.5D0,5.0D0,atol,tol)
          status = 0
          CALL dinvr(status,df,fx,qleft,qhi)
  310     IF (.NOT. (status.EQ.1)) GO TO 350
          CALL cumchi(x,df,cum,ccum)
          IF (.NOT. (qporq)) GO TO 320
          fx = cum - p
          GO TO 330

  320     fx = ccum - q
  330     IF (.NOT. ((fx+porq).GT.1.5D0)) GO TO 340
          status = 10
C          RETURN
          GO TO 666

  340     CALL dinvr(status,df,fx,qleft,qhi)
          GO TO 310

  350     IF (.NOT. (status.EQ.-1)) GO TO 370
          IF (.NOT. (qleft)) GO TO 360
          status = 1
          bound = zero
C          GO TO 370
          GO TO 666

  360     status = 2
          bound = inf
          GO TO 666
 370      CONTINUE
      END IF

      RETURN

  666 CONTINUE
      RETURN

      END
      SUBROUTINE cumchi(x,df,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE FUNCTION CUMCHI(X,DF,CUM,CCUM)
C             CUMulative of the CHi-square distribution
C
C
C                              Function
C
C
C     Calculates the cumulative chi-square distribution.
C
C
C                              Arguments
C
C
C     X       --> Upper limit of integration of the
C                 chi-square distribution.
C                                                 X is DOUBLE PRECISION
C
C     DF      --> Degrees of freedom of the
C                 chi-square distribution.
C                                                 DF is DOUBLE PRECISION
C
C     CUM <-- Cumulative chi-square distribution.
C                                                 CUM is DOUBLE PRECISIO
C
C     CCUM <-- Compliment of Cumulative chi-square distribution.
C                                                 CCUM is DOUBLE PRECISI
C
C
C                              Method
C
C
C     Calls incomplete gamma function (CUMGAM)
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION df,x,cum,ccum
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,xx
C     ..
C     .. External Subroutines ..
      EXTERNAL cumgam
C     ..
C     .. Executable Statements ..
      a = df*0.5D0
      xx = x*0.5D0
      CALL cumgam(xx,a,cum,ccum)
      RETURN

      END
      SUBROUTINE cumgam(x,a,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE CUMGAM(X,A,CUM,CCUM)
C           Double precision cUMulative incomplete GAMma distribution
C
C
C                              Function
C
C
C     Computes   the  cumulative        of    the     incomplete   gamma
C     distribution, i.e., the integral from 0 to X of
C          (1/GAM(A))*EXP(-T)*T**(A-1) DT
C     where GAM(A) is the complete gamma function of A, i.e.,
C          GAM(A) = integral from 0 to infinity of
C                    EXP(-T)*T**(A-1) DT
C
C
C                              Arguments
C
C
C     X --> The upper limit of integration of the incomplete gamma.
C                                                X is DOUBLE PRECISION
C
C     A --> The shape parameter of the incomplete gamma.
C                                                A is DOUBLE PRECISION
C
C     CUM <-- Cumulative incomplete gamma distribution.
C                                        CUM is DOUBLE PRECISION
C
C     CCUM <-- Compliment of Cumulative incomplete gamma distribution.
C                                                CCUM is DOUBLE PRECISIO
C
C
C                              Method
C
C
C     Calls the routine GRATIO.
C
C**********************************************************************
C
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,x,cum,ccum
C     ..
C     .. External Routines ..
      EXTERNAL gratio
C     ..
C     .. Executable Statements ..
      IF (.NOT. (x.LE.0.0D0)) GO TO 10
      cum = 0.0D0
      ccum = 1.0D0
      RETURN

   10 CALL gratio(a,x,cum,ccum,0)

C     Call gratio routine

      RETURN

      END
      SUBROUTINE dinvr(status,x,fx,qleft,qhi)
C**********************************************************************
C
C     SUBROUTINE DINVR(STATUS, X, FX, QLEFT, QHI)
C          Double precision
C          bounds the zero of the function and invokes zror
C                    Reverse Communication
C
C
C                              Function
C
C
C     Bounds the    function  and  invokes  ZROR   to perform the   zero
C     finding.  STINVR  must  have   been  called  before this   routine
C     in order to set its parameters.
C
C
C                              Arguments
C
C
C     STATUS <--> At the beginning of a zero finding problem, STATUS
C                 should be set to 0 and INVR invoked.  (The value
C                 of parameters other than X will be ignored on this cal
C
C                 When INVR needs the function evaluated, it will set
C                 STATUS to 1 and return.  The value of the function
C                 should be set in FX and INVR again called without
C                 changing any of its other parameters.
C
C                 When INVR has finished without error, it will return
C                 with STATUS 0.  In that case X is approximately a root
C                 of F(X).
C
C                 If INVR cannot bound the function, it returns status
C                 -1 and sets QLEFT and QHI.
C                         INTEGER STATUS
C
C     X <-- The value of X at which F(X) is to be evaluated.
C                         DOUBLE PRECISION X
C
C     FX --> The value of F(X) calculated when INVR returns with
C            STATUS = 1.
C                         DOUBLE PRECISION FX
C
C     QLEFT <-- Defined only if QMFINV returns .FALSE.  In that
C          case it is .TRUE. If the stepping search terminated
C          unsucessfully at SMALL.  If it is .FALSE. the search
C          terminated unsucessfully at BIG.
C                    QLEFT is LOGICAL
C
C     QHI <-- Defined only if QMFINV returns .FALSE.  In that
C          case it is .TRUE. if F(X) .GT. Y at the termination
C          of the search and .FALSE. if F(X) .LT. Y at the
C          termination of the search.
C                    QHI is LOGICAL

C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION fx,x,zabsst,zabsto,zbig,zrelst,zrelto,zsmall,
     +                 zstpmu
      INTEGER status
      LOGICAL qhi,qleft
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION absstp,abstol,big,fbig,fsmall,relstp,reltol,
     +                 small,step,stpmul,xhi,xlb,xlo,xsave,xub,yy,zx,zy,
     +                 zz
C     INTEGER i99999
      integer ijump
      LOGICAL qbdd,qcond,qdum1,qdum2,qincr,qlim,qok,qup
C     ..
C     .. External Subroutines ..
      EXTERNAL dstzr,dzror
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,min
C     ..
C     .. Statement Functions ..
      LOGICAL qxmon
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Statement Function definitions ..
      qxmon(zx,zy,zz) = zx .LE. zy .AND. zy .LE. zz
C     ..
C     .. Executable Statements ..

      IF (status.GT.0) GO TO 310

      qcond = .NOT. qxmon(small,x,big)
      IF (qcond) STOP ' SMALL, X, BIG not monotone in INVR'
      xsave = x
C
C     See that SMALL and BIG bound the zero and set QINCR
C
      x = small
C     GET-FUNCTION-VALUE
C     ASSIGN 10 TO i99999
      ijump = 10
      GO TO 300

   10 fsmall = fx
      x = big
C     GET-FUNCTION-VALUE
C     ASSIGN 20 TO i99999
      ijump = 20
      GO TO 300

   20 fbig = fx
      qincr = fbig .GT. fsmall
      IF (.NOT. (qincr)) GO TO 50
      IF (fsmall.LE.0.0D0) GO TO 30
      status = -1
      qleft = .TRUE.
      qhi = .TRUE.
      RETURN

   30 IF (fbig.GE.0.0D0) GO TO 40
      status = -1
      qleft = .FALSE.
      qhi = .FALSE.
      RETURN

   40 GO TO 80

   50 IF (fsmall.GE.0.0D0) GO TO 60
      status = -1
      qleft = .TRUE.
      qhi = .FALSE.
      RETURN

   60 IF (fbig.LE.0.0D0) GO TO 70
      status = -1
      qleft = .FALSE.
      qhi = .TRUE.
      RETURN

   70 CONTINUE
   80 x = xsave
      step = max(absstp,relstp*abs(x))
C      YY = F(X) - Y
C     GET-FUNCTION-VALUE
C     ASSIGN 90 TO i99999
      ijump = 90
      GO TO 300

   90 yy = fx
      IF (.NOT. (yy.EQ.0.0D0)) GO TO 100
      status = 0
      qok = .TRUE.
      RETURN

  100 qup = (qincr .AND. (yy.LT.0.0D0)) .OR.
     +      (.NOT.qincr .AND. (yy.GT.0.0D0))
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     HANDLE CASE IN WHICH WE MUST STEP HIGHER
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (.NOT. (qup)) GO TO 170
      xlb = xsave
      xub = min(xlb+step,big)
      GO TO 120

  110 IF (qcond) GO TO 150
C      YY = F(XUB) - Y
  120 x = xub
C     GET-FUNCTION-VALUE
C     ASSIGN 130 TO i99999
      ijump = 130
      GO TO 300

  130 yy = fx
      qbdd = (qincr .AND. (yy.GE.0.0D0)) .OR.
     +       (.NOT.qincr .AND. (yy.LE.0.0D0))
      qlim = xub .GE. big
      qcond = qbdd .OR. qlim
      IF (qcond) GO TO 140
      step = stpmul*step
      xlb = xub
      xub = min(xlb+step,big)
  140 GO TO 110

  150 IF (.NOT. (qlim.AND..NOT.qbdd)) GO TO 160
      status = -1
      qleft = .FALSE.
      qhi = .NOT. qincr
      x = big
      RETURN

  160 GO TO 240
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     HANDLE CASE IN WHICH WE MUST STEP LOWER
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  170 xub = xsave
      xlb = max(xub-step,small)
      GO TO 190

  180 IF (qcond) GO TO 220
C      YY = F(XLB) - Y
  190 x = xlb
C     GET-FUNCTION-VALUE
C     ASSIGN 200 TO i99999
      ijump = 200
      GO TO 300

  200 yy = fx
      qbdd = (qincr .AND. (yy.LE.0.0D0)) .OR.
     +       (.NOT.qincr .AND. (yy.GE.0.0D0))
      qlim = xlb .LE. small
      qcond = qbdd .OR. qlim
      IF (qcond) GO TO 210
      step = stpmul*step
      xub = xlb
      xlb = max(xub-step,small)
  210 GO TO 180

  220 IF (.NOT. (qlim.AND..NOT.qbdd)) GO TO 230
      status = -1
      qleft = .TRUE.
      qhi = qincr
      x = small
      RETURN

  230 CONTINUE
  240 CALL dstzr(xlb,xub,abstol,reltol)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      status = 0
      GO TO 260

  250 IF (.NOT. (status.EQ.1)) GO TO 290
  260 CALL dzror(status,x,fx,xlo,xhi,qdum1,qdum2)
      IF (.NOT. (status.EQ.1)) GO TO 280
C     GET-FUNCTION-VALUE
C     ASSIGN 270 TO i99999
      ijump = 270
      GO TO 300

  270 CONTINUE
  280 GO TO 250

  290 x = xlo
      status = 0
      RETURN

      ENTRY dstinv(zsmall,zbig,zabsst,zrelst,zstpmu,zabsto,zrelto)
C**********************************************************************
C
C      SUBROUTINE DSTINV( SMALL, BIG, ABSSTP, RELSTP, STPMUL,
C     +                   ABSTOL, RELTOL )
C      Double Precision - SeT INverse finder - Reverse Communication
C
C
C                              Function
C
C
C     Concise Description - Given a monotone function F finds X
C     such that F(X) = Y.  Uses Reverse communication -- see invr.
C     This routine sets quantities needed by INVR.
C
C          More Precise Description of INVR -
C
C     F must be a monotone function, the results of QMFINV are
C     otherwise undefined.  QINCR must be .TRUE. if F is non-
C     decreasing and .FALSE. if F is non-increasing.
C
C     QMFINV will return .TRUE. if and only if F(SMALL) and
C     F(BIG) bracket Y, i. e.,
C          QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
C          QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
C
C     if QMFINV returns .TRUE., then the X returned satisfies
C     the following condition.  let
C               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
C     then if QINCR is .TRUE.,
C          F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
C     and if QINCR is .FALSE.
C          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
C
C
C                              Arguments
C
C
C     SMALL --> The left endpoint of the interval to be
C          searched for a solution.
C                    SMALL is DOUBLE PRECISION
C
C     BIG --> The right endpoint of the interval to be
C          searched for a solution.
C                    BIG is DOUBLE PRECISION
C
C     ABSSTP, RELSTP --> The initial step size in the search
C          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
C                    ABSSTP is DOUBLE PRECISION
C                    RELSTP is DOUBLE PRECISION
C
C     STPMUL --> When a step doesn't bound the zero, the step
C                size is multiplied by STPMUL and another step
C                taken.  A popular value is 2.0
C                    DOUBLE PRECISION STPMUL
C
C     ABSTOL, RELTOL --> Two numbers that determine the accuracy
C          of the solution.  See function for a precise definition.
C                    ABSTOL is DOUBLE PRECISION
C                    RELTOL is DOUBLE PRECISION
C
C
C                              Method
C
C
C     Compares F(X) with Y for the input value of X then uses QINCR
C     to determine whether to step left or right to bound the
C     desired x.  the initial step size is
C          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
C     Iteratively steps right or left until it bounds X.
C     At each step which doesn't bound X, the step size is doubled.
C     The routine is careful never to step beyond SMALL or BIG.  If
C     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
C     after setting QLEFT and QHI.
C
C     If X is successfully bounded then Algorithm R of the paper
C     'Two Efficient Algorithms with Guaranteed Convergence for
C     Finding a Zero of a Function' by J. C. P. Bus and
C     T. J. Dekker in ACM Transactions on Mathematical
C     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
C     to find the zero of the function F(X)-Y. This is routine
C     QRZERO.
C
C**********************************************************************
      small = zsmall
      big = zbig
      absstp = zabsst
      relstp = zrelst
      stpmul = zstpmu
      abstol = zabsto
      reltol = zrelto
      RETURN

C     STOP '*** EXECUTION FLOWING INTO FLECS PROCEDURES ***'
C     TO GET-FUNCTION-VALUE
  300 status = 1
      RETURN

  310 CONTINUE
C      GO TO i99999
      if (ijump.eq.10) then
         go to 10
      else if (ijump.eq.20) then
         go to 20
      else if (ijump.eq.90) then
         go to 90
      else if (ijump.eq.130) then
         go to 130
      else if (ijump.eq.200) then
         go to 200
      else if (ijump.eq.270) then
         go to 270
      else
         print *, 'Programming error!  missing branch ',ijump
         call exit (-1)
      end if

      END
      SUBROUTINE dzror(status,x,fx,xlo,xhi,qleft,qhi)
C**********************************************************************
C
C     SUBROUTINE DZROR(STATUS, X, FX, XLO, XHI, QLEFT, QHI)
C     Double precision ZeRo of a function -- Reverse Communication
C
C
C                              Function
C
C
C     Performs the zero finding.  STZROR must have been called before
C     this routine in order to set its parameters.
C
C
C                              Arguments
C
C
C     STATUS <--> At the beginning of a zero finding problem, STATUS
C                 should be set to 0 and ZROR invoked.  (The value
C                 of other parameters will be ignored on this call.)
C
C                 When ZROR needs the function evaluated, it will set
C                 STATUS to 1 and return.  The value of the function
C                 should be set in FX and ZROR again called without
C                 changing any of its other parameters.
C
C                 When ZROR has finished without error, it will return
C                 with STATUS 0.  In that case (XLO,XHI) bound the answe
C
C                 If ZROR finds an error (which implies that F(XLO)-Y an
C                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
C                 this case, XLO and XHI are undefined.
C                         INTEGER STATUS
C
C     X <-- The value of X at which F(X) is to be evaluated.
C                         DOUBLE PRECISION X
C
C     FX --> The value of F(X) calculated when ZROR returns with
C            STATUS = 1.
C                         DOUBLE PRECISION FX
C
C     XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
C             inverval in X containing the solution below.
C                         DOUBLE PRECISION XLO
C
C     XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
C             inverval in X containing the solution above.
C                         DOUBLE PRECISION XHI
C
C     QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
C                at XLO.  If it is .FALSE. the search terminated
C                unsucessfully at XHI.
C                    QLEFT is LOGICAL
C
C     QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
C              search and .FALSE. if F(X) .LT. Y at the
C              termination of the search.
C                    QHI is LOGICAL
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION fx,x,xhi,xlo,zabstl,zreltl,zxhi,zxlo
      INTEGER status
      LOGICAL qhi,qleft
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,abstol,b,c,d,fa,fb,fc,fd,fda,fdb,m,mb,p,q,
     +                 reltol,tol,w,xxhi,xxlo,zx
      INTEGER ext
C     INTEGER i99999
      integer ijump
      LOGICAL first,qrzero
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,max,sign
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION ftol
C     ..
C     .. Statement Function definitions ..
      ftol(zx) = 0.5D0*max(abstol,reltol*abs(zx))
C     ..
C     .. Executable Statements ..

      IF (status.GT.0) GO TO 280
      xlo = xxlo
      xhi = xxhi
      b = xlo
      x = xlo
C     GET-FUNCTION-VALUE
C     ASSIGN 10 TO i99999
      ijump = 10
      GO TO 270

   10 fb = fx
      xlo = xhi
      a = xlo
      x = xlo
C     GET-FUNCTION-VALUE
C     ASSIGN 20 TO i99999
      ijump = 20
      GO TO 270
C
C     Check that F(ZXLO) < 0 < F(ZXHI)  or
C                F(ZXLO) > 0 > F(ZXHI)
C
   20 IF (.NOT. (fb.LT.0.0D0)) GO TO 40
      IF (.NOT. (fx.LT.0.0D0)) GO TO 30
      status = -1
      qleft = fx .LT. fb
      qhi = .FALSE.
      RETURN

   30 CONTINUE
   40 IF (.NOT. (fb.GT.0.0D0)) GO TO 60
      IF (.NOT. (fx.GT.0.0D0)) GO TO 50
      status = -1
      qleft = fx .GT. fb
      qhi = .TRUE.
      RETURN

   50 CONTINUE
   60 fa = fx
C
      first = .TRUE.
   70 c = a
      fc = fa
      ext = 0
   80 IF (.NOT. (abs(fc).LT.abs(fb))) GO TO 100
      IF (.NOT. (c.NE.a)) GO TO 90
      d = a
      fd = fa
   90 a = b
      fa = fb
      xlo = c
      b = xlo
      fb = fc
      c = a
      fc = fa
  100 tol = ftol(xlo)
      m = (c+b)*.5D0
      mb = m - b
      IF (.NOT. (abs(mb).GT.tol)) GO TO 240
      IF (.NOT. (ext.GT.3)) GO TO 110
      w = mb
      GO TO 190

  110 tol = sign(tol,mb)
      p = (b-a)*fb
      IF (.NOT. (first)) GO TO 120
      q = fa - fb
      first = .FALSE.
      GO TO 130

  120 fdb = (fd-fb)/ (d-b)
      fda = (fd-fa)/ (d-a)
      p = fda*p
      q = fdb*fa - fda*fb
  130 IF (.NOT. (p.LT.0.0D0)) GO TO 140
      p = -p
      q = -q
  140 IF (ext.EQ.3) p = p*2.0D0
      IF (.NOT. ((p*1.0D0).EQ.0.0D0.OR.p.LE. (q*tol))) GO TO 150
      w = tol
      GO TO 180

  150 IF (.NOT. (p.LT. (mb*q))) GO TO 160
      w = p/q
      GO TO 170

  160 w = mb
  170 CONTINUE
  180 CONTINUE
  190 d = a
      fd = fa
      a = b
      fa = fb
      b = b + w
      xlo = b
      x = xlo
C     GET-FUNCTION-VALUE
C     ASSIGN 200 TO i99999
      ijump = 200
      GO TO 270

  200 fb = fx
      IF (.NOT. ((fc*fb).GE.0.0D0)) GO TO 210
      GO TO 70

  210 IF (.NOT. (w.EQ.mb)) GO TO 220
      ext = 0
      GO TO 230

  220 ext = ext + 1
  230 GO TO 80

  240 xhi = c
      qrzero = (fc.GE.0.0D0 .AND. fb.LE.0.0D0) .OR.
     +         (fc.LT.0.0D0 .AND. fb.GE.0.0D0)
      IF (.NOT. (qrzero)) GO TO 250
      status = 0
      GO TO 260

  250 status = -1
  260 RETURN

      ENTRY dstzr(zxlo,zxhi,zabstl,zreltl)
C**********************************************************************
C
C     SUBROUTINE DSTZR( XLO, XHI, ABSTOL, RELTOL )
C     Double precision SeT ZeRo finder - Reverse communication version
C
C
C                              Function
C
C
C
C     Sets quantities needed by ZROR.  The function of ZROR
C     and the quantities set is given here.
C
C     Concise Description - Given a function F
C     find XLO such that F(XLO) = 0.
C
C          More Precise Description -
C
C     Input condition. F is a double precision function of a single
C     double precision argument and XLO and XHI are such that
C          F(XLO)*F(XHI)  .LE.  0.0
C
C     If the input condition is met, QRZERO returns .TRUE.
C     and output values of XLO and XHI satisfy the following
C          F(XLO)*F(XHI)  .LE. 0.
C          ABS(F(XLO)  .LE. ABS(F(XHI)
C          ABS(XLO-XHI)  .LE. TOL(X)
C     where
C          TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
C
C     If this algorithm does not find XLO and XHI satisfying
C     these conditions then QRZERO returns .FALSE.  This
C     implies that the input condition was not met.
C
C
C                              Arguments
C
C
C     XLO --> The left endpoint of the interval to be
C           searched for a solution.
C                    XLO is DOUBLE PRECISION
C
C     XHI --> The right endpoint of the interval to be
C           for a solution.
C                    XHI is DOUBLE PRECISION
C
C     ABSTOL, RELTOL --> Two numbers that determine the accuracy
C                      of the solution.  See function for a
C                      precise definition.
C                    ABSTOL is DOUBLE PRECISION
C                    RELTOL is DOUBLE PRECISION
C
C
C                              Method
C
C
C     Algorithm R of the paper 'Two Efficient Algorithms with
C     Guaranteed Convergence for Finding a Zero of a Function'
C     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
C     Mathematical Software, Volume 1, no. 4 page 330
C     (Dec. '75) is employed to find the zero of F(X)-Y.
C
C**********************************************************************
      xxlo = zxlo
      xxhi = zxhi
      abstol = zabstl
      reltol = zreltl
      RETURN

C     STOP '*** EXECUTION FLOWING INTO FLECS PROCEDURES ***'
C     TO GET-FUNCTION-VALUE
  270 status = 1
      RETURN

  280 CONTINUE
C     GO TO i99999
      if (ijump.eq.10) then
         go to 10
      else if (ijump.eq.20) then
         go to 20
      else if (ijump.eq.200) then
         go to 200
      else
         print *, 'Programming error!  missing branch ',ijump
     1,' in dstzr'
         call exit (-1)
      end if
      

      END
      DOUBLE PRECISION FUNCTION erf(x)
C-----------------------------------------------------------------------
C             EVALUATION OF THE REAL ERROR FUNCTION
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ax,bot,c,t,top,x2
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a(5),b(3),p(8),q(8),r(5),s(4)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,sign
C     ..
C     .. Data statements ..
C-------------------------
C-------------------------
C-------------------------
C-------------------------
      DATA c/.564189583547756D0/
      DATA a(1)/.771058495001320D-04/,a(2)/-.133733772997339D-02/,
     +     a(3)/.323076579225834D-01/,a(4)/.479137145607681D-01/,
     +     a(5)/.128379167095513D+00/
      DATA b(1)/.301048631703895D-02/,b(2)/.538971687740286D-01/,
     +     b(3)/.375795757275549D+00/
      DATA p(1)/-1.36864857382717D-07/,p(2)/5.64195517478974D-01/,
     +     p(3)/7.21175825088309D+00/,p(4)/4.31622272220567D+01/,
     +     p(5)/1.52989285046940D+02/,p(6)/3.39320816734344D+02/,
     +     p(7)/4.51918953711873D+02/,p(8)/3.00459261020162D+02/
      DATA q(1)/1.00000000000000D+00/,q(2)/1.27827273196294D+01/,
     +     q(3)/7.70001529352295D+01/,q(4)/2.77585444743988D+02/,
     +     q(5)/6.38980264465631D+02/,q(6)/9.31354094850610D+02/,
     +     q(7)/7.90950925327898D+02/,q(8)/3.00459260956983D+02/
      DATA r(1)/2.10144126479064D+00/,r(2)/2.62370141675169D+01/,
     +     r(3)/2.13688200555087D+01/,r(4)/4.65807828718470D+00/,
     +     r(5)/2.82094791773523D-01/
      DATA s(1)/9.41537750555460D+01/,s(2)/1.87114811799590D+02/,
     +     s(3)/9.90191814623914D+01/,s(4)/1.80124575948747D+01/
C     ..
C     .. Executable Statements ..
C-------------------------
      ax = abs(x)
      IF (ax.GT.0.5D0) GO TO 10
      t = x*x
      top = ((((a(1)*t+a(2))*t+a(3))*t+a(4))*t+a(5)) + 1.0D0
      bot = ((b(1)*t+b(2))*t+b(3))*t + 1.0D0
      erf = x* (top/bot)
      RETURN
C
   10 IF (ax.GT.4.0D0) GO TO 20
      top = ((((((p(1)*ax+p(2))*ax+p(3))*ax+p(4))*ax+p(5))*ax+p(6))*ax+
     +      p(7))*ax + p(8)
      bot = ((((((q(1)*ax+q(2))*ax+q(3))*ax+q(4))*ax+q(5))*ax+q(6))*ax+
     +      q(7))*ax + q(8)
      erf = 0.5D0 + (0.5D0-exp(-x*x)*top/bot)
      IF (x.LT.0.0D0) erf = -erf
      RETURN
C
   20 IF (ax.GE.5.8D0) GO TO 30
      x2 = x*x
      t = 1.0D0/x2
      top = (((r(1)*t+r(2))*t+r(3))*t+r(4))*t + r(5)
      bot = (((s(1)*t+s(2))*t+s(3))*t+s(4))*t + 1.0D0
      erf = (c-top/ (x2*bot))/ax
      erf = 0.5D0 + (0.5D0-exp(-x2)*erf)
      IF (x.LT.0.0D0) erf = -erf
      RETURN
C
   30 erf = sign(1.0D0,x)
      RETURN

      END
      DOUBLE PRECISION FUNCTION erfc1(ind,x)
C-----------------------------------------------------------------------
C         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
C
C          ERFC1(IND,X) = ERFC(X)            IF IND = 0
C          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER ind
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ax,bot,c,e,t,top,w
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a(5),b(3),p(8),q(8),r(5),s(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION exparg
      EXTERNAL exparg
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,exp
C     ..
C     .. Data statements ..
C-------------------------
C-------------------------
C-------------------------
C-------------------------
      DATA c/.564189583547756D0/
      DATA a(1)/.771058495001320D-04/,a(2)/-.133733772997339D-02/,
     +     a(3)/.323076579225834D-01/,a(4)/.479137145607681D-01/,
     +     a(5)/.128379167095513D+00/
      DATA b(1)/.301048631703895D-02/,b(2)/.538971687740286D-01/,
     +     b(3)/.375795757275549D+00/
      DATA p(1)/-1.36864857382717D-07/,p(2)/5.64195517478974D-01/,
     +     p(3)/7.21175825088309D+00/,p(4)/4.31622272220567D+01/,
     +     p(5)/1.52989285046940D+02/,p(6)/3.39320816734344D+02/,
     +     p(7)/4.51918953711873D+02/,p(8)/3.00459261020162D+02/
      DATA q(1)/1.00000000000000D+00/,q(2)/1.27827273196294D+01/,
     +     q(3)/7.70001529352295D+01/,q(4)/2.77585444743988D+02/,
     +     q(5)/6.38980264465631D+02/,q(6)/9.31354094850610D+02/,
     +     q(7)/7.90950925327898D+02/,q(8)/3.00459260956983D+02/
      DATA r(1)/2.10144126479064D+00/,r(2)/2.62370141675169D+01/,
     +     r(3)/2.13688200555087D+01/,r(4)/4.65807828718470D+00/,
     +     r(5)/2.82094791773523D-01/
      DATA s(1)/9.41537750555460D+01/,s(2)/1.87114811799590D+02/,
     +     s(3)/9.90191814623914D+01/,s(4)/1.80124575948747D+01/
C     ..
C     .. Executable Statements ..
C-------------------------
C
C                     ABS(X) .LE. 0.5
C
      ax = abs(x)
      IF (ax.GT.0.5D0) GO TO 10
      t = x*x
      top = ((((a(1)*t+a(2))*t+a(3))*t+a(4))*t+a(5)) + 1.0D0
      bot = ((b(1)*t+b(2))*t+b(3))*t + 1.0D0
      erfc1 = 0.5D0 + (0.5D0-x* (top/bot))
      IF (ind.NE.0) erfc1 = exp(t)*erfc1
      RETURN
C
C                  0.5 .LT. ABS(X) .LE. 4
C
   10 IF (ax.GT.4.0D0) GO TO 20
      top = ((((((p(1)*ax+p(2))*ax+p(3))*ax+p(4))*ax+p(5))*ax+p(6))*ax+
     +      p(7))*ax + p(8)
      bot = ((((((q(1)*ax+q(2))*ax+q(3))*ax+q(4))*ax+q(5))*ax+q(6))*ax+
     +      q(7))*ax + q(8)
      erfc1 = top/bot
      GO TO 40
C
C                      ABS(X) .GT. 4
C
   20 IF (x.LE.-5.6D0) GO TO 60
      IF (ind.NE.0) GO TO 30
      IF (x.GT.100.0D0) GO TO 70
      IF (x*x.GT.-exparg(1)) GO TO 70
C
   30 t = (1.0D0/x)**2
      top = (((r(1)*t+r(2))*t+r(3))*t+r(4))*t + r(5)
      bot = (((s(1)*t+s(2))*t+s(3))*t+s(4))*t + 1.0D0
      erfc1 = (c-t*top/bot)/ax
C
C                      FINAL ASSEMBLY
C
   40 IF (ind.EQ.0) GO TO 50
      IF (x.LT.0.0D0) erfc1 = 2.0D0*exp(x*x) - erfc1
      RETURN

   50 w = dble(x)*dble(x)
      t = w
      e = w - dble(t)
      erfc1 = ((0.5D0+ (0.5D0-e))*exp(-t))*erfc1
      IF (x.LT.0.0D0) erfc1 = 2.0D0 - erfc1
      RETURN
C
C             LIMIT VALUE FOR LARGE NEGATIVE X
C
   60 erfc1 = 2.0D0
      IF (ind.NE.0) erfc1 = 2.0D0*exp(x*x)
      RETURN
C
C             LIMIT VALUE FOR LARGE POSITIVE X
C                       WHEN IND = 0
C
   70 erfc1 = 0.0D0
      RETURN

      END
      DOUBLE PRECISION FUNCTION exparg(l)
C--------------------------------------------------------------------
C     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
C     EXP(W) CAN BE COMPUTED.
C
C     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
C     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
C
C     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
C--------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER l
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION lnb
      INTEGER b,m
C     ..
C     .. External Functions ..
      INTEGER ipmpar
      EXTERNAL ipmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,dlog
C     ..
C     .. Executable Statements ..
C
      b = ipmpar(4)
      IF (b.NE.2) GO TO 10
      lnb = .69314718055995D0
      GO TO 40

   10 IF (b.NE.8) GO TO 20
      lnb = 2.0794415416798D0
      GO TO 40

   20 IF (b.NE.16) GO TO 30
      lnb = 2.7725887222398D0
      GO TO 40

   30 lnb = dlog(dble(b))
C
   40 IF (l.EQ.0) GO TO 50
      m = ipmpar(9) - 1
      exparg = 0.99999D0* (m*lnb)
      RETURN

   50 m = ipmpar(10)
      exparg = 0.99999D0* (m*lnb)
      RETURN

      END
      DOUBLE PRECISION FUNCTION gam1(a)
C     ------------------------------------------------------------------
C     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5
C     ------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bot,d,s1,s2,t,top,w
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION p(7),q(5),r(9)
C     ..
C     .. Data statements ..
C     -------------------
C     -------------------
C     -------------------
C     -------------------
      DATA p(1)/.577215664901533D+00/,p(2)/-.409078193005776D+00/,
     +     p(3)/-.230975380857675D+00/,p(4)/.597275330452234D-01/,
     +     p(5)/.766968181649490D-02/,p(6)/-.514889771323592D-02/,
     +     p(7)/.589597428611429D-03/
      DATA q(1)/.100000000000000D+01/,q(2)/.427569613095214D+00/,
     +     q(3)/.158451672430138D+00/,q(4)/.261132021441447D-01/,
     +     q(5)/.423244297896961D-02/
      DATA r(1)/-.422784335098468D+00/,r(2)/-.771330383816272D+00/,
     +     r(3)/-.244757765222226D+00/,r(4)/.118378989872749D+00/,
     +     r(5)/.930357293360349D-03/,r(6)/-.118290993445146D-01/,
     +     r(7)/.223047661158249D-02/,r(8)/.266505979058923D-03/,
     +     r(9)/-.132674909766242D-03/
      DATA s1/.273076135303957D+00/,s2/.559398236957378D-01/
C     ..
C     .. Executable Statements ..
C     -------------------
      t = a
      d = a - 0.5D0
      IF (d.GT.0.0D0) t = d - 0.5D0

C     IF (t) 40,10,20

      IF (t.eq.1) GO TO 40
      IF (t.eq.2) GO TO 10
      IF (t.eq.3) GO TO 20
C
   10 gam1 = 0.0D0
      RETURN
C
   20 top = (((((p(7)*t+p(6))*t+p(5))*t+p(4))*t+p(3))*t+p(2))*t + p(1)
      bot = (((q(5)*t+q(4))*t+q(3))*t+q(2))*t + 1.0D0
      w = top/bot
      IF (d.GT.0.0D0) GO TO 30
      gam1 = a*w
      RETURN

   30 gam1 = (t/a)* ((w-0.5D0)-0.5D0)
      RETURN
C
   40 top = (((((((r(9)*t+r(8))*t+r(7))*t+r(6))*t+r(5))*t+r(4))*t+r(3))*
     +      t+r(2))*t + r(1)
      bot = (s2*t+s1)*t + 1.0D0
      w = top/bot
      IF (d.GT.0.0D0) GO TO 50
      gam1 = a* ((w+0.5D0)+0.5D0)
      RETURN

   50 gam1 = t*w/a
      RETURN

      END
      DOUBLE PRECISION FUNCTION gamma(a)
C-----------------------------------------------------------------------
C
C         EVALUATION OF THE GAMMA FUNCTION FOR REAL ARGUMENTS
C
C                           -----------
C
C     GAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
C     BE COMPUTED.
C
C-----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C          NAVAL SURFACE WEAPONS CENTER
C          DAHLGREN, VIRGINIA
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bot,d,g,lnx,pi,r1,r2,r3,r4,r5,s,t,top,w,x,z
      INTEGER i,j,m,n
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION p(7),q(7)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION exparg,spmpar
      EXTERNAL exparg,spmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,exp,int,mod,sin
C     ..
C     .. Data statements ..
C--------------------------
C     D = 0.5*(LN(2*PI) - 1)
C--------------------------
C--------------------------
C--------------------------
      DATA pi/3.1415926535898D0/
      DATA d/.41893853320467274178D0/
      DATA p(1)/.539637273585445D-03/,p(2)/.261939260042690D-02/,
     +     p(3)/.204493667594920D-01/,p(4)/.730981088720487D-01/,
     +     p(5)/.279648642639792D+00/,p(6)/.553413866010467D+00/,
     +     p(7)/1.0D0/
      DATA q(1)/-.832979206704073D-03/,q(2)/.470059485860584D-02/,
     +     q(3)/.225211131035340D-01/,q(4)/-.170458969313360D+00/,
     +     q(5)/-.567902761974940D-01/,q(6)/.113062953091122D+01/,
     +     q(7)/1.0D0/
      DATA r1/.820756370353826D-03/,r2/-.595156336428591D-03/,
     +     r3/.793650663183693D-03/,r4/-.277777777770481D-02/,
     +     r5/.833333333333333D-01/
C     ..
C     .. Executable Statements ..
C--------------------------
      gamma = 0.0D0
      x = a
      IF (abs(a).GE.15.0D0) GO TO 110
C-----------------------------------------------------------------------
C            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15
C-----------------------------------------------------------------------
      t = 1.0D0
      m = int(a) - 1
C
C     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2
C
C     IF (m) 40,30,10

      IF (m.eq.1) GO TO 40
      IF (m.eq.2) GO TO 30
      IF (m.eq.3) GO TO 10

   10 DO 20 j = 1,m
          x = x - 1.0D0
          t = x*t
   20 CONTINUE
   30 x = x - 1.0D0
      GO TO 80
C
C     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1
C
   40 t = a
      IF (a.GT.0.0D0) GO TO 70
      m = -m - 1
      IF (m.EQ.0) GO TO 60
      DO 50 j = 1,m
          x = x + 1.0D0
          t = x*t
   50 CONTINUE
   60 x = (x+0.5D0) + 0.5D0
      t = x*t
      IF (t.EQ.0.0D0) RETURN
C
   70 CONTINUE
C
C     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
C     CODE MAY BE OMITTED IF DESIRED.
C
      IF (abs(t).GE.1.D-30) GO TO 80
      IF (abs(t)*spmpar(3).LE.1.0001D0) RETURN
      gamma = 1.0D0/t
      RETURN
C
C     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1
C
   80 top = p(1)
      bot = q(1)
      DO 90 i = 2,7
          top = p(i) + x*top
          bot = q(i) + x*bot
   90 CONTINUE
      gamma = top/bot
C
C     TERMINATION
C
      IF (a.LT.1.0D0) GO TO 100
      gamma = gamma*t
      RETURN

  100 gamma = gamma/t
      RETURN
C-----------------------------------------------------------------------
C            EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
C-----------------------------------------------------------------------
  110 IF (abs(a).GE.1.D3) RETURN
      IF (a.GT.0.0D0) GO TO 120
      x = -a
      n = x
      t = x - n
      IF (t.GT.0.9D0) t = 1.0D0 - t
      s = sin(pi*t)/pi
      IF (mod(n,2).EQ.0) s = -s
      IF (s.EQ.0.0D0) RETURN
C
C     COMPUTE THE MODIFIED ASYMPTOTIC SUM
C
  120 t = 1.0D0/ (x*x)
      g = ((((r1*t+r2)*t+r3)*t+r4)*t+r5)/x
C
C     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
C     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.
C
      lnx = dlog(x)
C
C     FINAL ASSEMBLY
C
      z = x
      g = (d+g) + (z-0.5D0)* (lnx-1.D0)
      w = g
      t = g - dble(w)
      IF (w.GT.0.99999D0*exparg(0)) RETURN
      gamma = exp(w)* (1.0D0+t)
      IF (a.LT.0.0D0) gamma = (1.0D0/ (gamma*s))/x
      RETURN

      END
      SUBROUTINE gratio(a,x,ans,qans,ind)
C ----------------------------------------------------------------------
C        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
C                      P(A,X) AND Q(A,X)
C
C                        ----------
C
C     IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X
C     ARE NOT BOTH 0.
C
C     ANS AND QANS ARE VARIABLES. GRATIO ASSIGNS ANS THE VALUE
C     P(A,X) AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER.
C     IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS
C     POSSIBLE (UP TO 14 SIGNIFICANT DIGITS). OTHERWISE, IF
C     IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE
C     6-TH SIGNIFICANT DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY
C     IS REQUESTED TO WITHIN 1 UNIT OF THE 3RD SIGNIFICANT DIGIT.
C
C     ERROR RETURN ...
C        ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE,
C     WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT.
C     P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN
C     X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE.
C ----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WEAPONS CENTER
C        DAHLGREN, VIRGINIA
C     --------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,ans,qans,x
      INTEGER ind
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a2n,a2nm1,acc,alog10,am0,amn,an,an0,apn,b2n,
     +                 b2nm1,c,c0,c1,c2,c3,c4,c5,c6,cma,d10,d20,d30,d40,
     +                 d50,d60,d70,e,e0,g,h,j,l,r,rt2pin,rta,rtpi,rtx,s,
     +                 sum,t,t1,third,tol,twoa,u,w,x0,y,z
      INTEGER i,iop,m,max,n
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION acc0(3),big(3),d0(13),d1(12),d2(10),d3(8),d4(6),
     +                 d5(4),d6(2),e00(3),wk(20),x00(3)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION erf,erfc1,gam1,gamma,rexp,rlog,spmpar
      EXTERNAL erf,erfc1,gam1,gamma,rexp,rlog,spmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,dmax1,exp,int,sqrt
C     ..
C     .. Data statements ..
C     --------------------
C     --------------------
C     ALOG10 = LN(10)
C     RT2PIN = 1/SQRT(2*PI)
C     RTPI   = SQRT(PI)
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
C     --------------------
      DATA acc0(1)/5.D-15/,acc0(2)/5.D-7/,acc0(3)/5.D-4/
      DATA big(1)/20.0D0/,big(2)/14.0D0/,big(3)/10.0D0/
      DATA e00(1)/.25D-3/,e00(2)/.25D-1/,e00(3)/.14D0/
      DATA x00(1)/31.0D0/,x00(2)/17.0D0/,x00(3)/9.7D0/
      DATA alog10/2.30258509299405D0/
      DATA rt2pin/.398942280401433D0/
      DATA rtpi/1.77245385090552D0/
      DATA third/.333333333333333D0/
      DATA d0(1)/.833333333333333D-01/,d0(2)/-.148148148148148D-01/,
     +     d0(3)/.115740740740741D-02/,d0(4)/.352733686067019D-03/,
     +     d0(5)/-.178755144032922D-03/,d0(6)/.391926317852244D-04/,
     +     d0(7)/-.218544851067999D-05/,d0(8)/-.185406221071516D-05/,
     +     d0(9)/.829671134095309D-06/,d0(10)/-.176659527368261D-06/,
     +     d0(11)/.670785354340150D-08/,d0(12)/.102618097842403D-07/,
     +     d0(13)/-.438203601845335D-08/
      DATA d10/-.185185185185185D-02/,d1(1)/-.347222222222222D-02/,
     +     d1(2)/.264550264550265D-02/,d1(3)/-.990226337448560D-03/,
     +     d1(4)/.205761316872428D-03/,d1(5)/-.401877572016461D-06/,
     +     d1(6)/-.180985503344900D-04/,d1(7)/.764916091608111D-05/,
     +     d1(8)/-.161209008945634D-05/,d1(9)/.464712780280743D-08/,
     +     d1(10)/.137863344691572D-06/,d1(11)/-.575254560351770D-07/,
     +     d1(12)/.119516285997781D-07/
      DATA d20/.413359788359788D-02/,d2(1)/-.268132716049383D-02/,
     +     d2(2)/.771604938271605D-03/,d2(3)/.200938786008230D-05/,
     +     d2(4)/-.107366532263652D-03/,d2(5)/.529234488291201D-04/,
     +     d2(6)/-.127606351886187D-04/,d2(7)/.342357873409614D-07/,
     +     d2(8)/.137219573090629D-05/,d2(9)/-.629899213838006D-06/,
     +     d2(10)/.142806142060642D-06/
      DATA d30/.649434156378601D-03/,d3(1)/.229472093621399D-03/,
     +     d3(2)/-.469189494395256D-03/,d3(3)/.267720632062839D-03/,
     +     d3(4)/-.756180167188398D-04/,d3(5)/-.239650511386730D-06/,
     +     d3(6)/.110826541153473D-04/,d3(7)/-.567495282699160D-05/,
     +     d3(8)/.142309007324359D-05/
      DATA d40/-.861888290916712D-03/,d4(1)/.784039221720067D-03/,
     +     d4(2)/-.299072480303190D-03/,d4(3)/-.146384525788434D-05/,
     +     d4(4)/.664149821546512D-04/,d4(5)/-.396836504717943D-04/,
     +     d4(6)/.113757269706784D-04/
      DATA d50/-.336798553366358D-03/,d5(1)/-.697281375836586D-04/,
     +     d5(2)/.277275324495939D-03/,d5(3)/-.199325705161888D-03/,
     +     d5(4)/.679778047793721D-04/
      DATA d60/.531307936463992D-03/,d6(1)/-.592166437353694D-03/,
     +     d6(2)/.270878209671804D-03/
      DATA d70/.344367606892378D-03/
C     ..
C     .. Executable Statements ..
C     --------------------
C     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
C            FLOATING POINT NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
C
      e = spmpar(1)
C
C     --------------------
      IF (a.LT.0.0D0 .OR. x.LT.0.0D0) GO TO 430
      IF (a.EQ.0.0D0 .AND. x.EQ.0.0D0) GO TO 430
      IF (a*x.EQ.0.0D0) GO TO 420
C
      iop = ind + 1
      IF (iop.NE.1 .AND. iop.NE.2) iop = 3
      acc = dmax1(acc0(iop),e)
      e0 = e00(iop)
      x0 = x00(iop)
C
C            SELECT THE APPROPRIATE ALGORITHM
C
      IF (a.GE.1.0D0) GO TO 10
      IF (a.EQ.0.5D0) GO TO 390
      IF (x.LT.1.1D0) GO TO 160
      t1 = a*dlog(x) - x
      u = a*exp(t1)
      IF (u.EQ.0.0D0) GO TO 380
      r = u* (1.0D0+gam1(a))
      GO TO 250
C
   10 IF (a.GE.big(iop)) GO TO 30
      IF (a.GT.x .OR. x.GE.x0) GO TO 20
      twoa = a + a
      m = int(twoa)
      IF (twoa.NE.dble(m)) GO TO 20
      i = m/2
      IF (a.EQ.dble(i)) GO TO 210
      GO TO 220

   20 t1 = a*dlog(x) - x
      r = exp(t1)/gamma(a)
      GO TO 40
C
   30 l = x/a
      IF (l.EQ.0.0D0) GO TO 370
      s = 0.5D0 + (0.5D0-l)
      z = rlog(l)
      IF (z.GE.700.0D0/a) GO TO 410
      y = a*z
      rta = sqrt(a)
      IF (abs(s).LE.e0/rta) GO TO 330
      IF (abs(s).LE.0.4D0) GO TO 270
C
      t = (1.0D0/a)**2
      t1 = (((0.75D0*t-1.0D0)*t+3.5D0)*t-105.0D0)/ (a*1260.0D0)
      t1 = t1 - y
      r = rt2pin*rta*exp(t1)
C
   40 IF (r.EQ.0.0D0) GO TO 420
      IF (x.LE.dmax1(a,alog10)) GO TO 50
      IF (x.LT.x0) GO TO 250
      GO TO 100
C
C                 TAYLOR SERIES FOR P/R
C
   50 apn = a + 1.0D0
      t = x/apn
      wk(1) = t
      DO 60 n = 2,20
          apn = apn + 1.0D0
          t = t* (x/apn)
          IF (t.LE.1.D-3) GO TO 70
          wk(n) = t
   60 CONTINUE
      n = 20
C
   70 sum = t
      tol = 0.5D0*acc
   80 apn = apn + 1.0D0
      t = t* (x/apn)
      sum = sum + t
      IF (t.GT.tol) GO TO 80
C
      max = n - 1
      DO 90 m = 1,max
          n = n - 1
          sum = sum + wk(n)
   90 CONTINUE
      ans = (r/a)* (1.0D0+sum)
      qans = 0.5D0 + (0.5D0-ans)
      RETURN
C
C                 ASYMPTOTIC EXPANSION
C
  100 amn = a - 1.0D0
      t = amn/x
      wk(1) = t
      DO 110 n = 2,20
          amn = amn - 1.0D0
          t = t* (amn/x)
          IF (abs(t).LE.1.D-3) GO TO 120
          wk(n) = t
  110 CONTINUE
      n = 20
C
  120 sum = t
  130 IF (abs(t).LE.acc) GO TO 140
      amn = amn - 1.0D0
      t = t* (amn/x)
      sum = sum + t
      GO TO 130
C
  140 max = n - 1
      DO 150 m = 1,max
          n = n - 1
          sum = sum + wk(n)
  150 CONTINUE
      qans = (r/x)* (1.0D0+sum)
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
C             TAYLOR SERIES FOR P(A,X)/X**A
C
  160 an = 3.0D0
      c = x
      sum = x/ (a+3.0D0)
      tol = 3.0D0*acc/ (a+1.0D0)
  170 an = an + 1.0D0
      c = -c* (x/an)
      t = c/ (a+an)
      sum = sum + t
      IF (abs(t).GT.tol) GO TO 170
      j = a*x* ((sum/6.0D0-0.5D0/ (a+2.0D0))*x+1.0D0/ (a+1.0D0))
C
      z = a*dlog(x)
      h = gam1(a)
      g = 1.0D0 + h
      IF (x.LT.0.25D0) GO TO 180
      IF (a.LT.x/2.59D0) GO TO 200
      GO TO 190

  180 IF (z.GT.-.13394D0) GO TO 200
C
  190 w = exp(z)
      ans = w*g* (0.5D0+ (0.5D0-j))
      qans = 0.5D0 + (0.5D0-ans)
      RETURN
C
  200 l = rexp(z)
      w = 0.5D0 + (0.5D0+l)
      qans = (w*j-l)*g - h
      IF (qans.LT.0.0D0) GO TO 380
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
C             FINITE SUMS FOR Q WHEN A .GE. 1
C                 AND 2*A IS AN INTEGER
C
  210 sum = exp(-x)
      t = sum
      n = 1
      c = 0.0D0
      GO TO 230
C
  220 rtx = sqrt(x)
      sum = erfc1(0,rtx)
      t = exp(-x)/ (rtpi*rtx)
      n = 0
      c = -0.5D0
C
  230 IF (n.EQ.i) GO TO 240
      n = n + 1
      c = c + 1.0D0
      t = (x*t)/c
      sum = sum + t
      GO TO 230

  240 qans = sum
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
C              CONTINUED FRACTION EXPANSION
C
  250 tol = dmax1(5.0D0*e,acc)
      a2nm1 = 1.0D0
      a2n = 1.0D0
      b2nm1 = x
      b2n = x + (1.0D0-a)
      c = 1.0D0
  260 a2nm1 = x*a2n + c*a2nm1
      b2nm1 = x*b2n + c*b2nm1
      am0 = a2nm1/b2nm1
      c = c + 1.0D0
      cma = c - a
      a2n = a2nm1 + cma*a2n
      b2n = b2nm1 + cma*b2n
      an0 = a2n/b2n
      IF (abs(an0-am0).GE.tol*an0) GO TO 260
C
      qans = r*an0
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
C                GENERAL TEMME EXPANSION
C
  270 IF (abs(s).LE.2.0D0*e .AND. a*e*e.GT.3.28D-3) GO TO 430
      c = exp(-y)
      w = 0.5D0*erfc1(1,sqrt(y))
      u = 1.0D0/a
      z = sqrt(z+z)
      IF (l.LT.1.0D0) z = -z
C     IF (iop-2) 280,290,300

      IF (iop-2.eq.1) GO TO 280
      IF (iop-2.eq.2) GO TO 290
      IF (iop-2.eq.3) GO TO 300

C
  280 IF (abs(s).LE.1.D-3) GO TO 340
      c0 = ((((((((((((d0(13)*z+d0(12))*z+d0(11))*z+d0(10))*z+d0(9))*z+
     +     d0(8))*z+d0(7))*z+d0(6))*z+d0(5))*z+d0(4))*z+d0(3))*z+d0(2))*
     +     z+d0(1))*z - third
      c1 = (((((((((((d1(12)*z+d1(11))*z+d1(10))*z+d1(9))*z+d1(8))*z+
     +     d1(7))*z+d1(6))*z+d1(5))*z+d1(4))*z+d1(3))*z+d1(2))*z+d1(1))*
     +     z + d10
      c2 = (((((((((d2(10)*z+d2(9))*z+d2(8))*z+d2(7))*z+d2(6))*z+
     +     d2(5))*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20
      c3 = (((((((d3(8)*z+d3(7))*z+d3(6))*z+d3(5))*z+d3(4))*z+d3(3))*z+
     +     d3(2))*z+d3(1))*z + d30
      c4 = (((((d4(6)*z+d4(5))*z+d4(4))*z+d4(3))*z+d4(2))*z+d4(1))*z +
     +     d40
      c5 = (((d5(4)*z+d5(3))*z+d5(2))*z+d5(1))*z + d50
      c6 = (d6(2)*z+d6(1))*z + d60
      t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0
      GO TO 310
C
  290 c0 = (((((d0(6)*z+d0(5))*z+d0(4))*z+d0(3))*z+d0(2))*z+d0(1))*z -
     +     third
      c1 = (((d1(4)*z+d1(3))*z+d1(2))*z+d1(1))*z + d10
      c2 = d2(1)*z + d20
      t = (c2*u+c1)*u + c0
      GO TO 310
C
  300 t = ((d0(3)*z+d0(2))*z+d0(1))*z - third
C
  310 IF (l.LT.1.0D0) GO TO 320
      qans = c* (w+rt2pin*t/rta)
      ans = 0.5D0 + (0.5D0-qans)
      RETURN

  320 ans = c* (w-rt2pin*t/rta)
      qans = 0.5D0 + (0.5D0-ans)
      RETURN
C
C               TEMME EXPANSION FOR L = 1
C
  330 IF (a*e*e.GT.3.28D-3) GO TO 430
      c = 0.5D0 + (0.5D0-y)
      w = (0.5D0-sqrt(y)* (0.5D0+ (0.5D0-y/3.0D0))/rtpi)/c
      u = 1.0D0/a
      z = sqrt(z+z)
      IF (l.LT.1.0D0) z = -z
C     IF (iop-2) 340,350,360

      IF (iop-2.eq.1) GO TO 340
      IF (iop-2.eq.2) GO TO 350
      IF (iop-2.eq.3) GO TO 360

C
  340 c0 = ((((((d0(7)*z+d0(6))*z+d0(5))*z+d0(4))*z+d0(3))*z+d0(2))*z+
     +     d0(1))*z - third
      c1 = (((((d1(6)*z+d1(5))*z+d1(4))*z+d1(3))*z+d1(2))*z+d1(1))*z +
     +     d10
      c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20
      c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30
      c4 = (d4(2)*z+d4(1))*z + d40
      c5 = (d5(2)*z+d5(1))*z + d50
      c6 = d6(1)*z + d60
      t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0
      GO TO 310
C
  350 c0 = (d0(2)*z+d0(1))*z - third
      c1 = d1(1)*z + d10
      t = (d20*u+c1)*u + c0
      GO TO 310
C
  360 t = d0(1)*z - third
      GO TO 310
C
C                     SPECIAL CASES
C
  370 ans = 0.0D0
      qans = 1.0D0
      RETURN
C
  380 ans = 1.0D0
      qans = 0.0D0
      RETURN
C
  390 IF (x.GE.0.25D0) GO TO 400
      ans = erf(sqrt(x))
      qans = 0.5D0 + (0.5D0-ans)
      RETURN

  400 qans = erfc1(0,sqrt(x))
      ans = 0.5D0 + (0.5D0-qans)
      RETURN
C
  410 IF (abs(s).LE.2.0D0*e) GO TO 430
  420 IF (x.LE.a) GO TO 370
      GO TO 380
C
C                     ERROR RETURN
C
  430 ans = 2.0D0
      RETURN

      END
      DOUBLE PRECISION FUNCTION rexp(x)
C-----------------------------------------------------------------------
C            EVALUATION OF THE FUNCTION EXP(X) - 1
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION p1,p2,q1,q2,q3,q4,w
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
C     .. Data statements ..
      DATA p1/.914041914819518D-09/,p2/.238082361044469D-01/,
     +     q1/-.499999999085958D+00/,q2/.107141568980644D+00/,
     +     q3/-.119041179760821D-01/,q4/.595130811860248D-03/
C     ..
C     .. Executable Statements ..
C-----------------------
      IF (abs(x).GT.0.15D0) GO TO 10
      rexp = x* (((p2*x+p1)*x+1.0D0)/ ((((q4*x+q3)*x+q2)*x+q1)*x+1.0D0))
      RETURN
C
   10 w = exp(x)
      IF (x.GT.0.0D0) GO TO 20
      rexp = (w-0.5D0) - 0.5D0
      RETURN

   20 rexp = w* (0.5D0+ (0.5D0-1.0D0/w))
      RETURN

      END
      DOUBLE PRECISION FUNCTION rlog(x)
C     -------------------
C     COMPUTATION OF  X - 1 - LN(X)
C     -------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,b,p0,p1,p2,q1,q2,r,t,u,w,w1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,dlog
C     ..
C     .. Data statements ..
C     -------------------
      DATA a/.566749439387324D-01/
      DATA b/.456512608815524D-01/
      DATA p0/.333333333333333D+00/,p1/-.224696413112536D+00/,
     +     p2/.620886815375787D-02/
      DATA q1/-.127408923933623D+01/,q2/.354508718369557D+00/
C     ..
C     .. Executable Statements ..
C     -------------------
      IF (x.LT.0.61D0 .OR. x.GT.1.57D0) GO TO 40
      IF (x.LT.0.82D0) GO TO 10
      IF (x.GT.1.18D0) GO TO 20
C
C              ARGUMENT REDUCTION
C
      u = (x-0.5D0) - 0.5D0
      w1 = 0.0D0
      GO TO 30
C
   10 u = dble(x) - 0.7D0
      u = u/0.7D0
      w1 = a - u*0.3D0
      GO TO 30
C
   20 u = 0.75D0*dble(x) - 1.D0
      w1 = b + u/3.0D0
C
C               SERIES EXPANSION
C
   30 r = u/ (u+2.0D0)
      t = r*r
      w = ((p2*t+p1)*t+p0)/ ((q2*t+q1)*t+1.0D0)
      rlog = 2.0D0*t* (1.0D0/ (1.0D0-r)-r*w) + w1
      RETURN
C
C
   40 r = (x-0.5D0) - 0.5D0
      rlog = r - dlog(x)
      RETURN

      END
