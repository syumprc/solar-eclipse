C  FUNCTION DCOVAR is the discrete version of COVAR...
C    It bridges the gap between the standalone and "Fisher" versions

      DOUBLE PRECISION FUNCTION DCOVAR (PERI,PERJ,TRAITI,TRAITJ
     1,PAR,NPAR,VARDATA,NVAR,NTOT,MALEI,MALEJ)

      DOUBLE PRECISION PAR(NPAR),VARDATA(NVAR,NTOT)
      INTEGER PERI,PERJ,TRAITI,TRAITJ,NVAR,NTOT
      LOGICAL MALEI,MALEJ

      DOUBLE PRECISION KIN2,DELTA7,COVAR
      INTEGER NTRAIT,ITSE

C Currently, only one trait is supported

      NTRAIT = 1

C Intra trait SE measurement feature not currently supported

      ITSE = 0
      
C In the discrete environment, the phi2.gz matrix must be loaded
C   to supply the phi2 and delta7 matrices (superceding internal
C   values):  load matrix phi2.gz phi2 delta7
C   The following "KIN2" and "DELTA7" values are then not used.

      KIN2 = 0.0
      DELTA7 = 0.0
      
C In discrete environment VARI, VARJ are positions in the VARDATA array
C  (See COVAR invocation below.)

      DCOVAR = COVAR(PAR,VARDATA(1,PERI),VARDATA(1,PERJ),DELTA7,KIN2
     1,NPAR,NTRAIT,NVAR,PERI,PERJ,TRAITI,TRAITJ,ITSE,MALEI,MALEJ)

C     PRINT *,"  ",PERI,"  ",PERJ,"  ",DCOVAR

      END
