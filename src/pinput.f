C Modified 1996-1997 by Charles Peterson, SFBR.
C This is the New Fisher signature.  Do not remove.  fI6%.
C
      SUBROUTINE PINPUT(GRLIST,PARENT,PERSON,VAR,VMAX,VMEAN,VMIN,VVAR
     1,FATHER,GROUP,MOTHER,NTL,PERM,PTL,SEX,CARRAY,ABSENT,IERROR,IPED
     2,MXTWIN,MVAR,NMALES,NPBAND,NPTOT,NSAMP,NTWINS,NVAR,NV2READ
     3,UNIT1,UNIT7,UNIT3,FRMT2,ECHO,CONOUT,NTRAIT
     4,NCOVAR,NSAMPA,VBINARY,NSAMPTOT,FIRSTPID,NPEO,MVFATHER,MVMOTHER
     5,MVGROUP,MVVAR,VTRAITS,NTRASAMP,SAMPLEDAT,evdphase,mibdid
     6,ifam,nimatrix,EVDstarted)
C
C     THIS SUBROUTINE INPUTS THE PEDIGREE DATA AND ATTEMPTS TO
C     DETECT ERRORS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION GRLIST(MXTWIN),PARENT(2,NPTOT),PERSON(NPTOT)
     1,VAR(MVAR),VMAX(NVAR),VMEAN(NVAR),VMIN(NVAR),VVAR(NVAR),GETPHENO
     2,UNKNOWN,DSELECT,MVVAR(*),OLD,VOLD,SAMPLE,VSAMPLE,VV,UNBAL
     3,MEAN1,MEAN2,SD1,SD2,ZSCORE,RSAVE
      INTEGER FATHER(NPTOT),GROUP(NPTOT),MOTHER(NPTOT),NTL(MXTWIN)
     1,PERM(NPTOT),PTL(MXTWIN),SEX(NPTOT),UNIT1,UNIT7,UNIT3,CONOUT
     2,NV2READ,NTRAIT,NCOVAR,I,J,K,L,IPBINDX,ISELECT,ISEX,IGNORE
     3,M,NSAMPA(NVAR),FIRSTNP,NSAMPTOT,NPEO,MVINDEX,MVPEO,MVPROB
     4,IVARZ,ITIN,IGROUP,MVFATHER(*),MVMOTHER(*),MVGROUP(*)
     5,VTRAITS,JSTART,NSAMP,NTRASAMP(NTRAIT),ISAVE,LBASE
     6,TRAITNEEDS,OLDFIRST,evdphase,mibdid,ifam,nimatrix
      LOGICAL FORCEMISS(NTRAIT),EVDstarted
      CHARACTER*8 CARRAY(*)*(*),FRMT2*(*)
      LOGICAL ECHO,MISS,VBINARY(NVAR),IPBAND,ONETRAIT,ZSD,SAMPLEDAT
      CHARACTER*18 PERMID(NPTOT),FAMID(NPTOT),FIRSTPID,IDSAVE
      INTEGER UNITNO, ipacked(nptot+1),iun,uperson(nptot),IIBDID
     1,notfound
C
C **** end of data declarations
C
C     FORMAT STATEMENTS FOR OUTPUT OF DATA AND ERROR MESSAGES.
C
 500  FORMAT(I5,5(1X,A8),(T51,3(1X,A8),:))
 550  FORMAT('    *** Error!  MZ twin ',A,' has a different parent.')
 560  FORMAT('    *** Error!  MZ twin ',A,' has opposite sex.')
 570  FORMAT('    *** Error!  Pedigree with individual ',A
     1,/,' exceeds the maximum of ',I3,' twin sets.')
 590  FORMAT(' *** ERROR *** PERSON NUMBER',I4,' IN PEDIGREE NUMBER'
     1,I4,' IS AN MZ TWIN'/' WITH NO IDENTICAL SIBS.')
 620  FORMAT(' *** ERROR *** PERSON NUMBER',I4,' IN PEDIGREE NUMBER'
     1,I4,' HAS PARENTS OF'/' THE SAME SEX.')
 630  FORMAT(' *** ERROR *** PARENT NUMBER',I2,' OF PERSON NUMBER'
     1,I4,' IN PEDIGREE NUMBER',/,I4,' IS NOT IN THE PEDIGREE.')
 640  FORMAT(' *** ERROR *** PERSONS NUMBERED',I4,' AND',I4,' IN'
     1,' PEDIGREE NUMBER',I4,' HAVE'/' THE SAME ID.')
 650  FORMAT(' *** ERROR *** SOMEONE NEAR PERSON NUMBER',I4,' IN'
     1,' PEDIGREE NUMBER',/,I4,' IS HIS/HER OWN ANCESTOR.')
C
C See if unbalanced traits are permitted
C
C
      iun=0
      IGNORE=0
      CALL DOPTION ('UnbalancedTraits', UNBAL)
      IF (UNBAL.EQ.0D0) THEN
         IGNORE=9999
      ENDIF
      IF (EVDPHASE.GT.0) THEN
         IGNORE=9999
      ENDIF
      ZSD=.FALSE.
      CALL DOPTION ('zscore',ZSCORE)
      IF (ZSCORE.NE.0D0) THEN
         CALL DOPTION ('ZMEAN1',MEAN1)
         CALL DOPTION ('ZSD1',SD1)
         IF (SD1.EQ.0D0) ZSD=.TRUE.
         IF (NTRAIT.GT.1) THEN
            CALL DOPTION ('ZMEAN2',MEAN2)
            CALL DOPTION ('ZSD2',SD2)
            IF (SD2.EQ.0D0) ZSD=.TRUE.
         END IF
      END IF
      IF (ZSD) THEN
         CALL ZSDERROR
C        The above call never returns
         RETURN
      END IF

C
C We figure out how many probands there actually are here
C First phenotype after last trait is the proband var (0=not proband)
C
      UNKNOWN=0.0
      NPBAND=0
C
C     FOR EACH PERSON IN THE PEDIGREE, INPUT THE PERSONAL DATA AND
C     CHECK FOR IRREGULARITIES.
C
      LGR=0
      DO 10 I=1,NPTOT

      CALL NEXTPERSON (PERSON(I),PARENT(1,I),PARENT(2,I),SEX(I)
     1,GRP,PERMID(I),FAMID(I))

      IF (I.EQ.1) FIRSTPID = PERMID(I)
C
C INCREMENT MALE COUNT IF NECESSARY
C
      IF (SEX(I).EQ.1) THEN
         NMALES = NMALES + 1
      END IF
C
C     CHECK WHETHER THIS PERSON IS AN MZ TWIN.
C
      IF (GRP.NE.0) THEN
      NTWINS=NTWINS+1
      DO 30 J=1,LGR
      IF (GRP.EQ.GRLIST(J)) THEN
C
C     THE CURRENT TWIN IS A COTWIN.  CHECK THAT SEX AND PARENTS
C     ARE THE SAME AS THOSE OF THE PRIMARY TWIN.
C
      JPTL=PTL(J)
      IF (PARENT(1,I).EQ.PARENT(1,JPTL).AND.PARENT(2,I).EQ.
     1PARENT(2,JPTL)) GO TO 40
      IF (PARENT(1,I).EQ.PARENT(2,JPTL).AND.PARENT(2,I).EQ.
     1PARENT(1,JPTL)) GO TO 40
      IERROR=1
      WRITE(UNIT3,550) PERMID(I)
      WRITE(CONOUT,550) PERMID(I)
 40   IF (SEX(I).NE.SEX(JPTL)) THEN
      IERROR=1
      WRITE(UNIT3,560) PERMID(I)
      WRITE(CONOUT,560) PERMID(I)
      END IF
      GROUP(I)=J*100000+I
      NTL(J)=NTL(J)+1
      GO TO 50
      END IF
 30   CONTINUE
C
C     THE CURRENT TWIN IS A PRIMARY TWIN.  UPDATE THE PRIMARY TWIN LIST,
C     OR NOTE IF THE MAXIMUM TWIN COUNT HAS BEEN EXCEEDED.
C
      IF (LGR.LT.MXTWIN) THEN
      LGR=LGR+1
      PTL(LGR)=I
      NTL(LGR)=1
      GRLIST(LGR)=GRP
      GROUP(I)=LGR*100000+I
      ELSE
      GROUP(I)=I
      IERROR=1
      WRITE(UNIT3,570) FIRSTPID,MXTWIN
      WRITE(CONOUT,570) FIRSTPID,MXTWIN
      END IF
C
C     THE CURRENT PERSON IS NOT A TWIN.
C
      ELSE
      GROUP(I)=I
      END IF
C
C     CONVERT THE QUANTITATIVE CHARACTER VARIABLES TO NUMBERS AND CHECK
C     THAT NO PROBAND HAS MISSING TRAITS OR COVARS
C
 50   MISS = .FALSE.
      ONETRAIT = .FALSE.
      L = 0
C
C See if individual is blanked by matrix
      iibdid = person(i)
      call matcheck (iibdid,notfound)

C The first NTRAIT phenotypes are traits

      DO 55 J = 1, NTRAIT
         L = L + 1
         M = (I-1)*NVAR + L
         var(m) = getpheno (j)
C
C If restricting sample by matrix, do that now
C   by blanking all traits
C
         if (notfound.eq.1) then
            if (var(m).ne.ABSENT.and.J.eq.1) then
               nimatrix = nimatrix+1
            end if
            var(m) = ABSENT
         end if
C
         IF (VAR(M).EQ.ABSENT) THEN
            MISS = .TRUE.
         ELSE
            ONETRAIT = .TRUE.
            NSAMPA(J) = NSAMPA(J) + 1
            IF (ZSCORE.NE.0) THEN
               IF (J.EQ.1) THEN
                  VAR(M) = (VAR(M)-MEAN1)/SD1
               ELSE
                  VAR(M) = (VAR(M)-MEAN2)/SD2
               END IF
            END IF
         END IF
 55   CONTINUE
C
C Only one trait required if unbalanced allowed
C
      IF (ONETRAIT.AND.VTRAITS.GT.IGNORE) THEN
         MISS = .FALSE.
      END IF


C NTRAIT+1 is reserved for the IBDID value, which must be saved as a data
C element

      L = L + 1
      M = (I-1) * NVAR + L
      VAR(M) = PERSON(I)
      NSAMPA(NTRAIT+1) = NSAMPA(NTRAIT+1) + 1

C The NTRAIT+2 phenotype is the proband indicator (always present)
C We check for it here, then it is saved as if it were another
C covariate.

      L = L + 1
      M = (I-1) * NVAR + L
      IPBINDX = M
      NSAMPA(NTRAIT+2)=NSAMPA(NTRAIT+2)+1
      IPBAND=.FALSE.
      IF (0.NE.GETPHENO(NTRAIT+2)) THEN
         NPBAND=NPBAND+1
         IPBAND=.TRUE.
         VAR(IPBINDX) = 1
      ELSE
         VAR(IPBINDX) = 0
      ENDIF

C The remaining phenotypes are called covariates here, but might actually be
C   used for some other purposes

      DO 59 J = 1, VTRAITS
         FORCEMISS(J) = .FALSE.
 59   CONTINUE

      DO 60 J = 3, NCOVAR
         L = L + 1
         M = (I-1)*NVAR + L
         VAR(M) = GETPHENO (NTRAIT+J)
         IF (VAR(M).EQ.ABSENT) THEN
            IF (VTRAITS.LE.IGNORE) THEN
               MISS = .TRUE.
            ELSE
               DO 47 IT=1,NTRAIT
                  IF (TRAITNEEDS(IT,NTRAIT+J).NE.0) THEN
                     FORCEMISS(IT)=.TRUE.
                  ENDIF
 47            CONTINUE
            ENDIF
         ELSE
            NSAMPA(J+NTRAIT) = NSAMPA(J+NTRAIT) + 1
         END IF
 60   CONTINUE
C
C Now, make a second pass through traits, forcing those that are
C missing related covariates to themselves appear missing,
C unless we've already decided too many traits are missing
C
      IF (.NOT.MISS.AND.VTRAITS.GT.IGNORE) THEN
         ONETRAIT = .FALSE.
         DO 61 J = 1, NTRAIT
            M = (I-1)*NVAR + J
            IF (FORCEMISS(J).OR.VAR(M).EQ.ABSENT) THEN
               IF (VAR(M).NE.ABSENT) THEN
                  NSAMPA(J) = NSAMPA(J) - 1
                  VAR(M) = ABSENT
               ENDIF
               MISS = .TRUE.
            ELSE
               ONETRAIT = .TRUE.
            ENDIF
 61      CONTINUE

         IF (ONETRAIT) THEN
            MISS = .FALSE.
         ENDIF

      ENDIF




C
C PROBANDS ARE NOT ALLOWED TO HAVE MISSING VARIABLES
C   So, we simply make them non-probands
C   They won't contribute anything anyway
C
      IF (MISS.AND.IPBAND) THEN
         NPBAND = NPBAND - 1
         IPBAND = .FALSE.
         VAR(IPBINDX) = 0
      END IF
C
C     UPDATE THE DESCRIPTIVE STATISTICS FOR EACH VARIABLE.
C
      IF (.NOT.MISS) THEN
         NSAMPTOT = NSAMPTOT + 1
      ENDIF

      IF (.NOT.MISS.AND.(.NOT.IPBAND)) THEN
C
C Rolling counters for non-trait variables
C
         OLD=NSAMP
         NSAMP=NSAMP+1
         SAMPLE=NSAMP
         DO 70 J=1,NVAR
            VARJ=VAR((I-1)*NVAR+J)
            IF (J.LE.NTRAIT) THEN
               IF (VARJ.NE.ABSENT) THEN
C
C Rolling counters and stats for trait variables
C
                  VOLD=NTRASAMP(J)
                  NTRASAMP(J)=NTRASAMP(J)+1
                  VSAMPLE=NTRASAMP(J)
                  VMAX(J)=MAX(VMAX(J),VARJ)
                  VMIN(J)=MIN(VMIN(J),VARJ)
                  VMEAN(J)=(VOLD*VMEAN(J)+VARJ)/VSAMPLE
                  VV=(VARJ-VMEAN(J))**2
                  IF (VARJ.NE.VMIN(J).AND.VARJ.NE.VMAX(J)) THEN
                     VBINARY(J)=.FALSE.
                  ENDIF
                  IF (NTRASAMP(J).GT.1) VV=VV/VOLD
                  VVAR(J)=VOLD*VVAR(J)/VSAMPLE+VV
               END IF
            ELSE
C
C Stats for non-trait variables
C
               IF (VARJ.NE.ABSENT) THEN
                  VMAX(J)=MAX(VMAX(J),VARJ)
                  VMIN(J)=MIN(VMIN(J),VARJ)
                  VMEAN(J)=(OLD*VMEAN(J)+VARJ)/SAMPLE
                  VV=(VARJ-VMEAN(J))**2
                  IF (VARJ.NE.VMIN(J).AND.VARJ.NE.VMAX(J)) THEN
                     VBINARY(J)=.FALSE.
                  ENDIF
                  IF (NSAMP.GT.1) VV=VV/OLD
                  VVAR(J)=OLD*VVAR(J)/SAMPLE+VV
               ENDIF
            ENDIF
 70      CONTINUE
      END IF

 10   CONTINUE
C
C     CHECK THAT NO MZ TWIN LACKS A COTWIN OR HAS MORE THAN ONE COTWIN.
C
      DO 80 J=1,LGR
      IF (NTL(J).LE.1) THEN
      IERROR=1
      WRITE(UNIT3,590) PTL(J),IPED
      WRITE(CONOUT,590) PTL(J),IPED
      END IF
 80   CONTINUE
C
C     RETURN IF THERE HAVE BEEN ANY ERRORS DETECTED SO FAR.
C
      IF (IERROR.NE.0) RETURN
C
C     MOVE PROBANDS TO BEGINNING OF PEDIGREE
C     (PROBAND STATUS IS VARIABLE AFTER LAST TRAIT)
C
      IF (NPBAND.GT.0) THEN
       FIRSTNP=0
       DO 91 I=1,NPTOT
        IF (VAR((I-1)*NVAR+NTRAIT+2).EQ.0) THEN
C This is not a proband
         IF (FIRSTNP.EQ.0) FIRSTNP=I
        ELSE
C This is a proband
         IF (FIRSTNP.NE.0) THEN

C This is a proband, and there is a non-proband in front of it
C So, swap this person with the first non-proband in front of it

          RSAVE=PERSON(I)
          PERSON(I)=PERSON(FIRSTNP)
          PERSON(FIRSTNP)=RSAVE
          DO 92 J=1,2
           RSAVE=PARENT(J,I)
           PARENT(J,I)=PARENT(J,FIRSTNP)
           PARENT(J,FIRSTNP)=RSAVE
 92       CONTINUE
          ISAVE=SEX(I)
          SEX(I)=SEX(FIRSTNP)
          SEX(FIRSTNP)=ISAVE
          ISAVE=GROUP(I)
          GROUP(I)=GROUP(FIRSTNP)
          GROUP(FIRSTNP)=ISAVE
          IDSAVE=PERMID(I)
          PERMID(I)=PERMID(FIRSTNP)
          PERMID(FIRSTNP)=IDSAVE
          IDSAVE=FAMID(I)
          FAMID(I)=FAMID(FIRSTNP)
          FAMID(FIRSTNP)=IDSAVE
          DO 93 J=1,NVAR
           RSAVE=VAR((I-1)*NVAR+J)
           VAR((I-1)*NVAR+J)=VAR((FIRSTNP-1)*NVAR+J)
           VAR((FIRSTNP-1)*NVAR+J)=RSAVE
 93       CONTINUE
C
C Now we have to find if there are any non-probands
C between FIRSTNP and I
C 
          OLDFIRST=FIRSTNP
          FIRSTNP=0
          DO 94 J=OLDFIRST+1,I
           IF (VAR((J-1)*NVAR+NTRAIT+2).EQ.0) THEN
            FIRSTNP=J
            GO TO 95
           ENDIF
 94       CONTINUE
 95       CONTINUE
         ENDIF
        ENDIF
 91    CONTINUE
      ENDIF
C
C     MOVE PEOPLE WITH MISSING VALUES TO THE END OF THE PEDIGREE.
C
      if (evdphase.eq.1) then
         do 902 i=1,nptot
            ipacked(i) = 0
 902     continue
      end if

      NPEO=0
      IEND=NPTOT+1
      DO 103 I=1,NPTOT
C
C IF WE ARE DOING VIRTUAL TRAITS, ONLY ONE TRAIT NEEDED
C
         JSTART=1
         IF (VTRAITS.GT.IGNORE) THEN
            JSTART = VTRAITS+1
            DO 96 J=1,VTRAITS
               IF (VAR((I-1)*NVAR+J).NE.ABSENT) GOTO 98
 96         CONTINUE
            GOTO 110
         ENDIF
 98      CONTINUE
C FOR VIRTUAL TRAITS, MISSING COVAR HANDLED AS MISSING TRAIT
         IF (VTRAITS.LE.IGNORE) THEN
            DO 100 J=JSTART,NVAR
               IF (VAR((I-1)*NVAR+J).EQ.ABSENT) GO TO 110
 100        CONTINUE
         ENDIF
         if (evdphase.eq.1) then
            ipacked(npeo+1)=i
         end if
         GO TO 102
C
C MISING DATA!
C FIND SOMEONE AT THE END OF THE PEDIGREE WITH ALL VALUES PRESENT
C (OR WITH AT LEAST ONE TRAIT IF VIRTUAL TRAITS, AND ALL OTHERS)
C
 110  IEND=IEND-1
C
C for standard and evd1, we do swap
C
      if (evdphase.eq.0) then
      IF (IEND.LE.I) GO TO 120
      JSTART=1
      IF (VTRAITS.GT.IGNORE) THEN
         JSTART=VTRAITS+1
         DO 116 J=1,VTRAITS
            IF (VAR((IEND-1)*NVAR+J).NE.ABSENT) GO TO 118
 116     CONTINUE
         GOTO 110
      ENDIF
 118  CONTINUE
C FOR VIRTUAL TRAITS, MISSING COVAR HANDLED AS MISSING TRAIT
      IF (VTRAITS.LE.IGNORE) THEN
         DO 130 J=JSTART,NVAR
            IF (VAR((IEND-1)*NVAR+J).EQ.ABSENT) GO TO 110
 130     CONTINUE
      ENDIF
C
C     SWITCH THE CURRENT PERSON INTO THE LOCATION AT THE END OF THE
C     PEDIGREE.
C
      RSAVE=PERSON(I)
      PERSON(I)=PERSON(IEND)
      PERSON(IEND)=RSAVE
      DO 140 J=1,2
      RSAVE=PARENT(J,I)
      PARENT(J,I)=PARENT(J,IEND)
 140  PARENT(J,IEND)=RSAVE
      ISAVE=SEX(I)
      SEX(I)=SEX(IEND)
      SEX(IEND)=ISAVE
      ISAVE=GROUP(I)
      GROUP(I)=GROUP(IEND)
      GROUP(IEND)=ISAVE
      IDSAVE=PERMID(I)
      PERMID(I)=PERMID(IEND)
      PERMID(IEND)=IDSAVE
      IDSAVE=FAMID(I)
      FAMID(I)=FAMID(IEND)
      FAMID(IEND)=IDSAVE
      DO 150 J=1,NVAR
      RSAVE=VAR((I-1)*NVAR+J)
      VAR((I-1)*NVAR+J)=VAR((IEND-1)*NVAR+J)
 150  VAR((IEND-1)*NVAR+J)=RSAVE
      else
C
C for EVD2 we don't need any variable or pedigree data for the untyped.
C the pedigree data isn't needed because we actually use stored phi2.gz
C we'll just repeat fake pedigree info from first typed person, but we
C will save the person ID so that can still be tested
C
         iun = iun + 1
         uperson(iun) = person(i)
C         uparent(1,iun) = parent(1,i)
C         uparent(2,iun) = parent(2,i)
C         usex(iun) = sex(i)
C         ugroup(iun) = group(i)
C         upermid(iun) = permid(i)
C         ufamid(iun) = famid(i)
         goto 103
      end if

 102  NPEO=NPEO+1
 103  CONTINUE

C 103 ends loop, 120 allows break outside loop

 120  CONTINUE

C
C Now handle reordering for EVD2
C   This is done differently so we can retain increasing IBDID order
C   among typed individuals, but avoiding O(N*N) data moving.  Thus
C   we only move people after final positions are determined.  In future,
C   this method may be used in other cases, such as sampledata, where
C   a rational ordering has some value.  The method below minimizes
C   data movement within the constraint of preserving IBDID order
C
      if (iun.gt.0) then
         do 912 i = 1, npeo
            if (i.eq.ipacked(i)) goto 912
C
C copy person at ipacked up to this position
C
            ip = ipacked(i)
            person(i) = person(ip)
            parent(1,i) = parent(1,ip)
            parent(2,i) = parent(2,ip)
            sex(i) = sex(ip)
            group(i) = group(ip)
            permid(i) = permid(ip)
            famid(i) = famid(ip)
            do 911 j = 1, nvar
               var((i-1)*nvar+j) = var((ip-1)*nvar+j)
 911        continue
 912     continue
C
C use real IBDID but fake pedigree data and ABSENT data for untyped
C   OK, we don't need the ABSENT data either!
C
         do 914 i = 1, iun
            ip = npeo+i
            person(ip) = uperson(i)
            parent(1,ip) = parent(1,1)
            parent(2,ip) = parent(2,1)
            sex(ip) = sex(1)
            group(ip) = group(1)
            permid(ip) = permid(1)
            famid(ip) = famid(1)
C            do 913 j = 1,nvar
C               var((ip-1)*nvar+j) = ABSENT
 913           continue
 914        continue
      end if
C
C SEE IF THERE ARE NO PEOPLE WITH DATA IN THE PEDIGREE
C
      IF (NPEO.EQ.0) THEN
         IERROR=3
         RETURN
      ENDIF
C
C     REPORT IF THERE ARE NO NON-PROBANDS (with data) IN THE PEDIGREE.
C
      IF (NPEO.LE.NPBAND) THEN
         IERROR=3
         RETURN
      END IF
C
C Return now if not selecting this pedigree
C
      CALL PEDSELECT (IPED, ISELECT)
      IF (ISELECT.EQ.0) THEN
         IERROR=3
         RETURN
      END IF

C 22222

C
C     NUMBER THE PEOPLE IN THE PEDIGREE ACCORDING TO THEIR ORDER
C     ON DATA ENTRY.  FATHER STORES THE NUMBER CORRESPONDING TO THE
C     INDIVIDUAL'S FATHER AND SIMILARLY FOR MOTHER.
C
C EVD2 bypasses these tests, uses phi2.gz
      if (evdphase.eq.0) then
      DO 160 I=1,NPTOT
      FATHER(I)=0
      MOTHER(I)=0
      DO 160 J=1,2
      PAR=PARENT(J,I)
      IF (PAR.NE.UNKNOWN) THEN
      DO 170 K=1,NPTOT
      IF (PERSON(K).EQ.PAR) THEN
C
C     PARENT J IS IN THE PEDIGREE.
C
      IF (SEX(K).EQ.1.AND.FATHER(I).EQ.0) THEN
      FATHER(I)=K
      ELSE IF (SEX(K).EQ.2.AND.MOTHER(I).EQ.0) THEN
      MOTHER(I)=K
      ELSE
C
C     THE PARENTS ARE OF THE SAME SEX.
C
      IERROR=1
      II=MOD(GROUP(I),100000)
      WRITE(UNIT3,620) II,IPED
      WRITE(CONOUT,620) II,IPED
      END IF
      GO TO 160
      END IF
 170  CONTINUE
C
C     PARENT J IS NOT IN THE PEDIGREE.
C
      IERROR=1
      II=MOD(GROUP(I),100000)
      WRITE(UNIT3,630) J,II,IPED
      WRITE(CONOUT,630) J,II,IPED
      END IF
 160  CONTINUE  
      endif
C
C     CHECK WHETHER ANY PAIR OF PEOPLE SHARE AN INDIVIDUAL ID.
C
      if (evdphase.gt.0) then
C EVD2 uses a fast version of test because ID's are compacted and
C we can ignore untyped people
         if (npeo.gt.1) then
            ii = mod(group(1),100000)
            do 179 i=2,npeo
               jj= mod(group(i),100000)
               if (ii.eq.jj) then
                  ierror=1
                  WRITE(UNIT3,640) II,JJ,IPED
                  WRITE(CONOUT,640) II,JJ,IPED
               end if
               ii = jj
 179        continue
         end if
C
C Write evddata.out file, and return here, skipping the remaining tests
C
         call evdout (var,mvar,npeo,nvar,iped,vtraits,sex,mibdid
     *,ifam,conout,EVDstarted)
         return
      end if
C
C use old slow version of duplicate ID test
C
      DO 180 I=1,NPTOT
      PER=PERSON(I)
      DO 180 J=I+1,NPTOT
      IF (PERSON(J).EQ.PER) THEN
      IERROR=1
      II=MOD(GROUP(I),100000)
      JJ=MOD(GROUP(J),100000)
      WRITE(UNIT3,640) II,JJ,IPED
      WRITE(CONOUT,640) II,JJ,IPED
      END IF
 180  CONTINUE
C
C     CHECK WHETHER THE PEDIGREE HAS ANY DIRECTED CYCLES OR LOOPS
C     BY CALLING LOOP.
C
      IF (IERROR.NE.0) RETURN
      CALL LOOP(FATHER,MOTHER,PERM,IERROR,NPTOT,NPTOT)
      IF (IERROR.NE.0) THEN
      IERROR=MOD(GROUP(IERROR),100000)
      WRITE(UNIT3,650) IERROR,IPED
      WRITE(CONOUT,650) IERROR,IPED
      IERROR=1
      RETURN
      END IF
C
C     INCORPORATE THE SEX CODE AS PART OF GROUP.
C
      DO 190 I=1,NPTOT
 190  IF (SEX(I).EQ.1) GROUP(I)=-GROUP(I)
C
C If requested, write sample data to file
C
      IF (SAMPLEDAT) THEN
C        OPEN (UNIT=25, FILE='sampledata.out',ACCESS='APPEND')
C        OPEN (UNIT=25, FILE='sampledata.out',POSITION='APPEND')
C POSITION='APPEND' Broken in current gfortran
C Old fortran workaround follows
         UNITNO=25
         OPEN (UNITNO, FILE='sampledata.out')
         CALL ENDPOS (UNITNO)
C        print *,"NPEO is ",NPEO," out of ",NPTOT
         DO 210 I=1,NPEO
            LBASE = 1+(NVAR*(I-1))
            WRITE (25,211) PERMID(I)(1:LNBLNK(PERMID(I)))
     1           ,FAMID(I)(1:LNBLNK(FAMID(I))),(VAR(L),L=LBASE
     2           ,LBASE+NVAR-1)
 210     CONTINUE
 211     FORMAT (1X,A,',',A,20000(',',E25.16))
         CLOSE (25)
      END IF
C
C
C IF NOT VIRTUAL MULTIVARIATE, WRITE THE DATA TO THE SCRATCH FILE
C
      IF (VTRAITS.LE.0) THEN
         CALL ISCRAT (NPBAND,1,UNIT1,.TRUE.)
         CALL ISCRAT (NPEO,1,UNIT1,.TRUE.)
         CALL ISCRAT (NPTOT,1,UNIT1,.TRUE.)
         CALL ISCRAT(FATHER,NPTOT,UNIT1,.TRUE.)
         CALL ISCRAT(GROUP,NPTOT,UNIT1,.TRUE.)
         CALL ISCRAT(MOTHER,NPTOT,UNIT1,.TRUE.)
         CALL RSCRAT(VAR,MVAR,UNIT1,.TRUE.)
      ELSE
C
C If Multivariate, we inflate the data into person-traits
C and make everybody a "founder" (we use external phi2 matrix)
C
C The first "Trait" is really the trait number
C The real Trait is the second trait
C
         MVINDEX = 0
         MVPEO =  0
         MVPROB = 0

         DO 230 I = 1,NPEO
            IVARZ = (I-1)*NVAR
            DO 220 J = 1,VTRAITS
               ITIN = IVARZ + J
               IF (VAR(ITIN).NE.ABSENT.AND.(UNBAL.GT.0D0.OR.
     1              VAR(ITIN+(3-2*J)).NE.ABSENT)) THEN
C
C Include this person-trait
C
                  MVPEO = MVPEO + 1
                  IF (I.LE.NPBAND) THEN
                     MVPROB = MVPROB + 1
                  ENDIF
C
C Make everybody a "founder" (we use external phi2 matrix)
C Keep the upper part of group, recount the lower
C
                  MVFATHER(MVPEO) = 0
                  MVMOTHER(MVPEO) = 0
C
C Group is tricky.  Recover sex and upper part of group (the "group" part)
C Re-create group value using new person-trait index, MVPEO
C
                  IGROUP = GROUP(I)
                  ISEX = 0
                  IF (IGROUP.LT.0) THEN
                     ISEX = 1
                     IGROUP = -IGROUP
                  ENDIF
                  IGROUP = IGROUP - MOD(IGROUP,100000)
                  IGROUP = IGROUP + MVPEO
                  IF (ISEX.EQ.1) THEN
                     IGROUP = - IGROUP
                  ENDIF 
                  MVGROUP(MVPEO) = IGROUP
C
C Save trait number as first variable
C Save actual trait as second variable
C
                  MVINDEX = MVINDEX + 1
                  MVVAR(MVINDEX) = J
                  MVINDEX = MVINDEX + 1
                  MVVAR(MVINDEX) = VAR(ITIN)
C
C Fill the rest of the trait slots with zeros
C
                  IF (VTRAITS.GT.2) THEN
                     DO 223 K = 3,VTRAITS
                        MVINDEX = MVINDEX + 1
                        MVVAR(MVINDEX) = 0D0
 223                 CONTINUE
                  ENDIF
C
C Copy the remaining variables
C
                  DO 225 K = VTRAITS+1, NVAR
                     MVINDEX = MVINDEX + 1
                     MVVAR(MVINDEX) = VAR(IVARZ+K)
 225              CONTINUE
               ENDIF
 220        CONTINUE
 230     CONTINUE
C
C Copy the inflated data into the "scratchfile"
C
         CALL ISCRAT (MVPROB,1,UNIT1,.TRUE.)
         CALL ISCRAT (MVPEO, 1,UNIT1,.TRUE.)
         CALL ISCRAT (MVPEO, 1,UNIT1,.TRUE.)
C
C In MV everyone appears as a "founder." No additional individuals
C are needed to complete pedigrees.  *PEO (sample) and *TOT (total)
C are the same
C
         CALL ISCRAT (MVFATHER,MVPEO,  UNIT1,.TRUE.)
         CALL ISCRAT (MVGROUP, MVPEO,  UNIT1,.TRUE.)
         CALL ISCRAT (MVMOTHER,MVPEO,  UNIT1,.TRUE.)
         CALL RSCRAT (MVVAR,   MVINDEX,UNIT1,.TRUE.)
      ENDIF

      RETURN
      END

      SUBROUTINE ENDPOS (UNITNO)
      INTEGER UNITNO,NEXT
      CHARACTER*1 CHA
      DO 99992 NEXT = 1,999999999
         READ (UNITNO,99991,END=99993) CHA
99991    FORMAT (A1)
99992 CONTINUE
99993 BACKSPACE(UNITNO)
      RETURN
      END

