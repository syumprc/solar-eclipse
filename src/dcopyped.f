C
C File:       dcopyped.f
C 
      SUBROUTINE DCOPYPED (UNIT1,NPED,NVAR,NASCER,NIND,VARDATA,NTOT
     1,RARRAY,RLEN,NCUMPED,IARRAY,MALE,MAXPEO,VMEAN,VMIN,NTRAIT,IDISC
     2,TRTTYPE)

C This is the "discrete trait" version of copyped.
C Instead of copying data into a scratch file, it is copied into an array.
C Also, the order of probands is reversed (put at the end, instead of the
C beginning) since that is the way the discrete code wants it.

      INTEGER UNIT1,NPED,NVAR,NTOT
      INTEGER*8 RLEN
      INTEGER NASCER(NPED),NIND(NPED),NCUMPED(0:NPED)
      INTEGER IPED,IINDEX,NPBAND,NPEO,NPTOT,MALE(NTOT)
      INTEGER IARRAY(MAXPEO)
      INTEGER I, J, K, NREAD, NP, IENDZONE, MTEMP,NTRAIT,IDISC(NTRAIT)
      INTEGER ITRAIT,ITOFF,ISDISC,ITRT,JTRT,ITEMP,NWRITE,ID,IC
      DOUBLE PRECISION VARDATA(NVAR,NTOT),RARRAY(RLEN),AFFECT,VTEMP
      DOUBLE PRECISION VMEAN(NVAR),VMIN(NVAR),PREVALENCE,ORDER,FNP
      DOUBLE PRECISION SPREV(2*NTRAIT),TEMP
      INTEGER ST(2*NTRAIT),ITIN,IU,ITN,ITA,NSTARTING,NNSTARTING
      CHARACTER*1 TRTTYPE(NTOT),CTRT
      LOGICAL AFIRST, AREDISC, ARECONT

      CALL DOPTION ("DiscreteOrder", ORDER)
      IF (ORDER.GT.2.OR.ORDER.LT.-2) THEN
         OPEN (20, FILE='discrete.out')
      ENDIF
C
C Determine presence of discrete and continuous traits
C

      AREDISC = .FALSE.
      ARECONT = .FALSE.
      DO 3 ITRT = 1,NTRAIT
         IF (IDISC(ITRT).NE.0) THEN
            AREDISC = .TRUE.
         ELSE
            ARECONT = .TRUE.
         ENDIF
 3    CONTINUE

C
C Now, read through data and set up arrays
C            
      CALL RESCRATCH (UNIT1)
      ITOFF = 0
      IINDEX = 0
      NCUMPED(0) = 0
      DO 5000 IPED = 1, NPED
         NSTARTING = 1
         AFFECT=0.0
         CALL ISCRAT (NPBAND,1,UNIT1,.FALSE.)
         CALL ISCRAT (NPEO,1,UNIT1,.FALSE.)
         CALL ISCRAT (NPTOT,1,UNIT1,.FALSE.)

C NPEO is all typed (including probands), NPTOT is all typed or not
C Note: We do not copy non-typed individuals since optima does not
C need them

         NASCER(IPED) = NPBAND
         NIND(IPED) = NPEO
C Number of typed non-probands
         NP = NPEO - NPBAND
         NCUMPED(IPED) = NPEO + NCUMPED(IPED-1)

C The order of scratch file data is Father, Group, Mother, and Var

C We don't need anything from Father data (we're using stored phi2.gz)
         CALL ISCRAT(IARRAY,NPTOT,UNIT1,.FALSE.)

C Save sex from Group data into 'Male' array
C Put probands in the back, as will be required

         CALL ISCRAT(IARRAY,NPTOT,UNIT1,.FALSE.)
         IF (NP.GT.0) THEN
            DO 25 I = 1, NP
               IF (IARRAY(I+NPBAND).LT.0) THEN
                  MALE(IINDEX+I) = 1
               ELSE
                  MALE(IINDEX+I) = 0
               END IF
 25         CONTINUE
         END IF
         IF (NPBAND.GT.0) THEN
            DO 26 I = 1, NPBAND
               IF (IARRAY(I).LT.0) THEN
                  MALE(IINDEX+NP+I) = 1
               ELSE
                  MALE(IINDEX+NP+I) = 0
               END IF
 26         CONTINUE
         END IF

C We don't need anything from Mother data
         CALL ISCRAT(IARRAY,NPTOT,UNIT1,.FALSE.)

C Read all variable data in, then parcel out
C  (structure of data file requires this apparently)

         NREAD = NPTOT * NVAR
         CALL RSCRAT (RARRAY, NREAD, UNIT1, .FALSE.)

C Copy variable data for probands into back of array
C   To accomplish required ordering, copy continuous traits to the
C     beginning of the end, and copy discrete traits starting from the
C     end of the end.
         
         K = 1
         IC = 1
         ID = 0
         IF (NPBAND .GT. 0) THEN
            DO 35 I = 1, NPBAND

C Determine trait number
C   Then determine whether this trait is discrete

               ITRAIT = 1
               ITOFF = 0
               IF (NTRAIT.GT.1) THEN
                  ITOFF = 1
                  ITRAIT = RARRAY (K)
                  ISDISC = IDISC (ITRAIT)
                  IF (ISDISC.EQ.0) THEN
                     VARDATA(1, IINDEX+NP+IC) = RARRAY (K)
                  ELSE
                     VARDATA(1, IINDEX+NPEO-ID) = RARRAY (K)
                  END IF
               ENDIF
               ISDISC = IDISC (ITRAIT)

C Save trait type
C If discrete, adjust trait to minimum value (becomes 0/1)

               IF (ISDISC.EQ.0) THEN
                  TRTTYPE(NP+IINDEX+IC) = 'c'
                  VARDATA(1+ITOFF, NP+IINDEX+IC) = RARRAY(K+ITOFF)
                  IC = IC + 1
               ELSE
                  TRTTYPE(IINDEX+NPEO-ID) = 'd'
                  VARDATA(1+ITOFF,IINDEX+NPEO-ID)=RARRAY(K+ITOFF)-
     1                 VMIN(ITRAIT)
                  ID = ID + 1
                  AFFECT = AFFECT + RARRAY(K+ITOFF) - VMIN(ITRAIT)
               ENDIF

               K = K + 1 + ITOFF

               IF (2+ITOFF.LE.NVAR) THEN
                  DO 30 J = 2+ITOFF, NVAR
                     IF (ISDISC.EQ.0) THEN
                        VARDATA(J, NP+IINDEX+IC-1) = RARRAY(K)
                     ELSE
                        VARDATA(J, NPEO+IINDEX+1-ID) = RARRAY(K)
                     END IF
                     K = K + 1
 30               CONTINUE
               ENDIF

C Note: MALE status for probands was already copied to proband section

 35         CONTINUE
         END IF

C Copy variable data for non-probands into front of array

         IF (NP .GT. 0) THEN
            DO 45 I = 1, NP

C Determine trait number
C Then determine whether this trait is discrete

               ITRAIT = 1
               ITOFF = 0
               IF (NTRAIT.GT.1) THEN
                  ITOFF = 1
                  ITRAIT = RARRAY(K)
                  VARDATA(1, IINDEX+I) = RARRAY(K)
               ENDIF
               ISDISC = IDISC (ITRAIT)

C Save trait type
C If discrete, adjust trait to minimum value (becomes 0/1)

               IF (ISDISC.EQ.0) THEN
                  TRTTYPE(IINDEX+I) = 'c'
                  VARDATA(1+ITOFF, IINDEX+I) = RARRAY(K+ITOFF)
               ELSE
                  TRTTYPE(IINDEX+I) = 'd'
                  VARDATA(1+ITOFF,IINDEX+I)=RARRAY(K+ITOFF) -
     1                 VMIN(ITRAIT)
                  AFFECT = AFFECT + RARRAY(K+ITOFF) - VMIN(ITRAIT)
               ENDIF

C Copy remaining variables

               K = K + 1 + ITOFF

               IF (2+ITOFF.LE.NVAR) THEN
                  DO 40 J = 2+ITOFF, NVAR
                     VARDATA(J, IINDEX+I) = RARRAY(K)
                     K = K + 1
 40               CONTINUE
               ENDIF

 45         CONTINUE
         END IF

C If any trait is discrete, we must do ordering.
C
C Ordering puts the continuous trait(s), if any, first, since they
C   have exact likelihood calculations.  Then, for each discrete trait,
C   affected and unaffected groups are ranked by p value, with the highest
C   p value first.

C If there's only one trait, use original single-trait ordering code
C (Which has more options, now considered obsolescent)

         IF (0.0 .EQ. ORDER) GO TO 4990

         IF (NTRAIT.EQ.1.AND.ORDER.LE.5) THEN
            IF (IDISC(1).NE.0) THEN

C Determine ordering of non-probands: affecteds first or last
C   If prevalence < 0.5, put affecteds first
C   If prevalence >= 0.5, put unaffecteds first
C   This is the default for DiscreteOrder = 1
C Reverse above if DiscreteOrder < 0
C No sorting if DiscreteOrder = 0
C Per-pedigree sorting for -2 and +2


               PREVALENCE = VMEAN(1)-VMIN(1)

               IF (ORDER.EQ.2.0  .OR.  ORDER.EQ.-2.0) THEN
                  FNP = NP
                  PREVALENCE = AFFECT / FNP
               ENDIF

               AFIRST = .TRUE.
               IF (PREVALENCE .GE. 0.5) THEN
                  AFIRST = .FALSE.
               ENDIF

               IF (ORDER .LT. 0) THEN
                  IF (AFIRST) THEN
                     AFIRST = .FALSE.
                  ELSE
                     AFIRST = .TRUE.
                  ENDIF
               ENDIF

               IF (.FALSE.) THEN
                  IF (AFIRST) THEN
C                    PRINT *,"AFFECTED FIRST"
                  ELSE
C                    PRINT *,"UNAFFECTED FIRST"
                  ENDIF
               ENDIF

C Sort non-probands into optimal ordering for integration

               IENDZONE = NP
               DO 80 I = 1, NP

C If we're in the endzone already, quit
                  IF (I.GT.IENDZONE) GOTO 4990

C See if this record needs moving
                  IF ((VARDATA(1,IINDEX+I).EQ.0 .AND. AFIRST) .OR.
     1             (VARDATA(1,IINDEX+I).GT.0 .AND. (.NOT.AFIRST))) THEN

C Swap this record with first record from the end having
C   opposite status
                     DO 70 J = IENDZONE, I+1, -1
                        IF ((VARDATA(1,IINDEX+J).GT.0 .AND. AFIRST) .OR.
     1              (VARDATA(1,IINDEX+J).EQ.0 .AND. (.NOT.AFIRST))) THEN

C These records can be swapped
                           DO 60 K = 1, NVAR
                              VTEMP = VARDATA(K,IINDEX+I)
                              VARDATA(K,IINDEX+I) = VARDATA(K,IINDEX+J)
                              VARDATA(K,IINDEX+J) = VTEMP
 60                        CONTINUE

C Swap MALE (with group info) also
                           MTEMP = MALE(IINDEX+I)
                           MALE(IINDEX+I) = MALE(IINDEX+J)
                           MALE(IINDEX+J) = MTEMP

C Swap TRTTYPE
                           CTRT = TRTTYPE(IINDEX+I)
                           TRTTYPE(IINDEX+I) = TRTTYPE(IINDEX+J)
                           TRTTYPE(IINDEX+J) = CTRT

C Done! Update endzone pointer and return to main relocation loop
                           IENDZONE = J
                           GO TO 80

C Loop and try previous record for swapping
                        ENDIF
 70                  CONTINUE

C Fell through without finding record to swap
C This means we're done swapping!
                     GO TO 85

C End of "this individual needs to be moved" clause
                  ENDIF

C End of affected/unaffected relocation loop
 80            CONTINUE
 85            CONTINUE
               GOTO 4990
            ENDIF

         ELSE
C There is more than one trait

C Now consider other possibilities
C IF no discrete trait, then no sorting required either

         IF (.NOT.AREDISC) GOTO 4990
         
C Move continuous traits to the front by moving other traits to end

         IF (ARECONT) THEN
            DO 150 ITRT = 1, NTRAIT

               IENDZONE = NP
               NNSTARTING = NSTARTING

               IF (IDISC(ITRT).EQ.0) THEN

                  DO 180 I = NSTARTING, NP
                     IF (VARDATA(1,IINDEX+I).EQ.ITRT) THEN
                        NNSTARTING = I + 1
                     ELSE

C Swap this record with first record from end having desired trait

                        DO 170 J = IENDZONE, I+1, -1
                           IF (VARDATA(1,IINDEX+J).EQ.ITRT) THEN
                              DO 160 K = 1, NVAR
                                 VTEMP = VARDATA(K,IINDEX+I)
                                 VARDATA(K,IINDEX+I)=VARDATA(K,IINDEX+J)
                                 VARDATA(K,IINDEX+J) = VTEMP
 160                          CONTINUE

C Swap MALE and TRTTYPE too

                              MTEMP = MALE(IINDEX+I)
                              MALE(IINDEX+I) = MALE(IINDEX+J)
                              MALE(IINDEX+J) = MTEMP
                              CTRT = TRTTYPE(IINDEX+I)
                              TRTTYPE(IINDEX+I) = TRTTYPE(IINDEX+J)
                              TRTTYPE(IINDEX+J) = CTRT
                              IENDZONE = J
                              NNSTARTING = I + 1
                              GOTO 180
                           ENDIF
 170                    CONTINUE

C If we get here, we can't find any more candidates for swapping...
C So we have completely copied continuous trait Cn to the front

                        NSTARTING = NNSTARTING
                        GOTO 150
                        
                     ENDIF
 180              CONTINUE
                  NSTARTING=NNSTARTING
               ENDIF
 150        CONTINUE
         ENDIF

C Create unsorted compressed index of discrete trait "probabilities"
C    So that, for trait i, the nth discrete trait,
C     SPREV((n-1)*2+1) = probability of affect
C     SPREV((n-1)*2+2) = probability of not affected
C     ST((n-1)*2+1) = 2*i
C     ST((n-1)*2+2) = (i-1)*2+2
C The ST variable is used to remember the original uncompressed index
C   after sorting
    
         ITIN = 1
         DO 210 ITRT = 1, NTRAIT
            IF (IDISC(ITRT).NE.0) THEN
               SPREV(ITIN) = VMEAN(ITRT)-VMIN(ITRT)
               SPREV(ITIN+1) = 1.0 - SPREV(ITIN)
               ST(ITIN) = (2*ITRT)-1
               ST(ITIN+1) = 2*ITRT
               ITIN = ITIN + 2
            ENDIF
 210     CONTINUE
         ITIN = ITIN - 1

C Now, reverse sort this index so that st reflects trait indexing
C  Trait index is (ITIN-1)*2+AFF where ITIN is 1-based trait
C  index of discrete trait in array of discrete traits only,
C  and AFF=1 means the affected group.

         DO 230 ITRT = 1, ITIN-1
            DO 220 JTRT = ITRT+1, ITIN
               IF ((ORDER.GT.0.AND.SPREV(ITRT).GT.SPREV(JTRT)).OR.
     1             (ORDER.LT.0.AND.SPREV(ITRT).LE.SPREV(JTRT))) THEN
                   TEMP = SPREV(ITRT)
                   SPREV(ITRT) = SPREV(JTRT)
                   SPREV(JTRT) = TEMP
                   ITEMP = ST(ITRT)
                   ST(ITRT) = ST(JTRT)
                   ST(JTRT) = ITEMP
                END IF
 220         CONTINUE
 230      CONTINUE

C Now, procede to move indicated traits to front

          DO 300 ITRT = 1, ITIN-1
            ITN = (ST(ITRT)+1)/2
            IU = ST(ITRT)/2
            ITA = ST(ITRT)-2*IU

C Move to front: 
C    Trait ITN, Status ITA (1 affected, 0 unaffected--as stored)

            NNSTARTING = NSTARTING
            IENDZONE = NP
            DO 280 I = NSTARTING, NP
               IF (VARDATA(1,I+IINDEX).EQ.ITN.AND.
     1              VARDATA(2,I+IINDEX).EQ.ITA) THEN
                  NNSTARTING = I + 1
               ELSE
                  DO 270 J = IENDZONE , I+1, -1

                     IF (VARDATA(1,J+IINDEX).EQ.ITN.AND.
     1                    VARDATA(2,J+IINDEX).EQ.ITA) THEN

                        DO 260 K = 1, NVAR
                           VTEMP = VARDATA(K,IINDEX+I)
                           VARDATA(K,IINDEX+I)=VARDATA(K,IINDEX+J)
                           VARDATA(K,IINDEX+J)=VTEMP
 260                    CONTINUE
                        MTEMP = MALE(IINDEX+I)
                        MALE(IINDEX+I) = MALE(IINDEX+J)
                        MALE(IINDEX+J) = MTEMP
                        CTRT = TRTTYPE(IINDEX+I)
                        TRTTYPE(IINDEX+I)=TRTTYPE(IINDEX+J)
                        TRTTYPE(IINDEX+J)=CTRT
                        IENDZONE = J
                        NNSTARTING = I + 1
                        GOTO 280
                     ENDIF
 270              CONTINUE
C
C If we get here, it means that there was nothing in back to swap
C to the front, so we've done all we can do for this trait.
C 
                  NSTARTING = NNSTARTING
                  GOTO 300
               ENDIF
 280        CONTINUE
            NSTARTING=NNSTARTING
 300     CONTINUE

C NOW, sort discrete probands in the same way
C Remember, iindex+npeo+1-id is first discrete proband if id > 0

         IF (ID.GT.1) THEN
            IENDZONE = NPEO
            NSTARTING = NPEO+1-ID
            NNSTARTING = NSTARTING

            DO 400 ITRT = 1, ITIN-1
               ITN = (ST(ITRT)+1)/2
               IU = ST(ITRT)/2
               ITA = ST(ITRT)-2*IU

C Move to front: 
C    Trait ITN, Status ITA (1 affected, 0 unaffected--as stored)

               NNSTARTING = NSTARTING
               IENDZONE = NPEO
               DO 380 I = NSTARTING, NPEO
                  IF (VARDATA(1,I+IINDEX).EQ.ITN.AND.
     1                 VARDATA(2,I+IINDEX).EQ.ITA) THEN
                     NNSTARTING = I + 1
                  ELSE
                     DO 370 J = IENDZONE , I+1, -1

                        IF (VARDATA(1,J+IINDEX).EQ.ITN.AND.
     1                       VARDATA(2,J+IINDEX).EQ.ITA) THEN

                           DO 360 K = 1, NVAR
                              VTEMP = VARDATA(K,IINDEX+I)
                              VARDATA(K,IINDEX+I)=VARDATA(K,IINDEX+J)
                              VARDATA(K,IINDEX+J)=VTEMP
 360                       CONTINUE
                           MTEMP = MALE(IINDEX+I)
                           MALE(IINDEX+I) = MALE(IINDEX+J)
                           MALE(IINDEX+J) = MTEMP
                           CTRT = TRTTYPE(IINDEX+I)
                           TRTTYPE(IINDEX+I)=TRTTYPE(IINDEX+J)
                           TRTTYPE(IINDEX+J)=CTRT
                           IENDZONE = J
                           NNSTARTING = I + 1
                           GOTO 380
                        ENDIF
 370                 CONTINUE
C
C If we get here, it means that there was nothing in back to swap
C to the front, so we've done all we can do for this trait.
C 
                     NSTARTING = NNSTARTING
                     GOTO 400
                  ENDIF
 380           CONTINUE
               NSTARTING=NNSTARTING
 400        CONTINUE

C End if more than one discrete probands to sort
         END IF

C End handling more than one trait
      ENDIF

 4990 CONTINUE

C  Copy discrete data to a file for examination
C  Pedigrees are separated by two blank lines
C  Probands (if any) are separated by one blank line from rest of pedigree
C  trttype and male status begin each record

 1301 FORMAT (1X,A1,4X,I1,2X,(F10.3))
 1305 FORMAT (7(F10.3))
 1302 FORMAT (" ")
 1303 FORMAT ("TYPE SEX VARS")

      IF (ABS(ORDER).GT.2) THEN
         NWRITE = NVAR
C         IF (NWRITE.GT.7) THEN
C            NWRITE = 7
C         END IF
         IF (IPED.EQ.1.AND.ABS(ORDER).LT.4) THEN
            WRITE (20,1303)
         END IF
         DO 1310 I = IINDEX+1, IINDEX+NP
            IF (ABS(ORDER).LT.4) THEN
               WRITE (20, 1301) TRTTYPE(I),MALE(I),
     1              (VARDATA(K,I),K=1,NWRITE)
            ELSE
               WRITE (20, 1305) (VARDATA(K,I),K=1,NWRITE)
            END IF
 1310     CONTINUE

         IF (NPBAND.GT.0) THEN
            WRITE (20,1302)
            DO 1320 I = IINDEX+NP+1, IINDEX+NPEO
               IF (ABS(ORDER).LT.4) THEN
                  WRITE (20, 1301) TRTTYPE(I),MALE(I),
     1                 (VARDATA(K,I),K=1,NWRITE)
               ELSE
                  WRITE (20, 1305) (VARDATA(K,I),K=1,NWRITE)
               END IF
 1320        CONTINUE
         END IF
 
         WRITE (20,1302) 
         WRITE (20,1302) 
      END IF

C Procede to next pedigree

      IINDEX = IINDEX + NPEO
 5000 CONTINUE


C     READ VMEAN AND VVAR INTO THE FIRST TWO SECTIONS OF RARRAY.
C     THESE WILL BE PASSED TO THE REST OF THE PROGRAM AND REAPPEAR
C     WITH THEIR PROPER NAMES IN SUBROUTINE SEARCH.

      CALL RSCRAT (RARRAY, NVAR, UNIT1, .FALSE.)
      CALL RSCRAT (RARRAY(NVAR+1), NVAR, UNIT1, .FALSE.)

      IF (ORDER.GT.2.OR.ORDER.LT.-2) THEN
         CLOSE (20)
      END IF

      END
