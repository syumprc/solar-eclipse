C Modified 1996-1997 by Charles Peterson, SFBR.
C This is the New Fisher signature.  Do not remove.  fI6%.
C
      SUBROUTINE PREPED(RARRAY,VMAX,VMEAN,VMIN,VVAR,IARRAY,CARRAY
     1,ABSENT,LENI,LENR,MAXPEO,MXTWIN,NIPED,TSAMP,NVAR,UNIT1
     2,UNIT7,UNIT3,FRMT1,FRMT2,ECHO,UNIT4,TITLE
     3,PEDFIL,MODFIL,CONOUT,NV2READ,NTRAIT,NCOVAR
     4,VNAMES,TSAMPA,IQUIT,VBINARY,PEDFILE,PHENFILE
     5,TPROB,TSAMPIP,FIRSTPID,ALLPED,ERRMSG,MSGLNB,TPEOP,RAWSD,RAWMU
     6,WHO,VTRAITS,NTRASAMP,TTRASAMP,SAMPLEDAT)
C
C     THIS SUBROUTINE CONTROLS THE READING OF THE PEDIGREE DATA
C     FROM THE PEDIGREE INPUT FILE.  ONCE THE DATA IS READ AND
C     CHECKED, IT IS WRITTEN ON A SCRATCH FILE.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER*8 LENR
      DOUBLE PRECISION RARRAY(*),VMAX(NVAR),VMEAN(NVAR),VMIN(NVAR)
     1,VVAR(NVAR),VDELTA,RAWSD,RAWMU
      INTEGER IARRAY(LENI),UNIT1,UNIT7,UNIT3,UNIT4,CONOUT,NV2READ
     1,NTRAIT,NCOVAR,LTITLE,TTITLE,ALLPED,WHO,NTRASAMP(*),TTRASAMP(*)
     2,TSAMPA(NVAR),IQUIT,IBDID,VTRAITS
      CHARACTER*8 CARRAY(*)*(*),FRMT1*(*),FRMT2*(*)
     1,TITLE*40,PEDFIL*40,MODFIL*40,SDATE*24
      CHARACTER VNAMES(*)*(*)
      CHARACTER*18 FIRSTPID(ALLPED), THISPID
      LOGICAL ECHO,QUIT,VSAMPLE,VERBOSE,VBINARY(NVAR),NEXTPEDINFO
      LOGICAL IFPED,IFBINARY,PEDXO
      LOGICAL*4 SAMPLEDAT
      CHARACTER*256 PEDFILE
      CHARACTER*256 PEDTITLE,SUBTITLE
      CHARACTER*1024 PHENFILE,PHENTITLE
      CHARACTER*1 ASTER
      CHARACTER*96 VARNAM
      CHARACTER*96 SPACES
      CHARACTER*120 PEDSS,CPEDSS
      INTEGER VNL,PEDNUM,NIPED,NPEOP,TPEOP,NMALES,TMALES,NPROB,TPROB
      INTEGER NTWINS,TTWINS,NSAMP,TSAMP,NSAMPIP,TSAMPIP,SKIPPED
      INTEGER SKIPPEOP,SKIPMALE,SKIPTWIN,SKIPPROB,LSTR,LNBLNK
      INTEGER M2,M3,M4,M5,M6,M7,M8,M9,MEND,MSIZE,ISPACES,VARLEN
      INTEGER N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,NSIZE,LONGNAME,NEWNAME
      integer allsped,nimatrix
      CHARACTER*(*) ERRMSG
      INTEGER MSGLNB,evdphase,mibdid,ifam
      LOGICAL THERE,EVDstarted
C
C     FORMAT STATEMENTS FOR OUTPUT OF DATA AND ERROR MESSAGES.
C
 510  FORMAT(A)
 511  FORMAT(I7)
 520  FORMAT(' Formats for Input of Phenotype Data:',2(/,' ',A))
 540  FORMAT(' *** ERROR *** Pedigree Number',I4,' has',I6,' People.'
     1,'  The Allowed Range is'/' 1 TO 99999.')
 580  FORMAT(/,' Solar has stopped because of errors in the'
     1,' phenotype file.  These'/' are noted above.')
 581  FORMAT(/,'Errors in phenotype file.  File solar.out may have'
     1,' more details.')
 591  FORMAT(T31,'Pedigrees Included',/,T31,'--------- --------',/
     1,T12,'Pedigrees  People   Females   Males    MZTwins'
     2,' Probands',/,12X,I5,1X,5(4X,I5)/)
 592  FORMAT(/,T31,'Pedigrees Excluded',/,T31,'--------- --------',/
     1,T12,'Pedigrees  People   Females   Males    MZTwins'
     2,' Probands',/,12X,I5,1X,5(4X,I5)/)
 593  FORMAT (T27,'No pedigrees were skipped.',/)
 594  FORMAT (T8
     1,'Incomplete data in pedigree(s) which include these individual'
     2,' IDs:')
 595  FORMAT (T6,A)
 596  FORMAT (A,' ',A)
 597  FORMAT ('No individuals have complete data')
 598  FORMAT 
     1   ('All pedigrees excluded because no non-probands had data')
 599  FORMAT 
     1   (/'No complete sample with covariates and trait ',A)
 601  FORMAT ('   pedno')
 602  FORMAT (1X,I7)
 640  FORMAT(1X,A,A1,F12.6,1X,F12.6,2X,F12.6,1X,F12.6,1X,I6)
 650  FORMAT(T13
     1,'Descriptive Statistics for the Variables (Non-Probands)'
     1//'           Variable       Mean       Std Dev'
     1,'       Minimum      Maximum   Avail')
 651  FORMAT(T16,I7,' pedigrees skipped with',I7,' people',/)
C
C     INITIALIZE SOME VARIABLES.
C
      EVDstarted = .false.
      mibdid = 0
      PEDXO=.FALSE.
      PEDNUM=0
      IFBINARY=.FALSE.
      VSAMPLE = VERBOSE ("SAMPLE")
      DO 10 I=1,NVAR
         VMAX(I)=-1.0D20
         VMEAN(I)=0.0D0
         VMIN(I)=1.0D20
         VVAR(I)=0.0D0
         VBINARY(I)=.TRUE.
 10   CONTINUE

C Clear out total accumulators (for all families)

      nimatrix=0
      PEDSS=' '
      SKIPPED=0
      SKIPPEOP=0
      SKIPMALE=0
      SKIPTWIN=0
      SKIPPROB=0
      NIPED=0
      TPEOP=0
      TMALES=0
      TPROB=0
      TTWINS=0
      TSAMP=0
      allsped=0
      TSAMPIP=0
      DO 15 I=1,NVAR
         TSAMPA(I)=0
 15   CONTINUE
      DO 16 I=1,NTRAIT
         TTRASAMP(I)=0
 16   CONTINUE


      QUIT=.FALSE.
      FRMT1=' '
      FRMT2=' '
C
C Open evddata.out file
C
      call ioption ("EVDPhase",evdphase)
      if (evdphase.eq.1) then
         call unitactive (0,26)
         ifam = 1
C
C write evddata.out header in C++ then re-open to append here
C
         open (unit=26, file='evddata.out')
         call unitactive (1,26)
         call endpos (26)
      end if
C
C     READ THE PEDIGREE DATA FAMILY BY FAMILY.  CHECK FOR ERRORS AND
C     AND AN END OF FILE INDICATOR.
C
 20   IFPED = NEXTPEDINFO (NPEOP)
      IF (.NOT. IFPED) GO TO 60

C     PEDNUM is a pedigree INDEX (NIPED is the count of included pedigrees)
      PEDNUM=PEDNUM+1

C
C     CHECK THAT THE NUMBER OF PEOPLE IN EACH PEDIGREE IS IN BOUNDS.
C
      IF (NPEOP.LE.0.OR.NPEOP.GE.100000) THEN
      WRITE(UNIT3,540) PEDNUM,NPEOP
      WRITE(CONOUT,540) PEDNUM,NPEOP
      IQUIT = 1
      RETURN
      ENDIF
C
C     SET START POINTS FOR ARRAYS TO BE CARVED OUT OF RARRAY.
C
C  Consistency check
      IF (NPEOP.GT.MAXPEO) THEN
         WRITE (UNIT3,541) MAXPEO,NPEOP
         WRITE (CONOUT,541) MAXPEO,NPEOP
 541     FORMAT ('Maximum pedigree size not consistent: ',I8,',',I8)
         IQUIT = 1
         RETURN
      ENDIF

      MVAR=NPEOP*NVAR
      M2=MXTWIN+1
      M3=2*NPEOP+M2
      M4=NPEOP+M3
      M5=MVAR+M4
      M6=NVAR+M5
      M7=NVAR+M6
      M8=NVAR+M7
      M9=NVAR+M8

C 1 is added in case VTRAITS is 0

      MEND=M9+VTRAITS*NPEOP*NVAR+1

C Additional NVAR*4 is required because of VMAX,VMEAN,VMIN,VVAR
C which were actually carved out of beginning of RARRAY.
C 1 is subtracted because it is a 1-based array index.

      MSIZE=MEND+(NVAR*4)-1

C Undersize should not happen anymore
      IF (MSIZE.GT.LENR) THEN
         WRITE (UNIT3,5001)
         WRITE (CONOUT,5001)
 5001    FORMAT (
     1'Inadequate Real allocation; contact solar@txbiomedgenetics.org')
         IQUIT = 1
         RETURN
      ENDIF
C
C     SET START POINTS FOR ARRAYS TO BE CARVED OUT OF IARRAY.
C
      ISLICE = NPEOP
C
C Force 8-byte alignment for portability
C
      IF (MOD(NPEOP,2).EQ.1) THEN
         ISLICE = NPEOP + 1
      ENDIF
      N2=ISLICE+1
      N3=ISLICE+N2
      N4=ISLICE+N3
      N5=MXTWIN+N4
      N6=ISLICE+N5
      N7=MXTWIN+N6
      N8=ISLICE+N7
      N9=NVAR+N8
      N10=ISLICE*VTRAITS+N9
      N11=ISLICE*VTRAITS+N10
C 1 is subtracted because the starting point N's are 1-based
      NSIZE=ISLICE*VTRAITS+N11-1
      IF (NSIZE.GT.LENI) THEN
         WRITE (UNIT3,5002)
         WRITE (CONOUT,5002)
 5002    FORMAT (
     1'Inadequate int allocation; contact solar@txbiomedgenetics.org')
         IQUIT = 1
         RETURN
      ENDIF
C
C     CALL PINPUT TO READ AND CHECK THE DATA ON EACH INDIVIDUAL
C     IN THE PEDIGREE.
C
      IERROR=0
      NMALES=0
      NPROB=0
      NTWINS=0

C IARRAY(N8) is the 1-pedigree copy of TSAMPA (would be NSAMPA)
      DO 46 I=1,NVAR
         IARRAY(N8+I-1)=0
 46   CONTINUE

C NSAMP and NSAMPIP are working copies (not yet committed)
C They must not start from zero (to allow proper stat calculations)

      NSAMP=TSAMP
      NSAMPIP=TSAMPIP
      DO 465 I=1,NTRAIT
 465     NTRASAMP(I)=TTRASAMP(I)
         

C Save current values of VMAX, VMEAN, VMIN, and VVAR to allow restore
      DO 47 I=1,NVAR
         RARRAY(M5+I-1) = VMAX(I)
         RARRAY(M6+I-1) = VMEAN(I)
         RARRAY(M7+I-1) = VMIN(I)
         RARRAY(M8+I-1) = VVAR(I)
 47   CONTINUE

      CALL PINPUT(RARRAY,RARRAY(M2),RARRAY(M3),RARRAY(M4),VMAX,VMEAN
     1,VMIN,VVAR,IARRAY,IARRAY(N2),IARRAY(N3),IARRAY(N4),IARRAY(N5)
     2,IARRAY(N6),IARRAY(N7),CARRAY,ABSENT,IERROR,PEDNUM,MXTWIN
     3,MVAR,NMALES,NPROB,NPEOP,NSAMP,NTWINS,NVAR,NV2READ,UNIT1,UNIT7
     4,UNIT3,FRMT2,ECHO,CONOUT,NTRAIT,NCOVAR
     5,IARRAY(N8),VBINARY,NSAMPIP,FIRSTPID(NIPED+1),NDATA,IARRAY(N9)
     6,IARRAY(N10),IARRAY(N11),RARRAY(M9),VTRAITS,NTRASAMP,SAMPLEDAT
     7,evdphase,mibdid,ifam,nimatrix,EVDstarted)
C
C count all peds that contain at least one fully typed individual
C even if that ped cannot be included because no non-probands are typed.
C This is done for error reporting purposes.
C
      if (nsampip.gt.tsampip) then
         allsped = allsped + 1
      end if
C
C No non-probands with data
C Skip this pedigree silently...not included in running totals
C
      IF (IERROR.EQ.3) THEN
C
C Restore previous values of VMAX, VMEAN, VMIN, and VVAR, and TSAMP
C
         DO 48 I=1,NVAR
            VMAX(I) = RARRAY(M5+I-1)
            VMEAN(I) = RARRAY(M6+I-1)
            VMIN(I) = RARRAY(M7+I-1) 
            VVAR(I) = RARRAY(M8+I-1)
 48      CONTINUE
         SKIPPED = SKIPPED + 1
         SKIPPEOP = SKIPPEOP + NPEOP
         SKIPMALE = SKIPMALE + NMALES
         SKIPTWIN = SKIPTWIN + NTWINS
         SKIPPROB = SKIPPROB + NPROB
C
C Show which pedigrees have been skipped
C
         IF (SKIPPED.EQ.1) THEN
            WRITE (UNIT3,510) ' '
            WRITE (UNIT3,594)
            IF (VSAMPLE) WRITE (CONOUT,594)
         ENDIF

         LSTR = LNBLNK (PEDSS)
         IF (LSTR.GT.62) THEN
            WRITE (UNIT3,595) PEDSS(1:LSTR)
            IF (VSAMPLE) WRITE (CONOUT,595) PEDSS(1:LSTR)
            LSTR = 0
         ENDIF

         CPEDSS = PEDSS
         IF (LSTR.EQ.0) THEN
            PEDSS = ' '
            LSTR = 1
         ENDIF

         THISPID = FIRSTPID(NIPED+1)
         WRITE (PEDSS,596) CPEDSS(1:LSTR),THISPID(1:LNBLNK(THISPID))
C
C If pedlike in effect, write list of excluded pedigrees to file
C pedexclude.dat
C
         CALL DOPTION ('PEDLIKE',PEDLIKE)
         IF (PEDLIKE.EQ.1) THEN
            IF (.NOT.PEDXO) THEN
               OPEN (35,FILE='pedexclude.dat')
               WRITE (35,601)
               PEDXO=.TRUE.
            ELSE
C              OPEN (35,FILE='pedexclude.dat',ACCESS='APPEND')
C              OPEN (35,FILE='pedexclude.dat',POSITION='APPEND')
C gfortran bug workaround
               INQUIRE (FILE='pedexclude.dat',EXIST=THERE)
               OPEN (35,FILE='pedexclude.dat')
               IF (THERE) THEN
                  CALL ENDPOS (35)
               END IF

            END IF
            WRITE (35,602) PEDNUM
            CLOSE (35)
         END IF
         GOTO 20
      ENDIF
C
C NOW that we have skipped unusable pedigrees, add to cumulative counts
C for pedigrees, people, males, probands, mztwins
C
      NIPED  = NIPED  + 1
      TPEOP  = TPEOP  + NPEOP
      TMALES = TMALES + NMALES
      TPROB  = TPROB  + NPROB
      TTWINS = TTWINS + NTWINS
C
C Now commit updated sample sizes
C
      TSAMP  = NSAMP
      TSAMPIP= NSAMPIP
      DO 54 I=1,NTRAIT
         TTRASAMP(I)= NTRASAMP(I)
 54   CONTINUE
      DO 55 I=1,NVAR
         TSAMPA(I)=TSAMPA(I)+IARRAY(N8+I-1)
 55   CONTINUE
C
C IERROR >= 2, Quit now, setting IQUIT error status
C
      IF (IERROR.GE.2) THEN
         IQUIT = 1
         if (evdphase.eq.1) then
            call unitactive (0,26)
         end if
         RETURN
      ENDIF
C
C IERROR >= 2, continue reading pedigree to find more errors.
C Error status will be set on return
C
      IF (IERROR.NE.0) QUIT=.TRUE.
C
C If we are saving who is included, add to file here
C
      IF (WHO.NE.0) THEN
C        OPEN (UNIT=31,FILE='who.out',ACCESS='APPEND')
C        OPEN (UNIT=31,FILE='who.out',POSITION='APPEND')
C gfortran bug workaround
         INQUIRE (FILE='who.out', EXIST=THERE)
         OPEN (UNIT=31,FILE='who.out')
         IF (THERE) THEN
            CALL ENDPOS (31)
         END IF

         DO 56 I=0,NDATA-1
            IBDID = RARRAY(M3+I)
            WRITE (31,511) IBDID
 56      CONTINUE
         CLOSE (UNIT=31)
      ENDIF
C
C This is the end of the "get next family" loop
C
      GO TO 20
C
C Finished reading pedigree.  If QUIT flag was set, return with
C IQUIT error status
C
 60   LSTR = LNBLNK (PEDSS)
      IF (LSTR.GT.0) THEN
         WRITE (UNIT3,595) PEDSS(1:LSTR)
         IF (VSAMPLE) WRITE (CONOUT,595) PEDSS(1:LSTR)
      ENDIF

      IF (QUIT) THEN
      WRITE(UNIT3,580)
      WRITE(CONOUT,581)
      IQUIT = 1
      RETURN
      END IF
C
C     WRITE HEADER FOR OUTPUT FILES
C
 700  FORMAT(A/)

C  Write "title"


      WRITE(UNIT4,510) ' '
      CALL FDATE (SDATE)
      SUBTITLE = 'Solar  '//SDATE
      LTITLE = LNBLNK (SUBTITLE)
      TTITLE = ((80-LTITLE)/2)-1
C     CALL TABOVER (UNIT3, TTITLE)
      CALL TABOVER (UNIT4, TTITLE)
      WRITE(UNIT4,510) SUBTITLE (1:LTITLE)

C  Write pedigree filename

      WRITE (UNIT4,510) ' '
      PEDTITLE = 'Pedigree:  '//PEDFILE
      LTITLE = LNBLNK (PEDTITLE)
      IF (LTITLE.LT.78) THEN
         TTITLE = ((80-LTITLE)/2)-1
         CALL TABOVER (UNIT4, TTITLE)
      ENDIF
      WRITE (UNIT4,510) PEDTITLE(1:LTITLE)

C  Write phenotype filename

      PHENTITLE = 'Phenotypes:  '//PHENFILE
      LTITLE = LNBLNK (PHENTITLE)
      IF (LTITLE.LT.78) THEN
         TTITLE = ((80-LTITLE)/2)-1
         CALL TABOVER (UNIT4, TTITLE)
      ENDIF
      WRITE (UNIT4,510) PHENTITLE(1:LTITLE)
C
C     OUTPUT DESCRIPTIVE STATISTICS FOR THE PEDIGREE DATA.
C
      IF (SKIPPED.GT.0) THEN
         WRITE (UNIT3,592) SKIPPED,SKIPPEOP,SKIPPEOP-SKIPMALE
     *        ,SKIPMALE,SKIPTWIN,SKIPPROB
      ELSE
         WRITE(UNIT3,700) ' '
         WRITE (UNIT3,593)
      ENDIF
      WRITE(UNIT3,591)  NIPED,TPEOP,TPEOP-TMALES,TMALES,TTWINS,TPROB
      if (nimatrix.ne.0) then
         WRITE (UNIT3, 645) nimatrix
      end if
      WRITE(UNIT3,650)

      IF (VSAMPLE) THEN
      IF (SKIPPED.GT.0) THEN
         WRITE (CONOUT,592) SKIPPED,SKIPPEOP,SKIPPEOP-SKIPMALE
     *        ,SKIPMALE,SKIPTWIN,SKIPPROB
      ELSE
         WRITE (CONOUT,593)
      ENDIF
      WRITE(CONOUT,591) NIPED,TPEOP,TPEOP-TMALES,TMALES,TTWINS,TPROB
      if (nimatrix.ne.0) then
         if (VSAMPLE) then
            WRITE (CONOUT, 645) nimatrix
         end if
      end if
      WRITE(CONOUT,650)
      ENDIF

      if (allsped.eq.0) then
         WRITE (UNIT3, 597)
         ERRMSG = 'No individuals have complete data.  Check that ped an
     1d phen IDs match.'
         MSGLNB = LNBLNK (ERRMSG)
         IQUIT = 1
         RETURN
      ENDIF
         

      IF (TSAMP.EQ.0) THEN
         WRITE (UNIT3, 598)
         ERRMSG = 'No non-probands have complete data.  Every ped must h
     1ave one non-proband.'
         MSGLNB = LNBLNK (ERRMSG)
         IQUIT = 1
         RETURN
      ENDIF

      IF (NIPED.EQ.0) THEN
         WRITE (UNIT3, 598)
         ERRMSG = 'No non-probands have complete data.  Every ped must h
     1ave one non-proband.'
         MSGLNB = LNBLNK (ERRMSG)
         IQUIT = 1
         RETURN
      ENDIF

      DO 69 I=1,NTRAIT
         IF (TSAMPA(I).EQ.0) THEN
            WRITE (UNIT3,599) VNAMES(I)
            ERRMSG = 'No complete sample with covariates and trait '
     1//VNAMES(I)
            MSGLNB = LNBLNK (ERRMSG)
            IQUIT = 1
            RETURN
         ENDIF
 69   CONTINUE

      LONGNAME = LNBLNK(VNAMES(1))
      DO 70 I=1,NVAR
         NEWNAME = LNBLNK(VNAMES(I))
         IF (NEWNAME.GT.LONGNAME) LONGNAME = NEWNAME
 70   CONTINUE
      IF (LONGNAME.LT.18) LONGNAME = 18

      DO 71 I=1,NVAR
      IF (I.GT.NTRAIT.AND.I.LE.NTRAIT+2) GO TO 71
      SD=SQRT(VVAR(I))
      VDELTA=VMAX(I)-VMIN(I)
C VBINARY might also be set to false in pinput if there are 
C More than 2 values
      IF (VDELTA.NE.1.0.AND.VDELTA.NE.0) THEN
         VBINARY(I) = .FALSE.
      ENDIF
      IF (VBINARY(I)) THEN
         ASTER='*'
         IFBINARY=.TRUE.
      ELSE
         ASTER=' '
      ENDIF
C     '123456789012345678901234567890123456789012345678901234567890'
      SPACES=
     1'                                                            '//
     2'                                    '
      VNL = LNBLNK(VNAMES(I))
      ISPACES = LONGNAME - VNL
      IF (ISPACES.GT.96) ISPACES = 96
      IF (ISPACES.LT.1) ISPACES = 0
      VARNAM=SPACES(1:ISPACES)//VNAMES(I)(1:VNL)
      VARLEN = LNBLNK (VARNAM)
C      if ("_evd" .eq. VARNAM(VARLEN-3:VARLEN)) then
C         VARLEN = VARLEN - 4
C      end if
      WRITE(UNIT3,640) VARNAM(1:VARLEN),ASTER,VMEAN(I),SD,VMIN(I)
     1,VMAX(I),TSAMPA(I)
      IF (VSAMPLE) THEN
      WRITE(CONOUT,640) VARNAM(1:VARLEN),ASTER,VMEAN(I),SD,VMIN(I)
     1,VMAX(I),TSAMPA(I)
      ENDIF
      CALL SAVEMSD (VNAMES(I), VMEAN(I), SD)

 71   CONTINUE
      WRITE(UNIT3,641) TSAMP
      WRITE(UNIT3,642)
      WRITE(UNIT3,644) TSAMPIP
      IF (IFBINARY) WRITE(UNIT3,643)
      IF (VSAMPLE) THEN
         WRITE(CONOUT,641) TSAMP
         WRITE(CONOUT,642)
         WRITE(CONOUT,644) TSAMPIP
         IF (IFBINARY) WRITE(CONOUT,643)
      ENDIF

C trap here if evd2 phase 1
      if (evdphase.eq.1) then
        call unitactive (0,26)
         call evdtrap (mibdid)
      end if

 641  FORMAT (/T20,'The non-proband analysis sample size is',I6,'.')
 642  FORMAT (T14
     1,'(Note: THIS is the sample size for the means shown above.)')
 643  FORMAT (T10
     1,'(Variables marked with * are binary: covar is adj to low'
     1,' value.)')
 644  FORMAT (T17,'The sample size including probands used is',I6,'.')
 645  FORMAT (T8,'** Warning!',I6
     1,' Individual(s) skipped because not in matrix. **'/)
C
C     WRITE VMEAN AND VVAR ON THE FIRST SCRATCH FILE.
C
      CALL RSCRAT(VMEAN,NVAR,UNIT1,.TRUE.)
      CALL RSCRAT(VVAR,NVAR,UNIT1,.TRUE.)
      RAWMU = VMEAN(1)
      RAWSD = SQRT (VVAR(1))

      END

      SUBROUTINE TABOVER (UN, COUNT)
C     ruler                1234567890123456789012345678901234567890
      CHARACTER*40 BLANKS/'                                        '/
      INTEGER UN, COUNT, RCOUNT
      RCOUNT = COUNT
      IF (RCOUNT .LT. 1) RETURN
      IF (RCOUNT .GT. 40) RCOUNT = 40
      WRITE(UN,10) BLANKS(1:RCOUNT)
 10   FORMAT (A,$)
      END

