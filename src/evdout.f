C evdout.f
C purpose: write evddata.out file for evdphase 1
C written by: Charles Peterson, June 16, 2011
C called by: pinput.f
C
      subroutine evdout (var,mvar,npeo,nvar,iped,vtraits,sex,mibdid,
     * ifam,conout,EVDstarted)

      double precision var(mvar)
      integer npeo,nvar,iped,vtraits,sex(npeo),ifam,conout
      integer firstped
      logical EVDstarted

C allocated here

      double precision phi2(npeo,npeo),evec(npeo,npeo),
     *  tpheno(nvar,npeo),eval(npeo),evali(npeo),tmean(npeo),
     *  sex_evd(npeo)
      double precision phi2temp,pheno,vmin,vmax,value,tsum,traitmean
      integer ibdid(npeo),unitno,i,j,idoff,lbase,libdid,ivar,ierr,kvar,
     *  ntraits,id,ifevectors,ifevdmat,nevdcovs,matrixsize,midx,midy,
     *  lastpeo,ifmethod
      logical vsample,verbose
      character*80 efilename

c save needed for matrix output over multiple pedigrees
      save midx,midy,matrixsize,lastpeo
      
c procedure begins here
      call ioption ("Eigenvectors", ifevectors)
      call ioption ("EVDmat", ifevdmat)
      if (ifevdmat.ne.0) then
         call ioption ("FPHIMethod",ifmethod)
      else
         ifmethod = 1
      end if
      nevdcovs = 0
      
c determine ntraits and offset to ibdid
      vmin = 0
      vmax = 0
      if (vtraits.eq.0) then
         idoff = 2
         ntraits = 1
      else
         idoff = vtraits + 1
         ntraits = vtraits
      end if

c get ibdid vector

      do 301 i=1,npeo
         kvar = (i-1)*nvar+idoff
         ibdid(i) = var(kvar)
c        print *,"ibdid(",i,") = ",ibdid(i)
         if (ibdid(i).gt.mibdid) then
            mibdid = ibdid(i)
         end if
 301  continue

c get phi2 matrix from loaded matrix file through C++
c note: non-default model type forces loading phi2 matrix during maximization

c     open (unit=29,file='phi2copy.out')

      do 303 i=1,npeo
         do 303 j=i,npeo
            call ibdid2phi (ibdid(i),ibdid(j),phi2temp)
            if (phi2temp.ne.0) then
c               write (29,309) ibdid(i),ibdid(j),phi2temp
 309           format (1x,i5,1x,i5,3x,G12.6)
            end if
c           print *,"phi2temp(",i,",",j,") is ",phi2temp
            phi2(i,j) = phi2temp
            phi2(j,i) = phi2temp
 303  continue
c     close (29)

c If making Y matrix for FPHI, demean the trait variable
c but not if FPHI method 2

      if (ifevdmat.gt.0) then
         if (ntraits .gt. 1) then
            print *,"EVDOUT only supports 1 trait\n"
            call exit
         endif
         if (ifmethod.ne.2) then
c            print *,"DEMEANING TRAIT!"
            tsum=0
            do 304 j=1,npeo
               tsum = tsum + var( 1 + ((j-1)*nvar))
 304        continue
            traitmean = tsum/npeo
            do 305 j=1,npeo
               kvar = 1 + ((j-1) * nvar)
               var(kvar) = var(kvar) - traitmean
 305        continue
         endif
      endif

c do eigenvalue decomposition for this pedigree

      if (.not.EVDstarted) then
         vsample = verbose ("SAMPLE")
         if (vsample) then
            write (conout,311)
 311        format (" Calculating EVD...")
         endif
      end if
      ierr = 0
      call tred2 (npeo,npeo,phi2,eval,evali,evec)
      call tql2 (npeo,npeo,eval,evali,evec,ierr)
      if (ierr.ne.0) then
         print *,"tql2 error: ",ierr
         call exit
      end if 
c     print *,"   Calculated."

c compute tmean

      do 321 i=1,npeo
          tmean(i) = 0
          do 320 j=1,npeo
             tmean(i) = tmean(i) + evec(j,i)
 320      continue
 321   continue

C compute transformed trait values
C true matrix multiply with eigenmatrix inverse (flipped indexes)
C canonical form: c(i) = c(i) + A(k,i) * B(k)
       
       do 340 ivar=1,nvar
          if (ivar.eq.ntraits+1.or.ivar.eq.ntraits+2) then 
C ibdid or group indicator
             go to 340
          end if
          do 339 i=1,npeo
             tpheno(ivar,i)=0
             do 338 k=1,npeo
                kvar = (k-1)*nvar+ivar
                pheno = var(kvar)
                tpheno(ivar,i) = tpheno(ivar,i) + evec (k,i) * pheno
 338         continue
 339      continue
 340   continue

c compute transformed sex value
       do 350 i=1,npeo
          sex_evd(i) = 0
          do 349 k=1,npeo
C                sex_evd(i) = sex_evd(i) + evec (k,i) * (sex(k)-1)
                sex_evd(i) = sex_evd(i) + evec (k,i) * (2-sex(k))
 349      continue
 350   continue


c write out results (header written in ccsearch)

       if (ifevdmat.lt.3) then

       do 390 i=1,npeo
         kvar = (i-1)*nvar+idoff
         id = var(kvar)
         write (26,399) id,sex(i),sex_evd(i),tmean(i),eval(i),
     *     (tpheno(j,i),j=1,ntraits),(tpheno(j,i),j=idoff+2,nvar)
 390  continue
 399  format (I5,",0,0,",I1,9000(',',D18.12))

      else

C WARNING!  The evdmat code below has problems that have not yet
C been solved.  It has been decided not to support evdmat in version
C 8.0.2.

c Make X matrix directly
c  note that some variables may be unused base variables, so use actual
c  count of covariates

         call ioption ("evdcovs",nevdcovs)
         if (.not.EVDstarted) then
            print *,"making NEW new  X matrix"
            call fmatrix_start (npeo,nevdcovs+1,midx)
            matrixsize = npeo
            lastpeo = 0
         else
            lastpeo = matrixsize
            matrixsize = matrixsize + npeo
            call fmatrix_morerows (midx,matrixsize)
         end if
         value = 1;
         do 410 i=lastpeo+1,lastpeo+npeo
            call fmatrix_insert (midx,i,1,value)
            print *,"nevdcovs=",nevdcovs
            do 409 j=1,nevdcovs
               call fmatrix_insert (midx,i,j+1,tpheno(idoff+1+j,i))
 409        continue
 410     continue

c Make Y matrix directly

         if (.not.EVDstarted) then
            print *,"making new Y matrix"
            call fmatrix_start (npeo,1,midy)
         else
            call fmatrix_morerows (midy,matrixsize)
         end if
         do 420 i=lastpeo+1,lastpeo+npeo
            call fmatrix_insert (midy,i,1,tpheno(1,i))
 420     continue
      end if

C Write or Save Eigenvectors

      if (ifevectors.eq.1) then

c determine output filename this time
         if (ifevdmat.eq.2.or.ifevdmat.eq.4) then
c one file
            efilename = "evectors.mat.csv"
         else
            write (efilename,437) ifam
 437        format ('evectors.family',i5.5,'.mat.csv')
         end if
C        print *,"evectors filename is ",efilename

         open (unit=28, file=efilename)

 447     format (D18.12,1000000(',',D18.12))
         do 450 i=1,npeo
            write (28,447) (evec(i,j),j=1,npeo)
 450     continue
         close (28)
      end if
C
C Make eigenmatrix directly
C
      if (ifevectors.eq.2) then
         call fmatrix_start (npeo,npeo,mid)
         do 510 i=1,npeo
            do 500 j=1,npeo
               value = evec(i,j)
               if (value.gt.vmax) vmax = value
               if (value.lt.vmin) vmin = value
               call fmatrix_insert (mid,i,j,value)
 500        continue
 510     continue
      end if

      EVDstarted = .true.
      ifam = ifam + 1
      return
      end



C subroutine unitclof (below) MUST ONLY BE CALLED BY C++ routine unitactive
C this allows for proper closing of fortran unit under all circumstances

      subroutine unitclof(unit)
      integer unit

      close (unit)
      return
      end

      
