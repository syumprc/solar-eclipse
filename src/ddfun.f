      subroutine ddfun(df,obs,par,f,iter,kase,ncase,nobs,npar,
     &   npass,pass,differ,vardata,male,trttype,vtraits,nvar,ntot,
     &   nind,nascer,nped,ncumind,cov,mu,affect,maxpeo)

c  returns the likelihood function for a set of pedigrees

      implicit double precision(a-h,o-z)

      double precision df(npar),obs(nobs),par(npar),vardata(nvar,ntot)
      double precision cov(maxpeo,maxpeo),mu(maxpeo),affect(maxpeo)
      integer pass,nind(nped),nascer(nped),nped,ncumind(0:nped)
      integer maxpeo,trtnoi,trtnoj,vtraits,tindex
      logical differ(2),male(ntot)
      character*1 trttype(ntot)

      double precision lnl(nped),lnl_ascer(nped)
      double precision lnjac(nped),lnjac_ascer(nped)
      double precision dmean,dcovar,temp1,method,pedlike,ll
      double precision sqr2pi,zmin,udiag,phi2i
      logical evd,udiagset
      integer momega,omegatyp,ifirstper,maxj,ifdisc,ierr,evdphase
      character*1 disc(maxpeo)
      character*32 modeltype

      zmin = -1.d37
      sqr2pi = 0.3989422804014
      udiagset = .false.
      ierr = 0

c  check for pedlike option

C      print *, "ENTERING DDFUN"
      call doption ('PEDLIKE',pedlike)
      if (pedlike.eq.1) then
         open (34,file='pedlike.dat')
         write (34,9)
 9       format ('LOGLIKELIHOOD')
      end if

c  use Jeff's discrete code if requested

      call doption ("DiscreteMethod",method)
      if (method.eq.2.0) then
         call fun_mehd (par,f,npar,vardata,male,nvar,ntot,nind,
     &                  nascer,nped,ncumind,maxpeo)
         return
      endif

C  check ModelType option for EVD status

      evd = .false.
      evdphase = 0
      call soption ("ModelType",modeltype)

C  Only applicable to polygenic model type
C  But this feature not yet implemented so fake it
      momega = omegatyp()
      
c      print *, "Getting EVDPhase"
      call ioption ("EVDPhase",evdphase)
      call ifanydisc (ifdisc)
      if (modeltype(1:3).eq."Evd".and.(evdphase.eq.1.or.
     & (momega.ge.1.and.momega.le.2))
     & .and.ifdisc.eq.0) then
         evd = .true.
      endif

c      print *, "Got EVDPhase ",evdphase

      trtnoi = 1
      trtnoj = 1
      tindex = 1
      if (vtraits.ne.0) then
         tindex = 2
      endif
c
c  likelihood function
c
c**   zero vectors
c
      do i = 1, nped
         lnl(i) = 0.0
         lnl_ascer(i) = 0.0
         lnjac(i) = 0.0
         lnjac_ascer(i) = 0.0
      enddo
      biglik = 0.0
      biglik_ascer = 0.0
c
c**   calculate likelihood for each pedigree
c
      do 20 i = 1, nped
         itmp = 1
         ifirstper = ncumind(i-1)+1

         do 21 iind = ncumind(i-1)+1, ncumind(i)
            disc(itmp) = trttype(iind)
            affect(itmp) = vardata(tindex,iind)
            trtnoi = 1
            if (vtraits.ne.0) then
               trtnoi = vardata(1,iind)
            endif
            if (.not.evd) then
               temp1 = dcovar(iind,iind,trtnoi,trtnoi,par,npar,
     &                        vardata,nvar,ntot,male(iind),male(iind))
               temp1 = 1.0/dsqrt(temp1)
            endif
            mu(itmp) = dmean(iind,trtnoi,par,npar,vardata,nvar,ntot,
     &                       male(iind))

c           print *, "Got past mu"
            if (disc(itmp).eq.'c') then
               if (evd) then
                  mu(itmp) = affect(itmp) - mu(itmp)
               else
                  mu(itmp) = (affect(itmp) - mu(itmp))*temp1
                  lnjac(i) = lnjac(i) + dlog(temp1) - dlog(sqr2pi)
               endif
            elseif (disc(itmp).eq.'d') then
               if (.not.evd) then
                  mu(itmp) = mu(itmp)*temp1
               endif
            else
               print *, "invalid trait type!"
               stop
            endif
            jtmp = itmp

            maxj = ncumind(i)
            phi2i = 0
            if (evd) then
               maxj = iind
               call getphi2i (iind,vardata,nvar,phi2i)
c               print *,"Got phi2"
            end if

            do 22 jind = iind, maxj
               trtnoj = 1
               if (vtraits.ne.0) then
                  trtnoj = vardata(1,jind)
               endif
               
               if ((.not.udiagset).or.(phi2i.ne.1).or.(.not.evd)) then
C                 print *,"* computing omega for ",jind,iind
                  cov(itmp,jtmp) = dcovar(iind,jind,trtnoi,trtnoj,par,
     &                 npar,vardata,nvar,ntot,male(iind),male(jind))
                  if (phi2i.eq.1) then
                     udiag = cov(itmp,jtmp)
                     udiagset = .true.
                  end if
               else
C                  print *,"* NOT computing omega for ",jind,iind
                  cov(itmp,jtmp) = udiag
               endif
               if (.not.evd) then
                  temp2 = dcovar(jind,jind,trtnoj,trtnoj,par,npar,
     &                           vardata,nvar,ntot,male(jind),
     &                           male(jind))
                  temp2 = 1.0/dsqrt(temp2)
                  cov(itmp,jtmp) = cov(itmp,jtmp)*temp1*temp2
               endif
               cov(jtmp,itmp) = cov(itmp,jtmp)
               jtmp = jtmp + 1
22          continue
            itmp = itmp + 1
C           print *, "got to 21"
21       continue

C         print *, "Got to tred2"
         if (evd) then
C            call tred2(maxpeo,nind(i),tphi2,eval,evali,evec)
C            call tql2(maxpeo,nind(i),eval,evali,evec,ierr)
C            call evdlik(nind(i),maxpeo,mu,cov,lnl(i),par(4),eval,evec)

            call evdlikc (i,maxpeo,nind(i),mu,cov,lnl(i),par(4),
     &     vardata(1,ifirstper),nvar,evdphase,
     &     male,vtraits,ntot,nind,nascer,nped,ncumind,ierr)
            biglik = biglik + lnl(i)
C          print *,"Current likelihood is ",lnl(i)
         else
            call mvncdf(nind(i),maxpeo,affect,disc,mu,cov,lnl(i))
            biglik = biglik + lnl(i) + lnjac(i)
         endif
20    continue

c
c**   calculate ascertainment correction for each pedigree
c
      do 40 i = 1, nped
         if (nascer(i).eq.0) goto 40
         itmp = 1
         ifirstper = ncumind(i)+1-nascer(i)
         do 41 iind = ncumind(i)+1-nascer(i), ncumind(i)
            disc(itmp) = trttype(iind)
            affect(itmp) = vardata(tindex,iind)
            trtnoi = 1
            if (vtraits.ne.0) then
               trtnoi = vardata(1,iind)
            endif
            if (.not.evd) then
               temp1 = dcovar(iind,iind,trtnoi,trtnoi,par,npar,
     &                        vardata,nvar,ntot,male(iind),male(iind))
               temp1 = 1.0/dsqrt(temp1)
            endif
            mu(itmp) = dmean(iind,trtnoi,par,npar,vardata,nvar,ntot,
     &                       male(iind))
            if (disc(itmp).eq.'c') then
               if (evd) then
                  mu(itmp) = affect(itmp) - mu(itmp)
               else
                  mu(itmp) = (affect(itmp) - mu(itmp))*temp1
                  lnjac_ascer(i) =
     &               lnjac_ascer(i) + dlog(temp1) - dlog(sqr2pi)
               endif
            elseif (disc(itmp).eq.'d') then
               if (.not.evd) then
                  mu(itmp) = mu(itmp)*temp1
               endif
            else
               print *, "invalid trait type!"
               stop
            endif
            jtmp = itmp
            maxj = ncumind(i)
            if (evd) maxj = iind

            do 42 jind = iind, maxj
               trtnoj = 1
               if (vtraits.ne.0) then
                  trtnoj = vardata(1,jind)
               endif
               cov(itmp,jtmp) = dcovar(iind,jind,trtnoi,trtnoj,
     &                                 par,npar,vardata,nvar,ntot,
     &                                 male(iind),male(jind))
               if (.not.evd) then
                  temp2 = dcovar(jind,jind,trtnoj,trtnoj,par,npar,
     &                           vardata,nvar,ntot,male(jind),
     &                           male(jind))
                  temp2 = 1.0/dsqrt(temp2)
                  cov(itmp,jtmp) = cov(itmp,jtmp)*temp1*temp2
               endif
               cov(jtmp,itmp) = cov(itmp,jtmp)
               jtmp = jtmp + 1
42          continue
            itmp = itmp + 1
41       continue
        if (evd) then
C            call tred2(maxpeo,nind(i),tphi2,eval,evali,evec)
C            call tql2(maxpeo,nind(i),eval,evali,evec,ierr)
C            call evdlik(nascer(i),maxpeo,mu,cov,lnl_ascer(i),par(4),
C     &                  eval,evec)

C
C Proband matrices are stacked after non-proband ones
C
           call evdlikc (i+nped,maxpeo,nascer(i),mu,cov,
     &      lnl_ascer(i),par(4),vardata(1,ifirstper),nvar,evdphase,
     &      male,vtraits,ntot,nind,nascer,nped,ncumind,ierr)
            biglik_ascer = biglik_ascer + lnl_ascer(i)
         else
            call mvncdf(nascer(i),maxpeo,affect,disc,mu,cov,
     &                  lnl_ascer(i))
            biglik_ascer = biglik_ascer + lnl_ascer(i) + lnjac_ascer(i)
         endif
40    continue

666   f = -(biglik - biglik_ascer)
      if (evdphase.eq.1) then
C       evdtrap traps out of maximize and does not return here
         call evdtrap (0)
      end if

c**   write results to pedlike.dat if required

      if (pedlike.eq.1) then
         do 90 i = 1, nped
            ll = lnl(i) + lnjac(i) - (lnl_ascer(i) + lnjac_ascer(i))
            write (34,45) ll
45          format (1X,G26.18)
90       continue
         close (34)
      end if

      return
      end
