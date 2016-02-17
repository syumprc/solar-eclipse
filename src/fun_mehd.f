      SUBROUTINE FUN_MEHD(par,f,npar,vardata,male,nvar,ntot,nind,
     2 nascer,nped,ncumped,maxind_ped)
C                                                                             
C     Returns the likelihood function for a set of pedigrees
C                                                                             
      IMPLICIT NONE

      INTEGER npar,nvar,ntot,maxind_ped,nped
      DOUBLE PRECISION vardata(nvar,ntot),f
      LOGICAL male(maxind_ped)
      INTEGER nind(nped),nascer(nped),ncumped(0:nped)

      DOUBLE PRECISION ainf,biglik,biglik_ascer
      INTEGER i,itmp,iind,jtmp,jind

      DOUBLE PRECISION par(npar)                  
C     INCLUDE 'uvc.h'
      double precision cov(maxind_ped,maxind_ped),mu(maxind_ped)
      double precision lnl(nped),lnl_ascer(nped)
      double precision affect(maxind_ped)
      double precision MUCC,DCOVAR
      double precision a(maxind_ped),b(maxind_ped)

      data ainf /10.0d0/
      double precision w_mehd

      biglik=0.0
      biglik_ascer=0.0
C
C     CALCULATE LIKELIHOOD FOR EACH PEDIGREE  
C
      do 20 i=1,nped
        itmp=1
        do 21 iind=ncumped(i-1)+1,ncumped(i)
          affect(itmp)=vardata(1,iind)
          mu(itmp)=MUCC(male(iind),vardata(1,iind),par,1)

c     limits of integration
          if (affect(itmp) .eq. 1.0) then
c            a(itmp) = -mu(itmp)+mu(itmp)
            a(itmp) = 0.0
            b(itmp) = ainf+mu(itmp)
          else
            a(itmp) = -ainf+mu(itmp)
c            b(itmp) = -mu(itmp)+mu(itmp)
            b(itmp) = 0.0
          endif

          jtmp=itmp
          do 22 jind=iind,ncumped(i)
            cov(itmp,jtmp)=DCOVAR(iind,jind,1,1,par,npar,vardata,
     1                            nvar,ntot,male(iind),male(jind))
            cov(jtmp,itmp)=cov(itmp,jtmp)
            jtmp=jtmp+1
   22     continue
          itmp=itmp+1
   21   continue

        lnl(i)=w_mehd(maxind_ped,nind(i),mu,cov,a,b)
        biglik=biglik+lnl(i)

c        print *,"lnl(",i,")=",lnl(i)," biglik=",biglik
c        pause
   20 continue
c      print *, biglik
c      pause

C
C     CALCULATE ASCERTAINMENT CORRECTION FOR EACH PEDIGREE  
C
      do 40 I=1,nped
        if (nascer(i).eq.0) goto 40
        itmp=1
        do 41 iind=ncumped(i)+1-nascer(i),ncumped(i)
          affect(itmp)=vardata(1,iind)
          mu(itmp)=MUCC(male(iind),vardata(1,iind),par,1)
          jtmp=itmp
          do 42 jind=iind,ncumped(i)
            cov(itmp,jtmp)=DCOVAR(iind,jind,1,1,par,npar,vardata,
     1                            nvar,ntot,male(iind),male(jind))
            cov(jtmp,itmp)=cov(itmp,jtmp)
            jtmp=jtmp+1
   42     continue
          itmp=itmp+1
   41   continue
        lnl_ascer(i)=w_mehd(maxind_ped,nind(i),mu,cov,a,b)
        biglik_ascer=biglik_ascer+lnl_ascer(i)
   40 continue

  666 f = -(biglik-biglik_ascer)
      RETURN
      END
