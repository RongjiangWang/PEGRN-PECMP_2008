      subroutine pegspec(isp,nr1,nr2,nt,dt,accuracy)
      implicit none
c
      integer isp,nr1,nr2,nt
      double precision dt,accuracy
c
      include 'pegglob.h'
c
      integer i,ir,istp,l,lf,jf,istate,it,its,itga,itgb
      integer nfg,ntg,nkc
      double precision f,t,dfg,dtg,fc,fcut,rc,kc,dk0,dk
      double precision a,b,alpha,beta,thick,dkmin
      double precision fgrnabs,fgrnmax
      double precision dspabs(4),stsabs(4),stnabs(4),dcfabs(4)
      double precision dspmax(4),stsmax(4),stnmax(4),dcfmax(4)
      double precision tgrn(ntmax),obs(ntmax)
      double complex cs,ca,cb
      double complex fgrn(ntmax)
      logical again,lowpass
c
      integer nkmax
      double precision eps,pi,pi2
      double complex c1,c2
      data nkmax/100000/
      data eps,pi,pi2/0.05d0,3.14159265358979d0,6.28318530717959d0/
      data c1,c2/(1.d0,0.d0),(2.d0,0.d0)/      
c
      print *,' ==============================================='
      write(*,'(i3,a,F12.4,a,f12.4,a)')isp,'. sub-profile: ',
     &        r(nr1)/km2m,' -> ',r(nr2)/km2m,' km'
c
c     determine wavenumber sampling rate
c
      rc=rs(nr1)+0.5d0*(r(nr1)+r(nr2))
      kc=pi2/rc
      dk0=eps*pi2/dsqrt((zrec-zs)**2+rc**2)
      thick=0.d0
      do l=1,lp-1
        thick=thick+hp(l)
      enddo
      dkmin=eps*accuracy*pi2/dsqrt(thick**2+rc**2)
c
c     determine frequency sampling rate
c
      dfg=1.d0/(2.d0*dble(nt)*dt)
      alpha=pi2*dfg
c
c     istate = -1: elastic case
c
      istate=-1
      cs=(0.d0,0.d0)
      dk=dk0
      do istp=1,4
        do i=1,16
          do ir=nr1,nr2
            du(0,ir,i,istp)=(0.d0,0.d0)
          enddo
        enddo
      enddo
100   nkc=idint(2.d0*kc/dk)
      call pegwvint(istate,1,nr1,nr2,cs,nkc,nkmax,dk,accuracy,.false.)
      again=.false.
      do istp=1,4
        dspabs(istp)=0.d0
        stsabs(istp)=0.d0
        stnabs(istp)=0.d0
        dcfabs(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspabs(istp)=dspabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          do i=4,9
            stsabs(istp)=stsabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          stsabs(istp)=stsabs(istp)
     &                +geow(ir)*cdabs(du(1,ir,13,istp)-du(0,ir,13,istp))
          do i=10,12
            stnabs(istp)=stnabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          do i=14,16
            dcfabs(istp)=dcfabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
        enddo
        dspmax(istp)=0.d0
        stsmax(istp)=0.d0
        stnmax(istp)=0.d0
        dcfmax(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspmax(istp)=dspmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          do i=4,9
            stsmax(istp)=stsmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          stsmax(istp)=stsmax(istp)+geow(ir)*cdabs(du(1,ir,13,istp))
          do i=10,12
            stnmax(istp)=stnmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          do i=14,16
            dcfmax(istp)=dcfmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
        enddo
        again=again
     &        .or.dspabs(istp).gt.0.1d0*accuracy*dspmax(istp)
     &        .or.stsabs(istp).gt.0.1d0*accuracy*stsmax(istp)
     &        .or.stnabs(istp).gt.0.1d0*accuracy*stnmax(istp)
     &        .or.dcfabs(istp).gt.0.1d0*accuracy*dcfmax(istp)
      enddo
      if(again.and.dk.gt.dkmin)then
        do istp=1,4
          do i=1,16
            do ir=nr1,nr2
              du(0,ir,i,istp)=du(1,ir,i,istp)
            enddo
          enddo
        enddo
        dk=0.5d0*dk
        goto 100
      endif
      if(again)then
        print *,' Warning in pegspec: required acc. not satisfied!'
      endif
c
      istate=-1
      cs=(0.d0,0.d0)
      call pegwvint(istate,-1,nr1,nr2,cs,4*nkc,nkmax,
     &              0.25d0*dk,0.25*accuracy,.true.)
      if(nt.eq.1)then
        do istp=1,4
          do i=1,16
            do ir=nr1,nr2
              du(1,ir,i,istp)=dcmplx(dreal(du(-1,ir,i,istp)),0.d0)
            enddo
          enddo
        enddo
        return
      endif
c
      istate=1
      cs=dcmplx(alpha,dfg)
      dk=dk0
      do istp=1,4
        do i=1,16
          do ir=nr1,nr2
            du(0,ir,i,istp)=(0.d0,0.d0)
          enddo
        enddo
      enddo
200   nkc=idint(2.d0*kc/dk)
      call pegwvint(istate,1,nr1,nr2,cs,nkc,nkmax,dk,accuracy,.false.)
      again=.false.
      do istp=1,4
        dspabs(istp)=0.d0
        stsabs(istp)=0.d0
        stnabs(istp)=0.d0
        dcfabs(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspabs(istp)=dspabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          do i=4,9
            stsabs(istp)=stsabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          stsabs(istp)=stsabs(istp)
     &                +geow(ir)*cdabs(du(1,ir,13,istp)-du(0,ir,13,istp))
          do i=10,12
            stnabs(istp)=stnabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
          do i=14,16
            dcfabs(istp)=dcfabs(istp)
     &                  +geow(ir)*cdabs(du(1,ir,i,istp)-du(0,ir,i,istp))
          enddo
        enddo
        dspmax(istp)=0.d0
        stsmax(istp)=0.d0
        stnmax(istp)=0.d0
        dcfmax(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspmax(istp)=dspmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          do i=4,9
            stsmax(istp)=stsmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          stsmax(istp)=stsmax(istp)+geow(ir)*cdabs(du(1,ir,13,istp))
          do i=10,12
            stnmax(istp)=stnmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
          do i=14,16
            dcfmax(istp)=dcfmax(istp)+geow(ir)*cdabs(du(1,ir,i,istp))
          enddo
        enddo
        again=again
     &        .or.dspabs(istp).gt.accuracy*dspmax(istp)
     &        .or.stsabs(istp).gt.accuracy*stsmax(istp)
     &        .or.stnabs(istp).gt.accuracy*stnmax(istp)
     &        .or.dcfabs(istp).gt.accuracy*dcfmax(istp)
      enddo
      if(again.and.dk.gt.dkmin)then
        do istp=1,4
          do i=1,16
            do ir=nr1,nr2
              du(0,ir,i,istp)=du(1,ir,i,istp)
            enddo
          enddo
        enddo
        dk=0.5d0*dk
        goto 200
      endif
      if(again)then
        print *,' Warning in pegspec: required acc. not satisfied!'
      endif
c
      if(nt.eq.2)then
        istate=0
        cs=(0.d0,0.d0)
        call pegwvint(istate,0,nr1,nr2,cs,nkc,nkmax,dk,
     &                accuracy,.true.)
        do istp=1,4
          do i=1,16
            do ir=nr1,nr2
              du(1,ir,i,istp)=dcmplx(dreal(du(-1,ir,i,istp)),
     &                               dreal(du( 0,ir,i,istp)))
            enddo
          enddo
        enddo
        return
      endif
c
      print *,' Finding upper cut-off frequency'
      do istp=1,4
        dspmax(istp)=0.d0
        stsmax(istp)=0.d0
        stnmax(istp)=0.d0
        dcfmax(istp)=0.d0
      enddo
      istate=1
      lf=1
      nfg=1
500   fcut=dble(nfg-1)*dfg
      cs=dcmplx(alpha,pi2*fcut)
      call pegwvint(istate,lf,nr1,nr2,cs,nkc,nkmax,dk,accuracy,.true.)
      again=.false.
      do istp=1,4
        dspabs(istp)=0.d0
        stsabs(istp)=0.d0
        stnabs(istp)=0.d0
        dcfabs(istp)=0.d0
        do ir=nr1,nr2
          do i=1,3
            dspabs(istp)=dspabs(istp)+geow(ir)*cdabs(du(lf,ir,i,istp))
          enddo
          do i=4,9
            stsabs(istp)=stsabs(istp)+geow(ir)*cdabs(du(lf,ir,i,istp))
          enddo
          stsabs(istp)=stsabs(istp)+geow(ir)*cdabs(du(lf,ir,13,istp))
          do i=10,12
            stnabs(istp)=stnabs(istp)+geow(ir)*cdabs(du(lf,ir,i,istp))
          enddo
          do i=14,16
            dcfabs(istp)=dcfabs(istp)+geow(ir)*cdabs(du(lf,ir,i,istp))
          enddo
        enddo
        dspmax(istp)=dmax1(dspmax(istp),dspabs(istp))
        stsmax(istp)=dmax1(stsmax(istp),stsabs(istp))
        stnmax(istp)=dmax1(stnmax(istp),stnabs(istp))
        dcfmax(istp)=dmax1(dcfmax(istp),dcfabs(istp))
        again=again
     &        .or.dspabs(istp).gt.eps*dspmax(istp)
     &        .or.stsabs(istp).gt.eps*stsmax(istp)
     &        .or.stnabs(istp).gt.eps*stnmax(istp)
     &        .or.dcfabs(istp).gt.eps*dcfmax(istp)
c        write(*,'(4E12.4)')dspabs(istp),stsabs(istp),
c     &                     stnabs(istp),dcfabs(istp)
      enddo
c      pause
      if(again.and.nfg.le.nfmax.or.nfg.lt.max0(nfmin,nt/2))then
        lf=lf+1
        nfg=2*nfg
        goto 500
      endif
c
c     cut-off frequency found
c
      if(nfg.gt.nfmax)then
        write(*,'(a,i5,a)')'  Warning in pespectra: nfmax (=',
     &                   nfmax,') may be insufficient!'
        nfg=nfmax
        lowpass=.true.
      else if(nfg.lt.nfmin)then
        nfg=nfmin
        lowpass=.false.
      else
        lowpass=.false.
      endif
      fcut=dble(nfg-1)*dfg
      write(*,'(i6,a,E12.5,a)')nfg,' frequencies will be used,'
     &                             //' cut-off: ',fcut,'Hz'
      print *,' ==============================================='
c
      istate=1
      do lf=1,nfg
        f=dble(lf-1)*dfg
        cs=dcmplx(alpha,pi2*f)
        call pegwvint(istate,lf,nr1,nr2,cs,nkc,nkmax,dk,accuracy,.true.)
      enddo
c
c     convert to time domain by FFT
c
      dtg=1.d0/(dble(2*nfg)*dfg)
      ntg=2*nfg
c
      do istp=1,4
        do i=1,16
          if(select(i,istp))then
            do ir=nr1,nr2
              do lf=1,nfg
                fgrn(lf)=du(lf,ir,i,istp)
              enddo
              if(lowpass)then
c
c               low-pass filter has to be used
c
c               find a proper corner frequency as high as possible
c
                fc=fcut/eps
600             fgrnmax=0.d0
                beta=1.d0+alpha/(pi2*fc)
                do lf=1,nfg
                  fgrnabs=cdabs(fgrn(lf))/(beta**2
     &                   +(dble(lf-1)*dfg/fc)**2)**1.5d0
                  fgrnmax=dmax1(fgrnmax,fgrnabs)
                enddo
                if(fgrnabs.gt.eps*fgrnmax)then
                  fc=0.75d0*fc
                  goto 600
                endif
c
c               use low-pass filter with the corner frequency found
c
                do lf=1,nfg
                  fgrn(lf)=fgrn(lf)/(dcmplx(beta,dble(lf-1)*dfg/fc))**3
                enddo
              endif
	        jf=1
	        do lf=2*nfg,nfg+2,-1
	          jf=jf+1
	          fgrn(lf)=dconjg(fgrn(jf))
	        enddo
	        fgrn(nfg+1)=(0.d0,0.d0)
c
c	        convention for Fourier transform:
c	        f(t)=\int F(f) exp(i2\pi f t) df
c
	        call four1(fgrn,2*nfg,+1)
c
c             calculate response to the Heaviside source history
c
              tgrn(1)=0.d0
              do it=2,ntg
                tgrn(it)=tgrn(it-1)+dreal(fgrn(it-1))
     &                  *dtg*dfg*dexp(alpha*dble(it-1)*dtg)
              enddo
c
c             linear interpolation
c
              obs(1)=dreal(du(-1,ir,i,istp))
              do it=2,2*((1+nt)/2)
                t=dble(it-1)*dt
                b=dmod(t/dtg,1.d0)
                a=1.d0-b
                itga=min0(1+idint(t/dtg),ntg)
                itgb=min0(2+idint(t/dtg),ntg)
                obs(it)=tgrn(itga)*a+tgrn(itgb)*b
              enddo
c
              do it=1,(1+nt)/2
                du(it,ir,i,istp)=dcmplx(obs(2*it-1),obs(2*it))
              enddo
            enddo
          endif
        enddo
	enddo
      return
      end