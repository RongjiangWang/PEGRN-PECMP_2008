      subroutine pecout(neq,nrec,ntr,onlysc,nsc)
      implicit none
c
c     Last modified: Potsdam, April, 2003, by R. Wang
c
      integer neq,nrec,ntr,nsc
      logical onlysc
c
      include 'pecglob.h'
c
      integer i,iadd,nadd,it,ieq,ieq1,ieq2,isc
      integer itr,itr1,itr2,j,inx,irec,lend,m
      double precision dt,t,t1,t2,wl,wr,st1,di1,ra1,st2,di2,ra2
      double precision cmb0,s0xx,s0yy,s0zz,s0xy,s0yz,s0zx,p0
      double precision cmb,cmb1,cmb2,sig,sig1,sig2
      double precision sxx,syy,szz,sxy,syz,szx,p
      double precision stopt,diopt,raopt
      double precision obsout(16),obsadd(8)
      character*3 cmptxt(16)
      character*7 cmptxtadd(8)
      character*7 rtxt(NRECMAX)
      character*160 outfile
c
c     DATA OUTPUT
c     ===========
c
      do lend=80,1,-1
        if(outdir(lend:lend).ne.' ')goto 100
      enddo
100   continue
c
      if(lend.lt.1)then
        stop 'Error in edcmain: wrong for output dirtory!'
      endif
c
      cmptxt(1)=' Ux'
      cmptxt(2)=' Uy'
      cmptxt(3)=' Uz'
      cmptxt(4)='Sxx'
      cmptxt(5)='Syy'
      cmptxt(6)='Szz'
      cmptxt(7)='Sxy'
      cmptxt(8)='Syz'
      cmptxt(9)='Szx'
      cmptxt(10)=' Tx'
      cmptxt(11)=' Ty'
      cmptxt(12)='Rot'
      cmptxt(13)=' Pp'
      cmptxt(14)='Dcx'
      cmptxt(15)='Dcy'
      cmptxt(16)='Dcz'
      iadd=0
      if(insar.eq.1)then
        iadd=iadd+1
        cmptxtadd(iadd)='LOS_Dsp'
      endif
      if(icmb.eq.1)then
        iadd=iadd+1
        cmptxtadd(iadd)='CMB_Fix'
        iadd=iadd+1
        cmptxtadd(iadd)='Sig_Fix'
        iadd=iadd+1
        cmptxtadd(iadd)='CMB_Opt'
        iadd=iadd+1
        cmptxtadd(iadd)='Sig_Opt'
        iadd=iadd+1
        cmptxtadd(iadd)='Str_Opt'
        iadd=iadd+1
        cmptxtadd(iadd)='Dip_Opt'
        iadd=iadd+1
        cmptxtadd(iadd)='Slp_Opt'
c
c       p0 = -skempton*(sigma0(1)+sigma0(2)+sigma0(3))/3.d0
c       for undrained conditions
c       p0 = 0 for drained conditions (during preseismic period)
c
        p0=0.d0
        call prestress(sigma0(1),sigma0(2),sigma0(3),strike0,dip0,rake0,
     &                 p0,friction,cmb0,s0xx,s0yy,s0zz,s0xy,s0yz,s0zx)
      endif
      nadd=iadd
      inx=int(alog10(0.1+float(nrec)))+1
      do irec=1,nrec
        m=irec
        do j=1,inx
          i=m/10**(inx-j)
          rtxt(irec)(j:j)=char(ichar('0')+i)
          m=m-i*10**(inx-j)
        enddo
      enddo
c
      if(nt.gt.1)then
        dt=twindow/dble(nt-1)
      else
        dt=twindow
      endif
c
      do i=1,16
        if(itout(i).eq.1)then
          outfile=outdir(1:lend)//toutfile(i)
          open(30,file=outfile,status='unknown')
          write(30,'(a14,$)')'    Time_day'
          do irec=1,nrec-1
            write(30,'(a14,$)')cmptxt(i)//rtxt(irec)(1:inx)
          enddo
          write(30,'(a14)')cmptxt(i)//rtxt(nrec)(1:inx)
          do it=1,nt
            t=dble(it-1)*dt
            write(30,'(E14.6,$)')t/DAY2SEC
            do irec=1,nrec-1
              obsout(i)=obs(it,irec,i)
              do ieq=1,neq
                if(t.ge.eqtime(ieq))then
                  obsout(i)=obsout(i)+coobs(ieq,irec,i)
                endif
              enddo
              write(30,'(E14.6,$)')obsout(i)
            enddo
            obsout(i)=obs(it,nrec,i)
            do ieq=1,neq
              if(t.ge.eqtime(ieq))then
                obsout(i)=obsout(i)+coobs(ieq,nrec,i)
              endif
            enddo
            write(30,'(E14.6)')obsout(i)
          enddo
          close(30)
        endif
      enddo
c
c     scenario outputs
c
      do isc=1,nsc
        if(tsc(isc).gt.twindow)then
          nwarn=nwarn+1
          print *,' Warning in pecout: scenario outside time window,'
          print *,' no output for the ',isc,'. scenario!'
          goto 500
        endif
        outfile=outdir(1:lend)//scoutfile(isc)
        open(30,file=outfile,status='unknown')
        write(30,'(a24,$)')'    Lat[deg]    Lon[deg]'
        do i=1,16
          write(30,'(a14,$)')cmptxt(i)
        enddo
        write(30,'(8a14)')(cmptxtadd(i),i=1,nadd)
        if(onlysc)then
          it=min0(1+idint(tsc(isc)/dt),nt)
          itr1=0
          do itr=1,ntr
            if(it.eq.itsc(itr))itr1=itr
          enddo
          if(itr1.eq.0)then
            stop 'Error in pecout: snapshot time(-) not found!'
          endif
c
          it=min0(it+1,nt)
          itr2=0
          do itr=1,ntr
            if(it.eq.itsc(itr))itr2=itr
          enddo
          if(itr2.eq.0)then
            stop 'Error in pecout: snapshot time(+) not found!'
          endif
          t1=dble(itsc(itr1)-1)*dt
          t2=dble(itsc(itr2)-1)*dt
        else
          itr1=min0(1+idint(tsc(isc)/dt),nt)
          itr2=min0(itr1+1,nt)
          t1=dble(itr1-1)*dt
          t2=dble(itr2-1)*dt
        endif
        ieq1=0
        do ieq=1,neq
          if(t1.lt.eqtime(ieq).and.eqtime(ieq).le.tsc(isc))then
            ieq1=ieq
            t1=eqtime(ieq)
          endif
        enddo
        ieq2=0
        do ieq=neq,1,-1
          if(t2.gt.eqtime(ieq).and.eqtime(ieq).ge.tsc(isc))then
            ieq2=ieq
            t2=eqtime(ieq)
          endif
        enddo
        if(ieq1.eq.0)then
          do irec=1,nrec
            do i=1,16
              obs1(irec,i)=obs(itr1,irec,i)
            enddo
          enddo
        else
          do irec=1,nrec
            do i=1,16
              obs1(irec,i)=poobs(ieq1,irec,i)
            enddo
          enddo
        endif
        if(ieq2.eq.0)then
          do irec=1,nrec
            do i=1,16
              obs2(irec,i)=obs(itr2,irec,i)
            enddo
          enddo
        else
          do irec=1,nrec
            do i=1,16
              obs2(irec,i)=poobs(ieq2,irec,i)
            enddo
          enddo
        endif
        if(t2-t1.gt.0.1d-03*dt)then
          wr=(tsc(isc)-t1)/(t2-t1)
        else
          wr=0.d0
        endif
        wl=1.d0-wr
        do irec=1,nrec
          iadd=0
          if(insar.eq.1)then
            iadd=iadd+1
            obsadd(iadd)=xlos*(wl*obs1(irec,1)+wr*obs2(irec,1))
     &                  +ylos*(wl*obs1(irec,2)+wr*obs2(irec,2))
     &                  +zlos*(wl*obs1(irec,3)+wr*obs2(irec,3))
            do ieq=1,neq
              if(tsc(isc).ge.eqtime(ieq))then
                obsadd(iadd)=obsadd(iadd)+xlos*coobs(ieq,irec,1)
     &                 +ylos*coobs(ieq,irec,2)+zlos*coobs(ieq,irec,3)
              endif
            enddo
          endif
          if(icmb.eq.1)then
            sxx=wl*obs1(irec,4)+wr*obs2(irec,4)
            syy=wl*obs1(irec,5)+wr*obs2(irec,5)
            szz=wl*obs1(irec,6)+wr*obs2(irec,6)
            sxy=wl*obs1(irec,7)+wr*obs2(irec,7)
            syz=wl*obs1(irec,8)+wr*obs2(irec,8)
            szx=wl*obs1(irec,9)+wr*obs2(irec,9)
            p=wl*obs1(irec,13)+wr*obs2(irec,13)
            do ieq=1,neq
              if(tsc(isc).ge.eqtime(ieq))then
                sxx=sxx+coobs(ieq,irec,4)
                syy=syy+coobs(ieq,irec,5)
                szz=szz+coobs(ieq,irec,6)
                sxy=sxy+coobs(ieq,irec,7)
                syz=syz+coobs(ieq,irec,8)
                szx=szx+coobs(ieq,irec,9)
                p=p+coobs(ieq,irec,13)
              endif
            enddo
c
            call cmbfix(sxx,syy,szz,sxy,syz,szx,
     &                  p,friction,cmb,sig,strike0,dip0,rake0)
            iadd=iadd+1
            obsadd(iadd)=cmb
            iadd=iadd+1
            obsadd(iadd)=sig
            call cmbopt(s0xx+sxx,s0yy+syy,s0zz+szz,s0xy+sxy,
     &                  s0yz+syz,s0zx+szx,p,friction,1,
     &                  strike0,dip0,rake0,
     &                  cmb,sig,st1,di1,ra1,st2,di2,ra2)
            call cmbfix(s0xx,s0yy,s0zz,s0xy,s0yz,s0zx,
     &                  p0,friction,cmb1,sig1,st1,di1,ra1)
            call cmbfix(s0xx,s0yy,s0zz,s0xy,s0yz,s0zx,
     &                  p0,friction,cmb2,sig2,st2,di2,ra2)
            if(cmb1.le.cmb2)then
              call cmbfix(s0xx+sxx,s0yy+syy,s0zz+szz,s0xy+sxy,
     &            s0yz+syz,s0zx+szx,p,friction,cmb,sig,st1,di1,ra1)
              cmb=cmb-cmb1
              sig=sig-sig1
              stopt=st1
              diopt=di1
              raopt=ra1
            else
              call cmbfix(s0xx+sxx,s0yy+syy,s0zz+szz,s0xy+sxy,
     &            s0yz+syz,s0zx+szx,p,friction,cmb,sig,st2,di2,ra2)
              cmb=cmb-cmb2
              sig=sig-sig2
              stopt=st2
              diopt=di2
              raopt=ra2
            endif
            iadd=iadd+1
	      obsadd(iadd)=cmb
            iadd=iadd+1
            obsadd(iadd)=sig
            iadd=iadd+1
            obsadd(iadd)=stopt
            iadd=iadd+1
            obsadd(iadd)=diopt
            iadd=iadd+1
            obsadd(iadd)=raopt
          endif
          write(30,'(2f12.2,$)')xrec(irec),yrec(irec)
          do i=1,16
            obsout(i)=wl*obs1(irec,i)+wr*obs2(irec,i)
            do ieq=1,neq
              if(tsc(isc).ge.eqtime(ieq))then
                obsout(i)=obsout(i)+coobs(ieq,irec,i)
              endif
            enddo
          enddo
          t=tsc(isc)-t1
          if(isurfcon.eq.1.and.zrec.eq.0.d0.and.
     &       ieq1.gt.0.and.t.gt.0.d0)then
            obsout(16)=obsout(16)+wl*dc0(ieq1,irec)/dsqrt(t)
          endif
          write(30,'(24E14.6)')(obsout(i),i=1,16),(obsadd(i),i=1,iadd)
        enddo
        close(30)
500     continue
      enddo
c
1000  format(2f12.2,16E12.4)
      return
      end
