      subroutine pegwvint(istate,lf,nr1,nr2,cs,
     &                    nkc,nkmax,dk,accuracy,tty)
      implicit none
c
      integer istate,lf,nr1,nr2,nkc,nkmax
      double precision dk,accuracy
      double complex cs
      logical tty
c
      include 'pegglob.h'
c
c     u: 1=uz, 2=ur, 3=ut,
c        4=szz, 5=srr, 6=stt, 7=szr, 8=srt, 9=stz, 10=-dur/dz
c        11=-dut/dz, 12=rot(u)_z/2, 13=-xi*dp/dz, 14=-xi*dp/dr, 15=-xi*dp/dt/r
c     NOTE: uz, ur, szz, szr, srr, stt, dur/dz, dp/dz, dp/dr
c           have the same azimuth-factor as the poloidal mode (p-sv);
c           ut, srt, stz, dut/dz, rot(u)_z/2,dp/dt/r
c           have the same azimuth-factor as the
c           toroidal mode (sh);
c
      integer i,istp,ir,ik,nx,nnk
      double precision k,kgr,krt,pi,pi2,x,wl,wr,fac,y0abs,dyabs
      double precision bsj0(-1:3)
      double complex c0,c1,c2,c3,c4,c05,ck
      double complex clazr,cmuzr,xi,pe,evs
      double complex cics(4),cbsj(3),cms(4)
      double complex y0(8,4),cy(14,4)
      double complex urdr(4),utdr(4),urdr0(4)
      double complex obs(nrmax,16,4),obs0(nrmax,16,4)
      logical finish,tests,zerokern
c
      integer iret,ipot
      double precision xokada,yokada,dip
      double complex alpha,pot1,pot2,pot3,pot4,cs45,ss45,potfac
      double complex ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
      double complex uza,ura,uta,esfa,urza,utza,etta,erta
      double complex szza,srra,stta,szra,srta,szta
c
      double precision eps
      data eps/1.0d-03/
c
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
c
      c0=(0.d0,0.d0)
      c1=(1.d0,0.d0)
      c2=(2.d0,0.d0)
      c3=(3.d0,0.d0)
      c4=(4.d0,0.d0)
      c05=(0.5d0,0.d0)
      cs45=dcmplx(1.d0/dsqrt(2.d0),0.d0)
      ss45=dcmplx(1.d0/dsqrt(2.d0),0.d0)
c
c     ics = 1  when the azmuth-factor is cos(ms*theta) for poloidal mode
c             (psv) and sin(ms*theta) for the toroidal mode (sh);
c     ics = -1 otherwise.
c
      do istp=1,4
        cics(istp)=dcmplx(dble(ics(istp)),0.d0)
c
        cms(istp)=dcmplx(dble(ms(istp)),0.d0)
      enddo
c
      do istp=1,4
        urdr(istp)=(0.d0,0.d0)
        utdr(istp)=(0.d0,0.d0)
        urdr0(istp)=(0.d0,0.d0)
        do i=1,16
          do ir=nr1,nr2
            obs(ir,i,istp)=c0
            obs0(ir,i,istp)=c0
          enddo
        enddo
      enddo
c
      call pegmoduli(cs,istate)
      call pegsfct(1.d0)
      clazr=cla(nlr)
      cmuzr=cmu(nlr)
c
      nnk=1
c
      do ik=1,nkmax
        k=dble(ik)*dk
        ck=dcmplx(k,0.d0)
        call pegkern(y0,cs,k,istate,zerokern)
        if(zerokern)goto 100
        do istp=1,4
c
c         for displacement components
c
          cy(1,istp)=y0(1,istp)
          cy(2,istp)=c05*(y0(3,istp)+cics(istp)*y0(5,istp))
          cy(3,istp)=c05*(y0(3,istp)-cics(istp)*y0(5,istp))
c
c         for strain components
c
          cy(4,istp)=y0(2,istp)
          cy(5,istp)=ck*y0(3,istp)
          cy(6,istp)=c05*(y0(4,istp)+cics(istp)*y0(6,istp))
          cy(7,istp)=c05*(y0(4,istp)-cics(istp)*y0(6,istp))
          cy(8,istp)=ck*y0(5,istp)
          cy(9,istp)=ck*y0(1,istp)
c
c         for tilt components
c
          cy(10,istp)=ck*y0(1,istp)
c
c         for pore pressure
c
          cy(11,istp)=y0(7,istp)
c
c         for Darcy flux
c
          cy(12,istp)=y0(8,istp)
          cy(13,istp)=ck*y0(7,istp)
          cy(14,istp)=ck*y0(7,istp)
        enddo
c
c       obs1-3 are displacement components:
c       obs4 = normal stress: szz
c       obs5 = surface strain: err+ett
c       obs6 = ett for r > 0 will be derived later
c       obs7 = shear stress: szr
c       obs8 = 2*rot(u)_z = dut/dr - (dur/dt)/r + ut/r
c       obs9 = shear stress: szt
c       obs10 = tilt-r: duz/dr
c       obs11 = duz/dt/r (reserved for tilt-t)
c       obs12 = pore pressure p
c       obs13 = -xi*dp/dz
c       obs14 = -xi*dp/dr
c       obs15 = -xi*dp/dt/r
c       obs16 = obs8
c
        if(r(nr1).eq.0.d0)then
c
c          compute ur/r and ut/r for r -> 0
c
          fac=k*k*dk*dexp(-0.5d0*(k*rs(nr1))**2)
          do istp=1,4
            do i=1,3,2
              if(ms(istp)+i-2.eq.1)then
                cbsj(i)=dcmplx(0.5d0*fac,0.d0)
              else if(ms(istp)+i-2.eq.-1)then
                cbsj(i)=-dcmplx(0.5d0*fac,0.d0)
              else
                cbsj(i)=(0.d0,0.d0)
              endif
            enddo
            urdr(istp)=urdr(istp)
     &          +cy(2,istp)*cbsj(1)-cy(3,istp)*cbsj(3)
            utdr(istp)=utdr(istp)
     &          -cics(istp)*(cy(2,istp)*cbsj(1)+cy(3,istp)*cbsj(3))
          enddo
        endif  
        do ir=nr1,nr2
          fac=k*dk*dexp(-0.5d0*(k*rs(ir))**2)
          x=k*r(ir)
          nx=1+idint(x/dxbsj)
          if(nx.le.nbsjmax)then
            wl=(dble(nx)*dxbsj-x)/dxbsj
            wr=1.d0-wl
            do i=-1,3
              bsj0(i)=(wl*bsj(nx-1,i)+wr*bsj(nx,i))*fac
            enddo
          else
            do i=-1,3
              bsj0(i)=dcos(x-0.5d0*pi*(dble(i)+0.5d0))*fac
     &               /dsqrt(0.5d0*pi*x)
            enddo
          endif
          do istp=1,4
            do i=1,3
              cbsj(i)=dcmplx(bsj0(ms(istp)+i-2),0.d0)
            enddo
            obs(ir,1,istp)=obs(ir,1,istp)+cy(1,istp)*cbsj(2)
            obs(ir,2,istp)=obs(ir,2,istp)
     &         +cy(2,istp)*cbsj(1)-cy(3,istp)*cbsj(3)
            obs(ir,3,istp)=obs(ir,3,istp)
     &         -cics(istp)*(cy(2,istp)*cbsj(1)+cy(3,istp)*cbsj(3))
            obs(ir,4,istp)=obs(ir,4,istp)+cy(4,istp)*cbsj(2)
            obs(ir,5,istp)=obs(ir,5,istp)-cy(5,istp)*cbsj(2)
            obs(ir,7,istp)=obs(ir,7,istp)
     &          +cy(6,istp)*cbsj(1)-cy(7,istp)*cbsj(3)
            obs(ir,8,istp)=obs(ir,8,istp)+cy(8,istp)*cbsj(2)
            obs(ir,9,istp)=obs(ir,9,istp)
     &          -cics(istp)*(cy(6,istp)*cbsj(1)+cy(7,istp)*cbsj(3))
            obs(ir,10,istp)=obs(ir,10,istp)
     &          +c05*(cy(9,istp)*cbsj(1)-cy(10,istp)*cbsj(3))
            obs(ir,11,istp)=obs(ir,11,istp)-c05*cics(istp)
     &          *(cy(9,istp)*cbsj(1)+cy(10,istp)*cbsj(3))
            obs(ir,12,istp)=obs(ir,12,istp)+cy(11,istp)*cbsj(2)
            obs(ir,13,istp)=obs(ir,13,istp)+cy(12,istp)*cbsj(2)
            obs(ir,14,istp)=obs(ir,14,istp)
     &          +c05*(cy(13,istp)*cbsj(1)-cy(14,istp)*cbsj(3))
            obs(ir,15,istp)=obs(ir,15,istp)-c05*cics(istp)
     &          *(cy(13,istp)*cbsj(1)+cy(14,istp)*cbsj(3))        
          enddo
        enddo
        if(ik.eq.nnk*nkc)then
          nnk=2*nnk
          finish=.true.
          if(r(nr1).eq.0.d0)then
            do istp=1,4
              finish=finish.and.cdabs(urdr(istp)-urdr0(istp))
     &               .le.accuracy*cdabs(urdr(istp))
              urdr0(istp)=urdr(istp)
            enddo
          endif
          do istp=1,4
            do i=1,15
              y0abs=0.d0
              dyabs=0.d0
              do ir=nr1,nr2
                y0abs=y0abs+geow(ir)*cdabs(obs(ir,i,istp))
                dyabs=dyabs
     &               +geow(ir)*cdabs(obs(ir,i,istp)-obs0(ir,i,istp))
                obs0(ir,i,istp)=obs(ir,i,istp)
              enddo
              finish=finish.and.dyabs.le.accuracy*y0abs
            enddo
          enddo
          if(finish)goto 100
        endif
      enddo
      ik=ik-1
100   continue
c
      do ir=nr1,nr2
        do istp=1,4
c
c         obs6 is ett = ur/r + (dut/dt)/r
c
          if(r(ir).le.0.d0)then
            obs(ir,6,istp)=urdr(istp)
     &                    +cics(istp)*cms(istp)*utdr(istp)
          else
            obs(ir,6,istp)=(obs(ir,2,istp)+cics(istp)*cms(istp)
     &                    *obs(ir,3,istp))/dcmplx(r(ir),0.d0)
          endif
c 
c         obs5 now is err = obs5(before) - ett
c
          obs(ir,5,istp)=obs(ir,5,istp)-obs(ir,6,istp)
c
c         obs4 now is ezz
c
          obs(ir,4,istp)=(obs(ir,4,istp)
     &        -clazr*(obs(ir,5,istp)+obs(ir,6,istp)))/(clazr+c2*cmuzr)
c
c         obs7 now is ezr = erz
c
          obs(ir,7,istp)=obs(ir,7,istp)/(c2*cmuzr)
          obs(ir,16,istp)=c05*obs(ir,8,istp)
c
c         obs8 now is ert = etr
c                         = u16 + (dur/dt)/r - ut/r
c                         = 0.5 * (dut/dr + (dur/dt)/r - ut/r)
c
          if(r(ir).le.0.d0)then
            obs(ir,8,istp)=obs(ir,16,istp)
     &              -(cics(istp)*cms(istp)*urdr(istp)+utdr(istp))
          else
            obs(ir,8,istp)=obs(ir,16,istp)
     &                    -(cics(istp)*cms(istp)*obs(ir,2,istp)
     &                    +obs(ir,3,istp))/dcmplx(r(ir),0.d0)
          endif
c
c         obs9 now is ezt = etz
c
          obs(ir,9,istp)=obs(ir,9,istp)/(c2*cmuzr)
c
c         obs10 now is -dur/dz (vertical tilt)
c                      =  obs10(before) - 2 * ezr
c                      =  duz/dr - 2 * ezr
c
          obs(ir,10,istp)= obs(ir,10,istp)-c2*obs(ir,7,istp)
        enddo
      enddo
c
c     end of total wavenumber integral
c
      if(tty)then
        if(istate.eq.-1)then
          write(*,'(a,i7,a,f10.4)')'     Coseismic response:'
     &       //' samples = ',ik,' x_max = ',k*r(nr2)
        else if(istate.eq.0)then
          write(*,'(a,i7,a,f10.4)')'  Steady-state response:'
     &       //' wavenumber samples = ',ik,' x_max = ',k*r(nr2)
        else
          write(*,'(i7,a,E12.5,a,i7,a,f10.4)')lf,'.',dimag(cs)/pi2,
     &       ' Hz: samples = ',ik,' x_max = ',k*r(nr2)
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     transform to outputs                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      xi=dcmplx(alf(nlr)*dm(nlr)/qa(nlr),0.d0)
      do istp=1,4
        do ir=nr1,nr2
c
c         displacement components
c
          du(lf,ir,1,istp)=obs(ir,1,istp)
          du(lf,ir,2,istp)=obs(ir,2,istp)
          du(lf,ir,3,istp)=obs(ir,3,istp)
c
c         stress components
c
          if(isurfcon.eq.1.and.zrec.eq.0.d0)then
            pe=(0.d0,0.d0)
          else
            pe=dcmplx(alf(nlr),0.d0)*obs(ir,12,istp)
          endif
          evs=obs(ir,4,istp)+obs(ir,5,istp)+obs(ir,6,istp)
          du(lf,ir,4,istp)=clazr*evs+c2*cmuzr*obs(ir,4,istp)-pe
          du(lf,ir,5,istp)=clazr*evs+c2*cmuzr*obs(ir,5,istp)-pe
          du(lf,ir,6,istp)=clazr*evs+c2*cmuzr*obs(ir,6,istp)-pe
          du(lf,ir,7,istp)=c2*cmuzr*obs(ir,7,istp)
          du(lf,ir,8,istp)=c2*cmuzr*obs(ir,8,istp)
          du(lf,ir,9,istp)=c2*cmuzr*obs(ir,9,istp)
c
c         tilt components and rotation
c
          du(lf,ir,10,istp)=obs(ir,10,istp)
          du(lf,ir,11,istp)=obs(ir,11,istp)-c2*obs(ir,9,istp)
          du(lf,ir,12,istp)=obs(ir,16,istp)
c
c         pore pressure
c
          du(lf,ir,13,istp)=obs(ir,12,istp)
c
c         Darcy flux components
c
          du(lf,ir,14,istp)=obs(ir,13,istp)
          du(lf,ir,15,istp)=-xi*obs(ir,14,istp)
          du(lf,ir,16,istp)=-xi*obs(ir,15,istp)
        enddo
      enddo
c=============================================================================
c     end of wavenumber integration
c=============================================================================
      if(zrec.eq.0.d0.and.isurfcon.ne.0)then
        do ir=nr1,nr2
          do istp=1,4
            du(lf,ir,4,istp)=(0.d0,0.d0)
            du(lf,ir,7,istp)=(0.d0,0.d0)
            du(lf,ir,9,istp)=(0.d0,0.d0)
          enddo
        enddo
      endif
      return
      end
