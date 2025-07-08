      subroutine pegkern(y,cs,k,istate,zerokern)
      implicit none
c
c     calculation of response function in Laplace domain
c     y1(8,4): solution vector (complex)
c     k: wave number (input)
c
      integer istate
      double precision k
      double complex cs
      double complex y(8,4)
      logical zerokern
c
      include 'pegglob.h'
c
      integer i,n,istp
      double complex ck,cla0,cmu0,qa0,alf0,cdm0,et,etu
      double complex ca,cb,cd,p,dc,evs,devs,p0,dc0,evs0,devs0
      double complex y0(8,4),y1(8,4)
c
      double complex c0,c1,c2
      data c0,c1,c2/(0.d0,0.d0),(1.d0,0.d0),(2.d0,0.d0)/
c
      do istp=1,4
        do i=1,8
          y0(i,istp)=(0.d0,0.d0)
          y1(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
	ck=dcmplx(k,0.d0)
      alf0=dcmplx(alf(nlr),0.d0)
      qa0=dcmplx(qa(nlr),0.d0)
      cla0=dcmplx(la(nlr),0.d0)+alf0*qa0
      cmu0=dcmplx(mu(nlr),0.d0)
      cdm0=dcmplx(dm(nlr),0.d0)
      cd=cdm0/(c1+alf0*qa0/(cla0+c2*cmu0))
      ca=dcmplx(alf(nlr)*dm(nlr)/qa(nlr),0.d0)/cdsqrt(cd)
c
      if(istate.eq.-1)then
	  do n=1,n0
          ka(n)=ck
        enddo
        call pegpsv(y1,c0,k,-1)
        call peghskern(y0,k,cla0,cmu0)
      else
	  do n=1,n0
          et=cla(n)+c2*cmu(n)
          etu=et+dcmplx(alf(n)*qa(n),0.d0)
	    cb=et*dcmplx(dm(n),0.d0)/etu
          ka(n)=cdsqrt(ck*ck+cs/cb)
          if(cdabs(cs/cb)/(k*k).le.0.1d0.and.
     &       cdabs(cs/(cb*(ck+ka(n))))*h(n).le.10.d0)then
            smalls(n)=.true.
          else
            smalls(n)=.false.
          endif
        enddo
        call pegpsv(y1,cs,k,istate)
        call pegsh(y1,cs,k,istate)
c
	  do n=1,n0
          ka(n)=ck
        enddo
        call pegpsv(y0,c0,k,-1)
      endif
c
      do istp=1,4
        if(istate.eq.-1)then
          evs=(y1(2,istp)-c2*cmu0*ck*y1(3,istp))/(cla0+c2*cmu0)
          devs=ck*(c2*cmu0*ck*y1(1,istp)-y1(4,istp))/(cla0+c2*cmu0)
          p=-qa0*evs
          dc=cdm0*alf0*devs
          y1(7,istp)=p
          y1(8,istp)=dc
        endif
        evs0=(y0(2,istp)-c2*cmu0*ck*y0(3,istp))/(cla0+c2*cmu0)
        devs0=ck*(c2*cmu0*ck*y0(1,istp)-y0(4,istp))/(cla0+c2*cmu0)
        p0=-qa0*evs0
        dc0=cdm0*alf0*devs0
        y0(7,istp)=p0
        y0(8,istp)=dc0
c
c       subtract effective stress from normal stress
c
        y1(2,istp)=y1(2,istp)+alf0*y1(7,istp)
        y0(2,istp)=y0(2,istp)+alf0*y0(7,istp)
c
        if(zrec.eq.0.d0)then
          if(isurfcon.eq.1)then
            y1(7,istp)=(0.d0,0.d0)
            y0(7,istp)=(0.d0,0.d0)
          else if(isurfcon.eq.2)then
            y1(8,istp)=(0.d0,0.d0)
            y0(8,istp)=(0.d0,0.d0)
          else if(isurfcon.eq.3)then
            if(istate.eq.-1)y1(7,istp)=(0.d0,0.d0)
            y0(7,istp)=(0.d0,0.d0)
          endif
        endif
c
        if(istate.eq.1.and.nearsurf.eq.1)then
          y0(8,istp)=y0(8,istp)-ca*cdsqrt(cs)*p0
     &              *cdexp(-cdsqrt(cs/cd)*dcmplx(zrec,0.d0))
        endif
c
        do i=1,8
          y(i,istp)=y1(i,istp)-y0(i,istp)
        enddo
      enddo
c
      zerokern=.true.
      do istp=1,4
        do i=1,8
          if(cdabs(y(i,istp)).le.1.d-06*cdabs(y0(i,istp)))then
            y(i,istp)=(0.d0,0.d0)
          else
            zerokern=.false.
          endif
        enddo
      enddo
c
      return
      end	  
