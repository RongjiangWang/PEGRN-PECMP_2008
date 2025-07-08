      subroutine pegmatrix(a,cs,k,z,n,istate)
      implicit none
c
	include 'pegglob.h'
c
c       istate = -1, instantaneous response (s = infinite)
c                 0, stationary response (s = 0)
c                 1, transient response
c
      integer n,istate
      double precision k,z
	double complex cs
      double complex a(6,6)
c
      integer i,j
      double complex ck,cx,cz,c2bk2,cmu0,d1,d2,d3,dx
      double complex ca,cb,cqa,cdm,cp,cm,cqaa,et,et0,et1,xi,xi0,xi1
c
	double complex c1,c2,c3,c4,c5,c6
	data c1,c2,c3,c4,c5,c6/(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0),
     &                       (4.d0,0.d0),(5.d0,0.d0),(6.d0,0.d0)/
c
      ck=dcmplx(k,0.d0)
      cz=dcmplx(z,0.d0)
      cx=dcmplx(k*z,0.d0)
	ca=dcmplx(alf(n),0.d0)
	cqa=dcmplx(qa(n),0.d0)
	cdm=dcmplx(dm(n),0.d0)
	cp=c1+cx
	cm=c1-cx
	cqaa=cqa*ca
      et=cla(n)+cmu(n)
	et1=et+cqaa
      xi=cla(n)+c2*cmu(n)
	xi1=xi+cqaa
	cb=cdm*xi/xi1
      c2bk2=c2*cb*ck*ck
c
      if(istate.eq.1)then
c
        a(1,1)=c1
	  a(2,1)=c2*cmu(n)*ck
	  a(3,1)=c1
	  a(4,1)=c2*cmu(n)*ck
	  a(5,1)=(0.d0,0.d0)
	  a(6,1)=(0.d0,0.d0)
c
        a(1,2)=c1
	  a(2,2)=-c2*cmu(n)*ck
	  a(3,2)=-c1
	  a(4,2)=c2*cmu(n)*ck
	  a(5,2)=(0.d0,0.d0)
	  a(6,2)=(0.d0,0.d0)
c
        a(1,3)=c1+et1*cm/cmu(n)
	  a(2,3)=c2*et1*cm*ck
	  a(3,3)=-c1-et1*cx/cmu(n)
	  a(4,3)=-c2*et1*cx*ck
	  a(5,3)=-c2*cqa*ck
	  a(6,3)=c2*ca*cdm*ck*ck
c
        a(1,4)=c1+et1*cp/cmu(n)
	  a(2,4)=-c2*et1*cp*ck
	  a(3,4)=c1-et1*cx/cmu(n)
	  a(4,4)=c2*et1*cx*ck
	  a(5,4)=c2*cqa*ck
	  a(6,4)=c2*ca*cdm*ck*ck
c
        if(smalls(n))then
c
c         yi5 => (yi5 - alpha*yi1*exp(k-kp)*z))/s
c         yi6 => (yi6 - alpha*yi2*exp(kp-k)*z))/s
c
c         d1 = (k/kp - 1)/s = -1/b/(kp*(k+kp))
c         d2 = (1 - exp((k-kp)*z))/s = (1 - exp(kp*s*z*d1))/s
c
          d1=-c1/(cb*ka(n)*(ck+ka(n)))
          dx=ka(n)*cs*cz*d1
          if(cdabs(dx).gt.1.0d-02)then
            d2=(c1-cdexp( dx))/cs
            d3=(c1-cdexp(-dx))/cs
          else
            d2=-dx/cs*(c1+dx/c2*(c1+dx/c3*(c1
     &         +dx/c4*(c1+dx/c5*(c1+dx/c6)))))
            d3= dx/cs*(c1-dx/c2*(c1-dx/c3*(c1
     &         -dx/c4*(c1-dx/c5*(c1-dx/c6)))))
          endif
c
          a(1,5)=ca*d2
	    a(2,5)=c2*ca*cmu(n)*ck*(d1+d2)
	    a(3,5)=ca*(d1+d2)
	    a(4,5)=c2*ca*cmu(n)*ck*d2
	    a(5,5)=xi/(cb*ka(n))
	    a(6,5)=-ca*xi1/cqa
c
          a(1,6)=ca*d3
	    a(2,6)=-c2*ca*cmu(n)*ck*(d1+d3)
	    a(3,6)=-ca*(d1+d3)
	    a(4,6)=c2*ca*cmu(n)*ck*d3
	    a(5,6)=-xi/(cb*ka(n))
	    a(6,6)=-ca*xi1/cqa
        else
          a(1,5)=ca
	    a(2,5)=c2*ca*cmu(n)*ck*ck/ka(n)
	    a(3,5)=ca*ck/ka(n)
	    a(4,5)=c2*ca*cmu(n)*ck
	    a(5,5)=xi*cs/(cb*ka(n))
	    a(6,5)=-cs*ca*xi1/cqa
c
          a(1,6)=ca
	    a(2,6)=-c2*ca*cmu(n)*ck*ck/ka(n)
	    a(3,6)=-ca*ck/ka(n)
	    a(4,6)=c2*ca*cmu(n)*ck
	    a(5,6)=-xi*cs/(cb*ka(n))
	    a(6,6)=-cs*ca*xi1/cqa
        endif
c
c       transform y6 -> y6+a*s*y1
c
        do j=1,6
          a(6,j)=a(6,j)+ca*cs*a(1,j)
        enddo
      else
        cmu0=dcmplx(mu(n),0.d0)
        if(istate.eq.-1)then
c
c         for instantaneous responses (cs = infinite)
c
          xi0=dcmplx(la(n)+2.d0*mu(n)+alf(n)*qa(n),0.d0)
          et0=dcmplx(la(n)+mu(n)+alf(n)*qa(n),0.d0)
        else
c
c         for stationary responses (cs = 0)
c
          xi0=dcmplx(la(n)+2.d0*mu(n),0.d0)
          et0=dcmplx(la(n)+mu(n),0.d0)
        endif
c
        a(1,1)=c1
	  a(2,1)=c2*cmu0*ck
	  a(3,1)=c1
	  a(4,1)=a(2,1)
        a(5,1)=(0.d0,0.d0)
        a(6,1)=(0.d0,0.d0)
c
        a(1,2)=c1
	  a(2,2)=-a(2,1)
	  a(3,2)=-c1
	  a(4,2)=a(2,1)
        a(5,2)=(0.d0,0.d0)
        a(6,2)=(0.d0,0.d0)
c
        a(1,3)=c1+et0*cm/cmu0
	  a(2,3)=c2*et0*cm*ck
	  a(3,3)=-c1-et0*cx/cmu0
	  a(4,3)=-c2*et0*cx*ck
        a(5,3)=(0.d0,0.d0)
        a(6,3)=(0.d0,0.d0)
c
        a(1,4)=c1+et0*cp/cmu0
	  a(2,4)=-c2*et0*cp*ck
	  a(3,4)=c1-et0*cx/cmu0
	  a(4,4)=c2*et0*cx*ck
        a(5,4)=(0.d0,0.d0)
        a(6,4)=(0.d0,0.d0)
c
c
        a(1,5)=(0.d0,0.d0)
        a(2,5)=(0.d0,0.d0)
        a(3,5)=(0.d0,0.d0)
        a(4,5)=(0.d0,0.d0)
        a(5,5)=(1.d0,0.d0)
        a(6,5)=ck*cmu0
c
        a(1,6)=(0.d0,0.d0)
        a(2,6)=(0.d0,0.d0)
        a(3,6)=(0.d0,0.d0)
        a(4,6)=(0.d0,0.d0)
        a(5,6)=(1.d0,0.d0)
        a(6,6)=-ck*cmu0
	endif
c
      return
      end