      subroutine pegmatinv(a,cs,k,n,istate)
      implicit none
c
	include 'pegglob.h'
c
c     istate = -1, instantaneous response (s = infinite)
c               0, stationary response (s = 0)
c               1, transient response
c
      integer n,istate
      double precision k
	double complex cs
      double complex a(6,6)
c
      integer i,j
      double complex ck,ck2,cx,cz,c2bk2,cmu0
      double complex d1,cv1,cv2,cv3,cv2s,cv3s
      double complex ca,cb,cqa,cdm,cp,cm,cqaa,et,et0,et1,xi,xi0,xi1
c
	double complex c1,c2,c3,c4
	data c1,c2,c3,c4/(1.d0,0.d0),(2.d0,0.d0),
     &                 (3.d0,0.d0),(4.d0,0.d0)/
c
      ck=dcmplx(k,0.d0)
      ck2=dcmplx(k*k,0.d0)
	ca=dcmplx(alf(n),0.d0)
	cqa=dcmplx(qa(n),0.d0)
	cdm=dcmplx(dm(n),0.d0)
	cqaa=cqa*ca
      et=cla(n)+cmu(n)
	et1=et+cqaa
      xi=cla(n)+c2*cmu(n)
	xi1=xi+cqaa
	cb=cdm*xi/xi1
      c2bk2=c2*cb*ck2
	cv1=cqaa*cmu(n)*c2bk2
	cv2=c2*xi*xi1
	cv3=c2*cmu(n)*ck*cv2
	cv2s=cv2*cs
	cv3s=cv3*cs
c
	if(istate.eq.1)then
c
c       for transient responses
c
        a(3,1)=cmu(n)/(c2*xi1)
        a(3,2)=c1/(c4*xi1*ck)
        a(3,3)=-a(3,1)
        a(3,4)=-a(3,2)
        a(3,5)=(0.d0,0.d0)
        a(3,6)=(0.d0,0.d0)
c
        a(4,1)=a(3,1)
        a(4,2)=-a(3,2)
        a(4,3)=a(3,1)
        a(4,4)=-a(3,2)
        a(4,5)=(0.d0,0.d0)
        a(4,6)=(0.d0,0.d0)
c
        if(smalls(n))then
c
c         a(1,j)=a(1,j)+ca*a(5,j)
c         a(2,j)=a(2,j)+ca*a(6,j)
c         a(5,j)=a(5,j)*cs
c         a(6,j)=a(6,j)*cs
c
c         d1 = (k/kp - 1)/s = -1/b/(kp*(k+kp))
c
          d1=-c1/(cb*ka(n)*(ck+ka(n)))
c
          a(1,1)=(0.d0,0.d0)
          a(1,2)=(-cv1*ka(n)*d1/ck+xi*cmu(n))/cv3
          a(1,3)=(cv1*ka(n)*d1/ck+xi*et1)/cv2
          a(1,4)=xi*xi1/cv3
          a(1,5)=-ca*cb*ka(n)*d1/(c2*xi)
          a(1,6)=(0.d0,0.d0)
c
          a(2,1)=a(1,1)
          a(2,2)=-a(1,2)
          a(2,3)=-a(1,3)
          a(2,4)=a(1,4)
          a(2,5)=-a(1,5)
          a(2,6)=a(1,6)
c
          a(5,1)=cqa*cmu(n)*c2bk2/cv2
          a(5,2)=cb*cqa*ka(n)/cv2
          a(5,3)=-c2*cmu(n)*ck*a(5,2)
          a(5,4)=-cb*cqa*ck/cv2
          a(5,5)=cb*ka(n)/(c2*xi)
          a(5,6)=-cqa/(c2*ca*xi1)
c
          a(6,1)=a(5,1)
          a(6,2)=-a(5,2)
          a(6,3)=-a(5,3)
          a(6,4)=a(5,4)
          a(6,5)=-a(5,5)
          a(6,6)=a(5,6)
        else
          a(1,1)=-cv1/cv2s
          a(1,2)=(-cv1+cs*xi*cmu(n))/cv3s
          a(1,3)=(cv1+cs*xi*et1)/cv2s
          a(1,4)=(cv1+cs*xi*xi1)/cv3s
          a(1,5)=-ca*cb*ck/(c2*xi*cs)
          a(1,6)=cqa/(c2*xi1*cs)
c
          a(2,1)=a(1,1)
          a(2,2)=-a(1,2)
          a(2,3)=-a(1,3)
          a(2,4)=a(1,4)
          a(2,5)=-a(1,5)
          a(2,6)=a(1,6)
c
          a(5,1)=cqa*cmu(n)*c2bk2/cv2s
          a(5,2)=cb*cqa*ka(n)/cv2s
          a(5,3)=-c2*cmu(n)*ck*a(5,2)
          a(5,4)=-cb*cqa*ck/cv2s
          a(5,5)=cb*ka(n)/(c2*xi*cs)
          a(5,6)=-cqa/(c2*ca*xi1*cs)
c
          a(6,1)=a(5,1)
          a(6,2)=-a(5,2)
          a(6,3)=-a(5,3)
          a(6,4)=a(5,4)
          a(6,5)=-a(5,5)
          a(6,6)=a(5,6)
        endif
c
c       transform y6 -> y6+a*s*y1
c
        do i=1,6
          a(i,1)=a(i,1)-ca*cs*a(i,6)
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
        a(1,1)=(0.d0,0.d0)
        a(1,2)=c1/(c4*xi0*ck)
        a(1,3)=et0/(c2*xi0)
        a(1,4)=c1/(c4*cmu0*ck)
        a(1,5)=(0.d0,0.d0)
        a(1,6)=(0.d0,0.d0)
c
        a(2,1)=(0.d0,0.d0)
        a(2,2)=-a(1,2)
        a(2,3)=-a(1,3)
        a(2,4)=a(1,4)
        a(2,5)=(0.d0,0.d0)
        a(2,6)=(0.d0,0.d0)
c
        a(3,1)=cmu0/(c2*xi0)
        a(3,2)=c1/(c4*xi0*ck)
        a(3,3)=-a(3,1)
        a(3,4)=-a(3,2)
        a(3,5)=(0.d0,0.d0)
        a(3,6)=(0.d0,0.d0)
c
        a(4,1)=a(3,1)
        a(4,2)=-a(3,2)
        a(4,3)=a(3,1)
        a(4,4)=-a(3,2)
        a(4,5)=(0.d0,0.d0)
        a(4,6)=(0.d0,0.d0)
c
        a(5,1)=(0.d0,0.d0)
        a(5,2)=(0.d0,0.d0)
        a(5,3)=(0.d0,0.d0)
        a(5,4)=(0.d0,0.d0)
        a(5,5)=(0.5d0,0.d0)
        a(5,6)=(0.5d0,0.d0)/(ck*cmu0)
c
        a(6,1)=(0.d0,0.d0)
        a(6,2)=(0.d0,0.d0)
        a(6,3)=(0.d0,0.d0)
        a(6,4)=(0.d0,0.d0)
        a(6,5)=(0.5d0,0.d0)
        a(6,6)=-(0.5d0,0.d0)/(ck*cmu0)
      endif
c
      return
      end