      subroutine pegsfctud(dislocation)
      implicit none
c
      double precision dislocation
c
      include 'pegglob.h'
c
      integer i,istp
      double complex cdam,clas,cmus,kappa
c
      double precision pi
      double complex c2,c3,c4
      data pi/3.14159265358979d0/
      data c2,c3,c4/(2.d0,0.d0),(3.d0,0.d0),(4.d0,0.d0)/
c
      do istp=1,4
        do i=1,8
          sfctud0(i,istp)=(0.d0,0.d0)
          sfctud1(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      cdam=dcmplx(dislocation/(2.d0*pi),0.d0)
      cmus=dcmplx(mu(nls),0.d0)
      clas=dcmplx(la(nls),0.d0)+dcmplx(alf(nls)*qa(nls),0.d0)
      kappa=clas+cmus*c2/c3
c
c     explosion (m11=m22=m33=1*kappa)
c
      ms(1)=0
      ics(1)=1
      sfctud0(1,1)=-cdam*kappa/(clas+c2*cmus)
      sfctud1(4,1)=c2*cmus*sfctud0(1,1)
c
c     strike-slip (m12=m21=1*mue)
c
      ms(2)=2
      ics(2)=-1
      sfctud1(4,2)=cdam*cmus
      sfctud1(6,2)=-sfctud1(4,2)
c
c     dip-slip (m13=m31=1*mue)
c
      ms(3)=1
      ics(3)=1
      sfctud0(3,3)=-cdam
      sfctud0(5,3)=sfctud0(3,3)
c
c     compensated linear vector dipole (CLVD) (m11=m22=-1*mue/2, m33=1*mue)
c
      ms(4)=0
      ics(4)=1
      sfctud0(1,4)=-cdam*cmus/(clas+c2*cmus)
      sfctud1(4,4)=cdam*cmus*(c3-c4*cmus/(clas+c2*cmus))/c2
c
      return
      end
