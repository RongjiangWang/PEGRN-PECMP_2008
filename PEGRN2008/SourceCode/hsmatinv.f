      subroutine hsmatinv(a,k,clahs,cmuhs)
      implicit none
c
      double precision k
      double complex clahs,cmuhs,a(6,6)
c
      double complex ck,xihs,eths
c
      double complex c1,c2,c3,c4
      data c1,c2,c3,c4/(1.d0,0.d0),(2.d0,0.d0),
     &                 (3.d0,0.d0),(4.d0,0.d0)/
c
      ck=dcmplx(k,0.d0)
      xihs=clahs+(2.d0,0.d0)*cmuhs
      eths=clahs+cmuhs
c
      a(1,1)=(0.d0,0.d0)
      a(1,2)=c1/(c4*xihs*ck)
      a(1,3)=eths/(c2*xihs)
      a(1,4)=c1/(c4*cmuhs*ck)
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
      a(3,1)=cmuhs/(c2*xihs)
      a(3,2)=c1/(c4*xihs*ck)
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
      a(5,6)=(0.5d0,0.d0)/(ck*cmuhs)
c
      a(6,1)=(0.d0,0.d0)
      a(6,2)=(0.d0,0.d0)
      a(6,3)=(0.d0,0.d0)
      a(6,4)=(0.d0,0.d0)
      a(6,5)=(0.5d0,0.d0)
      a(6,6)=-(0.5d0,0.d0)/(ck*cmuhs)
c
      return
      end
