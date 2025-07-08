      subroutine pegmoduli(cs,istate)
      implicit none
c
      integer istate
      double complex cs
c
      include 'pegglob.h'
c
      integer n
c
c     only for elastic case
c
      do n=1,n0
        cmu(n)=dcmplx(mu(n),0.d0)
        cla(n)=dcmplx(la(n),0.d0)
      enddo
c
      return
      end
