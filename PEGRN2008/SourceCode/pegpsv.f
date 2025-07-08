      subroutine pegpsv(y,cs,k,istate)
      implicit none
c
c     calculation of response to p-sv source
c     y(8,4): solution vector (complex)
c     k: wave number
c
      integer istate
      double precision k
      double complex cs
      double complex y(8,4)
c
      include 'pegglob.h'
c
c     work space
c
      integer i,istp,j,l,n,lup,llw,lmd,ly,key
      double precision hpl
      double complex ck
      double complex b(6,4),cy(6,4)
      double complex yup(6,3),ylw(6,3)
      double complex ma0(6,6),mat(6,6),mai(6,6),c0(6,3),coef(6,6)
      double complex yup0(6,3),ylw0(6,3)
c
c     psv layer matrics
c
      double complex maup(6,6,nzmax),maiup(6,6,nzmax)
      double complex malw(6,6,nzmax),mailw(6,6,nzmax)
c
      ck=dcmplx(k,0.d0)
c
c===============================================================================
c
      lup=1
      llw=lp
      lmd=ls
c
      do l=lup,ls-1
        ly=l
        n=nno(ly)
        hpl=hp(ly)
        call pegmatrix(maup(1,1,ly),cs,k,hpl,n,istate)
        call pegmatinv(maiup(1,1,ly),cs,k,n,istate)
      enddo
      do l=ls,llw-1
        ly=l
        n=nno(ly)
        hpl=-hp(ly)
        call pegmatrix(malw(1,1,ly),cs,k,hpl,n,istate)
        call pegmatinv(mailw(1,1,ly),cs,k,n,istate)
      enddo
c
c     matrix propagation from surface to source
c
      do j=1,3
        do i=1,6
          yup(i,j)=(0.d0,0.d0)
        enddo
      enddo
c
	if(isurfcon.eq.1)then
c
c	  free surface: p = 0
c
	  yup(1,1)=(1.d0,0.d0)
        yup(3,2)=(1.d0,0.d0)
        if(istate.eq.1)then
	    yup(6,3)=(1.d0,0.d0)
        else
	    yup(5,3)=(1.d0,0.d0)
        endif
	else if(isurfcon.eq.2)then
c
c	  confined surface surface: dp/dz = 0
c
	  yup(1,1)=(1.d0,0.d0)
        if(istate.eq.1)then
	    yup(6,1)=dcmplx(alf(nno(lup)),0.d0)*cs*yup(1,1)
        endif
        yup(3,2)=(1.d0,0.d0)
	  yup(5,3)=(1.d0,0.d0)
	else
c
c	  unsaturated surface: dp/dt = -v_z*rho*g/n
c
	  yup(1,1)=(1.d0,0.d0)
        if(istate.eq.1)then
	    yup(6,1)=dcmplx(alf(nno(lup)),0.d0)*cs*yup(1,1)
        endif
        yup(3,2)=(1.d0,0.d0)
	  yup(5,3)=(1.d0,0.d0)
        if(istate.eq.1)then
          yup(6,3)=-yup(5,3)*dcmplx(porosity/rhog,0.d0)*cs
        endif
	endif
c
      if(lup.eq.lzrec)call cmemcpy(yup,yup0,18)
c
      call pegproppsv(maup,maiup,lup,lmd,yup,yup0,k)
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
c     coefficient vectors in the half-space
c
      do j=1,3
        do i=1,6
          c0(i,j)=(0.d0,0.d0)
        enddo
      enddo
      c0(2,1)=(1.d0,0.d0)
      c0(4,2)=(1.d0,0.d0)
      c0(6,3)=(1.d0,0.d0)
      n=nno(llw)
      call pegmatrix(ma0,cs,k,0.d0,n,istate)
      call caxcb(ma0,c0,6,6,3,ylw)
c
      if(llw.eq.lzrec)call cmemcpy(ylw,ylw0,18)
c
      call pegproppsv(malw,mailw,llw,lmd,ylw,ylw0,k)
c
c===============================================================================
c
c     conditions on the source surface
c
c
c     source function
c
      if(istate.eq.-1)then
        do istp=1,4
          do i=1,6
            b(i,istp)=sfctud0(i,istp)+sfctud1(i,istp)*ck
          enddo
        enddo
      else
        do istp=1,4
          do i=1,4
            b(i,istp)=sfct0(i,istp)+sfct1(i,istp)*ck
          enddo
          do i=5,6
            b(i,istp)=sfct0(i+2,istp)+sfct1(i+2,istp)*ck
          enddo
        enddo
      endif
      do i=1,6
        do j=1,3
          coef(i,j)=yup(i,j)
          coef(i,j+3)=-ylw(i,j)
        enddo
      enddo
      key=0
      call cdsvd500(coef,b,6,4,0.d0,key)
      if(key.eq.0)then
        print *,'warning in pegpsv: anormal exit from cdgemp!'
        return
      endif
      if(lzrec.lt.ls)then
        do istp=1,4
          do i=1,6
            cy(i,istp)=(0.d0,0.d0)
            do j=1,3
              cy(i,istp)=cy(i,istp)+b(j,istp)*yup0(i,j)
            enddo
          enddo
        enddo
      else if(lzrec.gt.ls)then
        do istp=1,4
          do i=1,6
            cy(i,istp)=(0.d0,0.d0)
            do j=1,3
              cy(i,istp)=cy(i,istp)+b(j+3,istp)*ylw0(i,j)
            enddo
          enddo
        enddo
      else
        do istp=1,4
          do i=1,6
            cy(i,istp)=(0.d0,0.d0)
            do j=1,3
              cy(i,istp)=cy(i,istp)+(0.5d0,0.d0)
     &             *(b(j,istp)*yup0(i,j)+b(j+3,istp)*ylw0(i,j))
            enddo
          enddo
        enddo
      endif
c
      if(istate.eq.1)then
        do istp=1,4
          do i=1,4
            y(i,istp)=cy(i,istp)
          enddo
          y(7,istp)=cy(5,istp)
          y(8,istp)=cy(6,istp)-dcmplx(alf(nlr),0.d0)*cs*cy(1,istp)
        enddo
      else
        do istp=1,4
          do i=1,6
            y(i,istp)=cy(i,istp)
          enddo
        enddo
      endif
c
      return
      end
