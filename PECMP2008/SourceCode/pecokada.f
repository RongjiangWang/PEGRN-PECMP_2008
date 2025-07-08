      subroutine pecokada(ns,nrec,tmax)
      implicit none
c
      integer ns,nrec
      double precision tmax
c
      include 'pecglob.h'
c
c     Last modified: Potsdam, Nov, 2003, by R. Wang
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     from Okada's subroutine DC3D0:
c
      INTEGER IRET
      REAL*4 ALPHA,X,Y,Z,DEPTH,DIPS,DX,DY,DZ,
     &       UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
c
c     more from Okada's subroutine DC3D:
c
      REAL*4 AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer ieq,is,iptch,irec
      double precision larecu,st,di,step_s,step_d,disn,dise
      double precision csst,ssst,csra,ssra,csdi,ssdi
      double precision cs2st,ss2st,ds,eii0,eii1,eii2,eii3
      double precision strain(6)
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     PROCESSING:
c     ===========
c     coobs(1)=' Ux'
c     coobs(2)=' Uy'
c     coobs(3)=' Uz'
c     coobs(4)='Sxx'
c     coobs(5)='Syy'
c     coobs(6)='Szz'
c     coobs(7)='Sxy'
c     coobs(8)='Syz'
c     coobs(9)='Szx'
c     coobs(10)=' Tx'
c     coobs(11)=' Ty'
c     coobs(12)='Rot'
c     coobs(13)=' Pp'
c     coobs(14)='Dcx'
c     coobs(15)='Dcy'
c     coobs(16)='Dcz'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     receiver and source independent variables
c
      larecu=larec+alfrec*qarec
      ALPHA=sngl((larecu+murec)/(larecu+2.d0*murec))
      Z=-sngl(zrec)
c
      do is=1,ns
        if(tstart(is).gt.tmax)goto 100
        ieq=ieqno(is)
        st=strike(is)*DEG2RAD
        csst=dcos(st)
        ssst=dsin(st)
        cs2st=dcos(2.d0*st)
        ss2st=dsin(2.d0*st)
c
	  di=dip(is)*DEG2RAD
        csdi=dcos(di)
        ssdi=dsin(di)
c
        DEPTH=sngl(zref(is))
        DIPS=sngl(dip(is))
c
        step_s=length(is)/dble(nptch_s(is))
        step_d=width(is)/dble(nptch_d(is))
c
        do iptch=1,nptch_s(is)*nptch_d(is)
          DISL1=sngl( slip_s(is,iptch))
          DISL2=sngl(-slip_d(is,iptch))
          DISL3=sngl(opening(is,iptch))
c
c         for extended source
c
c
          AL1=sngl(ptch_s(is,iptch)-0.5d0*step_s)
          AL2=sngl(ptch_s(is,iptch)+0.5d0*step_s)
          AW1=sngl(-ptch_d(is,iptch)-0.5d0*step_d)
          AW2=sngl(-ptch_d(is,iptch)+0.5d0*step_d)
          ds=0.001d0*step_s
          DZ=sngl(ds)
          do irec=1,nrec
c
c           transform from Aki's to Okada's system
c
            call disazi(REARTH,xref(is),yref(is),
     &                  xrec(irec),yrec(irec),disn,dise)
            X=sngl(disn*csst+dise*ssst)
            Y=sngl(disn*ssst-dise*csst)
            IRET=1
            call DC3D(ALPHA,X,Y,Z,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &            DISL1,DISL2,DISL3,UX,UY,UZ,
     &            UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
c
c           transform from Okada's to Aki's system
c
            coobs(ieq,irec,1)=coobs(ieq,irec,1)
     &                      +dble(UX)*csst+dble(UY)*ssst
            coobs(ieq,irec,2)=coobs(ieq,irec,2)
     &                      +dble(UX)*ssst-dble(UY)*csst
            coobs(ieq,irec,3)=coobs(ieq,irec,3)-dble(UZ)
c
            strain(1)=dble(UXX)*csst*csst+dble(UYY)*ssst*ssst
     &               +0.5d0*dble(UXY+UYX)*ss2st
            strain(2)=dble(UXX)*ssst*ssst+dble(UYY)*csst*csst
     &               -0.5d0*dble(UXY+UYX)*ss2st
            strain(3)=dble(UZZ)
            strain(4)=0.5d0*dble(UXX-UYY)*ss2st
     &               -0.5d0*dble(UXY+UYX)*cs2st
            strain(5)=-0.5d0*dble(UZX+UXZ)*ssst
     &               +0.5d0*dble(UYZ+UZY)*csst
            strain(6)=-0.5d0*dble(UZX+UXZ)*csst
     &               -0.5d0*dble(UYZ+UZY)*ssst
            eii0=dble(UXX+UYY+UZZ)
            coobs(ieq,irec,4)=coobs(ieq,irec,4)
     &                      +larec*eii0+2.d0*murec*strain(1)
            coobs(ieq,irec,5)=coobs(ieq,irec,5)
     &                      +larec*eii0+2.d0*murec*strain(2)
            coobs(ieq,irec,6)=coobs(ieq,irec,6)
     &                      +larec*eii0+2.d0*murec*strain(3)
            coobs(ieq,irec,7)=coobs(ieq,irec,7)+2.d0*murec*strain(4)
            coobs(ieq,irec,8)=coobs(ieq,irec,8)+2.d0*murec*strain(5)
            coobs(ieq,irec,9)=coobs(ieq,irec,9)+2.d0*murec*strain(6)
c
            coobs(ieq,irec,10)=coobs(ieq,irec,10)
     &                       -(dble(UXZ)*csst+dble(UYZ)*ssst)
            coobs(ieq,irec,11)=coobs(ieq,irec,11)
     &                       -(dble(UXZ)*ssst-dble(UYZ)*csst)
            coobs(ieq,irec,12)=coobs(ieq,irec,12)-0.5d0*dble(UYX-UXY)
            if((isurfcon.eq.1.or.isurfcon.eq.3).and.zrec.eq.0.d0)then
              coobs(ieq,irec,4)=coobs(ieq,irec,4)+alfrec*qarec*eii0
              coobs(ieq,irec,5)=coobs(ieq,irec,5)+alfrec*qarec*eii0
              coobs(ieq,irec,6)=coobs(ieq,irec,6)+alfrec*qarec*eii0
            else
              coobs(ieq,irec,13)=coobs(ieq,irec,13)-alfrec*qarec*eii0
            endif
c
            if(isurfcon.ne.1.and.isurfcon.ne.3.or.zrec.gt.0.d0)then
              DX=sngl(ds*csst)
              DY=sngl(ds*ssst)
              call DC3D(ALPHA,X+DX,Y+DY,Z,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &                  DISL1,DISL2,DISL3,UX,UY,UZ,
     &                  UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
              eii1=dble(UXX+UYY+UZZ)
              coobs(ieq,irec,14)=coobs(ieq,irec,14)
     &                          +alfrec*dmrec*(eii1-eii0)/ds
c
              DX=sngl(ds*ssst)
              DY=sngl(-ds*csst)
              call DC3D(ALPHA,X+DX,Y+DY,Z,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &                  DISL1,DISL2,DISL3,UX,UY,UZ,
     &                  UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
              eii2=dble(UXX+UYY+UZZ)
              coobs(ieq,irec,15)=coobs(ieq,irec,15)
     &                          +alfrec*dmrec*(eii2-eii0)/ds
            endif
c
            if(isurfcon.ne.2.or.zrec.gt.0.d0)then
              call DC3D(ALPHA,X,Y,Z-DZ,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &              DISL1,DISL2,DISL3,UX,UY,UZ,
     &              UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
              eii3=dble(UXX+UYY+UZZ)
              coobs(ieq,irec,16)=coobs(ieq,irec,16)
     &                          +alfrec*dmrec*(eii3-eii0)/ds
            endif
          enddo
        enddo
100     continue
      enddo
c
      return
      end
