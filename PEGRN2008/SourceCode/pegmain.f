      program pegmain
      implicit none
c
      include 'pegglob.h'
c
c     work space
c
      integer i,j,l,ir,nr,it,nt,izs,ierr
      integer istp,isp,nr1,nr2,nzs,nprf
      integer lend,lenf,leninp,iunit
      integer unit(16,4),idec(nrmax),nout(nrmax)
      double precision am,r1,r2,dr,dract,sampratio,rsmin,supalias
      double precision twindow,dt,pi,zs1,zs2,dzs,zrs2,beta,accuracy
      double precision swap,vp,vs,rho,lau,nu,nuu,skempton,diffus,gamma
      character*35 stype(4)
      character*35 comptxt(16)
      character*80 inputfile,fname(16),outdir
      character*163 green(16,4)
      character*180 dataline
c
      pi=4.d0*datan(1.d0)
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#                  Welcome to                        #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#     PPPP    EEEEE    GGGG   RRRR    N   N          #'
      print *,'#     P   P   E       G       R   R   NN  N          #'
      print *,'#     PPPP    EEEEE   G GGG   RRRR    N N N          #'
      print *,'#     P       E       G   G   R R     N  NN          #'
      print *,'#     P       EEEEE    GGGG   R  R    N   N          #'
      print *,'#                                                    #'
      print *,'#                  Version 2008                      #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#                     July 2008                      #'
      print *,'######################################################'
      print *,'                                                      '
c
      nwarn=0
      write(*,'(a,$)')' Please type the file name of input data: '
      read(*,'(a)')inputfile
      open(10,file=inputfile,status='old')
c
c     parameters for source-observation array
c     =======================================
c
      call getdata(10,dataline)
      read(dataline,*)zrec
      zrec=zrec*km2m
      call getdata(10,dataline)
      read(dataline,*)nr,r1,r2,sampratio
      if(sampratio.lt.1.d0)then
        stop 'Error: max. to min. sampling ratio < 1!'
      endif
      r1=r1*km2m
      r2=r2*km2m
      if(r1.gt.r2)then
        swap=r1
        r1=r2
        r2=swap
      endif
      if(r1.lt.0.d0.or.r2.lt.0.d0.or.nr.lt.1)then
        stop 'Error: wrong no of distance samples!'
      else if(nr.gt.nrmax)then
        stop 'Error: max. no of distance samples exceeded!'
      else if(nr.eq.1.or.r1.eq.r2)then
        r2=r1
        nr=1
        dr=0.d0
        r(1)=r1
      else if(nr.eq.2)then
        dr=r2-r1
        r(1)=r1
        r(2)=r2
      else
        dr=2.d0*(r2-r1)/dble(nr-1)/(1.d0+sampratio)
        r(1)=r1
        do i=2,nr
          dract=dr*(1.d0+(sampratio-1.d0)*dble(i-2)/dble(nr-2))
          r(i)=r(i-1)+dract
        enddo
      endif
c
      call getdata(10,dataline)
      read(dataline,*)nzs,zs1,zs2
      if(zs1.gt.zs2)then
        swap=zs1
        zs1=zs2
        zs2=swap
      endif
      zs1=zs1*km2m
      zs2=zs2*km2m
      if(zs1.le.0.d0.or.zs2.le.0.d0)then
        stop 'Error: source depths <= 0!'
      else if(nzs.lt.1)then
        stop 'Error: wrong no of source depths!'
      else if(nzs.gt.nzsmax)then
        stop 'Error: max. no of source depths exceeded!'
      else if(nzs.eq.1.or.zs1.eq.zs2)then
        nzs=1
        dzs=0.d0
      else
        dzs=(zs2-zs1)/dble(nzs-1)
        if(zs1.lt.zrec+0.5d0*dzs.and.zs2.gt.zrec-0.5d0*dzs)then
          zs=zrec+0.5d0*dzs
10        zs=zs-dzs
          if(zs.gt.zs1)goto 10
          if(zs.gt.0.d0)then
            zs1=zs
          else
            zs1=zs+dzs
          endif
        endif
        zs2=zs1+dble(nzs-1)*dzs
      endif
c
      call getdata(10,dataline)
      read(dataline,*)nt,twindow
      if(twindow.le.0.d0)then
        stop ' Error in input: wrong time window!'
      else if(nt.le.0)then
        stop ' Error in input: time sampling no <= 0!'
      endif
      twindow=twindow*day2sec
      if(nt.le.2)then
        dt=twindow
      else
        dt=twindow/dble(nt-1)
      endif
c
c     wavenumber integration parameters
c     =================================
c
      call getdata(10,dataline)
      read(dataline,*)accuracy
      if(accuracy.le.0.d0.or.accuracy.ge.1.d0)accuracy=0.1d0
c
c     parameters for output files
c     ===========================
c
      call getdata(10,dataline)
      read(dataline,*)outdir
c
      do lend=80,1,-1
        if(outdir(lend:lend).ne.' ')goto 100
      enddo
100   continue
c
      if(lend.lt.1)then
        stop 'Error: wrong format for output directory!'
      endif
c
      call getdata(10,dataline)
      read(dataline,*)(fname(i),i=1,3)
      call getdata(10,dataline)
      read(dataline,*)(fname(i),i=4,9)
      call getdata(10,dataline)
      read(dataline,*)(fname(i),i=10,12)
      call getdata(10,dataline)
      read(dataline,*)(fname(i),i=13,16)
      do i=1,16
        do lenf=80,1,-1
          if(fname(i)(lenf:lenf).ne.' ')goto 110
        enddo
110     continue
        green(i,1)=outdir(1:lend)//fname(i)(1:lenf)//'.ep'
        green(i,2)=outdir(1:lend)//fname(i)(1:lenf)//'.ss'
        green(i,3)=outdir(1:lend)//fname(i)(1:lenf)//'.ds'
        green(i,4)=outdir(1:lend)//fname(i)(1:lenf)//'.cl'
        do istp=1,4
          select(i,istp)=.true.
        enddo
      enddo
c
c     no tangential components for clvd sources
c
      select(3,1)=.false.
      select(8,1)=.false.
      select(9,1)=.false.
      select(11,1)=.false.
      select(12,1)=.false.
      select(16,1)=.false.
      select(3,4)=.false.
      select(8,4)=.false.
      select(9,4)=.false.
      select(11,4)=.false.
      select(12,4)=.false.
      select(16,4)=.false.
c
c     global model parameters
c     =======================
c
      call getdata(10,dataline)
      read(dataline,*)l,isurfcon
      if(l.gt.lmax)then
        stop ' Max. no of layers (lmax) too small defined!'
      else if(isurfcon.lt.1.or.isurfcon.gt.3)then
        stop ' Wrong selection for surface condition!'
      
      else if(isurfcon.eq.3)then
        read(dataline,*)l,isurfcon,porosity
        if(porosity.le.0.d0)then
          stop ' Wrong value for near-surface porosity!'
        endif
      endif
c
c      multilayered model parameters
c      =============================
c
      do i=1,l
        call getdata(10,dataline)
        read(dataline,*)j,h(i),vp,vs,rho,skempton,alf(i),diffus
        if(skempton.gt.1.d0.or.skempton.lt.0.d0)then
          stop 'Error in pegmain: wrong Skempton ratio!'
        endif
        if(alf(i).gt.1.d0.or.alf(i).lt.0.d0)then
          stop 'Error in pegmain: wrong eff. stress coefficient!'
        endif
        h(i)=h(i)*km2m
        vp=vp*km2m
        vs=vs*km2m
        mu(i)=rho*vs*vs
        lau=rho*vp*vp-2.d0*mu(i)
        if(lau.le.0.d0)then
          stop 'inconsistent Vp/Vs ratio!'
        endif
        if(skempton*alf(i).ge.1.d0.or.skempton*alf(i).le.0.d0)then
          stop 'inconsistent value for alpha*skempton!'
        endif
        if(diffus.le.0.d0)then
          stop 'inconsistent value for diffusivity!'
        endif
c
c       nu = drained Poisson ratio
c       nuu = undrained Poisson ratio
c       qa = alpha*Q (1/Q = Biot's compressibility)
c       dm = Q*xi (xi = Darcy conductivity)
c
        nuu=0.5d0*lau/(lau+mu(i))
        gamma=(1.d0-alf(i)*skempton)*(1.d0+nuu)/(1.d0-2.d0*nuu)
        nu=(gamma-1.d0)/(1.d0+2.d0*gamma)
        la(i)=mu(i)*2.d0*nu/(1.d0-2.d0*nu)
        qa(i)=(1.d0+nuu)*mu(i)*skempton/(1.5d0-3.d0*nuu)
        dm(i)=diffus*(1.d0-nuu)*(0.5d0-nu)/((1.d0-nu)*(0.5d0-nuu))
      enddo
      if(l.eq.1)h(l)=0.d0
c
c     end of inputs
c     =============
c
      close(10)
c
      comptxt(1)='Uz (vertical displacement)'
      comptxt(2)='Ur (radial displacement)'
      comptxt(3)='Ut (tangential displacement)'
      comptxt(4)='Szz (linear stress)'
      comptxt(5)='Srr (linear stress)'
      comptxt(6)='Stt (linear stress)'
      comptxt(7)='Szr (shear stress)'
      comptxt(8)='Srt (shear stress)'
      comptxt(9)='Stz (shear stress)'
      comptxt(10)='Tr (tilt -dUr/dz)'
      comptxt(11)='Tt (tilt -dUt/dz)'
      comptxt(12)='Rot (rotation ar. z-axis)'
	comptxt(13)='P (pore pressure)'
      comptxt(14)='Dcz (vertical Darcy flux)'
      comptxt(15)='Dcr (radial Darcy flux)'
      comptxt(16)='Dct (tangential Darcy flux)'
c
c     determine upper und lower parameter values of each layer
c
      l0=1
      z1(l0)=0.d0
      do i=2,l
        if(h(i).gt.h(i-1))then
          z1(l0)=h(i-1)
          la1(l0)=la(i-1)
          mu1(l0)=mu(i-1)
          alf1(l0)=alf(i-1)
          qa1(l0)=qa(i-1)
          dm1(l0)=dm(i-1)
c
          z2(l0)=h(i)
          la2(l0)=la(i)
          mu2(l0)=mu(i)
          alf2(l0)=alf(i)
          qa2(l0)=qa(i)
          dm2(l0)=dm(i)
          l0=l0+1
        else
          z1(l0)=h(i)
          la1(l0)=la(i)
          mu1(l0)=mu(i)
          alf1(l0)=alf(i)
          qa1(l0)=qa(i)
          dm1(l0)=dm(i)
        endif
      enddo
      z1(l0)=h(l)
      la1(l0)=la(l)
      mu1(l0)=mu(l)
      alf1(l0)=alf(l)
      qa1(l0)=qa(l)
      dm1(l0)=dm(l)
c
c     construction of sublayers
c
      write(*,*)'the multi-layered poroelastic model:'
c
      call pegsublay(ierr)
      if(ierr.eq.1)then
        stop 'the max. no of layers (lmax) too small defined!'
      endif
c
      zs=0.d0
      call peglayer(ierr)
      nlr=nno(lzrec)
c
      if((isurfcon.eq.1.or.isurfcon.eq.3).and.nlr.eq.1)then
        nearsurf=1
      else
        nearsurf=0
      endif
c
      leninp=index(inputfile,' ')-1
c
      stype(1)='explosion (M11=M22=M33=1*kappa)'
      stype(2)='strike-slip (M12=M21=1*mue)'
      stype(3)='dip-slip (M13=M31=1*mue)'
      stype(4)='clvd (M33=1*mue, M11=M22=-M33/2)'
c
      iunit=10
      do istp=1,4
        do i=1,16
          if(select(i,istp))then
            iunit=iunit+1
            unit(i,istp)=iunit
            open(unit(i,istp),file=green(i,istp),status='unknown')
            write(unit(i,istp),'(a)')'################################'
            write(unit(i,istp),'(a)')'# The input file used: '
     &                        //inputfile(1:leninp)
            write(unit(i,istp),'(a)')'################################'
            write(unit(i,istp),'(a)')'# Greens function component: '
     &                        //comptxt(i)
            write(unit(i,istp),'(a)')'#(Okada solutions subtracted)'
            write(unit(i,istp),'(a)')'# Source type: '//stype(istp)
            write(unit(i,istp),'(a)')'# Selections of surface condition'
     &                             //' and near surface correction:'
            write(unit(i,istp),'(2i10)')isurfcon,nearsurf
            write(unit(i,istp),'(a)')'# Observation distance sampling:'
            write(unit(i,istp),'(a)')'#    nr         r1[m]'
     &                             //'         r2[m] smp_ratio'
            write(unit(i,istp),'(i7,2d14.6,f10.4)')nr,r1,r2,sampratio
            write(unit(i,istp),'(a)')'# Uniform obs. site parameters:'
            write(unit(i,istp),'(a)')'#    depth[m]       la[Pa]       '
     &          //'mu[Pa]        alpha       qa[Pa]    dm[m^2/s]'
            write(unit(i,istp),'(6d14.6)')zrec,la(nlr),
     &           mu(nlr),alf(nlr),qa(nlr),dm(nlr)
            write(unit(i,istp),'(a)')'# Source depth sampling:'
            write(unit(i,istp),'(a)')'#   nzs       zs1[m]       zs2[m]'
            write(unit(i,istp),'(i7,2d14.6)')nzs,zs1,zs2
            write(unit(i,istp),'(a)')'# Time sampling:'
            write(unit(i,istp),'(a)')'#    nt        t-window[s]'
            write(unit(i,istp),'(i7,d14.6)')nt,twindow
            write(unit(i,istp),'(a)')'# Data in each source depth block'
            write(unit(i,istp),'(a)')'# ==============================='
            write(unit(i,istp),'(a)')'# Line 1: source layer parameters'
            write(unit(i,istp),'(a)')'#  s_depth, la, mu, alpha, qa, dm'
            write(unit(i,istp),'(a)')'# Line 2: coseismic responses '
     &                             //'(f(ir,it=1),ir=1,nr)'
            write(unit(i,istp),'(a)')'# Line 3: (idec(ir),ir=1,nr)'
     &                //'(decimal exponents for postseismic responses)'
            write(unit(i,istp),'(a)')'# Line 4: (f(ir,it=2),ir=1,nr)'
            write(unit(i,istp),'(a)')'# Line 5: (f(ir,it=3),ir=1,nr)'
            write(unit(i,istp),'(a)')'#  ...'
          endif
        enddo
      enddo
c
      call pegbsj(ierr)
      do izs=1,nzs
        zs=zs1+dble(izs-1)*dzs
        write(*,'(/,a,E13.4,a)')' Processing for the '
     &                 //'source at depth:',zs,' m.'
c
        call peglayer(ierr)
        nls=nno(ls)
        nlr=nno(lzrec)
        call pegsfcths(1.d0)
        call pegsfctud(1.d0)
c
        do istp=1,4
          do i=1,16
            do ir=1,nr
              do it=1,nfmax
                du(it,ir,i,istp)=(0.d0,0.d0)
              enddo
            enddo
          enddo
        enddo
c
        zrs2=(zrec-zs)**2
        rsmin=dmin1(0.5d0*dr,0.1d0*dsqrt(zrs2+r(nr)**2))
        do ir=1,nr
          rs(ir)=dmax1(rsmin,0.1d0*dsqrt(zrs2+r(ir)**2))
          geow(ir)=zrs2+(rs(ir)+r(ir))**2
        enddo
c
        swap=dsqrt(zrs2+(rs(nr)+r(nr))**2)/dsqrt(zrs2+(rs(1)+r(1))**2)
        nprf=1+idnint(dlog(swap)/dlog(2.5d0))
        if(nprf.gt.1)then
          beta=swap**(1.d0/dble(nprf-1))
        else
          beta=2.5d0
        endif
c
        isp=0
        nr2=0
200     isp=isp+1
        nr1=nr2+1
        nr2=nr1
        do ir=nr1+1,nr
          if(r(ir).le.beta*dsqrt(zrs2+(rs(nr1)+r(nr1))**2))nr2=ir
        enddo
        call pegspec(isp,nr1,nr2,nt,dt,accuracy)
        if(nr2.lt.nr)goto 200
c
        do istp=1,4
          do i=1,16
            if(.not.select(i,istp))goto 400
            write(unit(i,istp),'(a)')'#################################'
            write(unit(i,istp),'(a,i2,a)')'# the ',izs,'. source depth:'
            write(unit(i,istp),'(a)')'#################################'
            write(unit(i,istp),'(6E14.6)')zs,la(nls),
     &        mu(nls),alf(nls),qa(nls),dm(nls)
            do ir=1,nr-1
              write(unit(i,istp),'(E14.6,$)')dreal(du(1,ir,i,istp))
            enddo
            write(unit(i,istp),'(E14.6)')dreal(du(1,nr,i,istp))
            do ir=1,nr
              du(1,ir,i,istp)=dcmplx(0.d0,dimag(du(1,ir,i,istp)))
              am=0.d0
              do it=1,(nt+1)/2
                am=dmax1(am,dabs(dreal(du(it,ir,i,istp))),
     &                      dabs(dimag(du(it,ir,i,istp))))
              enddo
              if(am.le.0.d0)then
                idec(ir)=0
              else
                idec(ir)=idint(dlog10(am))-4
                do it=1,(nt+1)/2
                  du(it,ir,i,istp)=du(it,ir,i,istp)
     &                  *dcmplx(10.d0**(-idec(ir)),0.d0)
                enddo
              endif
            enddo
            call outint(unit(i,istp),idec,nr)
            do it=1,(nt+1)/2
              if(it.gt.1)then
                do ir=1,nr
                  nout(ir)=idnint(dreal(du(it,ir,i,istp)))
                enddo
                call outint(unit(i,istp),nout,nr)
              endif
              if(it*2.le.nt)then
                do ir=1,nr
                  nout(ir)=idnint(dimag(du(it,ir,i,istp)))
                enddo
                call outint(unit(i,istp),nout,nr)
              endif
            enddo
400         continue
          enddo
        enddo
      enddo
c
c     end of izs loop
c
      do istp=1,4
        do i=1,16
          if(select(i,istp))close(unit(i,istp))
        enddo
      enddo
      if(nwarn.eq.0)then
        print *,'####################################################'
        print *,'#                                                  #'
        print *,'#          End of computations with PEGRN          #'
        print *,'#                                                  #'
        print *,'####################################################'
      else
        print *,'####################################################'
        print *,'     Sorry, there have been',nwarn,' warnings.      '
        print *,'             Results may be inaccurate!             '
        print *,'####################################################'
      endif
      stop
      end
