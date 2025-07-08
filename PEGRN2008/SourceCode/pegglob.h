c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     nzmax: max. interface index;
c     lmax: max. no of total homogeneous layers (lmax <= nzmax-2);
c     nrmax: max. no of traces;
c     nfmax: max. no of frequency samples.
c     nzsmax: max. number of source depths
c
      integer lmax,nzmax,nrmax,nfmin,nfmax,ntmax,nzsmax
      parameter(lmax=100)
      parameter(nzmax=lmax+3)
      parameter(nrmax=151)
      parameter(nfmin=64)
	  parameter(nfmax=1024)
      parameter(ntmax=2*nfmax)
      parameter(nzsmax=50)
c
c     INDEX PARAMETERS FOR BESSEL FUNCTION TABLES
c     ===========================================
c
      integer dnx,nxmax,nbsjmax
      parameter(dnx=512)
      parameter(nxmax=1024)
      parameter(nbsjmax=dnx*nxmax)
c
c     GLOBAL CONSTANTS
c     ================
c
      double precision km2m,day2sec
      parameter(km2m=1.0d+03,day2sec=8.64d+04)
      double precision rhog
      parameter(rhog=9.8d+03)
c
      integer nwarn
      common /iwarning/ nwarn
c
c     DISCRETISATION ACCURACY FOR LAYERS WITH CONSTANT GRADIENT
c     =========================================================
c     reslm: for moduli
c     resld: for diffusivity
c
      double precision reslm,resld
      parameter(reslm=0.05d0,resld=0.05d0)
c
c     COMMON BLOCKS
c     =============
      integer lp,nls,nlr,nno(nzmax)
      double precision hp(nzmax)
      common /isublayer/ lp,nno
      common /dsublayer/ hp,nls,nlr
c
c     zrec: receiver depth
c     lzrec: sublayer no of receiver
c
c     isurfcon = 1: unconfined surface (p = 0)
c                2: confined surface (dp/dz = 0)
c                3: unsaturated surface (dp/dt = rho*g*dcf)
c
      integer isurfcon,nearsurf,lzrec
      double precision zrec,porosity
      common /ireceiver/ isurfcon,nearsurf,lzrec
      common /dreceiver/ zrec,porosity
c
c     original model parameters
c
      integer l0
      double precision z1(lmax),z2(lmax)
      double precision la1(lmax),la2(lmax),mu1(lmax),mu2(lmax)
      double precision alf1(lmax),alf2(lmax),qa1(lmax),qa2(lmax)
      double precision dm1(lmax),dm2(lmax)
      common /imodel0/ l0
      common /dmodel0/ z1,z2,la1,la2,mu1,mu2,alf1,alf2,qa1,qa2,
     &                 dm1,dm2
c       
c     model parameter:
c     n0: number of homogeneous layers
c
      integer n0
      double precision h(lmax),la(lmax),mu(lmax),alf(lmax)
      double precision qa(lmax),dm(lmax)
      logical smalls(lmax)
      common /lmodel/ smalls
      common /imodel/ n0
      common /dmodel/ h,la,mu,alf,qa,dm
c
      double complex cla(lmax),cmu(lmax),ka(lmax)
      common /dcpara/ cla,cmu,ka
c
c     source parameters
c
      integer ls,ms(4),ics(4)
      double precision zs
      double complex sfct0(8,4),sfct1(8,4)
      double complex sfctud0(8,4),sfctud1(8,4)
      common /isource/ ls,ms,ics
      common /dsource/ zs,sfct0,sfct1,sfctud0,sfctud1
c
c     half-space source parameters
c
      double complex sfcths0(8,4),sfcths1(8,4)
      common /dsourcehs/ sfcths0,sfcths1
c
c     table of J_n(x), n = -1, 0, 1, 2, 3
c
      double precision dxbsj,bsj(0:nbsjmax,-1:3)
      common /dbessels/ dxbsj,bsj
c
c     output data
c
      double precision r(nrmax),rs(nrmax),geow(nrmax)
      double complex du(-1:nfmax,nrmax,16,4)
      logical select(16,4)
      common /doutdata/ r,rs,geow,du
      common /loutdata/ select