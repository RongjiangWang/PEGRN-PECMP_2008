c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GLOBAL CONSTANTS
c     ================
c
c     NRECMAX = max. number of observation positions
c     NZSMAX = max. number of discrete source depths
c     NRMAX = max. number of discrete radial diatances
c     NSMAX = max. number of fault segments
c     NEQMAX = max. number of earthquakes
c     NPTCHMAX = max. number of patches at a source rectangle
c     NPSMAX = max. number of discrete point sources per source depth
c     NTMAX = max. number of time samples used for Green's functions
c     NTRMAX = max. number of time samples of the outputs
c     NSCMAX = max. number of scenario outputs (<= NTRMAX/2)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer NZSMAX,NRMAX,NEQMAX,NSMAX,NPTCHMAX
	  integer NPSMAX,NRECMAX,NTMAX,NTRMAX,NSCMAX
      parameter(NZSMAX=100,NRMAX=151)
      parameter(NEQMAX=20,NSMAX=1000,NPTCHMAX=1000)
	  parameter(NPSMAX=50000)
      parameter(NRECMAX=1000)
      parameter(NTMAX=4096,NTRMAX=4096)
      parameter(NSCMAX=NTRMAX/2)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     RECTANGULAR SOURCE PLANES
c     =========================
c
c     (xref,yref) = gegraphic coordinates of the reference point
c     zref = depth of the reference point.
c     all angles in degree.
c     NSMAX = the max. number of source rectangles
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer nptch_s(NSMAX),nptch_d(NSMAX),ieqno(NSMAX)
      double precision xref(NSMAX),yref(NSMAX),zref(NSMAX)
      double precision length(NSMAX),width(NSMAX)
      double precision strike(NSMAX),dip(NSMAX)
      double precision tstart(NSMAX),eqtime(NEQMAX)
      double precision ptch_s(NSMAX,NPTCHMAX),ptch_d(NSMAX,NPTCHMAX)
      double precision slip_s(NSMAX,NPTCHMAX)
	  double precision slip_d(NSMAX,NPTCHMAX)
      double precision opening(NSMAX,NPTCHMAX)
c
      common/irects/nptch_s,nptch_d,ieqno
      common/drects/xref,yref,zref,length,width,
     &              strike,dip,tstart,
     &              eqtime,ptch_s,ptch_d,slip_s,slip_d,opening
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     DISTRETE POINT SOURCES
c     ======================
c
c     (xs,ys,zs) = coordinates of the discrete point sources
c     with x = north, y = east, z = downward
c     angles in degree.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer nps(NZSMAX),isno(NPSMAX,NZSMAX)
      double precision px(NPSMAX,NZSMAX),py(NPSMAX,NZSMAX)
      double precision pz(NPSMAX,NZSMAX),pmwei(6,NPSMAX,NZSMAX)
c
      common/ipoints/nps,isno
      common/dpoints/px,py,pz,pmwei
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GREEN'S FUNCTION INFO
c     =====================
c
c     nzs,zs1,zs2 = number of depth samples, start and end depths used
c           in Green's functions
c     nr,r1,r2 = number of distance samples, start and end distances used
c           in Green's functions
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer isurfcon,nearsurf,nr,nzs,nt
      double precision r1,r2,sampratio,zs1,zs2
      double precision zrec,larec,murec,rhorec,alfrec,qarec,dmrec,
     +                 twindow
      common /igreeninfo/isurfcon,nearsurf,nr,nzs,nt
      common /greeninfo/r1,r2,sampratio,zs1,zs2,zrec,larec,murec,rhorec,
     +                  alfrec,qarec,dmrec,twindow
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     OBSERVATION POSITIONS AND OBSERVABLES
c     =====================================
c
c     (xrec(i),yrec(i))=coordinates of the observation positions
c     the 3 displcement/velocity/acceleration components: ux,uy,uz
c     NRECMAX = the max. number of observation positions
c     latlon: 1/0 = geographic/Cartesian
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer ntrec,itsc(2*NSCMAX),itout(16)
      double precision xrec(NRECMAX),yrec(NRECMAX)
      double precision tsc(NSCMAX),obs(NTRMAX,NRECMAX,16)
      double precision coobs(NEQMAX,NRECMAX,16)
      double precision poobs(NEQMAX,NRECMAX,16)
      double precision obs1(NRECMAX,16),obs2(NRECMAX,16)
      double precision dc0(NEQMAX,NRECMAX)
      character*80 grndir,green(16)
      character*80 outdir,toutfile(16),scoutfile(NSCMAX)
c
      common/iobsarray/ntrec,itsc,itout
      common/dobsarray/xrec,yrec,tsc,
     &                 obs,coobs,poobs,obs1,obs2,dc0
      common/cobsarray/grndir,green,outdir,toutfile,scoutfile
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     COSINES OF LOS TO INSAR ORBIT
c     =============================
c
c     insar = 1: output los displacements
c             0: not output los displacements
c     xlos, ylos, zlos = cosines of the los
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer insar,icmb
      double precision xlos,ylos,zlos
      double precision friction,strike0,dip0,rake0
      double precision sigma0(3)
c
      common/iothers/insar,icmb
      common/dothers/xlos,ylos,zlos,
     &      friction,strike0,dip0,rake0,sigma0
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     WARNING STATISTICS
c     ==================
c
c     nwarn = total number of warnings
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer nwarn
c
      common/warnings/nwarn
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL CONSTANTS
c     ==============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      double precision DEG2RAD,KM2M,DAY2SEC,REARTH,G0
      parameter(DEG2RAD=1.745329252d-02,KM2M=1.0d+03)
      parameter(DAY2SEC=8.64d+04,REARTH=6.371d+06,G0=9.82d+00)
