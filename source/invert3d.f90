!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! USE in any subroutine or function or other module.
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER :: pvi,pii,psi,nvpi,nipi,nspi,npi
INTEGER :: subdim,asds,invstep
INTEGER :: nvgi,nigi,nnfd,ntr
INTEGER :: ninp,nint
INTEGER, DIMENSION(:), ALLOCATABLE :: idvg,idvt,idig,ids
INTEGER, DIMENSION(:), ALLOCATABLE :: fcoln,cnfe,tfcoln,tcnfe
INTEGER, DIMENSION(:,:), ALLOCATABLE :: nvnr,nvnt,nvnp
REAL(KIND=i5), DIMENSION(:), ALLOCATABLE :: mo,mc,cm,ecmi,dm
REAL(KIND=i5), DIMENSION(:), ALLOCATABLE :: frech,tfrech
REAL(KIND=i5), DIMENSION(:), ALLOCATABLE :: dobs,dmod,cd
REAL(KIND=i5), DIMENSION(:), ALLOCATABLE :: dem,smv
REAL(KIND=i5) :: epsilon,eta
REAL(KIND=i5) :: etav,etai
!
! pvi = Perform velocity inversion (0=no, 1=yes)
! pii = Perform interface inversion (0=no, 1=yes)
! psi = Perform source inversion (0=no, 1=yes)
! epsv = Damping factor (epsilon) for velocity
! epsi = Damping factor (epsilon) for interfaces
! epss1 = Damping factor (epsilon) for source displacement
! epss2 = Damping factor (epsilon) for source time
! etav = Smoothing factor (eta) for velocity
! etai = Smoothing factor for interfaces
! subdim = Subspace dimension used in inversion
! asds = Apply second derivative smoothing? (0=no, 1=yes)
! epsilon = Global damping factor
! eta = Global smoothing factor
! invstep = Inversion step
! nvgi = Number of velocity grids for inversion
! nigi = Number of interface grids for inversion
! nspi = Number of source points for inversion
! idvg = Ids of velocity grid for inversion
! idvt = Associated velocity types (1 or 2) for inversion
! idig = Ids of interface grids for inversion
! ids = Ids of sources for inversion
! npi = Total number of parameters for inversion
! nvpi,nipi,nspi= No. of vel,int,source parameters for inversion
! mo = Reference model parameters
! mc = Current model parameters
! cm = Diagonal elements of a priori model covariance matrix
! ecmi = Inverse of cm with parameter pre-weighting
! dm = Model perturbation
! ntr = Number of traveltimes
! nnfd = Number of non-zero frechet derivatives
! frech = Array of frechet derivatives
! fcoln = Points to column number of frechet matrix (G)
! cnfe = Cumulative number of no-zero elements of G for row i
! dobs = Vector of observed data (traveltimes)
! dmod = Vector of model data (traveltimes)
! cd = a priori data covariance matrix (diagonal)
! tfrech,tcnfe,tfcoln = Same as frech,cnfe,fcoln but for transpose
! dem = Finite difference estimate of spatial derivative
! smv = Vector to which smoothing is applies
! nvnr,nvnt,nvnp = Number of velocity nodes in r,lat,long
! nint,ninp = Number of interface nodes in lat, long.
!
END MODULE globalp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program will perform a single iteration of an
! n-dimensional subspace inversion method using traveltime
! data from the 3-D fast marching code fm3d.
! LU decomposition is used to perform the matrix
! inversion and SVD is used to construct an orthonormal basis
! for the projection matrix A.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM invert
USE globalp
use lsmrModule, only:lsmr
IMPLICIT NONE
INTEGER :: i,j,k,l,m
INTEGER ::  minsd,comfd,ni,nvgt,vgid,mvnr,mvnp,mvnt,cstp
INTEGER :: ns,idm1,idm2,idm3,idm4,istep,jstep,kstep,iswt
INTEGER :: ntels,nrow,jup,rmtr,isw,isw2,istepo,isum,idm2o
INTEGER :: npgnr
INTEGER, DIMENSION(:), ALLOCATABLE :: tsid,istel,stpv
INTEGER, DIMENSION(:,:), ALLOCATABLE :: nsdf
INTEGER, DIMENSION(:), ALLOCATABLE :: lots,nps
INTEGER, PARAMETER :: maxsd=50,maxs=60,maxp=40
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: paths,patht
REAL(KIND=i10), PARAMETER :: pi=3.141592653589793
REAL(KIND=i5) :: epsv,epsi,epss1,epss2,rd1,dref
REAL(KIND=i10) :: gnsit,gnsip,goit,goip,earthr,mdbi,pgt,pgb,mpv
REAL(KIND=i5), DIMENSION(:,:,:,:,:), ALLOCATABLE :: veln,velnr,covvm
REAL(KIND=i5), DIMENSION(:,:,:), ALLOCATABLE :: intn,intnr,covim
REAL(KIND=i5), DIMENSION(:), ALLOCATABLE :: srad,slat,slon,stp
REAL(KIND=i5), DIMENSION(:), ALLOCATABLE :: sradr,slatr,slonr
REAL(KIND=i5), DIMENSION(:), ALLOCATABLE :: covsr,covst,covsp
REAL(KIND=i10), DIMENSION(:), ALLOCATABLE :: mtmean
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: gnsr,gnst,gnsp
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: gor,got,gop
CHARACTER (LEN=30) :: refvgfile,vgfile,refigfile,igfile
CHARACTER (LEN=30) :: refscfile,scfile,pgfile
CHARACTER (LEN=30) :: otfile,mtfile,rtfile,stfile
CHARACTER (LEN=30) :: frpfile,frdatfile,invifile
CHARACTER (LEN=10), DIMENSION(:), ALLOCATABLE :: tpath

! hongjian fang @ethz... adding surface wave data
!----------------------------------------------------------------
CHARACTER (LEN=30) :: frdatfilesurf
CHARACTER (LEN=30) :: otfilesurf,mtfilesurf
integer :: ntrsurf
integer :: surfjoint
REAL(KIND=i5) :: mean,std_surf

! variable definition for lsmr
integer leniw,lenrw
real(kind=i5) atol,btol
real(kind=i5) conlim
integer istop
integer itnlim
integer nout
integer itn
real(kind=i5) damp,acond,anorm,arnorm,rnorm,xnorm
integer localSize
integer,allocatable,dimension(:)::iw,col
real(kind=i5),allocatable,dimension(:)::rw
real(kind=i5),allocatable,dimension(:)::dtrav
real(kind=i5),allocatable,dimension(:)::dataweight
real,allocatable,dimension(:)::norm
!real,allocatable,dimension(:)::norm,norm_dws
!REAL(KIND=i5), DIMENSION(:,:,:,:,:), ALLOCATABLE :: velndws
integer checkstat
integer inversionScheme
integer is1,is2,is3
character*4 itnum
real(kind=i5) threshold
real(kind=i5) surfweight
integer nnode
integer vpvs
integer cnt,tmp
real(kind=i5) tmp1
integer jstep_tmp

character (len=40) cdum
!integer flex,choosedsrc
real(kind=i5) radall,latall,lonall
!----------------------------------------------------------------
!
! refvgfile = Reference velocity grid file
! vgfile = Current velocity grid file
! refigfile = Reference interface grid file
! igfile = Current interface grid file
! refscfile = Reference source coordinate file
! scfile = Current source coordinate file
! stfile = Source time perturbation file
! otfile = Observed traveltime file
! mtfile = Model traveltime file
! rtfile = Reference teleseismic traveltime file
! frpfile = Frechet parameter file
! frdatfile = Frechet data file
! invifile = Inversion iteration file
! minsd = Minimum permissable subspace dimension
! comfd = Have frechet derivatives been computed? (0=no,1=yes)
! idvg = Ids of velocity grid for inversion
! idvt = Associated velocity types (1 or 2) for inversion
! idig = Ids of interface grids for inversion
! ids = Ids of sources for inversion
! ni = Number of interfaces describing model
! nvgt = Number of velocity grid types (1 or 2)
! vgid = Velocity grids differ? (0=no,1=yes)
! veln = Value of velocity parameters
! velnr = Value of reference velocity parameters
! covvm = a priori model covariance matrix for velocity
! mvnr,mvnt,mvnp = Maximum dimensions of velocity grid
! intn = Value of interface parameters
! intnr = Value of reference interface parameters
! covim = a priori model covariance matrix for interface
! ns = Number of sources
! srad,slat,slon = Current source radius,latitude,longitude
! sradr,slatr,slonr = Reference source radius,latitude,longitude
! covsr,covst,covsp = a priori covariance of source locations
! stp = Source time perturbation
! epsv = Damping factor (epsilon) for velocity
! epsi = Damping factor (epsilon) for interfaces
! epss1 = Damping factor (epsilon) for source displacement
! epss2 = Damping factor (epsilon) for source time
! paths = Path sequence information
! ntels = Number of teleseismic sources
! tsid = Teleseismic source i.d.
! istel = Is traveltime teleseismic? (0=no, >0=source id.)
! nrow = Number of non-zero Frechet derivatives in row i
! istep,jstep,kstep = counting indices
! jup = jstep+cnfe(istep)-1
! dref = Reference teleseismic traveltimes
! stpv,cstp = variables for determining transpose of Frechet matrix
! gnsr,gnst,gnsp = Velocity grid node spacing in r,theta,phi
! gor,got,gop = Velocity grid origin in r,theta,phi
! gnsir,gnsit,gnsip = Interface grid node spacing in r,theta,phi
! goir,goit,goip = Interface grid origin in r,theta,phi
! maxs = Maximum number of path segments in ray
! maxp = Maximum number of path types for a source
! lots = local (0) or teleseismic (1) source
! nps = Number of path types from source
! tpath = Class of teleseismic path
! nsdf = Number of segments defining particular path 
! paths = Path sequence
! patht = Velocity type associated with path segment
! rmtr = Remove mean from predicted teleseismic residual (0=no, 1=yes)
! iswt = Switch for teleseismic arrival times
! earthr = Earth radius in km
! mtmean = Mean of model teleseismic traveltime residual
! pgfile = File containing propagation grid parameters
! mdbi = Minimum distance between interfaces (km)
! mpv = Minimum permitted velocity (km/s)
! pgt = Top of propagation grid
! pgb = Bottom of propagation grid
! npgnr = Number of propagation grid nodes in radius
!
! Read in the input parameters
!
OPEN(UNIT=10,FILE='invert3d.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)refvgfile
READ(10,1)vgfile
READ(10,1)refigfile
READ(10,1)igfile
READ(10,1)refscfile
READ(10,1)scfile
READ(10,1)stfile
READ(10,*)
READ(10,1)otfile
READ(10,1)mtfile
READ(10,1)rtfile
READ(10,1)frpfile
READ(10,1)frdatfile
READ(10,1)invifile
READ(10,1)pgfile
READ(10,*)mdbi
READ(10,*)mpv
READ(10,*)rmtr
READ(10,*)pvi,epsv,etav
READ(10,*)pii,epsi,etai
READ(10,*)psi,epss1,epss2
READ(10,*)subdim
READ(10,*)epsilon
READ(10,*)asds
READ(10,*)eta
READ(10,*)earthr
read(10,*)surfjoint
read(10,*)inversionScheme
read(10,*)threshold
read(10,*)surfweight
read(10,*)vpvs
!read(10,*)flex
1 FORMAT(a26)
CLOSE(10)
!
! Check legitimacy of requested subspace dimension
!
minsd=pvi+pii+2*psi
IF(minsd.LE.0)THEN
   WRITE(6,*)'According to the values in the inversion parameter'
   WRITE(6,*)'file, you do not want to invert for anything!!!'
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
IF(subdim.LT.minsd)THEN
   WRITE(6,*)'WARNING'
   WRITE(6,*)'For the selected parameter classes for inversion,'
   WRITE(6,*)'the minimum subspace dimension is',minsd
   WRITE(6,*)'Resetting subspace dimension to this value and'
   WRITE(6,*)'continuing.'
   subdim=minsd
ENDIF
IF(subdim.LE.0)THEN
   WRITE(6,*)'Requested subspace dimension is less than or'
   WRITE(6,*)'equal to zero!!!'
   WRITE(6,*)'Subspace dimension must be greater than zero!'
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
IF(subdim.GT.maxsd)THEN
   WRITE(6,*)'Requested subspace dimension is greater than the'
   WRITE(6,*)'maximum permitted, which is',maxsd
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
!
! Determine the iteration number
!
OPEN(UNIT=10,FILE=invifile,STATUS='old')
READ(10,*)invstep
CLOSE(10)
!
! Read in Frechet parameter file to see what partial derivatives
! are available.
!
OPEN(UNIT=10,FILE=frpfile,STATUS='old')
READ(10,*)comfd
IF(comfd.EQ.0)THEN
   WRITE(6,*)'Frechet derivatives were not computed during'
   WRITE(6,*)'the forward step with fm3d!!!!!!'
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
READ(10,*)nvgi
IF(nvgi.GT.0)THEN
   ALLOCATE(idvg(nvgi))
   ALLOCATE(idvt(nvgi))
   READ(10,*)idvg(1:nvgi)
   READ(10,*)idvt(1:nvgi)
ENDIF
READ(10,*)nigi
IF(nigi.GT.0)THEN
   ALLOCATE(idig(nigi))
   READ(10,*)idig(1:nigi)
ENDIF
READ(10,*)nspi
IF(nspi.GT.0)THEN
   ALLOCATE(ids(nspi))
   READ(10,*)ids(1:nspi)
ENDIF
CLOSE(10)
!
! Test that Frechet derivatives correlate with requested
! parameters for inversion.
!
IF(pvi.EQ.1.AND.nvgi.LE.0)THEN
   WRITE(6,*)'You want to invert for velocity, but no such'
   WRITE(6,*)'Frechet derivatives were computed by fm3d!!!'
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
IF(pii.EQ.1.AND.nigi.LE.0)THEN
   WRITE(6,*)'You want to invert for interface depth, but no such'
   WRITE(6,*)'Frechet derivatives were computed by fm3d!!!'
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
IF(psi.EQ.1.AND.nspi.LE.0)THEN
   WRITE(6,*)'You want to invert for interface depth, but no such'
   WRITE(6,*)'Frechet derivatives were computed by fm3d!!!'
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
!
! Start off by reading in complete velocity grid if velocity
! parameters are to be inverted for,
! 
mvnp=0
mvnr=0
mvnt=0
ntels=0
IF(pvi.EQ.1)THEN
!
!  Start with the reference velocity grid
!
   OPEN(UNIT=10,FILE=refvgfile,STATUS='old')
   READ(10,*)ni,nvgt
   ALLOCATE(nvnr(ni,nvgt),nvnt(ni,nvgt),nvnp(ni,nvgt))
   ALLOCATE(gnsr(ni,nvgt),gnst(ni,nvgt),gnsp(ni,nvgt))
   ALLOCATE(gor(ni,nvgt),got(ni,nvgt),gop(ni,nvgt))
   ni=ni+1
   vgid=0
   DO i=1,nvgt
!     IF(i.eq.2)DEALLOCATE(velnr,covvm)
      DO j=1,ni-1
         READ(10,*)nvnr(j,i),nvnt(j,i),nvnp(j,i)
         READ(10,*)gnsr(j,i),gnst(j,i),gnsp(j,i)
         READ(10,*)gor(j,i),got(j,i),gop(j,i)
         IF(j.GT.1)THEN
            IF(nvnr(j,i).NE.nvnr(j-1,i))vgid=1
            IF(nvnt(j,i).NE.nvnt(j-1,i))vgid=1
            IF(nvnp(j,i).NE.nvnp(j-1,i))vgid=1
         ENDIF
         IF(i.EQ.1.AND.j.EQ.1)THEN
            ALLOCATE(velnr(nvnr(j,i),nvnt(j,i),nvnp(j,i),ni-1,nvgt))
            ALLOCATE(covvm(nvnr(j,i),nvnt(j,i),nvnp(j,i),ni-1,nvgt))
            mvnr=nvnr(j,i)
            mvnt=nvnt(j,i)
            mvnp=nvnp(j,i)
         ENDIF
         IF(vgid.EQ.0)THEN
            DO k=1,nvnr(j,i)
               DO l=1,nvnt(j,i)
                  DO m=1,nvnp(j,i)
                     READ(10,*)velnr(k,l,m,j,i),covvm(k,l,m,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            !DO k=1,nvnr(j,i)*nvnt(j,i),nvnp(j,i)
            DO k=1,nvnr(j,i)*nvnt(j,i)*nvnp(j,i)
               !READ(10,*)idm1
               READ(10,*)
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   CLOSE(10)
!
!  In the event that node spacing is not identical on all velocity grids,
!  read in the grid again for memory allocation purposes.
!
   IF(vgid.NE.0)THEN
      DEALLOCATE(velnr,covvm)
      DO i=1,nvgt
         DO j=2,ni-1
            IF(nvnr(j,i).GT.mvnr)mvnr=nvnr(j,i)
            IF(nvnt(j,i).GT.mvnt)mvnt=nvnt(j,i)
            IF(nvnp(j,i).GT.mvnp)mvnp=nvnp(j,i)
         ENDDO
      ENDDO
      ALLOCATE(velnr(mvnr,mvnt,mvnp,ni-1,nvgt))
      ALLOCATE(covvm(mvnr,mvnt,mvnp,ni-1,nvgt))
      OPEN(UNIT=10,FILE=refvgfile,STATUS='old')
      READ(10,*)
      DO i=1,nvgt
         DO j=1,ni-1
            READ(10,*)idm1
            READ(10,*)rd1
            READ(10,*)rd1
            DO k=1,nvnr(j,i)
               DO l=1,nvnt(j,i)
                  DO m=1,nvnp(j,i)
                     READ(10,*)velnr(k,l,m,j,i),covvm(k,l,m,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CLOSE(10)
   ENDIF
!
!  Now read in current velocity grid values if iteration
!  number is greater than one.
!
   ALLOCATE(veln(mvnr,mvnt,mvnp,ni-1,nvgt))
!   ALLOCATE(velndws(mvnr,mvnt,mvnp,ni-1,nvgt))
   IF(invstep.GT.1)THEN
      OPEN(UNIT=10,FILE=vgfile,STATUS='old')
      READ(10,*)
      DO i=1,nvgt
         DO j=1,ni-1
            READ(10,*)idm1
            READ(10,*)rd1
            READ(10,*)rd1
            DO k=1,nvnr(j,i)
               DO l=1,nvnt(j,i)
                  DO m=1,nvnp(j,i)
                     READ(10,*)veln(k,l,m,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CLOSE(10)
   ENDIF
ENDIF
!
! Read in complete interface grid if interface depths are
! to be inverted for.
!
IF(pii.EQ.1)THEN
!
!  Start with reference interface depths and associated
!  covariance matrix
!
   OPEN(UNIT=10,FILE=refigfile,STATUS='old')
   READ(10,*)ni
   READ(10,*)nint,ninp
   READ(10,*)gnsit,gnsip
   READ(10,*)goit,goip
   ALLOCATE(intnr(nint,ninp,ni),covim(nint,ninp,ni))
   DO i=1,ni
      DO j=1,nint
         DO k=1,ninp
            READ(10,*)intnr(j,k,i),covim(j,k,i)
         ENDDO
      ENDDO
   ENDDO
   CLOSE(10)
!
!  Read in current interface grid depths if required.
!
   ALLOCATE(intn(nint,ninp,ni))
   IF(invstep.GT.1)THEN
      OPEN(UNIT=10,FILE=igfile,STATUS='old')
      READ(10,*)idm1
      READ(10,*)idm1
      READ(10,*)rd1
      READ(10,*)rd1
      DO i=1,ni
         DO j=1,nint
            DO k=1,ninp
               READ(10,*)intn(j,k,i)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(10)
   ENDIF
ENDIF
!
! Read in Source coordinates if source locations are
! to be inverted for.
!
IF(psi.EQ.1)THEN
!
!  Start off with reference source coordinates and
!  a priori model covariance.
!
   OPEN(UNIT=10,FILE=refscfile,STATUS='old')
   READ(10,*)ns
   ALLOCATE(sradr(ns),slatr(ns),slonr(ns))
   ALLOCATE(covsr(ns),covst(ns),covsp(ns))
   ALLOCATE(lots(ns),nps(ns),tpath(ns))
   ALLOCATE(nsdf(maxs,ns))
   ALLOCATE(stp(ns))
   ALLOCATE(paths(maxs,maxp,ns),patht(maxs,maxp,ns))
   stp=0.0
   DO i=1,ns
      READ(10,*)lots(i)
      IF(lots(i).EQ.1)READ(10,'(a10)')tpath(i)
      isw=0
      DO j=1,nspi
         IF(ids(j).EQ.i)isw=1
      ENDDO
      IF(isw.EQ.0)THEN
         READ(10,*)sradr(i),slatr(i),slonr(i)
      ELSE
         READ(10,*)sradr(i),slatr(i),slonr(i),covsr(i),covst(i),covsp(i)
      ENDIF
      READ(10,*)nps(i)
      DO j=1,nps(i)
         READ(10,*)nsdf(j,i)
         READ(10,*)paths(1:2*nsdf(j,i),j,i)
         READ(10,*)patht(1:nsdf(j,i),j,i)
      ENDDO
   ENDDO
   CLOSE(10)
!
!  Read in current source locations if required, including
!  source time perturbations.
!
   ALLOCATE(srad(ns),slat(ns),slon(ns))
   IF(invstep.GT.1)THEN
      OPEN(UNIT=10,FILE=scfile,STATUS='old')
      READ(10,*)ns
      DO i=1,ns
         READ(10,*)idm1
         IF(idm1.EQ.1)READ(10,*)
         READ(10,*)srad(i),slat(i),slon(i)
         READ(10,*)idm2
         DO j=1,idm2
            READ(10,*)idm3
            READ(10,*)paths(1:2*idm3,j,i)
            READ(10,*)patht(1:idm3,j,i)
         ENDDO
      ENDDO
      CLOSE(10)
      OPEN(UNIT=10,FILE=stfile,STATUS='old')
      READ(10,*)idm1
      DO i=1,nspi
         READ(10,*)stp(i)
      ENDDO
      CLOSE(10)
   ENDIF
ENDIF
!
! Using the above extracted information, create the
! Covariance matrix, the Model vector and the 
! reference model vector. First determine the dimension
! of these vectors.
!
nvpi=0
nipi=0
IF(pvi.EQ.1.AND.nvgi.GT.0)THEN
   DO i=1,nvgi
      idm1=nvnr(idvg(i),idvt(i))*nvnt(idvg(i),idvt(i))
      nvpi=nvpi+idm1*nvnp(idvg(i),idvt(i))
   ENDDO
ENDIF
IF(pii.EQ.1.AND.nigi.GT.0)THEN
   nipi=nipi+nigi*nint*ninp
ENDIF
IF(psi.NE.1)nspi=0
npi=nvpi+nipi+4*nspi
!
! Now create the vectors
!
ALLOCATE(mo(npi),mc(npi),cm(npi),ecmi(npi),dm(npi))
istep=1
IF(nvpi.GT.0)THEN
   DO i=1,nvgi
      DO k=1,nvnp(idvg(i),idvt(i))
         DO l=1,nvnt(idvg(i),idvt(i))
            DO m=1,nvnr(idvg(i),idvt(i))
               mo(istep)=velnr(m,l,k,idvg(i),idvt(i))
               cm(istep)=covvm(m,l,k,idvg(i),idvt(i))**2
               ecmi(istep)=epsv/cm(istep)
               istep=istep+1
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF
IF(nipi.GT.0)THEN
   DO i=1,nigi
      DO j=1,ninp
         DO k=1,nint
            mo(istep)=intnr(k,j,idig(i))
            cm(istep)=covim(k,j,idig(i))**2
            ecmi(istep)=epsi/cm(istep)
            istep=istep+1
         ENDDO
      ENDDO
   ENDDO
ENDIF
IF(nspi.GT.0)THEN
   DO i=1,nspi
      mo(istep)=sradr(ids(i))
      cm(istep)=covsr(ids(i))**2
      ecmi(istep)=epss1/cm(istep)
      istep=istep+1
   ENDDO
   DO i=1,nspi
      mo(istep)=slatr(ids(i))
      cm(istep)=covst(ids(i))**2
      ecmi(istep)=epss1/cm(istep)
      istep=istep+1
   ENDDO
   DO i=1,nspi
      mo(istep)=slonr(ids(i))
      cm(istep)=covsp(ids(i))**2
      ecmi(istep)=epss1/cm(istep)
      istep=istep+1
   ENDDO
   DO i=1,nspi
      mo(istep)=0.0
      cm(istep)=1.0
      ecmi(istep)=epss2/cm(istep)
      istep=istep+1
   ENDDO
ENDIF
!
! Repeat for current model if it exists
!
IF(invstep.GT.1)THEN
   istep=1
   IF(nvpi.GT.0)THEN
      DO i=1,nvgi
         DO k=1,nvnp(idvg(i),idvt(i))
            DO l=1,nvnt(idvg(i),idvt(i))
               DO m=1,nvnr(idvg(i),idvt(i))
                  mc(istep)=veln(m,l,k,idvg(i),idvt(i))
                  dm(istep)=mc(istep)-mo(istep)
                  istep=istep+1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   IF(nipi.GT.0)THEN
      DO i=1,nigi
         DO j=1,ninp
            DO k=1,nint
               mc(istep)=intn(k,j,idig(i))
               dm(istep)=mc(istep)-mo(istep)
               istep=istep+1
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   IF(nspi.GT.0)THEN
      DO i=1,nspi
         mc(istep)=srad(ids(i))
         dm(istep)=mc(istep)-mo(istep)
         istep=istep+1
      ENDDO
      DO i=1,nspi
         mc(istep)=slat(ids(i))
         dm(istep)=(mc(istep)-mo(istep))*pi*(earthr-srad(ids(i)))/180.0
         istep=istep+1
      ENDDO
      DO i=1,nspi
         mc(istep)=slon(ids(i))
         dm(istep)=(mc(istep)-mo(istep))*pi*(earthr-srad(ids(i)))
         dm(istep)=dm(istep)*cos(slat(ids(i))*pi/180.0)/180.0
         istep=istep+1
      ENDDO
      DO i=1,nspi
         mc(istep)=stp(i)
         dm(istep)=mc(istep)-mo(istep)
         istep=istep+1
      ENDDO
   ENDIF
ELSE
   mc=mo
   dm=0.0
ENDIF
!
! Deallocate unnecessary arrays
!
IF(pvi.EQ.1)THEN
   IF(invstep.GT.1)THEN
      DEALLOCATE(velnr,covvm)
   ELSE
      DEALLOCATE(covvm)
   ENDIF
ENDIF
IF(pii.EQ.1)THEN
   IF(invstep.GT.1)THEN
      DEALLOCATE(intnr,covim)
   ELSE
      DEALLOCATE(covim)
   ENDIF
ENDIF
IF(psi.EQ.1)THEN
   IF(invstep.GT.1)THEN
      DEALLOCATE(sradr,slatr,slonr,covsr,covst,covsp)
   ELSE
      DEALLOCATE(covsr,covst,covsp)
   ENDIF
ENDIF
!
! Now read in the observed and model traveltime residuals. In the
! case of teleseismic traveltimes, we also need to input reference
! traveltimes, and in the case of traveltimes from sources which
! are relocated, the model perturbations must be added to the
! model traveltimes. Model traveltimes for invalid rays are set to
! -1.0 by fm3d. Note that reflection matching, which is a feature
! of fm3d, is not currently supported by this code.
!
nnfd = 0
ntr = 0
nnode = nvpi/2
if(surfjoint==0.or.surfjoint==2) then
OPEN(UNIT=10,FILE=otfile,STATUS='old')
OPEN(UNIT=20,FILE=mtfile,STATUS='old')
READ(10,*)ntr
!
! Open Source file and determine the source IDs of all
! teleseismic sources
!
OPEN(UNIT=30,FILE=scfile,STATUS='old')
READ(30,*)ns
ALLOCATE(tsid(ns))
IF(psi.NE.1)ALLOCATE(paths(maxs,maxp,ns),patht(maxs,maxp,ns))
ntels=0
tsid=0
DO i=1,ns
   READ(30,*)idm1
   IF(idm1.EQ.1)THEN
      ntels=ntels+1
      tsid(i)=1
      READ(30,*)
   ENDIF
   READ(30,*)
   READ(30,*)idm2
   DO j=1,idm2
      READ(30,*)idm3
      READ(30,*)paths(1:2*idm3,j,i)
      READ(30,*)patht(1:idm3,j,i)
   ENDDO
ENDDO
CLOSE(30)
!
! Since the Frechet matrix is very large, we want to
! allocate the correct amount of memory to it. Thus, we
! will first read in the entire matrix to determine the
! number of non-zeo elements.
!
OPEN(UNIT=30,FILE=frdatfile,STATUS='old')
DO i=1,ntr
   READ(30,*)idm1,idm2,idm3,idm4,nrow
   nnfd=nnfd+nrow
   cnt = 0
   IF(nrow.GT.0)THEN
      DO j=1,nrow
         READ(30,*)tmp
         if(tmp>nnode.and.tmp<=nvpi)cnt=cnt+1 
      ENDDO
if(vpvs==1) then
nnfd = nnfd+cnt
endif
   ENDIF
ENDDO
CLOSE(30)
endif

! hongjian fang @ethz... adding surface wave data
!----------------------------------------------------------------
!surfjoint = 1
frdatfilesurf = 'frechetsurf.dat'
otfilesurf = 'otimessurf.dat'
mtfilesurf = 'mtimessurf.dat'
ntrsurf = 0
if (pvi>0.and.(surfjoint == 1.or.surfjoint==2)) then
OPEN(UNIT=30,FILE=frdatfilesurf,STATUS='old')
OPEN(UNIT=35,FILE=otfilesurf,STATUS='old')
OPEN(UNIT=34,FILE=mtfilesurf,STATUS='old')
read(35,*) ntrsurf


DO i=1,ntrsurf
   READ(30,*)idm1,nrow
   nnfd=nnfd+nrow
   cnt = 0
   IF(nrow.GT.0)THEN
      DO j=1,nrow
         READ(30,*)tmp
         if(tmp>nnode.and.tmp<=nvpi)cnt=cnt+1 
      ENDDO
if(vpvs==1) then
nnfd = nnfd+cnt
endif
   ENDIF
ENDDO
CLOSE(30)
endif
ntr = ntr+ntrsurf
!----------------------------------------------------------------

ALLOCATE(frech(nnfd),fcoln(nnfd),cnfe(0:ntr))
ALLOCATE(dobs(ntr),dmod(ntr),cd(ntr))

idm2o=0
jstep=0
kstep=0
if(surfjoint==0.or.surfjoint==2) then

OPEN(UNIT=30,FILE=frdatfile,STATUS='old')
!
! If teleseismic sources exist, we need to
! read in reference teleseismic traveltimes
!
IF(ntels.GT.0)THEN
   OPEN(UNIT=40,FILE=rtfile,STATUS='old')
   ALLOCATE(istel(ntr))
   istel=0
!
!  Set up for removing mean from teleseismic residual.
!
   ALLOCATE(mtmean(ns))
   mtmean=0.0
ENDIF
istep=1
istepo=0
jstep=0
cnfe(0)=0
kstep=0
isum=0
iswt=0
isw2=0

! hidden bug, invert for vp/vs,  only works for 1 layer...
!if (flex==1) stp = 0
DO i=1,ntr-ntrsurf
   READ(10,*)idm1,idm2,idm3,idm4,dobs(istep),cd(istep)
!
!  Apply source time correction if required
!
   IF(nspi.GT.0.AND.invstep.GT.1)THEN
      rd1=0.0
      DO j=1,nspi
         IF(idm2.EQ.ids(j))THEN
            rd1=stp(j)
         ENDIF
      ENDDO
   ENDIF
   READ(20,*)idm1,idm2,idm3,idm4,dmod(istep)
   IF(ntels.GT.0)THEN
      IF(tsid(idm2).EQ.1)istel(i)=idm2
   ENDIF
   IF(nspi.GT.0.AND.invstep.GT.1)dmod(istep)=dmod(istep)-rd1
!
!  Read in reference teleseismic traveltimes if required
!
   IF(ntels.GT.0)THEN
      READ(40,*)idm1,idm2,idm3,idm4,dref 
      IF(istel(i).GT.0)THEN
         IF(dref.LE.0.0)THEN
            dmod(istep)=-1.0
         ELSE
            IF(dmod(istep).GT.0.0)iswt=1
            dmod(istep)=dmod(istep)-dref
            IF(istepo.NE.istep)mtmean(idm2)=mtmean(idm2)+dmod(istep)
            isum=isum+1
            IF(isw2.GT.0)THEN
               IF(idm2o.NE.idm2)THEN
                  IF(isum-1.GT.0)mtmean(idm2o)=mtmean(idm2o)/REAL(isum-1)
                  isum=1
               ELSE IF(i.EQ.ntr)THEN
                  mtmean(idm2o)=mtmean(idm2o)/REAL(isum)
               ENDIF
            ENDIF
            istepo=istep
            kstep=kstep+1
            istel(istep)=idm2
         ENDIF
         idm2o=idm2
         isw2=1
      ELSE IF(i.EQ.ntr)THEN
         mtmean(idm2o)=mtmean(idm2o)/REAL(isum)
      ENDIF
   ENDIF
   IF(dobs(istep).LT.-50.0)iswt=0
   READ(30,*)idm1,idm2,idm3,idm4,nrow
   if(inversionScheme==1) then
   cd(istep)=cd(istep)**2
!   else
!   cd(istep)=cd(istep)/3.0
   endif
   cnfe(istep)=jstep+nrow
   IF(nrow.GT.0)THEN
      jstep=jstep+1
      jup=jstep+nrow-1
      cnt = 0
      DO j=jstep,jup
         READ(30,*)fcoln(j),frech(j)
if(vpvs==1.and.fcoln(j)>nnode.and.fcoln(j)<=nvpi) then
cnt = cnt+1
tmp1 = frech(j)
frech(j) = tmp1*(-mc(fcoln(j))**2/mc(fcoln(j)-nnode))
frech(jup+cnt) = tmp1*(mc(fcoln(j))/mc(fcoln(j)-nnode))
fcoln(jup+cnt) = fcoln(j)-nnode
endif
      ENDDO
      IF(dobs(istep).GT.0.0.AND.dmod(istep).GT.0.0)THEN
         jstep=jup
if(vpvs==1) then
jstep=jstep+cnt
cnfe(istep)=cnfe(istep)+cnt
endif
      ELSE IF(iswt.EQ.1)THEN
         jstep=jup
if(vpvs==1) then
jstep=jstep+cnt
cnfe(istep)=cnfe(istep)+cnt
endif
      ELSE
         jstep=jstep-1
      ENDIF
   ENDIF
   IF(dobs(istep).GT.0.0.AND.dmod(istep).GT.0.0)THEN
      istep=istep+1
   ELSE IF(iswt.EQ.1)THEN
      istep=istep+1
   ENDIF
   iswt=0
ENDDO
CLOSE(10)
CLOSE(20)
CLOSE(30)
IF(ntels.GT.0)CLOSE(40)
endif ! for data type (surfjoint==0 or 2)
! hongjian fang @ethz... adding surface wave data
!----------------------------------------------------------------
!print*,istep,jstep
if (surfjoint==1) then
istep=1
jstep=0
endif
if (pvi>0.and.(surfjoint == 1.or.surfjoint==2)) then
OPEN(UNIT=36,FILE=frdatfilesurf,STATUS='old')
do i=1,ntrsurf
read(34,*) dmod(istep)
read(35,*) dobs(istep),cd(istep)
READ(36,*)idm1,nrow
if(inversionScheme==1) then
cd(istep)=cd(istep)**2
!else
!cd(istep)=cd(istep)/3.0
endif
cnfe(istep)=jstep+nrow
IF(nrow.GT.0)THEN
   jstep=jstep+1
   jup=jstep+nrow-1
   cnt = 0
   DO j=jstep,jup
      READ(36,*)fcoln(j),frech(j)
if(vpvs==1.and.fcoln(j)>nnode.and.fcoln(j)<=nvpi) then
cnt = cnt+1
tmp1 = frech(j)
frech(j) = tmp1*(-mc(fcoln(j))**2/mc(fcoln(j)-nnode))
frech(jup+cnt) = tmp1*(mc(fcoln(j))/mc(fcoln(j)-nnode))
fcoln(jup+cnt) = fcoln(j)-nnode
endif
   ENDDO
   IF(dobs(istep).GT.0.0.AND.dmod(istep).GT.0.0)THEN
      jstep=jup
if(vpvs==1) then
jstep=jstep+cnt
cnfe(istep)=cnfe(istep)+cnt
endif
   ELSE
      jstep=jstep-1
   ENDIF
ENDIF
IF(dobs(istep).GT.0.0.AND.dmod(istep).GT.0.0)THEN
   istep=istep+1
ENDIF
enddo
CLOSE(34)
CLOSE(35)
CLOSE(36)

!mean = sum(dobs(istep-ntrsurf:istep-1)-dmod(istep-ntrsurf:istep-1))/ntrsurf
!std_surf = sqrt(sum((dobs(istep-ntrsurf:istep-1)-dmod(istep-ntrsurf:istep-1))**2)/ntrsurf-mean**2)
!write(*,'(a,f8.1,f8.2)'),'mean,std_devs and rms:', 1000*mean, 1000*std_surf

endif ! for data type (surfjoint==1 or 2)

! hongjian fang @ethz... adding surface wave data
!----------------------------------------------------------------
!if (surfjoint == 1) then
!ntr = ntr-ntrsurf
WRITE(6,*)'------------------------------------------------------------------'
write(6,*)'no. of non-zero in G and nodes', nnfd,nnode
write(6,*)'data for body wave & surface wave', ntr-ntrsurf, ntrsurf
!endif
!----------------------------------------------------------------

!----------------------------------------------------------------

ntr=istep-1
!print*,'traces number',ntr
!
! Rearrange Frechet derivatives to separate out source parameter
! classes only if sources are inverted for.
!
IF(nspi.GT.0)THEN
   DO i=1,ntr
      IF(cnfe(i).GT.cnfe(i-1))THEN
         DO j=cnfe(i-1)+1,cnfe(i)
            IF(fcoln(j).GT.nvpi+nipi)THEN
               idm1=fcoln(j)-(nvpi+nipi)
               idm2=MOD(idm1,4)
               IF(idm2.EQ.0)idm2=4
               idm3=(idm1-idm2)/4+1
               fcoln(j)=(idm2-1)*nspi+idm3+nvpi+nipi
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDIF
!
! Remove mean from model traveltimes if required
!
IF(rmtr.EQ.1.AND.kstep.GT.0)THEN
   DO i=1,ntr-ntrsurf
      IF(istel(i).GT.0)dmod(i)=dmod(i)-mtmean(istel(i))
   ENDDO
ENDIF
!
! Now construct the transpose of the Frechet matrix
!
!print*,nnfd,npi,ntr
ALLOCATE(tfrech(nnfd),tfcoln(nnfd),tcnfe(0:npi),stpv(ntr))
stpv=1
jstep=0
tcnfe(0)=0
DO i=1,npi
   DO j=1,ntr
      IF(stpv(j).LE.cnfe(j)-cnfe(j-1))THEN
         cstp=fcoln(cnfe(j-1)+stpv(j))
         IF(cstp.EQ.i)THEN
            jstep=jstep+1
            tfrech(jstep)=frech(cnfe(j-1)+stpv(j))
            tfcoln(jstep)=j
            stpv(j)=stpv(j)+1
         ENDIF
      ENDIF
   ENDDO
   tcnfe(i)=jstep
ENDDO
DEALLOCATE(stpv)
!print*,'nozeros in transpose',jstep
!print*,tfrech(jstep-5:jstep+5)
!print*,tfrech(nnfd-20:nnfd)
!
! We have now set up all the vectors and matrices we
! require for the subspace inversion scheme. Call a
! routine for performing the inversion
!
if (inversionScheme == 1) then
CALL subspace
if (pvi>0) then
write(*,*) 'no. of vel/interfaces/sources:', nvpi,nipi,nspi
write(*,*) 'min. and max. velocity variation', minval(dm(1:nvpi)),maxval(dm(1:nvpi))
endif

else 
! call lsmr instead of subspace to solve Ax = b
! first smooth velocity and interface grids if both are included in the inversion
allocate(iw(2*(nnfd+9*npi)+1),stat=checkstat)
if(checkstat>0) stop 'error allocating iw'
allocate(rw(nnfd+9*npi),stat=checkstat)
if(checkstat>0) stop 'error allocating rw'
allocate(col(nnfd+9*npi),stat=checkstat)
if(checkstat>0) stop 'error allocating rw'
allocate(dtrav(ntr+9*npi),stat=checkstat)
if(checkstat>0) stop 'error allocating dtrav'
allocate(dataweight(ntr+9*npi),stat=checkstat)
if(checkstat>0) stop 'error allocating dataweight'
allocate(norm(npi),stat=checkstat)
if(checkstat>0) stop 'error allocating norm'
!allocate(norm_dws(npi),stat=checkstat)
!if(checkstat>0) stop 'error allocating norm_dws'
iw = 0
rw = 0.
jstep = 0
do m = 1,ntr
  do j = cnfe(m-1)+1,cnfe(m)
    jstep = jstep + 1
    iw(1+jstep) = m
    rw(jstep) = frech(jstep)/cd(m)
    col(jstep) = fcoln(jstep)
  enddo
enddo

jstep_tmp = jstep

!aguement smooth matrix under the sensitivity matrix
damp = epsilon
is2=0
IF(nvpi.GT.0)THEN
  DO i=1,nvgi
    DO k=1,nvnp(idvg(i),idvt(i))
      DO l=1,nvnt(idvg(i),idvt(i))
        DO m=1,nvnr(idvg(i),idvt(i))
          is2=is2+1
          IF(m.NE.1.AND.m.NE.nvnr(idvg(i),idvt(i)))THEN
            is1=is2-1
            is3=is2+1
            iw(1+jstep+1) = istep
            iw(1+jstep+2) = istep
            iw(1+jstep+3) = istep
            rw(jstep+1) = 1.0*etav
            rw(jstep+2) = -2.0*etav
            rw(jstep+3) = 1.0*etav
            col(jstep+1) = is1 
            col(jstep+2) = is2
            col(jstep+3) = is3
            istep = istep + 1
            jstep = jstep + 3
          else
            rw(jstep+1) = 5*etav
            col(jstep+1) = is2
            iw(1+jstep+1) = istep
            jstep = jstep+1
            istep = istep+1
          ENDIF
          IF(l.NE.1.AND.l.NE.nvnt(idvg(i),idvt(i)))THEN
            is1=is2-nvnr(idvg(i),idvt(i))
            is3=is2+nvnr(idvg(i),idvt(i))
            iw(1+jstep+1) = istep
            iw(1+jstep+2) = istep
            iw(1+jstep+3) = istep
            rw(jstep+1) = 1.0*etav
            rw(jstep+2) = -2.0*etav
            rw(jstep+3) = 1.0*etav
            col(jstep+1) = is1 
            col(jstep+2) = is2
            col(jstep+3) = is3
            istep = istep + 1
            jstep = jstep + 3
          else
            rw(jstep+1) = 5*etav
            col(jstep+1) = is2
            iw(1+jstep+1) = istep
            jstep = jstep+1
            istep = istep+1
          ENDIF
          IF(k.NE.1.AND.k.NE.nvnp(idvg(i),idvt(i)))THEN
            is1=is2-nvnr(idvg(i),idvt(i))*nvnt(idvg(i),idvt(i))
            is3=is2+nvnr(idvg(i),idvt(i))*nvnt(idvg(i),idvt(i))
            iw(1+jstep+1) = istep
            iw(1+jstep+2) = istep
            iw(1+jstep+3) = istep
            rw(jstep+1) = 1.0*etav
            rw(jstep+2) = -2.0*etav
            rw(jstep+3) = 1.0*etav
            col(jstep+1) = is1 
            col(jstep+2) = is2
            col(jstep+3) = is3
            istep = istep + 1
            jstep = jstep + 3
          else
            rw(jstep+1) = 5*etav
            col(jstep+1) = is2
            iw(1+jstep+1) = istep
            jstep = jstep+1
            istep = istep+1
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF
!
! Now add interface parameters if they exist.
!
IF(nipi.GT.0)THEN
  DO i=1,nigi
    DO j=1,ninp
      DO k=1,nint
        is2=is2+1
        IF(k.NE.1.AND.k.NE.nint)THEN
          is1=is2-1
          is3=is2+1
          iw(1+jstep+1) = istep
          iw(1+jstep+2) = istep
          iw(1+jstep+3) = istep
          rw(jstep+1) = 1.0*etai
          rw(jstep+2) = -2.0*etai
          rw(jstep+3) = 1.0*etai
          col(jstep+1) = is1 
          col(jstep+2) = is2
          col(jstep+3) = is3
          istep = istep + 1
          jstep = jstep + 3
        else
          rw(jstep+1) = 5*etai
          col(jstep+1) = is2
          iw(1+jstep+1) = istep
          jstep = jstep+1
          istep = istep+1
        ENDIF
        IF(j.NE.1.AND.j.NE.ninp)THEN
          is1=is2-nint
          is3=is2+nint
          iw(1+jstep+1) = istep
          iw(1+jstep+2) = istep
          iw(1+jstep+3) = istep
          rw(jstep+1) = 1.0*etai
          rw(jstep+2) = -2.0*etai
          rw(jstep+3) = 1.0*etai
          col(jstep+1) = is1 
          col(jstep+2) = is2
          col(jstep+3) = is3
          istep = istep + 1
          jstep = jstep + 3
        else
          rw(jstep+1) = 5*etai
          col(jstep+1) = is2
          iw(1+jstep+1) = istep
          jstep = jstep+1
          istep = istep+1
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDIF

! regularization for sources position and origin times
if (nspi>0) then
  do i=1,nspi
    iw(1+jstep+i) = istep
    rw(jstep+i) = 3.0*epss1
    col(jstep+i) = nvpi+nipi+i
    istep = istep+1
  enddo
  do i=1,2*nspi
    iw(1+jstep+nspi+i) = istep
    rw(jstep+nspi+i) = 1.0*epss1
    col(jstep+nspi+i) = nvpi+nipi+nspi+i
    istep = istep+1
  enddo
  do i=1,nspi
    iw(1+jstep+3*nspi+i) = istep
    rw(jstep+3*nspi+i) = 1.0*epss2
    col(jstep+3*nspi+i) = nvpi+nipi+3*nspi+i
    istep = istep+1
  enddo
jstep = jstep+4*nspi
endif

m = istep -1
l = npi
iw(1) = jstep
do i = 1,jstep
  iw(1+jstep+i) = col(i)
enddo

!DO i=1,ntr
!   dtrav(i)=(dmod(i)-dobs(i))/cd(i)
!ENDDO
DO i=1,ntr
   dtrav(i)=(dmod(i)-dobs(i))
ENDDO
print*,'weighted rms',sum(dtrav(1:ntr)**2)/ntr
! downweight data with large residual
dataweight = 1.0
if(surfjoint==0 .or.surfjoint==2) then
mean = sum(dtrav(1:ntr-ntrsurf))/(ntr-ntrsurf)
std_surf = sqrt(sum((dtrav(1:ntr-ntrsurf))**2)/(ntr-ntrsurf)-mean**2)
DO i=1,ntr-ntrsurf
if (abs(dtrav(i))>threshold*std_surf) then
dataweight(i) = exp(-(abs(dtrav(i))/(threshold*std_surf)-1))
endif
dtrav(i)=dtrav(i)*dataweight(i)
ENDDO
mean = sum(dtrav(1:ntr-ntrsurf))/(ntr-ntrsurf)
std_surf = sqrt(sum((dtrav(1:ntr-ntrsurf))**2)/(ntr-ntrsurf)-mean**2)
write(*,'(a,f10.1,f10.1)'),'mean,std_devs for body waves:', 1000*mean, 1000*std_surf
endif

if(surfjoint==1) then
mean = sum(dtrav(ntr-ntrsurf+1:ntr))/(ntrsurf)
std_surf = sqrt(sum((dtrav(ntr-ntrsurf+1:ntr))**2)/(ntrsurf)-mean**2)
!print*,mean,std_surf
DO i=ntr-ntrsurf+1,ntr
if (abs(dtrav(i))>threshold*std_surf) then
dataweight(i) = exp(-(abs(dtrav(i))/(threshold*std_surf)-1))
endif
dtrav(i)=dtrav(i)*dataweight(i)
ENDDO
mean = sum(dtrav(ntr-ntrsurf+1:ntr))/(ntrsurf)
std_surf = sqrt(sum((dtrav(ntr-ntrsurf+1:ntr))**2)/(ntrsurf)-mean**2)
write(*,'(a,f10.1,f10.1)'),'mean,std_devs for surface waves:', 1000*mean, 1000*std_surf
endif


if(surfjoint==2) then
mean = sum(dtrav(ntr-ntrsurf+1:ntr))/(ntrsurf)
std_surf = sqrt(sum((dtrav(ntr-ntrsurf+1:ntr))**2)/(ntrsurf)-mean**2)
DO i=ntr-ntrsurf+1,ntr
if (abs(dtrav(i))>threshold*std_surf) then
dataweight(i) = exp(-(abs(dtrav(i))/(threshold*std_surf)-1))* &
sqrt(real(ntr-ntrsurf)/ntrsurf)*surfweight
else  
dataweight(i) = surfweight*sqrt(real(ntr-ntrsurf)/ntrsurf)
endif
!dtrav(i)=dtrav(i)*dataweight(i)
ENDDO
mean = sum(dtrav(ntr-ntrsurf+1:ntr))/(ntrsurf)
std_surf = sqrt(sum((dtrav(ntr-ntrsurf+1:ntr))**2)/(ntrsurf)-mean**2)
write(*,'(a,f10.1,f10.1)'),'mean,std_devs for surface waves:', 1000*mean, 1000*std_surf
DO i=ntr-ntrsurf+1,ntr
dtrav(i)=dtrav(i)*dataweight(i)
ENDDO
endif

DO i=1,ntr
   dtrav(i)=dtrav(i)/cd(i)
ENDDO

do i = 1,jstep
rw(i) = rw(i)*dataweight(iw(1+i))
enddo
do i=ntr+1,m
  dtrav(i)=0.
enddo

!calculate dws
!norm_dws = 0
!do j = 1, jstep_tmp
!norm_dws(iw(1+jstep+j)) = norm_dws(iw(1+jstep+j))+abs(rw(j))
!enddo


norm = 0
do i=1,jstep
  norm(iw(1+jstep+i)) = norm(iw(1+jstep+i)) + rw(i)**2
enddo

do i=1,npi
  norm(i) = sqrt(norm(i)/m)
enddo

! normilize each column to use a single damping
do i = 1,jstep
  rw(i) = rw(i)/norm(iw(1+jstep+i))
enddo


    leniw = 2*jstep+1
    lenrw = jstep 
    dm = 0
    atol = 1e-5
    btol = 1e-5
    conlim = 100
    itnlim = 500
    istop = 0
    anorm = 0.0
    acond = 0.0
    arnorm = 0.0
    xnorm = 0.0
    localSize = l/4



        nout = 63
        write (itnum,'(I0)') invstep
        open(nout,file='lsmrout'//trim(itnum)//'.txt')
print*,'-----------------------------------------------------'
print*,'iteration:',invstep
!print*,'min. and max. dws',minval(norm_dws),maxval(norm_dws)
print*,'min. and max. dws',minval(norm),maxval(norm)
    call LSMR(m, l, leniw, lenrw,iw,rw,dtrav, damp,&
      atol, btol, conlim, itnlim, localSize, nout,&
      dm, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
    dm = -dm
    do i = 1,npi
      dm(i) = dm(i)/norm(i)
    enddo
    !if(istop==3) print*,'istop = 3, large condition number'
    deallocate(iw,col)
    deallocate(rw,dtrav,norm,dataweight)
    close(nout)
if (pvi>0) then
write(*,*) 'no. of vel/interfaces/sources:', nvpi,nipi,nspi
write(*,*) 'min. and max. velocity variation', minval(dm(1:nvpi)),maxval(dm(1:nvpi))
endif
if(psi>0) then
write(*,*) 'min. and max. srcs location variation: rad', minval(dm(nvpi+nipi+1:nvpi+nipi+nspi)),&
                maxval(dm(nvpi+nipi+1:nvpi+nipi+nspi))
write(*,*) 'min. and max. srcs location variation: lat', 0.009*minval(dm(nvpi+nipi+nspi+1:nvpi+nipi+2*nspi)),&
                0.009*maxval(dm(nvpi+nipi+nspi+1:nvpi+nipi+2*nspi))
write(*,*) 'min. and max. srcs location variation: lon', 0.009*minval(dm(nvpi+nipi+2*nspi+1:nvpi+nipi+3*nspi)),& 
                0.009*maxval(dm(nvpi+nipi+2*nspi+1:nvpi+nipi+3*nspi))
write(*,*) 'min. and max. srcs location variation: stp', minval(dm(nvpi+nipi+3*nspi+1:nvpi+nipi+4*nspi)),&
                maxval(dm(nvpi+nipi+3*nspi+1:nvpi+nipi+4*nspi))
print*,'-----------------------------------------------------'
endif

endif


!open(64,file='dm.dat')
!do i=1,npi
!write(64,*) dm(i),mc(i)
!enddo
!close(64)
!
! Now write new model to file
!
istep=0
IF(pvi.EQ.1)THEN
   IF(invstep.EQ.1)THEN
      veln=velnr
      DEALLOCATE(velnr)
   ENDIF
   OPEN(UNIT=10,FILE=vgfile,STATUS='unknown')
   WRITE(10,*)ni-1,nvgt
   OPEN(UNIT=11,FILE=trim(vgfile)//trim('vpvs'),STATUS='unknown')
   WRITE(11,*)ni-1,nvgt
!   OPEN(UNIT=12,FILE=trim(vgfile)//trim('dws'),STATUS='unknown')
!   WRITE(12,*)ni-1,nvgt

  ! print*,npi
  ! print*
  ! do j=1,npi
  ! print*,dm(j)
  ! enddo
      DO j=1,ni-1
   DO i=1,nvgt
    !     WRITE(10,*)nvnr(j,i),nvnt(j,i),nvnp(j,i)
    !     WRITE(10,*)gnsr(j,i),gnst(j,i),gnsp(j,i)
    !     WRITE(10,*)gor(j,i),got(j,i),gop(j,i)
         !idm1=0
         !DO k=1,nvgi
         !   IF(idvg(k).EQ.j)THEN
         !      IF(idvt(k).EQ.i)idm1=1
         !   ENDIF
         !ENDDO
         !IF(idm1.EQ.1)THEN
            DO k=1,nvnp(j,i)
               DO l=1,nvnt(j,i)
                  DO m=1,nvnr(j,i)
                     istep=istep+1
                     if (abs(dm(istep))>0.3) dm(istep) = dm(istep)/abs(dm(istep))*0.3
                     veln(m,l,k,j,i)=mc(istep)+dm(istep)
!                     velndws(m,l,k,j,i)=norm_dws(istep)
                     if(vpvs==1) then
                     if(istep<=nnode) then
                     veln(m,l,k,j,i)=mc(istep)+dm(istep)
                     else
                     veln(m,l,k,j,i)=mc(istep-nnode)/mc(istep)+dm(istep)
                     !veln(m,l,k,j,i) = mc(istep-nnode)/veln(m,l,k,j,i)
                     veln(m,l,k,j,i) = veln(m,l,k,j,i-1)/veln(m,l,k,j,i)
                     endif
                     endif
                     IF(veln(m,l,k,j,i).LT.mpv)veln(m,l,k,j,i)=mpv
                  ENDDO
               ENDDO
            ENDDO
         !ENDIF
         enddo
         enddo

  DO i=1,nvgt
     DO j=1,ni-1
         WRITE(10,*)nvnr(j,i),nvnt(j,i),nvnp(j,i)
         WRITE(10,*)gnsr(j,i),gnst(j,i),gnsp(j,i)
         WRITE(10,*)gor(j,i),got(j,i),gop(j,i)
         DO k=1,nvnr(j,i)
            DO l=1,nvnt(j,i)
               DO m=1,nvnp(j,i)
                  WRITE(10,*)veln(k,l,m,j,i)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   CLOSE(10)
! write out vp/vs
   DO i=1,nvgt
     DO j=1,ni-1
         WRITE(11,*)nvnr(j,i),nvnt(j,i),nvnp(j,i)
         WRITE(11,*)gnsr(j,i),gnst(j,i),gnsp(j,i)
         WRITE(11,*)gor(j,i),got(j,i),gop(j,i)
         DO k=1,nvnr(j,i)
            DO l=1,nvnt(j,i)
               DO m=1,nvnp(j,i)
                if(i==1) then
                  WRITE(11,*)veln(k,l,m,j,i)
                else
                  write(11,*)veln(k,l,m,j,i-1)/veln(k,l,m,j,i)
                endif
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   CLOSE(11)

! write out dws
!   DO i=1,nvgt
!     DO j=1,ni-1
!         WRITE(12,*)nvnr(j,i),nvnt(j,i),nvnp(j,i)
!         WRITE(12,*)gnsr(j,i),gnst(j,i),gnsp(j,i)
!         WRITE(12,*)gor(j,i),got(j,i),gop(j,i)
!         DO k=1,nvnr(j,i)
!            DO l=1,nvnt(j,i)
!               DO m=1,nvnp(j,i)
!                write(12,*) velndws(k,l,m,j,i)
!               ENDDO
!            ENDDO
!         ENDDO
!      ENDDO
!   ENDDO
!   CLOSE(12)


ENDIF
IF(pii.EQ.1)THEN
   IF(invstep.EQ.1)THEN
      intn=intnr
      DEALLOCATE(intnr)
   ENDIF
   OPEN(UNIT=10,FILE=igfile,STATUS='unknown')
   WRITE(10,*)ni
   WRITE(10,*)nint,ninp
   WRITE(10,*)gnsit,gnsip
   WRITE(10,*)goit,goip
!
!  Locate the top and bottom radius of the propagation grid.
!
   OPEN(UNIT=20,FILE=pgfile,STATUS='old')
   READ(20,*)npgnr
   READ(20,*)pgb
   READ(20,*)pgt
   CLOSE(20)
   pgt=pgt+earthr
   pgb=pgt-(npgnr-1)*pgb
   DO i=1,ni
      idm1=0
      DO k=1,nigi
         IF(idig(k).EQ.i)idm1=1
      ENDDO
      IF(idm1.EQ.1)THEN
         DO j=1,ninp
            DO k=1,nint
               istep=istep+1
               intn(k,j,i)=mc(istep)+dm(istep)
            ENDDO
         ENDDO
      ENDIF
!
!     Now test to see whether any interfaces intersect each other or
!     the top and bottom boundaries of the propagation grid. Apply corrections
!     where necessary.
!
      IF(i.EQ.1)THEN
!
!        Test top surface
!
         DO j=1,nint
            DO k=1,ninp  
               IF(intn(j,k,i).GT.pgt-mdbi)THEN
                  intn(j,k,i)=pgt-mdbi
               ELSE IF(intn(j,k,i).LT.pgb+ni*mdbi)THEN
                  intn(j,k,i)=pgb+ni*mdbi
               ENDIF
               WRITE(10,*)intn(j,k,i)
            ENDDO
         ENDDO
      ELSE
!
!        Test remaining surfaces
!
         DO j=1,nint
            DO k=1,ninp
               IF(intn(j,k,i).GT.intn(j,k,i-1)-mdbi)THEN
                  intn(j,k,i)=intn(j,k,i-1)-mdbi
               ELSE IF(intn(j,k,i).LT.pgb+(ni-i+1)*mdbi)THEN
                  intn(j,k,i)=pgb+(ni-i+1)*mdbi
               ENDIF
               WRITE(10,*)intn(j,k,i)
            ENDDO
         ENDDO
      ENDIF
   ENDDO
   CLOSE(10)
ENDIF
IF(psi.EQ.1)THEN
   IF(invstep.EQ.1)THEN
      srad=sradr
      slat=slatr
      slon=slonr
      DEALLOCATE(sradr,slatr,slonr)
   ENDIF
   if(pvi>0) then
      !pgt = earthr - (gor(1,1)+(nvnr(1,1)-3)*gnsr(1,1))
      !pgb = earthr - gor(1,1)-2*gnsr(1,1)
      pgt = earthr - (gor(1,1)+(nvnr(1,1)-1)*gnsr(1,1))
      pgb = earthr - gor(1,1)
    endif
   DO j=1,nspi
      istep=istep+1
      if (abs(dm(istep))>1.0) dm(istep) = dm(istep)/abs(dm(istep))*1.0
      srad(ids(j))=mc(istep)-dm(istep)
   if(pvi>0) then
      if (srad(ids(j)) < pgt) then
        print*,'warning: rad outside (up)',srad(ids(j)),mc(istep),dm(istep)
        srad(ids(j)) = pgt+gnsr(1,1)
      elseif(srad(ids(j)) > pgb) then
        print*,'warning: rad outside (down)',srad(ids(j)),mc(istep),dm(istep)
        srad(ids(j)) = pgb-gnsr(1,1)
      endif
    endif
   ENDDO
   if(pvi>0) then
      !pgt = (got(1,1)+(nvnt(1,1)-3)*gnst(1,1))*180.0/pi
      !pgb = (got(1,1)+2*gnst(1,1))*180.0/pi
      pgt = (got(1,1)+(nvnt(1,1)-1)*gnst(1,1))*180.0/pi
      pgb = got(1,1)*180.0/pi
    endif
   DO j=1,nspi
      istep=istep+1
! dis/R*180/pi  rad-->degree
      dm(istep)=dm(istep)*180.0/(pi*(earthr-srad(ids(j))))
      if (abs(dm(istep))>0.05) dm(istep) = dm(istep)/abs(dm(istep))*0.05
      slat(ids(j))=mc(istep)+dm(istep)
   if(pvi>0) then
      if (slat(ids(j)) > pgt) then
        print*,'warning: lat outside(north)',slat(ids(j)),mc(istep),dm(istep)
        slat(ids(j)) = pgt-gnst(1,1)
      elseif(slat(ids(j)) < pgb) then
        print*,'warning: lat outside(south)',slat(ids(j)),mc(istep),dm(istep)
        slat(ids(j)) = pgb+gnst(1,1)
      endif
    endif
   ENDDO
   if(pvi>0) then
      !pgt = (gop(1,1)+(nvnp(1,1)-3)*gnsp(1,1))*180.0/pi
      !pgb = (gop(1,1)+2*gnsp(1,1))*180.0/pi
      pgt = (gop(1,1)+(nvnp(1,1)-1)*gnsp(1,1))*180.0/pi
      pgb = gop(1,1)*180.0/pi
    endif
   DO j=1,nspi
      istep=istep+1
      dm(istep)=dm(istep)*180.0/(pi*(earthr-srad(ids(j))))
      ! bug here, seems very important bug 180*pi-->pi/180
      ! degree-->rad for cos
      dm(istep)=dm(istep)/cos(slat(ids(j))*pi/180)
      if (abs(dm(istep))>0.05) dm(istep) = dm(istep)/abs(dm(istep))*0.05
      slon(ids(j))=mc(istep)+dm(istep)
   if(pvi>0) then
      if (slon(ids(j)) > pgt) then
        print*,'warning: lon outside(right)',slon(ids(j)),mc(istep),dm(istep)
        slon(ids(j)) = pgt-gnsp(1,1)
      elseif(slon(ids(j)) < pgb) then
        print*,'warning: lon outside(left)',slon(ids(j)),mc(istep),dm(istep)
        slon(ids(j)) = pgb+gnsp(1,1)
      endif
    endif
   ENDDO
   DO j=1,nspi
      istep=istep+1
      stp(j)=mc(istep)+dm(istep)
   ENDDO
   OPEN(UNIT=10,FILE=scfile,STATUS='unknown')
   WRITE(10,*)ns
   DO i=1,ns
      WRITE(10,*)lots(i)
      IF(lots(i).EQ.1)WRITE(10,'(a10)')tpath(i)
      idm1=0
      WRITE(10,*)srad(i),slat(i),slon(i)
      WRITE(10,*)nps(i)
      DO j=1,nps(i)
         WRITE(10,*)nsdf(j,i)
         WRITE(10,*)paths(1:2*nsdf(j,i),j,i)
         WRITE(10,*)patht(1:nsdf(j,i),j,i)
      ENDDO
   ENDDO
   CLOSE(10)

!  if(psi==1) then
!  DEALLOCATE(paths,patht,tpath,nsdf,lots,nps,stp)
!   OPEN(UNIT=20,FILE='sourcesallref.in',STATUS='unknown')
!   OPEN(UNIT=10,FILE='sourcesall.in',STATUS='unknown')
!   read(20,*)ns
!   !ALLOCATE(sradr(ns),slatr(ns),slonr(ns))
!   !ALLOCATE(covsr(ns),covst(ns),covsp(ns))
!   ALLOCATE(lots(ns),nps(ns),tpath(ns))
!   ALLOCATE(nsdf(maxs,ns))
!   ALLOCATE(stp(ns))
!   ALLOCATE(paths(maxs,maxp,ns),patht(maxs,maxp,ns))
!   WRITE(10,*)ns
!   DO i=1,ns
!      read(20,*)lots(i),cdum,choosedsrc
!      WRITE(10,*)lots(i),cdum,choosedsrc
!      IF(lots(i).EQ.1)WRITE(10,'(a10)')tpath(i)
!      idm1=0
!      read(20,*)radall,latall,lonall
!      if(choosedsrc>0) then
!      WRITE(10,*)srad(ids(choosedsrc)),slat(ids(choosedsrc)),slon(ids(choosedsrc))
!      else
!      write(10,*)radall,latall,lonall
!      endif
!      read(20,*)nps(i)
!      WRITE(10,*)nps(i)
!      DO j=1,nps(i)
!         read(20,*)nsdf(j,i)
!         read(20,*)paths(1:2*nsdf(j,i),j,i)
!         read(20,*)patht(1:nsdf(j,i),j,i)
!         WRITE(10,*)nsdf(j,i)
!         WRITE(10,*)paths(1:2*nsdf(j,i),j,i)
!         WRITE(10,*)patht(1:nsdf(j,i),j,i)
!      ENDDO
!   ENDDO
!   CLOSE(10)
!   CLOSE(20)
!  endif


   OPEN(UNIT=10,FILE=stfile,STATUS='unknown')
   WRITE(10,*)nspi
   DO i=1,nspi
      WRITE(10,*)stp(i)
   ENDDO
   CLOSE(10)
ENDIF
!
! Final deallocation
!
DEALLOCATE(frech,fcoln,cnfe)
DEALLOCATE(tfrech,tfcoln,tcnfe)
DEALLOCATE(dobs,dmod,cd)
DEALLOCATE(mo,mc,cm,ecmi,dm)
!DEALLOCATE(paths,patht)
!DEALLOCATE(tsid)
!if(inversionScheme /= 1) deallocate(norm_dws)
IF(pvi.EQ.1)THEN
   !DEALLOCATE(veln,velndws)
   DEALLOCATE(veln)
   DEALLOCATE(nvnr,nvnt,nvnp)
   DEALLOCATE(gnsr,gnst,gnsp)
   DEALLOCATE(gor,got,gop)
ENDIF
IF(pii.EQ.1)DEALLOCATE(intn)
IF(psi.EQ.1)DEALLOCATE(srad,slat,slon,tpath,nsdf,lots,nps,stp)
IF(nvgi.GT.0)DEALLOCATE(idvg,idvt)
IF(nigi.GT.0)DEALLOCATE(idig)
IF(nspi.GT.0)DEALLOCATE(ids)
IF(ntels.GT.0)DEALLOCATE(istel,mtmean)
END PROGRAM invert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates the model perturbation using
! the subspace inversion method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE subspace
USE globalp
IMPLICIT NONE
INTEGER :: ii,i,jj,j,k,is1,is2,npc,nits
INTEGER, DIMENSION (4) :: inds,inde
REAL(KIND=i5) :: rsum1
REAL(KIND=i5), DIMENSION (:,:), ALLOCATABLE :: a,mati
REAL(KIND=i5), DIMENSION (:), ALLOCATABLE :: gamma,dtrav,ga
REAL(KIND=i5), DIMENSION (subdim) :: r1mat,r2mat
!
! a = Projection matrix
! gamma = gradient vector
! dtrav = weighted differential traveltime residuals
! ga = Ga, the product of Frechet matrix and vector of a
! mati = matrix for inversion
! r1mat,r2mat = RHS vectors of matrix equation for model perturbation
! npc = number of parameter classes
! inds = Index start for parameter class
! inde = Index end for parameter class
! nits = Number of iterations for subspace vector generation
!
! The first step is to compute the projection matrix a. Allocate
! memory to this array and set to zero.
!
ALLOCATE(a(npi,subdim))
a=0.0
!
! The first subspace dimension is given by the gradient vector in
! model space. Allocate memory to this vector.
!
ALLOCATE(gamma(npi))
!
! Allocate memory to weighted differential traveltimes
!
ALLOCATE(dtrav(ntr))
DO i=1,ntr
   dtrav(i)=(dmod(i)-dobs(i))/cd(i)
ENDDO
!
! Now compute gamma and the first subspace vector. There are
! two options here, depending on whether smoothing is
! applied or not.
!
DO i=1,npi
   rsum1=0.0
   IF(tcnfe(i).GT.tcnfe(i-1))THEN
      DO j=tcnfe(i-1)+1,tcnfe(i)
         rsum1=rsum1+tfrech(j)*dtrav(tfcoln(j))
      ENDDO
   ENDIF
   gamma(i)=rsum1+epsilon*dm(i)*ecmi(i)
ENDDO
!
! Apply second derivative smoothing if required. Note that
! only velocity and interface parameters can be smoothed.
!
IF(asds.EQ.1)THEN
   ALLOCATE(dem(npi),smv(npi))
   smv=dm
   CALL smoothing
   DO i=1,npi
      gamma(i)=gamma(i)+eta*smv(i)
   ENDDO
ENDIF
!
! Determine the number of parameter classes and the index
! of the start and end of each class
!
npc=0
is1=0
IF(nvpi.GT.0)THEN
   npc=npc+1
   is1=is1+1
   inds(npc)=is1
   is1=is1+nvpi-1
   inde(npc)=is1
ENDIF
IF(nipi.GT.0)THEN
   npc=npc+1
   is1=is1+1
   inds(npc)=is1
   is1=is1+nipi-1
   inde(npc)=is1
ENDIF
IF(nspi.GT.0)THEN
   npc=npc+1
   is1=is1+1
   inds(npc)=is1
   is1=is1+3*nspi-1
   inde(npc)=is1
   npc=npc+1
   is1=is1+1
   inds(npc)=is1
   is1=is1+nspi-1
   inde(npc)=is1
ENDIF
!
! Now compute the first npc basis vectors
!
DO i=1,npc
   DO j=inds(i),inde(i)
      a(j,i)=gamma(j)*cm(j)
   ENDDO
ENDDO
!
! Normalize the subspace vectors
!
DO i=1,npc
   rsum1=0.0
   DO j=inds(i),inde(i)
      rsum1=rsum1+a(j,i)**2
   ENDDO
   DO j=inds(i),inde(i)
      a(j,i)=a(j,i)/SQRT(rsum1)
   ENDDO
ENDDO
!
! Now compute the remaining subspace vectors if
! required
!
ALLOCATE(ga(ntr))
IF(subdim.GT.npc)THEN
   nits=npc
   is1=npc
   is2=npc
   DO ii=1,subdim
      DO i=1,nits
!
!        First, calculate Ga weighted by Cd
!
         DO j=1,ntr
            rsum1=0.0
            IF(cnfe(j).GT.cnfe(j-1))THEN
               DO k=cnfe(j-1)+1,cnfe(j)
                  rsum1=rsum1+frech(k)*a(fcoln(k),is2-nits+i)
               ENDDO
            ENDIF
            ga(j)=rsum1/cd(j)
         ENDDO
!
!        Compute the contribution from smoothing if applied
!
         IF(asds.EQ.1)THEN
            smv=a(1:npi,is2-nits+i)
            CALL smoothing
         ENDIF
!
!        Now calculate the new subspace vector
!
         DO jj=1,npc
            is1=is1+1
            IF(is1.GT.subdim)EXIT
            DO j=inds(jj),inde(jj)
               rsum1=0.0
               IF(tcnfe(j).GT.tcnfe(j-1))THEN
                  DO k=tcnfe(j-1)+1,tcnfe(j)
                     rsum1=rsum1+tfrech(k)*ga(tfcoln(k))
                  ENDDO
               ENDIF
               IF(asds.EQ.1)THEN
                  a(j,is1)=rsum1*cm(j)+epsilon*a(j,is2-nits+i)+eta*smv(j)
               ELSE
                  a(j,is1)=rsum1*cm(j)+epsilon*a(j,is2-nits+i)
               ENDIF
            ENDDO
            !
            ! Normalize the subspace vector
            !
            rsum1=0.0
            DO j=inds(jj),inde(jj)
               rsum1=rsum1+a(j,is1)**2
            ENDDO
            DO j=inds(jj),inde(jj)
               a(j,is1)=a(j,is1)/SQRT(rsum1)
            ENDDO
         ENDDO
         IF(is1.GT.subdim)EXIT
      ENDDO
      IF(is1.GT.subdim)EXIT
      nits=nits*npc
      is2=is2+nits
   ENDDO
ENDIF
!
! If subdim is greater than 1, then we need to orthogonalize
! the vectors that make up the projection matrix A to
! avoid interdependence. This can be done using SVD.
!
CALL svdcmp(a)
!
! Now we have our projection matrix A. Construct the matrix that
! is to be inverted.
!
ALLOCATE(mati(subdim,subdim))
mati=0.0
DO i=1,subdim
   DO j=1,ntr
      rsum1=0.0
      IF(cnfe(j).GT.cnfe(j-1))THEN
         DO k=cnfe(j-1)+1,cnfe(j)
            rsum1=rsum1+frech(k)*a(fcoln(k),i)
         ENDDO
      ENDIF
      ga(j)=rsum1/cd(j)
   ENDDO
!
!  Compute contribution from smoothing
!  if required.
! 
   IF(asds.EQ.1)THEN
       smv=a(1:npi,i)
       CALL smoothing
   ENDIF
   DO j=1,npi
      rsum1=0.0
      IF(tcnfe(j).GT.tcnfe(j-1))THEN
         DO k=tcnfe(j-1)+1,tcnfe(j)
            rsum1=rsum1+tfrech(k)*ga(tfcoln(k))
         ENDDO
      ENDIF
      rsum1=rsum1+epsilon*a(j,i)*ecmi(j)
      IF(asds.EQ.1)rsum1=rsum1+eta*smv(j)
      DO k=1,subdim
         mati(k,i)=mati(k,i)+rsum1*a(j,k)
      ENDDO
   ENDDO
ENDDO
!
! Now perform the matrix inversion using LU decomposition
!
CALL luinv(subdim,mati)
!
! Calculate the model perturbation
!
r1mat=0.0
DO i=1,subdim
   DO j=1,npi
      r1mat(i)=r1mat(i)+a(j,i)*gamma(j)
   ENDDO
ENDDO
r2mat=0.0
DO i=1,subdim
   DO j=1,subdim
      r2mat(i)=r2mat(i)+mati(i,j)*r1mat(j)
   ENDDO
ENDDO
dm=0.0
DO i=1,npi
   DO j=1,subdim
      dm(i)=dm(i)+a(i,j)*r2mat(j)
   ENDDO
   dm(i)=-dm(i)
ENDDO
DEALLOCATE(a,gamma,dtrav)
DEALLOCATE(ga,mati)
IF(asds.EQ.1)THEN
   DEALLOCATE(dem,smv)
ENDIF
END SUBROUTINE subspace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine applies the smoothing operator to a vector
! call "smv" and returns the result as "smv".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE smoothing
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,l,m,is1,is2,is3
REAL(KIND=i5) :: rsum1,ddem
!
! is1,is2,is3 = counters for smoothness operators
! ddem = dem premultipled by transform of D
!
dem=0.0
!
! Start with velocity parameters if they exist.
!
is2=0
IF(nvpi.GT.0)THEN
   DO i=1,nvgi
      DO k=1,nvnp(idvg(i),idvt(i))
         DO l=1,nvnt(idvg(i),idvt(i))
            DO m=1,nvnr(idvg(i),idvt(i))
               is2=is2+1
               IF(m.NE.1.AND.m.NE.nvnr(idvg(i),idvt(i)))THEN
                  is1=is2-1
                  is3=is2+1
                  rsum1=smv(is1)-2.0*smv(is2)+smv(is3)
                  dem(is2)=dem(is2)+rsum1
               ENDIF
               IF(l.NE.1.AND.l.NE.nvnt(idvg(i),idvt(i)))THEN
                  is1=is2-nvnr(idvg(i),idvt(i))
                  is3=is2+nvnr(idvg(i),idvt(i))
                  rsum1=smv(is1)-2.0*smv(is2)+smv(is3)
                  dem(is2)=dem(is2)+rsum1
               ENDIF
               IF(k.NE.1.AND.k.NE.nvnp(idvg(i),idvt(i)))THEN
                  is1=is2-nvnr(idvg(i),idvt(i))*nvnt(idvg(i),idvt(i))
                  is3=is2+nvnr(idvg(i),idvt(i))*nvnt(idvg(i),idvt(i))
                  rsum1=smv(is1)-2.0*smv(is2)+smv(is3)
                  dem(is2)=dem(is2)+rsum1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
! Now add interface parameters if they exist.
!
IF(nipi.GT.0)THEN
   DO i=1,nigi
      DO j=1,ninp
         DO k=1,nint
            is2=is2+1
            IF(k.NE.1.AND.k.NE.nint)THEN
               is1=is2-1
               is3=is2+1
               rsum1=smv(is1)-2.0*smv(is2)+smv(is3)
               dem(is2)=dem(is2)+rsum1
            ENDIF
            IF(j.NE.1.AND.j.NE.ninp)THEN
               is1=is2-nint
               is3=is2+nint
               rsum1=smv(is1)-2.0*smv(is2)+smv(is3)
               dem(is2)=dem(is2)+rsum1
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF
smv=0.0
!
!  Now repeat the loop to premultiply by the transform of D
!
is2=0
IF(nvpi.GT.0)THEN
   DO i=1,nvgi
      DO k=1,nvnp(idvg(i),idvt(i))
         DO l=1,nvnt(idvg(i),idvt(i))
            DO m=1,nvnr(idvg(i),idvt(i))
               ddem=0.0
               is2=is2+1
               IF(m.NE.1.AND.m.NE.nvnr(idvg(i),idvt(i)))THEN
                  ddem=ddem-2.0*dem(is2)
               ENDIF
               IF(m-2.GE.1)THEN
                  is1=is2-1
                  ddem=ddem+dem(is1)
               ENDIF
               IF(m+2.LE.nvnr(idvg(i),idvt(i)))THEN
                  is3=is2+1
                  ddem=ddem+dem(is3)
               ENDIF
               IF(l.NE.1.AND.l.NE.nvnt(idvg(i),idvt(i)))THEN
                  ddem=ddem-2.0*dem(is2)
               ENDIF
               IF(l-2.GE.1)THEN
                  is1=is2-nvnr(idvg(i),idvt(i))
                  ddem=ddem+dem(is1)
               ENDIF
               IF(l+2.LE.nvnt(idvg(i),idvt(i)))THEN
                  is3=is2+nvnr(idvg(i),idvt(i))
                  ddem=ddem+dem(is3)
               ENDIF
               IF(k.NE.1.AND.k.NE.nvnp(idvg(i),idvt(i)))THEN
                  ddem=ddem-2.0*dem(is2)
               ENDIF
               IF(k-2.GE.1)THEN
                  is1=is2-nvnr(idvg(i),idvt(i))*nvnt(idvg(i),idvt(i))
                  ddem=ddem+dem(is1)
               ENDIF
               IF(k+2.LE.nvnp(idvg(i),idvt(i)))THEN
                  is3=is2+nvnr(idvg(i),idvt(i))*nvnt(idvg(i),idvt(i))
                  ddem=ddem+dem(is3)
               ENDIF
               smv(is2)=etav*ddem
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF
IF(nipi.GT.0)THEN
   DO i=1,nigi
      DO j=1,ninp
         DO k=1,nint
            ddem=0.0
            is2=is2+1
            IF(k.NE.1.AND.k.NE.nint)THEN
               ddem=ddem-2.0*dem(is2)
            ENDIF
            IF(k-2.GE.1)THEN
               is1=is2-1
               ddem=ddem+dem(is1)
            ENDIF
            IF(k+2.LE.nint)THEN
               is3=is2+1
               ddem=ddem+dem(is3)
            ENDIF
            IF(j.NE.1.AND.j.NE.ninp)THEN
               ddem=ddem-2.0*dem(is2)
            ENDIF
            IF(j-2.GE.1)THEN
               is1=is2-nint
               ddem=ddem+dem(is1)
            ENDIF
            IF(j+2.LE.ninp)THEN
               is3=is2+nint
               ddem=ddem+dem(is3)
            ENDIF
            smv(is2)=etai*ddem
         ENDDO
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE smoothing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine (adapted from Numerical Recipes) computes
! the singular value decomposition of the projection
! matrix a to produce an orthonormal basis. The dimension of
! the subspace is automatically reduced if the subspace does
! not span subdim.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE svdcmp(a)
USE globalp
IMPLICIT NONE
INTEGER :: m,n
INTEGER :: i,its,j,jj,k,l,nm,isw,isum
INTEGER, DIMENSION (subdim) :: wnz
INTEGER, PARAMETER :: nmax=500
REAL(KIND=i5), DIMENSION (npi,subdim) :: a
REAL(KIND=i5), DIMENSION (subdim,subdim) :: v
REAL(KIND=i5), DIMENSION (subdim) :: w
REAL(KIND=i5), DIMENSION (nmax) :: rv1
REAL(KIND=i5) :: anorm,c,f,g,h,s,scale,x,y,z
REAL(KIND=i5), EXTERNAL :: pythag
REAL(KIND=i5), PARAMETER :: wtol=1.0e-7
m=npi
n=subdim
g=0.0
scale=0.0
anorm=0.0
nm = 0
DO i=1,n
   l=i+1
   rv1(i)=scale*g
   g=0.0
   s=0.0
   scale=0.0
   If(i.LE.m)THEN
      DO k=i,m
         scale=scale+ABS(a(k,i))
      ENDDO
      If(scale.NE.0.0)THEN
         DO k=i,m
            a(k,i)=a(k,i)/scale
            s=s+a(k,i)*a(k,i)
         ENDDO
         f=a(i,i)
         g=-SIGN(SQRT(s),f)
         h=f*g-s
         a(i,i)=f-g
         DO j=l,n
            s=0.0
            DO k=i,m
               s=s+a(k,i)*a(k,j)
            ENDDO
            f=s/h
            DO k=i,m
               a(k,j)=a(k,j)+f*a(k,i)
            ENDDO
         ENDDO
         DO k=i,m
            a(k,i)=scale*a(k,i)
         ENDDO
      ENDIF
   ENDIF
   w(i)=scale*g
   g=0.0
   s=0.0
   scale=0.0
   IF((i.LE.m).AND.(i.NE.n))THEN
      DO k=l,n
         scale=scale+ABS(a(i,k))
      ENDDO
      IF(scale.NE.0.0)THEN
         DO k=l,n
            a(i,k)=a(i,k)/scale
            s=s+a(i,k)*a(i,k)
         ENDDO
         f=a(i,l)
         g=-SIGN(SQRT(s),f)
         h=f*g-s
         a(i,l)=f-g
         DO k=l,n
            rv1(k)=a(i,k)/h
         ENDDO
         DO j=l,m
            s=0.0
            DO k=l,n
               s=s+a(j,k)*a(i,k)
            ENDDO
            DO k=l,n
               a(j,k)=a(j,k)+s*rv1(k)
            ENDDO
         ENDDO
         DO k=l,n
            a(i,k)=scale*a(i,k)
         ENDDO
      ENDIF
   ENDIF
   anorm=MAX(anorm,(ABS(w(i))+abs(rv1(i))))
ENDDO
DO i=n,1,-1
   IF(i.LT.n)THEN
      IF(g.NE.0.0)THEN
         DO j=l,n
            v(j,i)=(a(i,j)/a(i,l))/g
         ENDDO
         DO j=l,n
            s=0.0
            DO k=l,n
               s=s+a(i,k)*v(k,j)
            ENDDO
            DO k=l,n
               v(k,j)=v(k,j)+s*v(k,i)
            ENDDO
         ENDDO
      ENDIF
      DO j=l,n
         v(i,j)=0.0
         v(j,i)=0.0
      ENDDO
   ENDIF
   v(i,i)=1.0
   g=rv1(i)
   l=i
ENDDO
DO i=MIN(m,n),1,-1
   l=i+1
   g=w(i)
   DO j=l,n
      a(i,j)=0.0
   ENDDO
   IF(g.NE.0.0)THEN
      g=1.0/g
      DO j=l,n
         s=0.0
         DO k=l,m
            s=s+a(k,i)*a(k,j)
         ENDDO
         f=(s/a(i,i))*g
         DO k=i,m
            a(k,j)=a(k,j)+f*a(k,i)
         ENDDO
      ENDDO
      DO j=i,m
         a(j,i)=a(j,i)*g
      ENDDO
   ELSE
      DO j=i,m
         a(j,i)=0.0
      ENDDO
   ENDIF
   a(i,i)=a(i,i)+1.0
ENDDO
DO k=n,1,-1
   DO its=1,30
      isw=0
      DO l=k,1,-1
         nm=l-1
         IF((ABS(rv1(l))+anorm).EQ.anorm)THEN
            isw=1
            EXIT
         ENDIF
         IF((ABS(w(nm))+anorm).EQ.anorm)EXIT
      ENDDO
      IF(isw.EQ.0)THEN
         c=0.0
         s=1.0
         DO i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            IF((ABS(f)+anorm).EQ.anorm)EXIT
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c=g*h
            s=-f*h
            DO j=1,m
               y=a(j,nm)
               z=a(j,i)
               a(j,nm)=y*c+z*s
               a(j,i)=-y*s+z*c
            ENDDO
         ENDDO
      ENDIF
      z=w(k)
      IF(l.EQ.k)THEN
         IF(z.LT.0.0)THEN
            w(k)=-z
            DO j=1,n
               v(j,k)=-v(j,k)
            ENDDO
         ENDIF
         EXIT
      ENDIF
      IF(its.EQ.30)print*, 'No convergence in svdcmp!'
      x=w(l)
      nm=k-1
      y=w(nm)
      g=rv1(nm)
      h=rv1(k)
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
      g=pythag(f,1.0)
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x
      c=1.0
      s=1.0
      DO j=l,nm
         i=j+1
         g=rv1(i)
         y=w(i)
         h=s*g
         g=c*g
         z=pythag(f,h)
         rv1(j)=z
         c=f/z
         s=h/z
         f=x*c+g*s
         g=-x*s+g*c
         h=y*s
         y=y*c
         DO jj=1,n
            x=v(jj,j)
            z=v(jj,i)
            v(jj,j)=x*c+z*s
            v(jj,i)=-x*s+z*c
         ENDDO
         z=pythag(f,h)
         w(j)=z
         IF(z.NE.0.0)THEN
            z=1.0/z
            c=f*z
            s=h*z
         ENDIF
         f=c*g+s*y
         x=-s*g+c*y
         DO jj=1,m
            y=a(jj,j)
            z=a(jj,i)
            a(jj,j)=y*c+z*s
            a(jj,i)=-y*s+z*c
         ENDDO
      ENDDO
      rv1(l)=0.0
      rv1(k)=f
      w(k)=x
   ENDDO
ENDDO
!
! Finally, eliminate any vectors with very small w's.
!
isum=0
DO i=1,n
   IF(ABS(w(i)).GE.wtol)THEN
      isum=isum+1
      wnz(isum)=i
   ENDIF
ENDDO
IF(isum.LT.subdim)THEN
   WRITE(6,*)'subspace dimension decreased from ',subdim,' to ',isum
   WRITE(6,*)'PROGRAM invert successfully completed!!'
   WRITE(6,*)'------------------------------------------------------------------'
   DO i=1,isum
      DO j=1,npi
         a(j,i)=a(j,wnz(i))
      ENDDO
   ENDDO
   subdim=isum
ENDIF
END SUBROUTINE svdcmp

REAL FUNCTION pythag(a,b)
USE globalp
IMPLICIT NONE
REAL(KIND=i5), INTENT(IN) :: a,b
REAL(KIND=i5) :: absa,absb
absa=ABS(a)
absb=ABS(b)
IF(absa.GT.absb)THEN
   pythag=absa*SQRT(1.0+(absb/absa)**2)
ELSE
   IF(absb.EQ.0.0)THEN
      pythag=0.0
   ELSE
      pythag=absb*SQRT(1.0+(absa/absb)**2)
   ENDIF
ENDIF
END FUNCTION pythag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine (adapted from Numerical Recipes) inverts
! the given matrix a using LU decomposition.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE luinv(n,a)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,n
INTEGER, DIMENSION(n) :: indx
REAL(KIND=i5), DIMENSION (n,n) :: a,y
DO i=1,n
   DO j=1,n
      y(i,j)=0.0
   ENDDO
   y(i,i)=1.0
ENDDO
CALL ludcmp(a,n,indx)
DO j=1,n
   CALL lubksb(a,n,indx,y(1,j))
ENDDO
DO i=1,n
   DO j=1,n
      a(i,j)=y(i,j)
   ENDDO
ENDDO
END SUBROUTINE luinv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine (adapted from Numerical Recipes) replaces
! a given matrix by its LU decomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ludcmp(a,n,indx)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,n,imax
INTEGER, DIMENSION(n) :: indx
REAL(KIND=i5), DIMENSION (n) :: vv
REAL(KIND=i5), DIMENSION (n,n) :: a
REAL(KIND=i5) :: aamax,sum,dum
REAL(KIND=i5), PARAMETER :: tiny=1.0e-20
imax = 0
DO i=1,n
   aamax=0.0
   DO j=1,n
      IF(ABS(a(i,j)).GT.aamax)aamax=ABS(a(i,j))
   enddo
   IF(aamax.EQ.0.0)THEN
      WRITE(6,*)'Singular matrix. Try C.G. method'
      STOP
   ENDIF
   vv(i)=1.0/aamax
ENDDO
DO j=1,n
   DO i=1,j-1
      sum=a(i,j)
      DO k=1,i-1
         sum=sum-a(i,k)*a(k,j)
      ENDDO
      a(i,j)=sum
   ENDDO
   aamax=0.0
   DO i=j,n
      sum=a(i,j)
      DO k=1,j-1
         sum=sum-a(i,k)*a(k,j)
      ENDDO
      a(i,j)=sum
      dum=vv(i)*ABS(sum)
      IF(dum.GE.aamax)THEN
         imax=i
         aamax=dum
      ENDIF
   ENDDO
   IF(j.NE.imax)THEN
      DO k=1,n
         dum=a(imax,k)
         a(imax,k)=a(j,k)
         a(j,k)=dum
      ENDDO
      vv(imax)=vv(j)
   ENDIF
   indx(j)=imax
   IF(a(j,j).EQ.0.0)a(j,j)=tiny
   IF(j.NE.n)THEN
      dum=1.0/a(j,j)
      DO i=j+1,n
         a(i,j)=a(i,j)*dum
      ENDDO
   ENDIF
ENDDO
END SUBROUTINE ludcmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine (adapted from Numerical Recipes) solves the
! set of n linear equations AX=B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lubksb(a,n,indx,b)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,n,ii,ll
INTEGER, DIMENSION (n) :: indx
REAL(KIND=i5), DIMENSION (n) :: b
REAL(KIND=i5), DIMENSION (n,n) :: a
REAL(KIND=i5) :: sum
ii=0
DO i=1,n
   ll=indx(i)
   sum=b(ll)
   b(ll)=b(i)
   IF(ii.NE.0)THEN
      DO j=ii,i-1
         sum=sum-a(i,j)*b(j)
      ENDDO
   ELSE IF(sum.NE.0.0)THEN
      ii=1
   ENDIF
   b(i)=sum
ENDDO
DO i=n,1,-1
   sum=b(i)
   DO j=i+1,n
      sum=sum-a(i,j)*b(j)
   ENDDO
   b(i)=sum/a(i,i)
ENDDO
END SUBROUTINE lubksb
