program subspaceproj
  use lsmrModule, only:lsmr
  use netcdf
  !use omp_lib
  implicit none
  ! variable definition for lsmr
  real,parameter::pi=3.141592654,earthr=6371
  integer*8 leniw,lenrw
  integer, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
  integer:: weightornot
  real(kind=i5)::spfrac,ftol
  integer:: nloc
  real(kind=i5) atol,btol
  real(kind=i5) conlim
  integer istop
  integer itnlim
  integer nout
  integer itn
  real(kind=i5) damp,acond,anorm,arnorm,rnorm,xnorm
  integer localSize
  real(kind=i5),allocatable,dimension(:)::nzero
  integer,allocatable,dimension(:)::nzero_id,nrow
  integer,allocatable,dimension(:)::col,iw
  real*4,allocatable,dimension(:)::rw
  real(kind=i5),allocatable,dimension(:)::dres,weight,dobs,dsyn,phase
  integer*4 linsize(2)
  integer :: ncid,nzeroid,nrowid,nonid
  integer :: nzerodimid,nrowdimid
  integer nnfd

  integer nd,start
  integer nzid
  integer ii
  real,dimension(:),allocatable::xunknown
  character(len=130) filename
  integer nnzero
  integer reclen
  integer jj,kk,mm
  integer checkstat
  integer*8 jstep
  real damploc
  real rad,lat,lon,drad,dlat,dlon


integer idm,idm1,idm2,idm3,idm4,nloc_sub,ndata,sloc
real fdm6,fdm5
integer,dimension(:),allocatable::sourceidx
real,dimension(:),allocatable::loc
integer irow

  open(10,file='VoroTomo.in')
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) nloc
  read(10,*) 
  read(10,*) 
  read(10,*) damploc
  read(10,*) 
  read(10,*) weightornot
  close(10)
  nloc = nloc*4
  !weightornot = 0
  call check(nf90_open('frechet.nc',nf90_nowrite,ncid))
  call check(nf90_inq_dimid(ncid,"nzero",nzerodimid))
  call check(nf90_inq_dimid(ncid,'nrow',nrowdimid))
  call check(nf90_inquire_dimension(ncid,nzerodimid,len = nnfd ))
  call check(nf90_inquire_dimension(ncid,nrowdimid,len = nd ))
  call check(nf90_inq_varid(ncid,"Non_value",nzeroid))
  call check(nf90_inq_varid(ncid,"Non_row",nrowid))
  call check(nf90_inq_varid(ncid,"Non_id",nonid))
  allocate(nzero(nnfd),nrow(nd),nzero_id(nnfd),stat=checkstat)
  if(checkstat>0) stop 'error allocating iw'
  call check(nf90_get_var(ncid,nzeroid,nzero))
  call check(nf90_get_var(ncid,nrowid,nrow))
  call check(nf90_get_var(ncid,nonid,nzero_id))
  call check(nf90_close(ncid))

  allocate(rw(nnfd+nloc),stat=checkstat)
  if(checkstat>0) stop 'error allocating rw'
  allocate(col(nnfd+nloc),iw((nnfd+nloc)*2),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'
  allocate(dres(nd+nloc),dobs(nd),dsyn(nd),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'
  allocate(xunknown(nloc))

  nnzero = int(nd*4)

  !read residual
  dres = 0
  open(10,file='otimes.dat')
  read(10,*)nd
  allocate(weight(nd),phase(nd),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'
  do ii = 1,nd
    read(10,*) idm,idm,phase(ii),idm,dobs(ii),weight(ii)
  enddo
  close(10)
  open(10,file='mtimes.dat')
  do ii = 1,nd
    read(10,*) idm,idm,idm,idm,dsyn(ii)
  enddo
  close(10)

  !iw = 0
  print*,'reading body wave data begin...'
  rw = 0.
  jstep = 0
  irow = 0
  do kk = 1,nd
    start = sum(nrow(1:kk-1))
    do ii = 1,nrow(kk)
      jstep = jstep + 1
      rw(jstep) = nzero(start+ii)*weight(kk)
      iw(jstep) = kk
      col(jstep) = nzero_id(start+ii)
    enddo
  enddo
  print*,'finishing reading. no of nonzeros in the original G for body data: ',jstep
  nzid = jstep

 do ii = 1,nd
    !dres(ii)=dobs(ii)-dsyn(ii)
    dres(ii)=(dsyn(ii)-dobs(ii))*weight(ii)
  enddo
  deallocate(weight)
  print*,'data residual before weighting: ',sum(dres(1:nd)**2)**0.5

    !if(weightornot==1) then
    !allocate(weight(nd),stat=checkstat)
  !if(checkstat>0) stop 'error allocating col'
  !weight(1:nd) = 1.0/(1+0.05*exp(dres(1:nd)**2*0.1))
  !  do jj = 1,nd!jstep
  !    start = sum(nrow(1:jj-1))
  !    !do ii = 1,nrow(jj)
  !    !rw(start+ii) = rw(start+ii)*weight(jj)
  !    !enddo
  !    if(phase(jj)==2) then
  !            rw(start+1:start+nrow(jj))=0
  !    endif
  !  enddo
  !  do ii = 1,nd
  !    !dres(ii) = dres(ii)*weight(ii) 
  !    if(phase(ii)==2) then
  !            dres(ii)=0
  !    endif
  !  enddo
    !deallocate(weight,phase)
    print*,'finishing weighting,residual norm: ',sum(dres(1:nd)**2)**0.5
  !endif

   do ii = 1,nloc
   nzid = nzid+1
   iw(nzid) = nd+ii
   rw(nzid) = damploc
   if (mod(ii,4) == 1) rw(nzid) = damploc*5.0
   col(nzid) = ii
   dres(nd+ii) = 0
   enddo

    print*,nzid,size(col)
    iw(nzid+1:2*nzid) = col
    leniw = 2*nzid
    lenrw = nzid

    xunknown = 0
    atol = 1e-3
    btol = 1e-4
    conlim = 50
    itnlim = 100
    istop = 0
    anorm = 0.0
    acond = 0.0
    arnorm = 0.0
    xnorm = 0.0
    localSize = 0
    damp = 0!damploc
    ! using lsmr to solve for the projection coefficients
    print*, 'LSMR beginning ...'

    nout = 63
    open(nout,file='lsmrout_reloc.txt')
    print*,nd,nloc

    call LSMR(nd+nloc, nloc, leniw, lenrw,iw,rw,dres,damp,&
      atol, btol, conlim, itnlim, localSize,nout,&
      xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)

    write(*,*) 'min and max location perturbation: ',minval(xunknown(1:nloc)),&
            maxval(xunknown(1:nloc))

open(10,file='otimes.dat')
read(10,*)ndata
allocate(sourceidx(ndata))
do ii = 1,ndata
read(10,*) idm1,idm2,idm3,idm4,fdm5,fdm6,sourceidx(ii)
enddo
close(10)
open(10,file='receiversref.in')
open(11,file='receivers.in')
read(10,*) ndata
write(11,*) ndata
do ii = 1,ndata
read(10,*) rad,lat,lon
drad = xunknown(4*sourceidx(ii)+1)
if (abs(drad) > 1.0) drad = drad/abs(drad)*1.0
rad = rad-drad
dlat = xunknown(4*sourceidx(ii)+2)*180/(pi*(earthr-rad))
if (abs(dlat) > 0.05) dlat = dlat/abs(dlat)*0.05
lat = lat+dlat
!lat = lat+xunknown(ns+sourceidx(ii))*180/(pi*(earthr-rad))
dlon = xunknown(4*sourceidx(ii)+3)*180/(pi*(earthr-rad))
dlon=dlon/cos(lat*pi/180)
if (abs(dlon) > 0.05) dlon = dlon/abs(dlon)*0.05
lon = lon+dlon
!lon = lon+xunknown(2*ns+sourceidx(ii))*180/(pi*(earthr-rad))
write(11,*) rad,lat,lon
read(10,*) idm1
read(10,*) idm2
read(10,*) idm3
write(11,*) idm1
write(11,*) idm2
write(11,*) idm3
enddo
close(10)
close(11)
!print*,'max rad,lat,loc',maxval(abs(xunknown(1:ns))),maxval(abs(xunknown(ns+1:2*ns))),maxval(abs(xunknown(2*ns+1:3*ns)))


  deallocate(rw,dres)
  deallocate(xunknown,sourceidx)
  deallocate(nzero_id,nzero)
contains
subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
end subroutine check  

  end
