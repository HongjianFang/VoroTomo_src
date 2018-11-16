program subspaceproj
  use lsmrModule, only:lsmr
  use netcdf
  !use omp_lib
  implicit none
  ! variable definition for lsmr
  integer*8 leniw,lenrw
  integer, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
  integer:: weightornot
  real(kind=i5)::spfrac,ftol
  integer:: mdim,mdim2,nlat,nlon,nrad,nloc
  integer:: joint,veltype
  real::damploc
  integer nchoose,nchoose_surf
  real(kind=i5) weight_surf
  integer,allocatable,dimension(:):: choose,eqidx,choose_surf,phasebody,phaseboth
  real(kind=i5),allocatable,dimension(:,:)::Gh 
  integer irow
  real(kind=i5) atol,btol
  real(kind=i5) conlim
  integer istop
  integer itnlim
  integer nout
  integer itn
  real(kind=i5) damp,acond,anorm,arnorm,rnorm,xnorm
  integer localSize
  real(kind=i5),allocatable,dimension(:)::nzero,nzerosurf
  integer,allocatable,dimension(:)::nzero_id,nzero_idsurf,nrow,nrowsurf,nrow_choose
  !integer,allocatable,dimension(:)::iw
  integer,allocatable,dimension(:)::col
  real*4,allocatable,dimension(:)::rw
  integer,allocatable,dimension(:)::iw_p
  real(kind=i5),allocatable,dimension(:)::rw_p
  integer,allocatable,dimension(:)::iwgp,colgp
  real(kind=i5),allocatable,dimension(:)::rwgp
  real(kind=i5),allocatable,dimension(:)::dres,weight,dobs,dsyn
  !real(kind=i5),allocatable,dimension(:,:)::yy,yys
  integer*4 linsize(2)

  integer nd,start
  integer ii,inet
  integer nzid,nz
  integer subrow,subrow2
  real,dimension(:),allocatable::norm,xunknown
  real,allocatable,dimension(:):: x,rrow,rrow_vel!,xaug
  character(len=130) filename
  integer nnzero
  integer reclen
  !character(len=8) numthreads
  !integer numth,startidx
  !integer thread_num,max_thread
  integer jj,kk,mm
  integer coltmp
  integer iseed(4)
  integer :: ncid,nzeroid,nrowid,nonid
  integer :: nzerodimid,nrowdimid
  integer :: ndbody,ndsurf,nnfdsurf,nnfdbody,nnfd,npara
  integer checkstat
  integer*8 jstep
  integer idm

  character(len=32) arg

  open(10,file='VoroTomo.in')
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) nlat,nlon,nrad
  read(10,*) subrow
  read(10,*) 
  read(10,*) nloc
  read(10,*) veltype
  read(10,*) joint
  read(10,*) damploc
  read(10,*) spfrac,ftol
  read(10,*) weightornot
  close(10)
  !weightornot = 0
  call check(nf90_open('frechet.nc',nf90_nowrite,ncid))
  call check(nf90_inq_dimid(ncid,"nzero",nzerodimid))
  call check(nf90_inq_dimid(ncid,'nrow',nrowdimid))
  call check(nf90_inquire_dimension(ncid,nzerodimid,len = nnfdbody ))
  call check(nf90_inquire_dimension(ncid,nrowdimid,len = ndbody ))
  call check(nf90_inq_varid(ncid,"Non_value",nzeroid))
  call check(nf90_inq_varid(ncid,"Non_row",nrowid))
  call check(nf90_inq_varid(ncid,"Non_id",nonid))
  allocate(nzero(nnfdbody),nrow(ndbody),nzero_id(nnfdbody),stat=checkstat)
  if(checkstat>0) stop 'error allocating iw'
  call check(nf90_get_var(ncid,nzeroid,nzero))
  call check(nf90_get_var(ncid,nrowid,nrow))
  call check(nf90_get_var(ncid,nonid,nzero_id))
  call check(nf90_close(ncid))

  nnfd = nnfdbody
  nd = ndbody
  !allocate(iw(2*nnfd),stat=checkstat)
  !if(checkstat>0) stop 'error allocating iw'
  allocate(rw(nnfd),stat=checkstat)
  if(checkstat>0) stop 'error allocating rw'
  allocate(col(nnfd),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'
  allocate(dres(nd),dobs(nd),dsyn(nd),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'

    nnzero = int(nd*4)

  !read residual
  dres = 0
  open(10,file='otimes.dat')
  read(10,*)ndbody
  do ii = 1,ndbody
    read(10,*) idm,idm,idm,idm,dobs(ii)
  enddo
  close(10)
  open(10,file='mtimes.dat')
  do ii = 1,ndbody
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
    do ii = 1,nrow(jj)
      jstep = jstep + 1
      rw(jstep) = nzero(start+ii)
      col(jstep) = nzero_id(start+ii)
    enddo
  enddo
  print*,'finishing reading. no of nonzeros in the original G for body data: ',jstep

 do ii = 1,nchoose
    dres(ii)=dobs(choose(ii))-dsyn(choose(ii))
  enddo

    if(weightornot==1) then
    allocate(weight(nd),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'
  weight(1:nd) = 1.0/(1+0.05*exp(dres(1:nd)**2*0.1))
    do jj = 1,nd!jstep
      start = sum(nrow(1:jj-1))
      do ii = 1,nrow(jj)
      rw(start+ii) = rw(start+ii)*weight(jj)
    enddo
    enddo
    do ii = 1,nd
      dres(ii) = dres(ii)*weight(ii) 
    enddo
    deallocate(weight)
    print*,'finishing weighting'
  endif

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
    damp = 0.0
    ! using lsmr to solve for the projection coefficients
    print*, 'LSMR beginning ...'

    nout = 63
    open(nout,file='lsmrout_reloc.txt')

    call LSMR(nd, nloc, leniw, lenrw,iw,rw,dres,damp,&
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
drad = locmean(sourceidx(ii)+1)
if (abs(drad) > 1.0) drad = drad/abs(drad)*1.0
rad = rad-drad!locmean(sourceidx(ii))
dlat = locmean(ns+sourceidx(ii)+1)*180/(pi*(earthr-rad))
if (abs(dlat) > 0.05) dlat = dlat/abs(dlat)*0.05
lat = lat+dlat
!lat = lat+locmean(ns+sourceidx(ii))*180/(pi*(earthr-rad))
dlon = locmean(2*ns+sourceidx(ii)+1)*180/(pi*(earthr-rad))
dlon=dlon/cos(lat*pi/180)
if (abs(dlon) > 0.05) dlon = dlon/abs(dlon)*0.05
lon = lon+dlon
!lon = lon+locmean(2*ns+sourceidx(ii))*180/(pi*(earthr-rad))
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
print*,'max rad,lat,loc',maxval(abs(locmean(1:ns))),maxval(abs(locmean(ns+1:2*ns))),maxval(abs(locmean(2*ns+1:3*ns)))


    ! write out results

    write(filename,'("./tempdata/xloc",i0,"fmtomo.bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*nlocnew)
    write(10,rec=1) xunknown(subrow2+1:subrow2+nlocnew)
    close(10)


  deallocate(rw,dres)
  deallocate(xunknown)
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
