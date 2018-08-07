program subspaceproj
  use lsmrModule, only:lsmr
  use netcdf
  use omp_lib
  implicit none
  ! variable definition for lsmr
  integer*8 leniw,lenrw
  integer leniw_p,lenrw_p
  integer*8 leniwgp,lenrwgp
  integer, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
  integer,parameter:: subrow = 7200,numdir = 1
  integer,parameter:: weightornot=1,synthetic = 0
  real(kind=i5),parameter::spfrac= 0.1,ftol=0.1
  integer,parameter:: mdim = 1944320 ,nloc = 7664
  real(kind=i5),parameter:: unweightlim = 10
  real(kind=i5) atol,btol
  real(kind=i5) conlim
  integer istop
  integer itnlim
  integer nout
  integer itn
  real(kind=i5) damp,acond,anorm,arnorm,rnorm,xnorm
  integer localSize
  real(kind=i5),allocatable,dimension(:)::nzero,nzerosurf
  integer,allocatable,dimension(:)::nzero_id,nzero_idsurf,nrow,nrowsurf
  integer,allocatable,dimension(:)::iw
  integer,allocatable,dimension(:)::col
  real*4,allocatable,dimension(:)::rw
  integer,allocatable,dimension(:)::iw_p
  real(kind=i5),allocatable,dimension(:)::rw_p
  integer,allocatable,dimension(:)::iwgp,colgp
  real(kind=i5),allocatable,dimension(:)::rwgp
  real(kind=i5),allocatable,dimension(:)::dres,y,weight,dobs,dsyn
  real(kind=i5),allocatable,dimension(:,:)::yy
  integer*4 linsize(2)

  integer nd,start
  integer ii,inet
  integer,parameter:: nnets = 20 
  integer nzid,nz
  real pkrow(subrow)
  real xunknown(subrow+nloc)
  real,allocatable,dimension(:):: x,xaug,xtrue
  character(len=130) filename
  integer nnzero
  integer reclen
  character(len=8) numthreads
  integer numth,startidx
  integer thread_num,max_thread
  integer jj,kk,mm
  integer coltmp,dntmp,dnold,ntmp
  integer iseed(4)
  integer :: ncid,nzeroid,nrowid,nonid
  integer :: nzerodimid,nrowdimid
  integer :: ndbody,ndsurf,nnfdsurf,nnfdbody,nnfd,npara
  integer checkstat
  integer*8 jstep
  integer idm

  CALL get_environment_variable("OMP_NUM_THREADS",numthreads)
  if (numthreads=='') then
    numth = 1
    print*,' number of threads: ', numth
  else
    read(numthreads,'(i2)') numth
    print*,' number of threads: ', numth
  endif
  if (mod(subrow,numth)/=0) then 
    print*, 'num of threads must be divisible for subrow'
    stop
  endif

! read netcdf files

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

  nnfdsurf = 0
  ndsurf = 0
  ! read surf
!  call check(nf90_open("frechetsurf.nc",nf90_nowrite,ncid))
!  call check(nf90_inq_dimid(ncid,"nzero",nzerodimid))
!  call check(nf90_inq_dimid(ncid,'nrow',nrowdimid))
!  call check(nf90_inquire_dimension(ncid,nzerodimid,len = nnfdsurf ))
!  call check(nf90_inquire_dimension(ncid,nrowdimid,len = ndsurf ))
!  call check(nf90_inq_varid(ncid,"Non_value",nzeroid))
!  call check(nf90_inq_varid(ncid,"Non_row",nrowid))
!  call check(nf90_inq_varid(ncid,"Non_id",nonid))
!  allocate(nzerosurf(nnfdsurf),nrowsurf(ndsurf),nzero_idsurf(nnfdsurf),stat=checkstat)
!  if(checkstat>0) stop 'error allocating memnetcdf'
!  call check(nf90_get_var(ncid,nzeroid,nzerosurf))
!  call check(nf90_get_var(ncid,nrowid,nrowsurf))
!  call check(nf90_get_var(ncid,nonid,nzero_idsurf))
!  call check(nf90_close(ncid))
!
!
  ! read surf over
  nnfd = nnfdbody+nnfdsurf
  nd = ndbody+ndsurf
  allocate(iw(2*nnfd),stat=checkstat)
  if(checkstat>0) stop 'error allocating iw'
  allocate(rw(nnfd),stat=checkstat)
  if(checkstat>0) stop 'error allocating rw'
  allocate(col(nnfd),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'

  iw = 0
  rw = 0.
  jstep = 0
  do jj = 1,ndbody
    do ii = 1,nrow(jj)
      jstep = jstep + 1
      iw(jstep) = jj
      rw(jstep) = nzero(jstep)
      col(jstep) = nzero_id(jstep)
    enddo
  enddo
  do jj = 1,ndsurf
    do ii = 1,nrowsurf(jj)
      jstep = jstep + 1
      iw(jstep) = ndbody+jj
      rw(jstep) = nzerosurf(jstep-nnfdbody)
      col(jstep) = nzero_idsurf(jstep-nnfdbody)
    enddo
  enddo
  do ii = 1,jstep
    iw(jstep+ii) = col(ii)
  enddo
  lenrw = jstep
  leniw = 2*jstep
  deallocate(nzero_id,nzero_idsurf,nzerosurf,nzero)

  !read residual
  allocate(dres(nd),dobs(nd),dsyn(nd))
  dres = 0
  open(10,file='otimes.dat')
  read(10,*)ndbody
  do ii = 1,nd
    read(10,*) idm,idm,idm,idm,dobs(ii)
  enddo
  close(10)
  open(10,file='mtimes.dat')
  do ii = 1,nd
    read(10,*) idm,idm,idm,idm,dsyn(ii)
  enddo
  close(10)
  do ii = 1,nd
    dres(ii)=dobs(ii)-dsyn(ii)
  enddo

  if(weightornot==1) then
    allocate(weight(nd))
  weight = 1.0/(1+0.05*exp(dres**2*0.1))
      do ii = 1,nnzero
      rw(ii) = rw(ii)*weight(iw(ii))
    enddo
    do ii = 1,nd
      dres(ii) = dres(ii)*weight(ii) 
    enddo
    deallocate(weight)
  endif
 
  allocate(yy(nd,numth))
  allocate(x(mdim),xtrue(mdim),xaug(mdim+nloc))
  if (synthetic == 1) then
    open(10,file='truemod.bin',form='unformatted',access='direct',recl=4*mdim)
    read(10,rec=1) xtrue 
    close(10)
    xaug = 0
    xaug(1:mdim) = xtrue
    y = 0
    call aprod(1,nd,mdim+nloc,xaug,y,leniw,lenrw,iw,rw)
    dres(1:nd) = y(1:nd)
    do  ii=1,nd
      dres(ii)=dres(ii)*weight(ii)
      if(isnan(dres(ii))) dres(ii)=0
    enddo
    iseed(1:4) = (/33,64,364,55/)
    call slarnv(1,iseed,nd,y)
    dres = dres+y
    open(10,file='Dsyn.bin',form='unformatted',access='direct',recl=4*nd,status='replace')
    write(10,rec=1) dres(1:nd)
    close(10)
  print*,'max and min of dres for synthetic'
  print*, maxval(dres),minval(dres)
  endif

  ! for projection matrix
  do inet = 1,1+nnets!nnets 
    if (synthetic == 0) then
    write(filename,'("linsize",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=8)
    read(10,rec=1) linsize
    close(10)
    allocate(rw_p(linsize(2)),iw_p(2*linsize(2)+1))
    write(filename,'("S_p",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(2))
    read(10,rec=1) rw_p
    close(10)
    write(filename,'("C_p",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*(2*linsize(2)+1))
    read(10,rec=1) iw_p
    close(10)


    nnzero = int(spfrac*nd*subrow)
    print*, 'construct G*Ps beginning ...,max nnzero: ',nnzero
    allocate(iwgp(2*nnzero+1),rwgp(nnzero),colgp(nnzero))
    nnzero = linsize(2)
    print*,'Projection matrix: row,nnzeros'
    print*, subrow,nnzero
    start = 0
    leniw_p = nnzero*2+1
    lenrw_p = nnzero

    ! matrix multiplication 
    ! extract one row from P by multiplying delta
    nzid = 0
    xaug = 0
    startidx = 0

    ! G has mdim+4*loc cols
    ! ifort and gfortran diff in 'reclen'
    print*,inet,'th net',subrow/numth,' patchs in all'
      !$omp parallel &
      !$omp default(private) &
      !$omp shared(startidx,nzid,rwgp,iwgp,colgp,numth,nd,iw,rw,iw_p,rw_p,leniw_p,lenrw_p,leniw,lenrw,yy)  
    do jj = 1,subrow/numth
      !$omp do 
      do ii = (jj-1)*numth+1,jj*numth!subrow
        pkrow = 0
        pkrow(ii) = 1
        x = 0
        call aprod(2,subrow,mdim,x,pkrow,leniw_p,lenrw_p,iw_p,rw_p)
        xaug(1:mdim) = x
        y = 0
        call aprod(1,nd,mdim+nloc,xaug,y,leniw,lenrw,iw,rw)
        thread_num = omp_get_thread_num()+1
        yy(:,thread_num) = y
      enddo
      !$omp end do

!$omp master
    startidx = startidx +1
    do mm = 1,numth
      do kk = 1,nd
        if(abs(yy(kk,mm))>ftol) then
          nzid =nzid+1
          rwgp(nzid) = yy(kk,mm)
          iwgp(1+nzid) = kk 
          colgp(nzid) = mm+(startidx-1)*numth
          !colgp(nzid) = mm+(jj-1)*numth
        endif
      enddo
    enddo
    if (mod(startidx,50)==0) write(*,'(f4.1,a)'),real(startidx)*numth/subrow*100,' percent done'
!$omp end master
    enddo
!$omp end parallel 

    deallocate(iw_p,rw_p)


    ! add relocation,leniw-->lenrw
    do kk = 1,iw(1)
      coltmp = iw(1+lenrw+kk)
      if(coltmp>mdim) then
        nzid = nzid +1
        iwgp(1+nzid) = iw(1+kk)
        colgp(nzid) = coltmp-mdim+subrow 
        rwgp(nzid) = rw(kk)
      endif
    enddo


    iwgp(1) = nzid
    iwgp(nzid+2:2*nzid+1) = colgp
    leniwgp = 2*nzid+1
    lenrwgp = nzid

    linsize(1) = nd
    linsize(2) = nzid
    write(filename,'("linsize_gp",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=8,status='replace')
    write(10,rec=1) linsize
    close(10)
    write(filename,'("S_gp",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(2),status='replace')
    write(10,rec=1) rwgp(1:nzid)
    close(10)
    write(filename,'("C_gp",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*(2*linsize(2)+1),status='replace')
    write(10,rec=1) iwgp(1:2*nzid+1)
    close(10)

    deallocate(colgp)
    print*,'row,column and nnzero of Gp'
    print*,nd,subrow,nzid
  else

    write(filename,'("linsize_gp",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=8)
    read(10,rec=1) linsize
    close(10)
    allocate(rwgp(linsize(2)),iwgp(2*linsize(2)+1))
    write(filename,'("S_gp",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(2))
    read(10,rec=1) rwgp(1:linsize(2))
    close(10)
    write(filename,'("C_gp",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*(2*linsize(2)+1))
    read(10,rec=1) iwgp(1:2*linsize(2)+1)
    close(10)
    nzid=iwgp(1)
    leniwgp = 2*iwgp(1)+1
    lenrwgp = iwgp(1)
  endif


    xunknown = 0
    atol = 1e-4
    btol = 1e-5
    conlim = 100
    itnlim = 400
    istop = 0
    anorm = 0.0
    acond = 0.0
    arnorm = 0.0
    xnorm = 0.0
    localSize = 10
    damp = 0.001
    ! using lsmr to solve for the projection coefficients
    print*, 'LSMR beginning ...'

    nout = 63
    open(nout,file='lsmrout_sub.txt')

    call LSMR(nd, subrow+nloc, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
      atol, btol, conlim, itnlim, localSize,nout,&
      xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)
    write(*,'(a,f7.1,a,f7.1)')' damp is:',damp,' lsmr finished with condition number: ',acond
    write(*,'(a,f7.3,f7.3)') 'min and max of x:',minval(xunknown(1:subrow)),maxval(xunknown(1:subrow))

    ! write out results

    print*, 'writing out results beginning ...'
    write(filename,'("xcoef",i0,"fmtomo.bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*subrow)
    write(10,rec=1) xunknown(1:subrow)
    close(10)

    deallocate(iwgp,rwgp)
  enddo ! inet

  deallocate(rw,iw,dres,y)
  deallocate(yy)
  deallocate(x,xtrue,xaug)
contains
subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
end subroutine check  

  end
