program subspaceproj
  use lsmrModule, only:lsmr
  use netcdf
  use omp_lib
  implicit none
  ! variable definition for lsmr
  integer*8 leniw,lenrw
  integer*8 leniw_p,lenrw_p
  integer*8 leniwgp,lenrwgp
  integer, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
  integer,parameter:: numdir = 1
  integer,parameter:: weightornot=1,synthetic = 0
  real(kind=i5),parameter::spfrac= 0.1,ftol=0.1
  integer,parameter:: mdim = 129024*2,nloc=5490*4,nlocnew = 400*4
  integer nchoose,nchoose_surf
  real(kind=i5) weight_surf
  integer,allocatable,dimension(:):: choose,eqidx,choose_surf
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
  integer,allocatable,dimension(:)::nzero_id,nzero_idsurf,nrow,nrowsurf
  integer,allocatable,dimension(:)::iw
  integer,allocatable,dimension(:)::col
  real*4,allocatable,dimension(:)::rw
  integer,allocatable,dimension(:)::iw_p
  real(kind=i5),allocatable,dimension(:)::rw_p
  integer,allocatable,dimension(:)::iwgp,colgp
  real(kind=i5),allocatable,dimension(:)::rwgp
  real(kind=i5),allocatable,dimension(:)::dres,y,weight,dobs,dsyn
  real(kind=i5),allocatable,dimension(:,:)::yy,yys
  integer*4 linsize(2)

  integer nd,start
  integer ii,inet
  integer,parameter:: nnets = 100 
  integer nzid,nz
  integer subrow
  real,dimension(:),allocatable::norm,xunknown, pkrow
  real,allocatable,dimension(:):: x,xaug
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

  character(len=32) arg

  call getarg(1,arg)
  read(arg,'(i4)') subrow
  call getarg(2,arg)
  read(arg,'(i4)') inet
  allocate(pkrow(subrow),xunknown(subrow*2+nloc),stat=checkstat)
  if(checkstat>0) stop 'error allocating 1'
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
  call check(nf90_open("frechetsurf.nc",nf90_nowrite,ncid))
  call check(nf90_inq_dimid(ncid,"nzero",nzerodimid))
  call check(nf90_inq_dimid(ncid,'nrow',nrowdimid))
  call check(nf90_inquire_dimension(ncid,nzerodimid,len = nnfdsurf ))
  call check(nf90_inquire_dimension(ncid,nrowdimid,len = ndsurf ))
  call check(nf90_inq_varid(ncid,"Non_value",nzeroid))
  call check(nf90_inq_varid(ncid,"Non_row",nrowid))
  call check(nf90_inq_varid(ncid,"Non_id",nonid))
  allocate(nzerosurf(nnfdsurf),nrowsurf(ndsurf),nzero_idsurf(nnfdsurf),stat=checkstat)
  if(checkstat>0) stop 'error allocating memnetcdf'
  call check(nf90_get_var(ncid,nzeroid,nzerosurf))
  call check(nf90_get_var(ncid,nrowid,nrowsurf))
  call check(nf90_get_var(ncid,nonid,nzero_idsurf))
  call check(nf90_close(ncid))
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
  allocate(dres(nd),dobs(nd),dsyn(nd),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'

  allocate(x(mdim/2),xaug(mdim+nloc),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'

    nnzero = int(spfrac*nd*subrow)
    print*, 'construct G*Ps beginning ...,max nnzero: ',nnzero
    allocate(iwgp(2*nnzero),rwgp(nnzero),colgp(nnzero),stat=checkstat)
    if(checkstat>0) stop 'error allocating col'

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

  ! for surface wave
  open(10,file='otimessurf.dat')
  read(10,*)ndsurf
  do ii = 1,ndsurf
    read(10,*) dobs(ndbody+ii)
  enddo
  close(10)
  open(10,file='mtimessurf.dat')
  do ii = 1,ndsurf
    read(10,*) dsyn(ndbody+ii)
  enddo
  close(10)



  !do inet = 2,3!nnets!nnets 
  iw = 0
  rw = 0.
  jstep = 0
  write(filename,'("./tempdata/randomchoose",i0,".dat")'),inet
  open(10,file=filename)
  read(10,*)nchoose
  allocate(choose(nchoose),eqidx(nchoose),Gh(nchoose,4),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'
  do ii = 1, nchoose
    read(10,*) choose(ii),eqidx(ii)
  enddo
  close(10)
  if(checkstat>0) stop 'error allocating col'
  !do jj = 1,ndbody
  irow = 0
  do kk = 1,nchoose
    jj = choose(kk)
    irow = irow+1
    do ii = 1,nrow(jj)
      start = sum(nrow(1:jj-1))
      jstep = jstep + 1
      iw(jstep) = irow
      rw(jstep) = nzero(start+ii)
      col(jstep) = nzero_id(start+ii)
    enddo
  enddo

  write(filename,'("./tempdata/randomchoose_surf",i0,".dat")'),inet
  open(10,file=filename)
  read(10,*)nchoose_surf
  allocate(choose_surf(nchoose_surf),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'
  do ii = 1, nchoose_surf
    read(10,*) choose_surf(ii)
  enddo
  close(10)
  allocate(yy(nchoose+nchoose_surf,numth),yys(nchoose+nchoose_surf,numth),y(nchoose+nchoose_surf),stat=checkstat)

  weight_surf = sqrt(nchoose/float(nchoose_surf)*0.7)
  do kk = 1,nchoose_surf
    jj = choose_surf(kk)
    irow = irow+1
    do ii = 1,nrowsurf(jj)
      start = sum(nrowsurf(1:jj-1))
      jstep = jstep + 1
      iw(jstep) = irow!ndbody+jj
      rw(jstep) = nzerosurf(start+ii)*weight_surf
      col(jstep) = nzero_idsurf(start+ii)
    enddo
  enddo

  do ii = 1,jstep
    iw(jstep+ii) = col(ii)
  enddo
  !print*,maxval(nzero_id),maxval(choose)
  print*,'nnzero in G',jstep,minval(iw(1:jstep)),maxval(iw(1:jstep)),minval(iw(jstep+1:2*jstep)),maxval(iw(jstep+1:2*jstep))

  lenrw = jstep
  leniw = 2*jstep

  !do ii = 1,nd
  do ii = 1,nchoose
    dres(ii)=dobs(choose(ii))-dsyn(choose(ii))
  enddo
  do ii = 1,nchoose_surf
    dres(nchoose+ii)=(dobs(ndbody+choose_surf(ii))-dsyn(ndbody+choose_surf(ii)))*weight_surf
  enddo
  nd = nchoose+nchoose_surf

  if(weightornot==1) then
    allocate(weight(nd),stat=checkstat)
  if(checkstat>0) stop 'error allocating col'
  weight = 1.0/(1+0.05*exp(dres(1:nd)**2*1.0))
      do ii = 1,jstep
      rw(ii) = rw(ii)*weight(iw(ii))
    enddo
    do ii = 1,nd
      dres(ii) = dres(ii)*weight(ii) 
    enddo
    deallocate(weight)
  endif
 
  !if (synthetic == 1) then
  !  open(10,file='truemod.bin',form='unformatted',access='direct',recl=4*mdim)
  !  read(10,rec=1) xtrue 
  !  close(10)
  !  xaug = 0
  !  xaug(1:mdim) = xtrue
  !  y = 0
  !  call aprod(1,nd,mdim+nloc,xaug,y,leniw,lenrw,iw,rw)
  !  dres(1:nd) = y(1:nd)
  !  !do  ii=1,nd
  !  !  dres(ii)=dres(ii)*weight(ii)
  !  !  if(isnan(dres(ii))) dres(ii)=0
  !  !enddo
  !  iseed(1:4) = (/33,64,364,55/)
  !  call slarnv(1,iseed,nd,y)
  !  dres = dres+y
  !  open(10,file='Dsyn.bin',form='unformatted',access='direct',recl=4*nd,status='replace')
  !  write(10,rec=1) dres(1:nd)
  !  close(10)
  !print*,'max and min of dres for synthetic'
  !print*, maxval(dres),minval(dres)
  !endif

  ! for projection matrix
  !  if (synthetic == 0) then
    write(filename,'("./tempdata/linsize",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=8)
    read(10,rec=1) linsize
    close(10)
    allocate(rw_p(linsize(2)),iw_p(2*linsize(2)),stat=checkstat)
    if(checkstat>0) stop 'error allocating col'
    write(filename,'("./tempdata/S_p",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(2))
    read(10,rec=1) rw_p
    close(10)
    write(filename,'("./tempdata/C_p",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*(2*linsize(2)))
    read(10,rec=1) iw_p
    close(10)


    nnzero = linsize(2)
    print*,'Projection matrix: row,nnzeros'
    print*, subrow,nnzero
    start = 0
    leniw_p = nnzero*2
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
      !$omp shared(subrow,startidx,nzid,rwgp,iwgp,colgp,numth,nd,iw,rw,iw_p,rw_p,leniw_p,lenrw_p,leniw,lenrw,yy,yys)  
    do jj = 1,subrow/numth
      !$omp do 
      do ii = (jj-1)*numth+1,jj*numth!subrow
        pkrow = 0
        pkrow(ii) = 1
        x = 0
        xaug = 0
        call aprod(2,subrow,mdim,x,pkrow,leniw_p,lenrw_p,iw_p,rw_p)
        xaug(1:mdim/2) = x
        y = 0
        call aprod(1,nd,mdim+nloc,xaug,y,leniw,lenrw,iw,rw)
        thread_num = omp_get_thread_num()+1
        yy(:,thread_num) = y

        xaug = 0
        xaug(mdim/2+1:mdim) = x
        y = 0
        call aprod(1,nd,mdim+nloc,xaug,y,leniw,lenrw,iw,rw)
        thread_num = omp_get_thread_num()+1
        yys(:,thread_num) = y

      enddo
      !$omp end do

!$omp master
    startidx = startidx +1
    do mm = 1,numth
      do kk = 1,nd
        if(abs(yy(kk,mm))>ftol) then
          nzid =nzid+1
          rwgp(nzid) = yy(kk,mm)
          iwgp(nzid) = kk 
          colgp(nzid) = mm+(startidx-1)*numth
          !colgp(nzid) = mm+(jj-1)*numth
        endif
      enddo
    enddo

    do mm = 1,numth
      do kk = 1,nd
        if(abs(yys(kk,mm))>ftol) then
          nzid =nzid+1
          rwgp(nzid) = yys(kk,mm)
          iwgp(nzid) = kk 
          colgp(nzid) = subrow+mm+(startidx-1)*numth
          !colgp(nzid) = mm+(jj-1)*numth
        endif
      enddo
    enddo

    if (mod(startidx,50)==0) write(*,'(f5.1,a)'),real(startidx)*numth/subrow*100,' percent done'
!$omp end master
    enddo
!$omp end parallel 

    deallocate(iw_p,rw_p)


    ! add relocation,leniw-->lenrw
    !do kk = 1,lenrw
    !  coltmp = iw(lenrw+kk)
    !  if(coltmp>mdim) then
    !    nzid = nzid +1
    !    iwgp(nzid) = iw(kk)
    !    colgp(nzid) = coltmp-mdim+subrow 
    !    rwgp(nzid) = rw(kk)
    !  endif
    !enddo

    ! add relocation,leniw-->lenrw
    Gh = 0
    do kk = 1,lenrw
      coltmp = iw(lenrw+kk)
      if(coltmp>mdim) then
        Gh(iw(kk),(coltmp-1-mdim)/(nloc/4)+1)=rw(kk)
        !print*,iw(kk),(coltmp-1-mdim)/(nloc/4)+1,rw(kk)
      endif
    enddo
      do ii = 1,nchoose
        !choose(ii),eqidx(ii)
        iwgp(nzid+1:nzid+4) = ii
        colgp(nzid+1) = subrow*2+(eqidx(ii))*4+1
        colgp(nzid+2) = subrow*2+(eqidx(ii))*4+2
        colgp(nzid+3) = subrow*2+(eqidx(ii))*4+3
        colgp(nzid+4) = subrow*2+(eqidx(ii))*4+4
        rwgp(nzid+1:nzid+4) = Gh(ii,:)
        nzid = nzid+4
      enddo


    iwgp(nzid+1:2*nzid) = colgp
    leniwgp = 2*nzid
    lenrwgp = nzid

    !linsize(1) = nd
    !linsize(2) = nzid
    !write(filename,'("./tempdata/linsize_gp",i0,".bin")'),inet
    !open(10,file=filename,form='unformatted',access='direct',recl=8,status='replace')
    !write(10,rec=1) linsize
    !close(10)
    !write(filename,'("./tempdata/S_gp",i0,".bin")'),inet
    !open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(2),status='replace')
    !write(10,rec=1) rwgp(1:nzid)
    !close(10)
    !write(filename,'("./tempdata/C_gp",i0,".bin")'),inet
    !open(10,file=filename,form='unformatted',access='direct',recl=4*(2*linsize(2)),status='replace')
    !write(10,rec=1) iwgp(1:2*nzid)
    !close(10)

    !deallocate(colgp)
    !print*,'row,column and nnzero of Gp'
    !print*,nd,subrow,nzid
  !else

  !  write(filename,'("./tempdata/linsize_gp",i0,".bin")'),inet
  !  open(10,file=filename,form='unformatted',access='direct',recl=8)
  !  read(10,rec=1) linsize
  !  close(10)
  !  allocate(rwgp(linsize(2)),iwgp(2*linsize(2)))
  !  write(filename,'("./tempdata/S_gp",i0,".bin")'),inet
  !  open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(2))
  !  read(10,rec=1) rwgp(1:linsize(2))
  !  close(10)
  !  write(filename,'("./tempdata/C_gp",i0,".bin")'),inet
  !  open(10,file=filename,form='unformatted',access='direct',recl=4*(2*linsize(2)+1))
  !  read(10,rec=1) iwgp(1:2*linsize(2))
  !  close(10)
  !  nzid=linsize(2)
  !  leniwgp = 2*nzid
  !  lenrwgp = nzid
  !endif

	!scaling
!	norm = 0
!    do ii=1,nzid
!	norm(iwgp(ii+nzid)) = norm(iwgp(ii+nzid))+rwgp(ii)
!    enddo
!
!    do ii =1,subrow
!	norm(ii) = sqrt(norm(ii)/nd+0.01)
!    enddo
!
!    do ii =1,nzid
!	rwgp(ii) = rwgp(ii)/norm(iwgp(ii+nzid))
!    enddo

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
    open(nout,file='lsmrout_sub.txt')

    call LSMR(nd, subrow*2+nlocnew, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
      atol, btol, conlim, itnlim, localSize,nout,&
      xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)
!   do ii = 1,subrow
!   	xunknown(ii) = xunknown(ii)/norm(ii)
!   enddo
    !write(*,'(a,f7.1,f7.1,a,f7.1)')' norm range is:',minval(norm),maxval(norm),' lsmr finished with condition number: ',acond
    write(*,'(a,f7.3,f7.3)') 'min and max of xp:',minval(xunknown(1:subrow)),maxval(xunknown(1:subrow))
    write(*,'(a,f7.3,f7.3)') 'min and max of xs:',minval(xunknown(subrow+1:2*subrow)),maxval(xunknown(subrow+1:2*subrow))

    ! write out results

    print*, 'writing out results beginning ...'
    write(filename,'("./tempdata/xcoef",i0,"fmtomo_joint.bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*subrow*2)
    write(10,rec=1) xunknown(1:subrow*2)
    close(10)

!    write(filename,'("./tempdata/xloc",i0,"fmtomo.bin")'),inet
!    open(10,file=filename,form='unformatted',access='direct',recl=4*nloc)
!    write(10,rec=1) xunknown(subrow+1:subrow+nloc)
!    close(10)

    deallocate(yy,y)
  deallocate(Gh,eqidx)
  deallocate(choose,choose_surf)
!    if(inet==nnets) stop
!  enddo ! inet

  deallocate(iwgp,rwgp,colgp)
  deallocate(rw,iw,dres)
  !deallocate(yy)
  deallocate(x,xaug)
  deallocate(pkrow,xunknown)
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
