program tv_backproj
  implicit none

  ! variable difinition
  INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
  integer,allocatable,dimension(:)::iw
  real(kind=i5),allocatable,dimension(:)::rw
  integer*8 leniw,lenrw
  integer:: mdim,nlat,nlon,nrad,subrow2,mdim2,veltype
  integer:: subrow
  integer mrow,ncol,nd
  integer nnzero
  integer inet
  character(len=130) filename
  integer*4 linsize(2)
  integer,dimension(:),allocatable:: C_p
  real*4,dimension(:),allocatable:: S_p,xcoef
  real(kind=i5),allocatable,dimension(:)::x,dres,xp,xs
  integer jj
  character(len=32) arg!,arg2
  integer picknet

  call getarg(1,arg)
  read(arg,'(i4)') picknet
  print*,'pick net:',picknet
  !call getarg(2,arg2)
  !read(arg2,'(i4)') subrow

  open(10,file='VoroTomo.in')
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) nlat,nlon,nrad
  read(10,*) subrow
  read(10,*) 
  read(10,*) 
  read(10,*) veltype
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  close(10)
  mdim = nlat*nlon*nrad
  nnzero = mdim
  subrow2 = subrow
  mdim2 = mdim
  if (veltype == 1) then
     mdim2 = mdim*2
     subrow2 = subrow*2
     !nnzero = mdim/2
  endif

  allocate(iw(2*nnzero),rw(nnzero))
  nd = subrow
  allocate(x(mdim2),dres(nd),xp(nnzero))
  ! read all projection matrices into one
  nnzero = 0
  print*,'reading projection matrix'
    inet = picknet
    write(filename,'("./tempdata/linsize",i0,".bin")'),inet
    print*,filename
    open(10,file=filename,form='unformatted',access='direct',recl=8)
    read(10,rec=1) linsize
    close(10)
    allocate(S_p(linsize(2)),xcoef(subrow2),C_p(2*linsize(2)))
    write(filename,'("./tempdata/S_p",i0,".bin")'),inet
    !write(filename,'("./tempdata/S_invp",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(2))
    read(10,rec=1) S_p
    close(10)
    write(filename,'("./tempdata/C_p",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*(2*linsize(2)))
    read(10,rec=1) C_p
    close(10)
    !write(filename,'("xcoef",i0,"jointbs.bin")'),inet
    write(filename,'("./tempdata/xcoef",i0,"fmtomo_joint.bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*subrow2)
    read(10,rec=1) xcoef
    close(10)
    print*,'min and max values of model coefficients:',minval(xcoef),maxval(xcoef)

    nnzero = linsize(2)
    print*,linsize(1),linsize(2),nnzero
    iw(1:2*nnzero) = C_p(1:2*nnzero)
    rw(1:nnzero) = S_p(1:nnzero)
  print*,'finish reading with nnzero and rows: ', nnzero,subrow

  leniw = 2*nnzero
  lenrw = nnzero
  !nd = subrow
      print*,'finishing regularization matrix with nd,mdim,nnzero:',nd,mdim,nnzero
      if (veltype==0) then
      dres(1:subrow) = xcoef(1:subrow)
      xp = 0
      call aprod(2,nd,mdim,xp,dres,leniw,lenrw,iw,rw)
      x(1:mdim) = xp
      else
      dres(1:subrow) = xcoef(1:subrow)
      xp = 0
      call aprod(2,nd,mdim,xp,dres,leniw,lenrw,iw,rw)
      x(1:mdim) = xp
      xp = 0
      dres(1:subrow) = xcoef(subrow+1:subrow2)
      call aprod(2,nd,mdim,xp,dres,leniw,lenrw,iw,rw)
      x(mdim+1:mdim2) = xp
      endif

    ! write out results

    print*, 'writing out results beginning ...'
    write(filename,'("./tempdata/sjfz",i0,"joint.bin")'),picknet
    open(10,file=filename,form='unformatted',access='direct',recl=4*mdim2)
    write(10,rec=1) x
    close(10)
    print*,'min and max values of model:',minval(x(1:mdim)),maxval(x(1:mdim))
    if (veltype==1) then
    print*,'min and max values of model:',minval(x(mdim+1:2*mdim)),maxval(x(mdim+1:2*mdim))
    endif

  deallocate(S_p,C_p,xcoef)
  deallocate(iw,rw)
  deallocate(x,dres)
  end program


