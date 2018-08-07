program tv_backproj
  implicit none

  ! variable difinition
  INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
  integer,allocatable,dimension(:)::iw
  real(kind=i5),allocatable,dimension(:)::rw
  integer*8 leniw,lenrw
  integer,parameter:: mdim = 1944320 
  integer,parameter:: subrow = 7200
  integer mrow,ncol,nd
  integer nnzero
  integer inet
  character(len=30) filename
  integer*4 linsize(2)
  integer,dimension(:),allocatable:: C_p
  real*4,dimension(:),allocatable:: S_p,xcoef
  real(kind=i5),allocatable,dimension(:)::x,dres
  integer jj
  character(len=32) arg
  integer picknet

  call getarg(1,arg)
  read(arg,'(i4)') picknet
  print*,'pick net:',picknet


  nnzero = mdim
  allocate(iw(2*nnzero),rw(nnzero))
  nd = subrow
  allocate(x(mdim),dres(nd))
  ! read all projection matrices into one
  nnzero = 0
  print*,'reading projection matrix'
    inet = picknet
    write(filename,'("linsize",i0,".bin")'),inet
    print*,filename
    open(10,file=filename,form='unformatted',access='direct',recl=8)
    read(10,rec=1) linsize
    close(10)
    allocate(S_p(linsize(2)),xcoef(linsize(1)),C_p(2*linsize(2)+1))
    write(filename,'("S_p",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(2))
    read(10,rec=1) S_p
    close(10)
    write(filename,'("C_p",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*(2*linsize(2)+1))
    read(10,rec=1) C_p
    close(10)
    !write(filename,'("xcoef",i0,"jointbs.bin")'),inet
    write(filename,'("xcoef",i0,"joint_syn.bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(1))
    read(10,rec=1) xcoef
    close(10)
    print*,'min and max values of model coefficients:',minval(xcoef),maxval(xcoef)

    nnzero = C_p(1)
    print*,linsize(1),linsize(2),nnzero
    iw(1:2*nnzero) = C_p(2:2*nnzero+1)
    rw(1:nnzero) = S_p(1:nnzero)
    dres(1:subrow) = xcoef(1:subrow)
  print*,'finish reading with nnzero and rows: ', nnzero,subrow

  leniw = 2*nnzero+1
  lenrw = nnzero
  nd = subrow
      print*,'finishing regularization matrix with nd,mdim,nnzero:',nd,mdim,nnzero
      x = 0
      call aprod(2,nd,mdim,x,dres,leniw,lenrw,iw,rw)

    ! write out results

    print*, 'writing out results beginning ...'
    write(filename,'("mit_s",i0,"jointbs.bin")'),picknet
    open(10,file=filename,form='unformatted',access='direct',recl=4*mdim)
    write(10,rec=1) x
    close(10)
    print*,'min and max values of model:',minval(x),maxval(x)

  deallocate(S_p,C_p,xcoef)
  deallocate(iw,rw)
  deallocate(x,dres)
  end program


