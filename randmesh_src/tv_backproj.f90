program tv_backproj
  use lsmrModule, only:lsmr
  implicit none

  ! variable difinition
  INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
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
  integer*8 leniw,lenrw
  real(kind=i5),parameter::tolr = 1e-6
  integer,parameter:: initer=1 
  integer,parameter:: mdim = 145800 ,nloc = 10112
  integer,parameter::nnets = 1 
  integer,parameter:: subrow = 2000
  integer mrow,ncol,nd
  integer start,nnzero
  integer inet
  character(len=30) filename
  integer*4 linsize(2)
  integer,dimension(:),allocatable:: C_p
  real*4,dimension(:),allocatable:: S_p,xcoef
  integer ii,zero_tmp,col_tmp,nar,nar_tmp
  real*4,parameter:: weight = 0.00001
  integer iiter
  real(kind=i5),allocatable,dimension(:)::x,dres
  integer,parameter:: nlat = 45,nlon = 81,nrad = 40 
  integer colidx,kk,jj


  mrow = nnets*subrow
  ncol = mdim
  nnzero = nnets*ncol!+mdim*7
  nd = nnets*subrow+mdim
  print*,'matrix size and maximum nnzeo',nd,ncol,nnzero
  allocate(iw(2*nnzero),rw(nnzero),col(nnzero))
  allocate(x(mdim),dres(nd))
  ! read all projection matrices into one
  zero_tmp = 0
  col_tmp = 0
  nnzero = 0
  print*,'reading projection matrix'
  do inet = 1,nnets 
    write(filename,'("linsize",i0,".bin")'),inet
    print*,filename
    open(10,file=filename,form='unformatted',access='direct',recl=8)
    read(10,rec=1) linsize
    close(10)
    print*,linsize(1),linsize(2)
    allocate(S_p(linsize(2)),xcoef(linsize(1)),C_p(2*linsize(2)))
    write(filename,'("S_invp",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(2))
    read(10,rec=1) S_p
    close(10)
    write(filename,'("C_p",i0,".bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*(2*linsize(2)))
    read(10,rec=1) C_p
    close(10)
  !  write(filename,'("L_p",i0,".bin")'),inet
  !  open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(1))
  !  read(10,rec=1) L_p
  !  close(10)
    !write(filename,'("xcoef",i0,"joint.bin")'),inet
    write(filename,'("xcoef",i0,"joint_syn.bin")'),inet
    open(10,file=filename,form='unformatted',access='direct',recl=4*linsize(1))
    read(10,rec=1) xcoef
    close(10)
    print*,'min and max values of model coefficients:',minval(xcoef),maxval(xcoef)

    print*,'reading projection matrix'
    !print*,linsize(1),linsize(2),subrow
    !start = 0
    !do ii = 1,subrow
    !  iw(zero_tmp+start+2:zero_tmp+L_p(ii)+1) = col_tmp+ii
    !  start = L_p(ii)
    !enddo
    !col(zero_tmp+1:zero_tmp+L_p(subrow)) = C_p
    !rw(zero_tmp+1:zero_tmp+L_p(subrow)) = S_p
    !dres(col_tmp+1:col_tmp+subrow) = xcoef
    !zero_tmp = zero_tmp+L_p(subrow)
    !col_tmp = col_tmp+subrow
    nnzero = linsize(2)
    iw(zero_tmp+1:zero_tmp+nnzero) = C_p(1:nnzero)+col_tmp
    col(zero_tmp+1:zero_tmp+nnzero) = C_p(nnzero+1:2*nnzero)
    rw(zero_tmp+1:zero_tmp+nnzero) = S_p
    dres(col_tmp+1:col_tmp+subrow) = xcoef
    zero_tmp = zero_tmp+nnzero
    col_tmp = col_tmp+subrow
    deallocate(S_p,C_p,xcoef)
  enddo
  !col_tmp = col_tmp+subrow
  nar_tmp = zero_tmp
  nnzero = nar_tmp
  nar = nar_tmp
  print*,'finish reading with nnzero and rows: ',nnzero,col_tmp,maxval(iw(1:nar)),maxval(iw(nar+1:2*nar))

  !augument identity matrix 
  !do ii = 1,ncol
  !  rw(zero_tmp+ii) = 1.0
  !  iw(2+zero_tmp+ii) = col_tmp+ii
  !  col(zero_tmp+ii) = ii
  !enddo
  !nnzero = nnzero+ncol
  colidx = 0
  if(weight>1e-4) then 
  do kk = 1,nrad
    do jj = 1,nlon
      do ii = 1,nlat
        if(ii==1.or.ii==nlat.or.jj==1.or.jj==nlon.or.kk==1.or.kk==nrad) then
          colidx = colidx+1
          col(nar+1) = (kk-1)*nlat*nlon+(jj-1)*nlat+ii
          rw(nar+1) = 2.0*weight
          iw(1+nar) = col_tmp+colidx
          dres(col_tmp+colidx) = 0.0
          nar = nar+1
        else
             colidx=colidx+1
              col(nar+1)=(kk-1)*nlat*nlon+(jj-1)*nlat+ii
              rw(nar+1)=6.0*weight
              iw(nar+1)=col_tmp+colidx
              rw(nar+2)=-1.0*weight
              iw(nar+2)=col_tmp+colidx
              col(nar+2)=(kk-1)*nlat*nlon+(jj-1)*nlat+ii-1
              rw(nar+3)=-1.0*weight
              iw(nar+3)=col_tmp+colidx
              col(nar+3)=(kk-1)*nlat*nlon+(jj-1)*nlat+ii+1
              rw(nar+4)=-1.0*weight
              iw(nar+4)=col_tmp+colidx
              col(nar+4)=(kk-1)*nlat*nlon+(jj-2)*nlat+ii
              rw(nar+5)=-1.0*weight
              iw(nar+5)=col_tmp+colidx
              col(nar+5)=(kk-1)*nlat*nlon+jj*nlat+ii
              rw(nar+6)=-1.0*weight
              iw(nar+6)=col_tmp+colidx
              col(nar+6)=(kk-2)*nlat*nlon+(jj-1)*nlat+ii
              rw(nar+7)=-1.0*weight
              iw(nar+7)=col_tmp+colidx
              col(nar+7)=kk*nlat*nlon+(jj-1)*nlat+ii
              dres(col_tmp+colidx)=0
              nar=nar+7
            endif
          enddo
        enddo
      enddo
    endif
      nnzero = nar
  iw(nnzero+1:2*nnzero) = col(1:nnzero)
  leniw = 2*nnzero
  lenrw = nnzero
  nd = col_tmp+colidx
      print*,'finishing regularization matrix with nd,mdim,nnzero:',nd,mdim,nnzero

    !using IRLS to obtain final results
    print*, 'LSMR beginning ...'
    x = 0
    atol = 1e-6
    btol = 1e-6
    conlim = 100
    itnlim = 400
    istop = 0
    anorm = 0.0
    acond = 0.0
    arnorm = 0.0
    xnorm = 0.0
    localSize = 10
    damp = 0.2
    ! using lsmr to solve for the projection
    ! coefficients
    nout = 63
    dres = dres*1000
    open(nout,file='lsmrout_tv.txt')
    call LSMR(nd, mdim, leniw,lenrw,iw,rw,dres,damp,&
      atol, btol, conlim, itnlim,localSize,nout,&
      x, istop, itn, anorm, acond,rnorm,arnorm, xnorm)
    write(*,'(a,f7.1,a,f7.1)')' damp is:',damp,' lsmr finished with condition number: ',acond

    if(initer>1) then
    do iiter = 1, initer-1
      do ii=nar_tmp+1,nar
        if (abs(x(iw(ii))).lt.tolr) then
          rw(ii)=1.0/sqrt(tolr)*weight
        else
          rw(ii)=sqrt(1.0/abs(x(iw(ii))))*weight
        endif
      enddo

      x = 0
      istop = 0
      anorm = 0.0
      acond = 0.0
      arnorm = 0.0
      xnorm = 0.0

      call LSMR(nd, mdim, leniw,lenrw,iw,rw,dres,damp,&
      atol, btol, conlim, itnlim,localSize,nout,&
      x, istop, itn, anorm, acond,rnorm,arnorm, xnorm)
    write(*,'(a,f7.1,a,f7.1)')' damp is:',damp,' lsmr finished with condition number: ',acond
    enddo
  endif

    ! write out results

    x = x/1000
    print*, 'writing out results beginning ...'
    write(filename,'("mit_s",i0,"_tv_syn.bin")'),nnets
    open(10,file=filename,form='unformatted',access='direct',recl=4*mdim)
    write(10,rec=1) x
    close(10)
    print*,'min and max values of model:',minval(x),maxval(x)

  deallocate(iw,rw,col)
  deallocate(x,dres)
  end program


