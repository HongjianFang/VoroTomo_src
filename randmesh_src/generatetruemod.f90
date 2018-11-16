program truemod
  implicit none

  integer,parameter::nlat=512,nlon=1024,ndep=64
  integer ii,jj,kk
  integer ndim,reclen
  integer idx

  real,allocatable,dimension(:):: dv
  ndim=nlat*nlon*ndep
  allocate(dv(ndim))
  ! hawaii low,500 km *500 km*400 km
  ! lat 19.2N,-155.25
  ! lat:300:320 lon:686:706 dep 0:15
  do kk = 1,5
    do jj = 572,592
      do ii = 300,320
        idx = (kk-1)*nlon*nlat+(jj-1)*nlat+ii
        dv(idx) = -0.05
      enddo
    enddo
  enddo

  ! japan subduction
  ! lat 28N, lon 145
  ! lat 335, lon 412
  do kk = 1,10
    do jj = 412-kk*2-5,412-kk*2+5
      do ii = 335-20,335+20
        idx = (kk-1)*nlon*nlat+(jj-1)*nlat+ii
        dv(idx) = 0.025
      enddo
    enddo
  enddo

  ! south america sub
  !lat 28S, lon -70
  !lat 182, lon 288
  do kk = 1,35
    do jj = 825+kk*1-5,825+kk*1+5
      do ii = 176-10,176+25
        idx = (kk-1)*nlon*nlat+(jj-1)*nlat+ii
        dv(idx) = 0.035 
      enddo
    enddo
  enddo

  !north america craton
  !lat 32.4 -111,31.2,-81, 48.4 -78,49,-114
  do kk = 1,8
    do jj = 708,802
      do ii = 344,395
        idx = (kk-1)*nlon*nlat+(jj-1)*nlat+ii
        dv(idx) = 0.04
      enddo
    enddo
  enddo
  
  print*,ndim
  inquire(iolength=reclen) dv(1:ndim)
  print*,reclen
  open(34,file='truemod.bin',form='unformatted',access='direct',recl=reclen,status='replace')
  write(34,rec=1)dv(1:ndim)
  close(34)
  deallocate(dv)
  end program


