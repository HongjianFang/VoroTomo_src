program generatemesh
  !use omp_lib
  implicit none

  !integer,parameter:: npts = 500
  integer:: npts 
  !real,parameter:: radial = 6371.0,hvratio = 1.0
  !real,parameter:: cmb = 6339.345
  !real,parameter :: lat_s = 0.564426323505318,lon_w = 4.22002412744190
  !real,parameter :: dlat = 7.113040168330356E-004,dlon = 7.057201084883245E-004, drad = 1.15517241379310
  !integer,parameter:: nlat = 56,nlon = 72,nrad = 32 
  !integer,parameter:: nnets = 100
  real :: radial,hvratio
  real :: cmb 
  real :: lat_s,lon_w 
  real :: dlat,dlon,drad
  integer :: nlat,nlon,nrad
  integer idpts,ii,iseed(4)
!  real,dimension(:),allocatable::theta_surf,theta_cmb,theta_between
!  real,dimension(:),allocatable::phi_surf,phi_cmb,phi_between
!  real,dimension(:),allocatable::rad_surf,rad_cmb,rad_between
  real,dimension(:),allocatable::rad,phi,theta
  real,dimension(:),allocatable::lat,lon,radd
  real,dimension(:),allocatable::xpts,ypts,zpts,dis
  integer,dimension(:),allocatable::extractp,tsum!,indices
  real,dimension(:),allocatable::S,S_inv
  integer*4,dimension(:),allocatable::C
  integer*4 linsize(2)
  real,parameter::pi=3.141592654
  real xs,ys,zs
  integer nzall,nzero_sub
  integer irad,ilat,ilon
  integer idx,mdim,mloc,idx1
  integer reclen
  integer inet,picknet
  character(len=50) filename
  character(len=32) arg
  integer randseed
  real cmb0

  call getarg(1,arg)
  read(arg,'(i4)') picknet

!  npts_between = npts_upm+npts_lowm
!  allocate(theta_surf(npts_surf),phi_surf(npts_surf),rad_surf(npts_surf))
!  allocate(theta_cmb(npts_cmb),phi_cmb(npts_cmb),rad_cmb(npts_cmb))
!  allocate(theta_between(npts),phi_between(npts),rad_between(npts))
  open (10,file='VoroTomo.in')
  read(10,*) radial,cmb,hvratio
  read(10,*) lat_s,dlat
  read(10,*) lon_w,dlon
  read(10,*) drad
  read(10,*) nlat,nlon,nrad
  read(10,*) npts
  read(10,*) randseed
  close(10)

  !real,parameter:: radial = 6371.0,hvratio = 1.0
  !real,parameter:: cmb = 6339.345
  !real,parameter :: lat_s = 0.564426323505318,lon_w = 4.22002412744190
  !real,parameter :: dlat = 7.113040168330356E-004,dlon = 7.057201084883245E-004, drad = 1.15517241379310
  !integer,parameter:: nlat = 56,nlon = 72,nrad = 32 
  !integer,parameter:: nnets = 100
  allocate(theta(npts),phi(npts),rad(npts))
  allocate(xpts(npts),ypts(npts),zpts(npts))
  mdim = nlon*nlat*nrad
  allocate(extractp(mdim),dis(npts))
  allocate(tsum(npts))
  allocate(S(mdim),S_inv(mdim),C(2*mdim))
  allocate(lat(nlat),lon(nlon),radd(nrad))
  idpts = 0
  iseed(1:3) = (/38,62,346/)
  iseed(4) = randseed
  !do inet = 1,1+nnets
  !do inet = nnets,nnets
  do inet = 1,picknet
    print*, 'the',inet,'th subnet' 
 ! call slarnv(1,iseed,npts_surf,theta_surf)
 ! theta_surf = theta_surf*2*pi
 ! call slarnv(2,iseed,npts_surf,phi_surf)
 ! phi_surf = acos(phi_surf)
 ! rad_surf = radial

 ! call slarnv(1,iseed,npts_cmb,theta_cmb)
 ! theta_cmb = theta_cmb*2*pi
 ! call slarnv(2,iseed,npts_cmb,phi_cmb)
 ! phi_cmb = acos(phi_cmb)
 ! rad_cmb = cmb 

  call slarnv(1,iseed,npts,theta)
  theta = lon_w+theta*(nlon-1)*dlon
  call slarnv(1,iseed,npts,phi)
  phi = pi/2-(lat_s+phi*(nlat-1)*dlat)
  call slarnv(1,iseed,npts,rad)
  cmb0 = radial - (radial-cmb)*hvratio
  rad = rad*((nrad-1)*drad*hvratio)+cmb0

!  rad(1:npts_surf) = rad_surf
!  rad(npts_surf+1:npts_cmb+npts_surf) = rad_cmb
!  rad(npts_cmb+npts_surf+1:npts) = rad_between
!  phi(1:npts_surf) = phi_surf
!  phi(npts_surf+1:npts_cmb+npts_surf) = phi_cmb
!  phi(npts_cmb+npts_surf+1:npts) = phi_between
!  theta(1:npts_surf) = theta_surf
!  theta(npts_surf+1:npts_cmb+npts_surf) = theta_cmb
!  theta(npts_cmb+npts_surf+1:npts) = theta_between

  xpts = rad*sin(phi)*cos(theta)
  ypts = rad*sin(phi)*sin(theta)
  zpts = rad*cos(phi)

!to radial,:w

  do ilon = 1,nlon
    lon(ilon) = (lon_w+(ilon-1)*dlon)
  enddo

  do ilat = 1,nlat
    lat(ilat) = (lat_s+(ilat-1)*dlat)
  enddo

  do irad = 1,nrad
    radd(irad) = cmb0+(irad-1)*drad*hvratio
    !radd(irad) = radial-cmb+(irad-1)*drad
  enddo

  !idx = 1
  nzall = mdim
  tsum = 0
  idx1 = 0
 ! !$omp parallel &
 ! !$omp default(private) &
 ! !$omp shared(xpts,radd,lat,lon,C,nzall) 
 ! !$omp do
 ! open(10,file='test2.dat')
    do ilon = 1,nlon
      do ilat = 1,nlat
        do irad = 1,nrad
      !do ilat = 256,257
        xs = radd(irad)*sin(pi/2-lat(ilat))*cos(lon(ilon)) 
        ys = radd(irad)*sin(pi/2-lat(ilat))*sin(lon(ilon)) 
        zs = radd(irad)*cos(pi/2-lat(ilat))
        dis = (xpts-xs)**2+(ypts-ys)**2+(zpts-zs)**2
        mloc = minloc(dis,1)
        !idx = (irad-1)*nlat*nlon+(ilon-1)*nlat+ilat
        idx = (ilon-1)*nlat*nrad+(ilat-1)*nrad+irad
        idx1 = idx1+1
        extractp(idx1) = mloc!minloc(dis,1)
        C(idx1) = mloc!minloc(dis,1)
        C(nzall+idx1) = idx
        tsum(mloc) = tsum(mloc)+1
 !       write(10,*)extractp(idx),lon(ilon),radd(irad)
      enddo
    enddo
  enddo
  !!$omp end do
  !!$omp end parallel
  !C(1) = nzall
 !do ii=1,nzall
 !  idx = count(extractp==C(ii))
 !  S_inv(ii) = 1.0/sqrt(real(idx))
 !  !S_inv(ii) = 1.0/sqrt(real(tsum(C(ii))))
 !enddo
 !do ii=1,npts
 !print*,tsum(ii),count(extractp==ii)
 !enddo
 !print*,sum(tsum)
  S = 1.0
 ! print*,'stop here'
 ! stop

  !nzall = 0
  !print*,'begin sorting'
  !call quicksort(extractp,1,mdim)
  !print*,'removing repeated indices'
  !call unique_indices(extractp,mdim,indices)
  !!print*,npts,len(indices)
  !do ii = 1,npts
  !  nzero_sub = indices(ii+1)-indices(ii)
  !  S(indices(ii):indices(ii+1)) = 1.0/nzero_sub
  !  C(indices(ii):indices(ii+1)) = extractp(indices(ii))
  !  L(ii) = indices(ii+1)
  !enddo
  !nzall = indices(ii+1)
  !print*,'begin packing'
  !!$omp parallel &
  !!$omp default(private) &
  !!$omp shared(nrad,pi,nlat,nlon,xpts,radd,lat,lon,extractp) &
  !!$omp do
  !do ii = 1,npts
  !  nzero_sub = count(extractp==ii)
  !  !S(nzall:nzall+nzero_sub) = 1.0/nzero_sub
  !  S(nzall:nzall+nzero_sub) = 1.0
  !  C(nzall:nzall+nzero_sub) = pack(extractp,extractp==ii)
  !  nzall = nzall+nzero_sub
  !  L(ii) = nzall
  !enddo
!!$omp end do
!!$omp end parallel
  linsize(1) = npts
  linsize(2) = nzall

  !write out projection matrix
  print*,'writing out matrix'
 ! inquire(iolength=reclen) S_inv(1:nzall)
 ! write(filename,'("./tempdata/S_invp",i0,".bin")'),inet
 ! open(34,file=filename,form='unformatted',access='direct',recl=reclen,status='replace')
 ! write(34,rec=1)S_inv(1:nzall)
 ! close(34)
  inquire(iolength=reclen) S(1:nzall)
  write(filename,'("./tempdata/S_p",i0,".bin")'),inet
  open(34,file=filename,form='unformatted',access='direct',recl=reclen,status='replace')
  write(34,rec=1)S(1:nzall)
  close(34)
  inquire(iolength=reclen) C(1:2*nzall)
  write(filename,'("./tempdata/C_p",i0,".bin")'),inet
  open(34,file=filename,form='unformatted',access='direct',recl=reclen,status='replace')
  write(34,rec=1)C(1:2*nzall)
  close(34)
!  inquire(iolength=reclen) L(1:npts)
!  write(filename,'("L_p",i0,".bin")'),inet
!  open(34,file=filename,form='unformatted',access='direct',recl=reclen,status='replace')
!  write(34,rec=1)L(1:npts)
!  close(34)
  inquire(iolength=reclen) linsize(1:2)
  write(filename,'("./tempdata/linsize",i0,".bin")'),inet
  open(34,file=filename,form='unformatted',access='direct',recl=reclen,status='replace')
  write(34,rec=1)linsize(1:2)
  close(34)

enddo !inet

!  print*,'deallocate memory'
!  deallocate(theta_surf,theta_cmb,theta_between)
!  print*,'deallocate memory'
!  deallocate(phi_surf,phi_cmb,phi_between)
!  print*,'deallocate memory'
  deallocate(lat,lon,radd,S,dis,S_inv)
  print*,'deallocate memory'
  deallocate(rad,phi,theta,xpts,ypts,zpts,extractp,C)
  print*,'deallocate memory'
!  deallocate(xpts,ypts,zpts)
!  print*,'deallocate memory'
!  deallocate(extractp,C,L)
!  print*,'deallocate memory'
!  deallocate(lat,lon,radd)
!  print*,'deallocate memory'
!  deallocate(S)
!  print*,'deallocate memory'
!  deallocate(dis)
!  print*,'deallocate memory'

  end program

! subroutine findlocx(vec,n,val,nr,idx)
!   implicit none
!   integer vec(*)
