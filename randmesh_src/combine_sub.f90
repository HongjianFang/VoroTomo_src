program combine_proj
implicit none

integer nvnr,nvnt,nvnp
real gnsr,gnst,gnsp
real gor,got,gop
integer nvgt,ni

!integer,parameter:: numfiles=100
integer numfiles
integer :: nrad,nlat,nlon
integer :: mdim,mdim2!=nrad*nlat*nlon*2
real,dimension(:),allocatable:: velmean,stdvel,vel,velref
!real velmean3d(nrad,nlat,nlon),stdvel3d(ndep,nlat,nlon)
real,dimension(:,:,:),allocatable:: velmean3d,stdvel3d
integer inet
character(len=130) filename
integer k,l,m,idx
integer veltype

!for sources relocation
integer ns
real rad,lat,lon,drad,dlat,dlon
integer::nloc!=5493*4
real,dimension(:),allocatable:: locmean,locall,stack!,loc(nloc),stdloc(nloc)
!INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: paths,patht
real,parameter::pi=3.141592654,earthr=6371
integer ii,j
integer velloc

integer idm1,idm2,idm3,idm4,nloc_sub,ndata,sloc
real fdm6,fdm5
integer,dimension(:),allocatable::sourceidx,eloc
real,dimension(:),allocatable::loc
character(len=32) arg

open(10,file='VoroTomo.in')
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) nlat,nlon,nrad
  read(10,*) 
  read(10,*) 
  read(10,*) nloc
  read(10,*) veltype
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) 
  read(10,*) velloc 
  close(10)
  print*,velloc

  call getarg(1,arg)
  read(arg,'(i4)') numfiles

  nloc = nloc*4
  mdim = nrad*nlat*nlon
  mdim2 = mdim
  if (veltype==1) mdim2 = mdim*2
  allocate(velmean(mdim2),stdvel(mdim2),vel(mdim2),velref(mdim2))
  allocate(velmean3d(nrad,nlat,nlon),stdvel3d(nrad,nlat,nlon))

velmean = 0
stdvel = 0
do inet = 1,numfiles
write(filename,'("./tempdata/sjfz",i0,"joint.bin")'),inet
open(10,file=filename,form='unformatted',access='direct',recl=4*mdim2)
read(10,rec=1) vel
close(10)
velmean = velmean+vel
stdvel = stdvel+vel**2
enddo
velmean = velmean/numfiles
stdvel = sqrt(stdvel/numfiles-velmean**2)
print*,maxval(velmean),minval(velmean)

idx = 0
open(10,file='vgrids.in')
read(10,*) nvgt,ni
read(10,*) nvnr,nvnt,nvnp
read(10,*) gnsr,gnst,gnsp
read(10,*) gor,got,gop
do k = 1,nvnr
do l = 1,nvnt
do m = 1,nvnp
idx = idx+1
read(10,*) velref(idx)
enddo
enddo
enddo
if (veltype==1) then
read(10,*) nvnr,nvnt,nvnp
read(10,*) gnsr,gnst,gnsp
read(10,*) gor,got,gop
do k = 1,nvnr
do l = 1,nvnt
do m = 1,nvnp
idx = idx+1
read(10,*) velref(idx)
enddo
enddo
enddo
endif

close(10)

velmean3d = reshape(velmean(1:mdim),(/nrad,nlat,nlon/))
print*,'mean variation for Vp:',maxval(abs(velmean(1:mdim)))!,maxval(abs(velmean(mdim/2+1:mdim)))
idx = 0
open(10,file='vgrids.in')
write(10,*) nvgt,ni
write(10,*) nvnr,nvnt,nvnp
write(10,*) gnsr,gnst,gnsp
write(10,*) gor,got,gop
do k = 1,nvnr
do l = 1,nvnt
do m = 1,nvnp
idx = idx+1
write(10,*) velref(idx)+velmean3d(k,l,m)
!write(10,*) velmean3d(k,l,m)
enddo
enddo
enddo
if (veltype==1) then
velmean3d = reshape(velmean(mdim+1:mdim2),(/nrad,nlat,nlon/))
print*,'mean variation for Vs:',maxval(abs(velmean(mdim+1:mdim2)))!,maxval(abs(velmean(mdim/2+1:mdim)))
write(10,*) nvnr,nvnt,nvnp
write(10,*) gnsr,gnst,gnsp
write(10,*) gor,got,gop
do k = 1,nvnr
do l = 1,nvnt
do m = 1,nvnp
idx = idx+1
write(10,*) velref(idx)+velmean3d(k,l,m)
enddo
enddo
enddo
endif
close(10)

idx = 0
stdvel3d = reshape(stdvel(1:mdim),(/nrad,nlat,nlon/))
print*,'mean variation for std Vp:',maxval(abs(stdvel(1:mdim)))!,maxval(abs(velmean(mdim/2+1:mdim)))
open(10,file='vgrids_std.in')
write(10,*) nvgt,ni
write(10,*) nvnr,nvnt,nvnp
write(10,*) gnsr,gnst,gnsp
write(10,*) gor,got,gop
do k = 1,nvnr
do l = 1,nvnt
do m = 1,nvnp
idx = idx+1
write(10,*) stdvel3d(k,l,m)
enddo
enddo
enddo
if (veltype==1) then
stdvel3d = reshape(stdvel(mdim+1:mdim2),(/nrad,nlat,nlon/))
print*,'mean variation for std Vs:',maxval(abs(stdvel(mdim+1:mdim2)))!,maxval(abs(velmean(mdim/2+1:mdim)))
do k = 1,nvnr
do l = 1,nvnt
do m = 1,nvnp
idx = idx+1
write(10,*) stdvel3d(k,l,m)
enddo
enddo
enddo
endif
close(10)

print*,'finishing combing subproj'
stop

!!print*,nloc
if (velloc==1) then
allocate(locmean(nloc),locall(nloc),stack(nloc))
locmean = 0
ns = nloc/4
!print*,ns
stack = 0
!stdloc = 0
do inet = 1,numfiles
write(filename,'("./tempdata/updatesrc",i0,".dat")'),inet
open(10,file=filename)
read(10,*) nloc_sub
allocate(loc(nloc_sub*4),eloc(nloc_sub))
do ii = 1,nloc_sub
read(10,*) eloc(ii)
enddo
close(10)
write(filename,'("./tempdata/xloc",i0,"fmtomo.bin")'),inet
open(10,file=filename,form='unformatted',access='direct',recl=4*nloc_sub*4)
read(10,rec=1) loc
close(10)
!print*, nloc_sub
locall = 0
do ii =1,nloc_sub
sloc = eloc(ii)
stack(4*sloc+1:4*sloc+4)=stack(4*sloc+1:4*sloc+4)+1
locall(4*sloc+1:4*sloc+4) = loc((ii-1)*4+1:ii*4) 
enddo
locmean = locmean+locall
deallocate(loc,eloc)
enddo

do ii = 1,nloc
if (stack(ii)==0) then
        locmean(ii) = 0
else
locmean(ii) = locmean(ii)/stack(ii)
endif
enddo
!update receivers.in
!print*,'bug1'
!read otimes
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
drad = locmean(4*sourceidx(ii)+1)
if (abs(drad) > 1.0) drad = drad/abs(drad)*1.0
rad = rad-drad!locmean(sourceidx(ii))
dlat = locmean(4*sourceidx(ii)+2)*180/(pi*(earthr-rad))
if (abs(dlat) > 0.05) dlat = dlat/abs(dlat)*0.05
lat = lat+dlat
!lat = lat+locmean(ns+sourceidx(ii))*180/(pi*(earthr-rad))
dlon = locmean(4*sourceidx(ii)+3)*180/(pi*(earthr-rad))
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
deallocate(sourceidx)
deallocate(velmean,stdvel,vel,velref)
deallocate(velmean3d,stdvel3d)
deallocate(locmean,locall,stack)!,loc(nloc),stdloc(nloc)
print*,'finishing all'
endif

end program
