program combine_proj
implicit none

integer nvnr,nvnt,nvnp
real gnsr,gnst,gnsp
real gor,got,gop
integer nvgt,ni

integer,parameter:: numfiles=100
integer,parameter:: ndep = 32,nlat = 56, nlon = 72
integer,parameter:: mdim=ndep*nlat*nlon
real velmean(mdim),stdvel(mdim),vel(mdim),velref(mdim)
real velmean3d(ndep,nlat,nlon),stdvel3d(ndep,nlat,nlon)
integer inet
character(len=130) filename
integer k,l,m,idx

!for sources relocation
integer ns
REAL, DIMENSION(:), ALLOCATABLE :: srad,slat,slon,stp
REAL, DIMENSION(:), ALLOCATABLE :: sradnew,slatnew,slonnew
INTEGER, DIMENSION(:), ALLOCATABLE :: lots,nps
CHARACTER (LEN=10), DIMENSION(:), ALLOCATABLE :: tpath
INTEGER, PARAMETER :: maxsd=50,maxs=60,maxp=40
integer idm1,idm2,idm3
real stpnew
!integer,parameter::nloc=0
!real locmean(nloc),loc(nloc),stdloc(nloc)
!INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: paths,patht
!real,parameter::pi=3.141592654,earthr=6371
!integer i,j



velmean = 0
stdvel = 0
do inet = 1,numfiles
write(filename,'("./tempdata800/sjfz",i0,"body.bin")'),inet
open(10,file=filename,form='unformatted',access='direct',recl=4*mdim)
read(10,rec=1) vel
close(10)
velmean = velmean+vel
stdvel = stdvel+vel**2
enddo
velmean = velmean/numfiles
stdvel = sqrt(stdvel/numfiles-velmean**2)
velmean3d = reshape(velmean,(/ndep,nlat,nlon/))
stdvel3d = reshape(stdvel,(/ndep,nlat,nlon/))

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
close(10)

idx = 0
open(10,file='vgrids_subnetmean.in')
write(10,*) nvgt,ni
write(10,*) nvnr,nvnt,nvnp
write(10,*) gnsr,gnst,gnsp
write(10,*) gor,got,gop
do k = 1,nvnr
do l = 1,nvnt
do m = 1,nvnp
idx = idx+1
!write(10,*) velref(idx)+velmean3d(k,l,m)
write(10,*) velmean3d(k,l,m)
enddo
enddo
enddo
close(10)

idx = 0
open(10,file='vgrids_subnetstd.in')
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
close(10)

print*,'finishing combing subproj'
!locmean = 0
!stdloc = 0
!do inet = 1,numfiles
!write(filename,'("./tempdata/xloc",i0,"fmtomo.bin")'),inet
!open(10,file=filename,form='unformatted',access='direct',recl=4*nloc)
!read(10,rec=1) loc
!close(10)
!locmean = locmean+loc
!stdloc = stdloc+loc**2
!enddo
!locmean = locmean/numfiles
!stdloc = sqrt(stdloc/numfiles-locmean**2)
!
!
!   OPEN(UNIT=10,FILE='sources.in',STATUS='old')
!   OPEN(UNIT=11,FILE='sourcesnew.in')
!   OPEN(UNIT=12,FILE='sourcesstd.in')
!   READ(10,*)ns
!   write(11,*)ns
!   write(12,*)ns
!   ALLOCATE(paths(maxs,maxp,ns),patht(maxs,maxp,ns))
!   ALLOCATE(srad(ns),slat(ns),slon(ns),stp(ns))
!   ALLOCATE(sradnew(ns),slatnew(ns),slonnew(ns))
!	do i=1,ns
!	sradnew(i) = srad(i)-locmean(i)
!	slatnew(i) = slat(i)+locmean(ns+i)*180/(pi*(earthr-sradnew(i)))
!	slonnew(i) = slon(i)+locmean(2*ns+i)*180/(pi*(earthr-sradnew(i)))
!	enddo
!      DO i=1,ns
!         READ(10,*)idm1
!         write(11,*)idm1
!         write(12,*)idm1
!         IF(idm1.EQ.1)READ(10,*)
!         READ(10,*)srad(i),slat(i),slon(i)
!         write(11,*)sradnew(i),slatnew(i),slonnew(i)
!         write(12,*)stdloc(i),stdloc(ns+i),stdloc(2*ns+i)
!         READ(10,*)idm2
!         write(11,*)idm2
!         write(12,*)idm2
!         DO j=1,idm2
!            READ(10,*)idm3
!            write(11,*)idm3
!            write(12,*)idm3
!            READ(10,*)paths(1:2*idm3,j,i)
!            write(11,*)paths(1:2*idm3,j,i)
!            write(12,*)paths(1:2*idm3,j,i)
!            READ(10,*)patht(1:idm3,j,i)
!            write(11,*)patht(1:idm3,j,i)
!            write(12,*)patht(1:idm3,j,i)
!         ENDDO
!      ENDDO
!      CLOSE(10)
!      close(11)
!      close(12)
!
!	open(10,file='stimes.dat',status='old')
!	open(11,file='stimesnew.dat')
!	read(10,*)ns
!	write(11,*)ns
!	do i=1,ns
!	read(10,*)stp(i)
!	stpnew = stp(i)+locmean(3*ns+i)
!	write(11,*)stpnew
!	enddo
!	close(10)
!	close(11)
!	
!deallocate(srad,slat,slon,stp)
!deallocate(sradnew,slatnew,slonnew)
!deallocate(patht,paths)
end program
