!___________________________________________________________
!variables define for surface wave part
        program fmsurf
        use mod_surf
        use variable_def
        implicit none
        integer :: n,m,i,j,k
        real :: earthr
        integer numgrid2,numgrid1
        integer nz1,nz2
        integer invstep
        integer syn
        real noisethresh
        integer sgs,sgdl
         real sparsefrac

!----------------------------------------------------------
  write(unit=*,fmt='(a30)',advance='no')' initializing velocity grids '
  call initialize_velocity_grids
  print *,'......finished'

  write(unit=*,fmt='(a30)',advance='no')' initializing interfaces'
  call initialize_interfaces
  print *,'......finished'

	goxd=(90-vgrid(1,1)%lat(vgrid(1,1)%nlat-1)*180/pi)*pi/180
	!goxd=vgrid(1,1)%lat(vgrid(1,1)%nlat-1)
        !print*,vgrid(1,1)%nlat-1,vgrid(1,1)%lat(vgrid(1,1)%nlat-1)
	gozd=vgrid(1,1)%long(2)
        print*, 90-goxd*180/pi,gozd*180/pi
        nx = vgrid(1,1)%nlong
        ny = vgrid(1,1)%nlat
        dvxd = vgrid(1,1)%dlat0
        dvzd = vgrid(1,1)%dlong0
        earthr = vgrid(1,1)%r(vgrid(1,1)%nr-1)
	print*,'goxd and gozd',goxd,gozd
	OPEN(UNIT=64,FILE='tomo_surf.in',STATUS='old')
	read(64,*) minthk,noiselevel
        read(64,*) nsrcsurf
        nrc = nsrcsurf
	read(64,*)kmaxRc
	if(kmaxRc.gt.0)then
	write(*,*) 'number of periods for Rayleigh wave phase velocity'
	write(*,'(i4)') kmaxRc
	allocate(tRc(NP), stat=checkstat)
	if(checkstat > 0)then
	write(6,*)'error with allocate: program fmmin2d: real t(periods)'
	endif
	read(64,*)(tRc(i),i=1,kmaxRc)
	write(*,*) 'periods range(seconds)'
	write(*,'(60f6.2)') (tRc(i),i=1,kmaxRc)
	endif
	read(64,*)kmaxRg
	if(kmaxRg.gt.0)then
	write(*,*) 'number of periods for Rayleigh wave group velocity'
	write(*,'(i4)') kmaxRg
	allocate(tRg(NP), stat=checkstat)
	if(checkstat > 0)then
	write(6,*)'error with allocate: program fmmin2d: real t(periods)'
	endif
	read(64,*)(tRg(i),i=1,kmaxRg)
	write(*,*) 'periods range(seconds)'
	write(*,'(30f5.2)') (tRg(i),i=1,kmaxRg)
	endif
	read(64,*)kmaxLc
	if(kmaxLc.gt.0)then
	write(*,*) 'number of periods for Love wave phase velocity'
	write(*,'(i4)') kmaxLc
	allocate(tLc(NP), stat=checkstat)
	if(checkstat > 0)then
	write(6,*)'error with allocate: program fmmin2d: real t(periods)'
	endif
	read(64,*)(tLc(i),i=1,kmaxLc)
	write(*,*) 'periods range(seconds)'
	write(*,'(30f5.2)') (tLc(i),i=1,kmaxLc)
	endif
	read(64,*)kmaxLg
	if(kmaxLg.gt.0)then
	write(*,*) 'number of periods for Love wave group velocity'
	write(*,'(i4)') kmaxLg
	allocate(tLg(NP), stat=checkstat)
	if(checkstat > 0)then
	write(6,*)'error with allocate: program fmmin2d: real t(periods)'
	endif
	read(64,*)(tLg(i),i=1,kmaxLg)
	write(*,*) 'periods range(seconds)'
	write(*,'(30f5.2)') (tLg(i),i=1,kmaxLg)
	endif
        read(64,*)syn
        read(64,*)noisethresh
!       hidden bugs, read some parameters that control the forward oporator from outside 
        read(64,*) sgs,sgdl
        read(64,*) sparsefrac
	close(64)
        kmax=kmaxRc+kmaxRg+kmaxLc+kmaxLg

!___________________________________________________________
	dall=0
        open(unit=87,file='surfdata.dat',status='old')
        allocate(scxf(nsrcsurf,kmax),sczf(nsrcsurf,kmax),rcxf(nrc,nsrcsurf,kmax),rczf(nrc,nsrcsurf,kmax),stat=checkstat)
        IF(checkstat > 0)THEN
        WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL srcx,rczf'
        ENDIF
        allocate(periods(nsrcsurf,kmax),wavetype(nsrcsurf,kmax),&
        nrc1(nsrcsurf,kmax),nsrcsurf1(kmax),&
        igrt(nsrcsurf,kmax),stat=checkstat)
        IF(checkstat > 0)THEN
        WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL wavetype'
        ENDIF
        allocate(obst(nrc*nsrcsurf*kmax),noises(nrc*nsrcsurf*kmax),dsurf(nrc*nsrcsurf*kmax),stat=checkstat)
        IF(checkstat > 0)THEN
           WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL obst'
        ENDIF
        istep=0
        istep2=0
        dall=0
        knum=0
        knumo=12345
        do 
        read(87,'(a)',iostat=err) line
        if(err.eq.0) then
        if(line(1:1).eq.'#') then
        read(line,*) strf,sta1_lat,sta1_lon,period,wavetp,veltp
        if(sta1_lon<0) sta1_lon=sta1_lon+360
        if(wavetp.eq.2.and.veltp.eq.0) knum=period
        if(wavetp.eq.2.and.veltp.eq.1) knum=kmaxRc+period
        if(wavetp.eq.1.and.veltp.eq.0) knum=kmaxRg+kmaxRc+period
        if(wavetp.eq.1.and.veltp.eq.1) knum=kmaxLc+kmaxRg+kmaxRc+period
        if(knum.ne.knumo) then
        istep=0
        istep2=istep2+1
        endif
        istep=istep+1
        istep1=0
        !south hemosphere: negetive latitude
        sta1_lat=(90.0-sta1_lat)*pi/180.0
        sta1_lon=sta1_lon*pi/180.0
        scxf(istep,knum)=sta1_lat
        sczf(istep,knum)=sta1_lon
        periods(istep,knum)=period
        wavetype(istep,knum)=wavetp
        igrt(istep,knum)=veltp
        nsrcsurf1(knum)=istep
        knumo=knum
        else
        read(line,*) sta2_lat,sta2_lon,velvalue
        if(sta2_lon<0) sta2_lon=sta2_lon+360
        istep1=istep1+1
        dall=dall+1
        sta2_lat=(90.0-sta2_lat)*pi/180.0
        sta2_lon=sta2_lon*pi/180.0
        rcxf(istep1,istep,knum)=sta2_lat
        rczf(istep1,istep,knum)=sta2_lon
        call delsph(sta1_lat,sta1_lon,sta2_lat,sta2_lon,dist)
        !print *,sta1_lat,sta1_lon,sta2_lat,sta2_lon,dist
        obst(dall)=dist/velvalue
        noises(dall) = dist*noiselevel/velvalue
        if(noises(dall)<noisethresh) noises(dall) = noisethresh
        nrc1(istep,knum)=istep1
        endif
        else
        exit
        endif
        enddo
        close(87)
        !noises(1:dall) = 1.0/noises(1:dall)
        !noises(1:dall) = noises(1:dall)/maxval(noises(1:dall))
        open(44,file='inviter.in')
        read(44,*) invstep
        close(44)
        if(invstep==1.and.syn==0) then
          open(88,file='otimessurf.dat')
          write(88,'(i6)') dall
          do i=1,dall
            write(88,'(f8.4,4x,f7.2)') obst(i),noises(i)
          enddo
          close(88)
        endif

        numgrid1=0
        numgrid2=0
        ngrid1sep=0
        ngrid2sep=0
        nz1=0
        nz2=0

    if(n_interfaces==2) then
      allocate(vel(vgrid(1,1)%nlong,vgrid(1,1)%nlat,vgrid(1,1)%nr*2))
      allocate(depz(vgrid(1,1)%nr-1))
    do m=1,n_vtypes
     do n=1,n_vgrids
      do i=1,vgrid(n,m)%nr
         do j=1,vgrid(n,m)%nlat
            do k=1,vgrid(n,m)%nlong
                vel(k,j,(m-1)*vgrid(n,m)%nr+i)=vgrid(n,m)%velocity(vgrid(n,m)%nr-i+1,j,k)
            end do
         end do
      end do
      enddo
      enddo

        depz(1) = 0
      do m=2,vgrid(1,1)%nr-1
      depz(m) = earthr-vgrid(1,1)%r(vgrid(1,1)%nr-m)
      enddo
      nz = vgrid(1,1)%nr
      !print*,depz
      !print*,vel(5,5,1:nz+2)
      !print*,vel(5,5,nz+2+1:2*(nz+2))
      numgrid1 = vgrid(1,1)%nnode

    elseif(n_interfaces==3) then
      mface = sum(intrface(2)%r)/size(intrface(2)%r) 
      ngrid1sep=nint((vgrid(1,1)%r(vgrid(1,1)%nr-1)-mface)/vgrid(1,1)%dr0)
      ngrid2sep=nint((mface-vgrid(2,1)%r(2))/vgrid(2,1)%dr0)
      nz = ngrid1sep+ngrid2sep+2
      nz1 = vgrid(1,1)%nr
      nz2 = vgrid(2,1)%nr
      allocate(vel(vgrid(1,1)%nlong,vgrid(1,1)%nlat,(nz+2)*2))
      allocate(depz(nz+1))
     ! print*,ngrid1stop,ngrid2start
     ! print*,mface,vgrid(1,1)%r
     ! print*,vgrid(2,1)%r
      depz(ngrid1sep+1:1:-1)=earthr-vgrid(1,1)%r(vgrid(1,1)%nr-1-ngrid1sep:vgrid(1,1)%nr-1)
      depz(nz+1:nz-ngrid2sep:-1)=earthr-vgrid(2,1)%r(1:ngrid2sep+2)
    numgrid2 = vgrid(2,1)%nnode
    numgrid1 = vgrid(1,1)%nnode
        !print*,numgrid1,numgrid2,ngrid1sep,ngrid2sep,nz1,nz2
    !print*, ngrid1sep+2,nz+1-ngrid2sep
    !print*,vel(5,5,1:nz+2)
    do m=1,n_vtypes
         do j=1,vgrid(1,m)%nlat
            do k=1,vgrid(1,m)%nlong
           vel(k,j,(m-1)*(nz+2)+ngrid1sep+2:(m-1)*(nz+2)+1:-1)=vgrid(1,m)%velocity(vgrid(1,1)%nr-1-ngrid1sep:vgrid(1,1)%nr,j,k)
           vel(k,j,(m-1)*(nz+2)+(nz+2):(m-1)*(nz+2)+(nz+2)-ngrid2sep-1:-1)=vgrid(1,m)%velocity(1:ngrid2sep+2,j,k)
            end do
         end do
      enddo
    else
      print*, 'not implemented yet, only for 2 or 3 interfaces'
      stop
    endif
!-----------------------------------------------------------------------------------
       !if (vgrid(1,m).r(i)<mface) then
       !  n=2
       !else
       !  n=1
       !endif
! compute dispersion
! ray tracing
! kernel (b-spline in vertical direction)
! frechet direvitives for surface wave data
!print*,depz
!print*,nx,ny,nz,goxd,gozd,dvxd,dvzd,minthk
!print*,vel(5,5,1:nz)
        !open(34,file='frechetsurf.dat')
        open(35,file='ttimessurf.dat')
        call CalSurfG(nx,ny,nz,nz1,nz2,vel,dsurf, &
              goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg, &
              tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk, &
              scxf,sczf,rcxf,rczf,nrc1,nsrcsurf1,kmax,nsrcsurf,nrc, &
              ngrid1sep,ngrid2sep,n_interfaces,numgrid1,numgrid2,dall,& 
              sgs,sgdl,sparsefrac)
         !close(34)
         close(35)
         print*,'surface wave frechet direvitives done'

        deallocate(obst,noises,depz,dsurf)
 end program
!-----------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------
