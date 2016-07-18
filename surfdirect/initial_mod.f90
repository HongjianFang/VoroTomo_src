!*************************************************************
! This subroutine reads the initial parametrization of the velocity
! fields from which the values on the propagation grid are to be interpolated
! from file into the appropriate structures

subroutine initialize_velocity_grids
use mod_surf
use variable_def
implicit none

integer :: n,m,i,j,k

open(10,file='vgrids.in')

! read the velocity mode from the file

! read the number of regions in the input velocity structure
read (10,*) n_vgrids,n_vtypes

! allocate space for these regions
allocate(vgrid(n_vgrids,n_vtypes))
do m=1,n_vtypes
   do n=1,n_vgrids 
      call vgrid_defaults(vgrid(n,m)) 
   end do
end do
! read the grid properties and velocity values to be interpolated for each region

do m=1,n_vtypes

   do n=1,n_vgrids

      ! grid parameters
      read(10,*) vgrid(n,m)%nr,vgrid(n,m)%nlat,vgrid(n,m)%nlong
      read(10,*) vgrid(n,m)%dr0,vgrid(n,m)%dlat0,vgrid(n,m)%dlong0
      read(10,*) vgrid(n,m)%r0,vgrid(n,m)%lat0,vgrid(n,m)%long0


! initialize the grid

      allocate(vgrid(n,m)%r(vgrid(n,m)%nr),vgrid(n,m)%lat(vgrid(n,m)%nlat), &
           vgrid(n,m)%long(vgrid(n,m)%nlong))

      do i=1,vgrid(n,m)%nr
         vgrid(n,m)%r(i)=vgrid(n,m)%r0 + (i-1)*vgrid(n,m)%dr0
      end do

      do i=1,vgrid(n,m)%nlat
         vgrid(n,m)%lat(i)=vgrid(n,m)%lat0 + (i-1)*vgrid(n,m)%dlat0
      end do

      do i=1,vgrid(n,m)%nlong
         vgrid(n,m)%long(i)=vgrid(n,m)%long0 + (i-1)*vgrid(n,m)%dlong0
      end do

! read in the velocity values on the interpolation grid

      allocate(vgrid(n,m)%velocity(vgrid(n,m)%nr,vgrid(n,m)%nlat,vgrid(n,m)%nlong))

      do i=1,vgrid(n,m)%nr
         do j=1,vgrid(n,m)%nlat
            do k=1,vgrid(n,m)%nlong
               read (10,*) vgrid(n,m)%velocity(i,j,k)
            end do
         end do
      end do

      if (count(vgrid(n,m)%velocity > 20.0_dp) > 0 ) &
           print *,'*** WARNING *** : velocity grid contains values larger than 20 km/sec'
      if (count(vgrid(n,m)%velocity < 1.0_dp) > 0 ) &
           print *,'*** WARNING *** : velocity grid contains values less than 1 km/sec'

      vgrid(n,m)%nnode=vgrid(n,m)%nr*vgrid(n,m)%nlat*vgrid(n,m)%nlong

 ! allocate and initialize the activity flag

      allocate(vgrid(n,m)%active(vgrid(n,m)%nr,vgrid(n,m)%nlat,vgrid(n,m)%nlong))
      vgrid(n,m)%active = .false.


   end do  ! loop over interpolation regions

end do  ! vtypes

close(10)

end subroutine initialize_velocity_grids

!*************************************************************
! This subroutine reads the initial parametrization of the interface positions
! from which the values on the propagation grid are to be interpolated
! from file into the appropriate structures

subroutine initialize_interfaces
use mod_surf
use variable_def
implicit none

integer :: n,i,j
integer :: nlat,nlong
real(kind=dp) :: dlat0,dlong0,lat0,long0,h,hb

open(10,file='interfaces.in')


! read the number of interfaces
read (10,*) n_interfaces

! allocate space for these interfaces (and the associated intersections and regions for future use)
allocate(intrface(n_interfaces))
do n=1,n_interfaces ; call interface_defaults(intrface(n)) ; intrface(n)%id=n ; end do


! read the grid properties and radius values to be interpolated for the internal interfaces

! grid parameters
   read(10,*) nlat,nlong
   read(10,*) dlat0,dlong0
   read(10,*) lat0,long0


do n=1,n_interfaces

   intrface(n)%nlat = nlat
   intrface(n)%nlong = nlong
   intrface(n)%dlat0 = dlat0
   intrface(n)%dlong0 = dlong0
   intrface(n)%lat0 = lat0
   intrface(n)%long0 = long0


! initialize the grid

   allocate(intrface(n)%lat(intrface(n)%nlat),intrface(n)%long(intrface(n)%nlong))

   do i=1,intrface(n)%nlat
      intrface(n)%lat(i)=intrface(n)%lat0 + (i-1)*intrface(n)%dlat0
   end do

   do i=1,intrface(n)%nlong
      intrface(n)%long(i)=intrface(n)%long0 + (i-1)*intrface(n)%dlong0
   end do


! read in the radius values on the interpolation grid

   allocate(intrface(n)%r(intrface(n)%nlat,intrface(n)%nlong))

   do i=1,intrface(n)%nlat
      do j=1,intrface(n)%nlong
         read (10,*) intrface(n)%r(i,j)
      end do
   end do

   intrface(n)%nnode=intrface(n)%nlat*intrface(n)%nlong

end do  ! loop over interfaces

close(10)

! check top and bottom interfaces are not outside the propagation grid

     !if (count(intrface(1)%r > pgrid%r(pgrid%nr)) > 0) stop ' ERROR: surface above propagation grid'
     !if (count(intrface(n_interfaces)%r < pgrid%r(1)) > 0) stop ' ERROR: bottom below propagation grid'


! correct for intersecting interfaces, higher takes priority

  do n=2,n_interfaces

     do j=1,intrface(n)%nlong
        do i=1,intrface(n)%nlat

    ! higher has priority, EXCEPT bottom interface

           if (n < n_interfaces) then

              hb=intrface(n_interfaces)%r(i,j) 

              if (intrface(n)%r(i,j) < hb) then
                 intrface(n)%r(i,j) = hb
                 intrface(n)%pinched = .true.
                 intrface(n_interfaces)%pinched = .true.
              endif

           endif

    ! check if interface above is crossed

           h=intrface(n-1)%r(i,j)

           if (intrface(n)%r(i,j) > h) then
              intrface(n)%r(i,j) = h
              intrface(n)%pinched = .true.
              intrface(n-1)%pinched = .true.
           endif

        end do
     end do

  end do


end subroutine initialize_interfaces


