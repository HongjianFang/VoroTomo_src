!-----------------------------------------------------------------------------------
module mod_surf
INTEGER, PARAMETER     :: sp = SELECTED_REAL_KIND(6,37)
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)
type Tvelocity_grid   ! regular grid defining the velocity as a function of position
                      ! Note: usually each velocity grid is defined over the entire propagation
                      ! grid, but only some nodes actually influence the corresponding region  

    integer       :: nr,nlong,nlat     ! # of grid cells in each direction
    REAL(KIND=dp) :: dr0,dlat0,dlong0  ! grid step sizes
    REAL(KIND=dp) :: r0,lat0,long0     ! position of grid origin


    integer       :: nnode             ! total # of nodes

    integer       :: start_index       ! for use in inversion if the grid is a velocity grid
    logical       :: to_be_inverted ! for use in inversion if the grid is a velocity grid


    REAL(KIND=dp), DIMENSION (:), pointer :: r 
    REAL(KIND=dp), DIMENSION (:), pointer :: lat 
    REAL(KIND=dp), DIMENSION (:), pointer :: long 


    REAL(KIND=dp), DIMENSION (:,:,:), pointer     :: velocity 
    logical,dimension(:,:,:),pointer              :: active   ! set to true for the nodes that actually
                                                              ! influence the region to which the grid belongs 


end type Tvelocity_grid

type Tinterface
! this type defines the position of an interface
!  note: type Tinterface defines of the position of the interface (used as input by the cubic spline interpolation)
!        type Tintersection contains the actual nodes of the interface and everything related to them

! grid definitions
    integer :: nlat,nlong,id                        ! # of points in lat,long
    REAL(KIND=dp) :: dlat0,dlong0                   ! size of intervals
    REAL(KIND=dp) :: lat0,long0                     ! position of grid origin
    logical       :: pinched                        ! true if the interface touches another


! parametrs of the inversion
    integer       :: nnode                          ! # of parameters describing the interface position    
    integer       :: start_index                    ! start index of the interface parameters in the global parametr list
    logical       :: to_be_inverted        ! is the position of this interface is to be inverted for?


    REAL(KIND=dp), DIMENSION (:), pointer :: lat      
    REAL(KIND=dp), DIMENSION (:), pointer :: long 

    REAL(KIND=dp), DIMENSION (:,:), pointer :: r    ! the actual radius values for the interface at the nodes

end type Tinterface


                                                              
contains
subroutine vgrid_defaults(grid)

    type(Tvelocity_grid)  :: grid

    grid%to_be_inverted  = .false.
    grid%nnode = 0      
    nullify(grid%r) 
    nullify(grid%lat) 
    nullify(grid%long) 
    nullify(grid%velocity)
    nullify(grid%active)

 end subroutine vgrid_defaults

  subroutine interface_defaults(iface)
    type(Tinterface)  :: iface

    iface%pinched = .false.       
    iface%nnode = 0                 
    iface%to_be_inverted = .false. 
    nullify(iface%lat)
    nullify(iface%long)
    nullify(iface%r)

  end subroutine interface_defaults

end module mod_surf



module variable_def
        use mod_surf
        integer,parameter:: iflsph=1
        integer,parameter:: mode=1
	real goxd,gozd
	real dvxd,dvzd
	integer kmaxRc,kmaxRg,kmaxLc,kmaxLg,kmax
	real minthk
	real*8,dimension(:),allocatable::tRc,tRg,tLc,tLg
	REAL, DIMENSION (:), ALLOCATABLE :: dsurf
	real,dimension(:),allocatable::depz
	real,dimension(:),allocatable::obst
	real,dimension(:),allocatable::noises
	integer,parameter::NP=60
	real sta1_lat,sta1_lon,sta2_lat,sta2_lon
	real dist,dcal
	integer dall
	integer istep
	real,parameter::pi=3.1415926535898
	integer nxf,nyf,nzf,checkstat
	real,parameter::ftol=1e-5
	integer igr,iwave
	integer ii,jj,kk
	integer nvx,nvz
        REAL, DIMENSION (:,:), ALLOCATABLE :: scxf,sczf
        REAL, DIMENSION (:,:,:), ALLOCATABLE :: rcxf,rczf
	integer,dimension(:,:),allocatable::wavetype,igrt,nrc1
	integer,dimension(:),allocatable::nsrcsurf1,knum1
	integer,dimension(:,:),allocatable::periods
	character strf
	integer veltp,wavetp
	real velvalue
	integer knum,knumo,err
	integer istep1,istep2
	integer period
	integer knumi,srcnum,count1
        integer nsrcsurf,nrc
        real,dimension(:,:,:),allocatable::vel
        integer ngrid1sep,ngrid2sep
        character line*200
        real mface
        real noiselevel
        integer nx,ny,nz

!.........................................................................
! this array of type Tinterface defines  the location of the interfaces on the propagation grids
! by B-spline interpolation 
  integer                                                  :: n_interfaces
  type(Tinterface),dimension(:),pointer                    :: intrface 


!.........................................................................
! these are the velocity grids used to define the velocity on the propagation grids
! by B-spline interpolation
  integer                                                  :: n_vgrids
  integer                                                  :: n_vtypes
  type(Tvelocity_grid),dimension(:,:),pointer          :: vgrid 
end module variable_def 
