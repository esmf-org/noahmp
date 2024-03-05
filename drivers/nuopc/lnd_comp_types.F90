module lnd_comp_types

  use ESMF         , only: ESMF_Grid, ESMF_Mesh

  use lnd_comp_kind, only: cl => shr_kind_cl
  use lnd_comp_kind, only: r4 => shr_kind_r4
  use lnd_comp_kind, only: r8 => shr_kind_r8

  use Machine
  use NoahmpVarType, only: noahmp_type

  !----------------------------------------------------------------------------
  ! Land NUOPC cap specific data types
  !----------------------------------------------------------------------------

  public

  ! constants
  integer, parameter :: iglobal   = 1 ! mosaic grid, global
  integer, parameter :: iregional = 2 ! mosaic grid, regional
  integer, parameter :: imesh     = 3 ! generic mesh, defined with ESMF mesh file

  ! data type for domain related variables
  type domain_type
     type(ESMF_Grid)            :: grid      ! ESMF grid object, used only for mosaic case and converted to mesh
     type(ESMF_Mesh)            :: mesh      ! ESMF mesh object
     integer                    :: ntiles    ! number of tiles, just for mosaic grid
     integer                    :: layout(2) ! layout for domain decomposition, just for mosaic grid
     integer                    :: dims(2)   ! domain size in both direction, just for regional grid
     integer                    :: gtype     ! grid type: iglb, ireg, imsh
     integer                    :: begl      ! starting index of size in PE
     integer                    :: endl      ! ending index of size in PE
     integer                    :: im        ! number of grid points in each PE
     real(kind=r8), allocatable :: hgt  (:)  ! topography
     integer      , allocatable :: mask (:)  ! mesh land mask: 1 = land, 0 = ocean
     real(kind=r8), allocatable :: frac (:)  ! mesh fractional land
     real(kind=r8), allocatable :: lats (:)  ! mesh latitude
     real(kind=r8), allocatable :: lons (:)  ! mesh longitude
     real(kind=r8), allocatable :: garea(:)  ! mesh cell area
  end type domain_type

  ! data type for coupling level namelist options
  type namelist_type
     character(len=cl) :: mesh_file   ! name of ESMF mesh file
     character(len=cl) :: domain_file ! name of domain file, incl. mask, area and also land fraction
     character(len=cl) :: mosaic_file ! name of mosaic file
     character(len=cl) :: input_dir   ! input directory for tiled files
     integer           :: debug_level ! debug level
  end type namelist_type

  ! data type for forcing
  type forcing_type
    real(kind=kind_noahmp), allocatable :: q1     (:) ! mixing ratio/specific humidty at lowest model layer (kg/kg)
    real(kind=kind_noahmp), allocatable :: t1     (:) ! air temperature (K)
    real(kind=kind_noahmp), allocatable :: u1     (:) ! u-component of wind (m/s)
    real(kind=kind_noahmp), allocatable :: v1     (:) ! v-component of wind (m/s)
    real(kind=kind_noahmp), allocatable :: dswsfc (:) ! downward shortwave radiation (W/m2)
    real(kind=kind_noahmp), allocatable :: dlwflx (:) ! downward longwave radiation (W/m2)
    real(kind=kind_noahmp), allocatable :: pbot   (:) ! bottom layer pressure (Pa)
    real(kind=kind_noahmp), allocatable :: ps     (:) ! surface pressure (Pa)
    real(kind=kind_noahmp), allocatable :: tprcpc (:) ! convective component of precipitation (mm/s)
    real(kind=kind_noahmp), allocatable :: tprcpl (:) ! large-scale component of precipitation (mm/s)
    real(kind=kind_noahmp), allocatable :: snowc  (:) ! convective component of snow fall (mm/s)
    real(kind=kind_noahmp), allocatable :: snowl  (:) ! large-scale component of snow fall (mm/s)
    real(kind=kind_noahmp), allocatable :: graupel(:) ! graupel rate (mm/s)
    real(kind=kind_noahmp), allocatable :: hail   (:) ! hail rate (mm/s)

    !real(kind=kind_noahmp), allocatable :: tskin  (:) ! skin temperature (K)
    !real(kind=kind_noahmp), allocatable :: snet   (:) ! net shortwave radiation (W/m2)
    !real(kind=kind_noahmp), allocatable :: wind   (:) ! wind speed (m/s)
    !real(kind=kind_noahmp), allocatable :: tprcp  (:) ! total precipitation (mm/s)
    !real(kind=kind_noahmp), allocatable :: snow   (:) ! snow fall (mm/s)
    !real(kind=kind_noahmp), allocatable :: vegfrac(:) ! vegetation fraction (unitless, 0-1)
    !real(kind=kind_noahmp), allocatable :: hgt    (:) ! forcing height (m)
    !real(kind=kind_noahmp), allocatable :: prslk1 (:) ! dimensionless Exner function at the lowest model layer
    !real(kind=kind_noahmp), allocatable :: ustar1 (:) ! friction velocity (m/s)
    !real(kind=kind_noahmp), allocatable :: zorl   (:) ! surface roughness (m)
  end type forcing_type

  ! higher level data type for model 
  type model_type
     type(domain_type)   :: domain
     type(namelist_type) :: nmlist
     type(forcing_type)  :: forcing
     type(noahmp_type)   :: noahmp
  end type model_type

  ! field type for I/O
  type field_type
     real(r4), pointer :: ptr1r4(:) => null()   ! data pointer for 1d r4
     real(r8), pointer :: ptr1r8(:) => null()   ! data pointer for 1d r8
     integer , pointer :: ptr1i4(:) => null()   ! data pointer for 1d i4
     real(r4), pointer :: ptr2r4(:,:) => null() ! data pointer for 2d r4
     real(r8), pointer :: ptr2r8(:,:) => null() ! data pointer for 2d r8
     integer , pointer :: ptr2i4(:,:) => null() ! data pointer for 2d i4
     integer           :: id = -999             ! field id
     character(len=cl) :: short_name = ""       ! variable short name
     character(len=cl) :: units = ""            ! variable unit
     character(len=cl) :: long_name = ""        ! variable long name
     character(len=cl) :: zaxis = ""            ! name of z-axis
     integer           :: nlev                  ! number of layers in z-axis
     integer           :: nrec                  ! number of record in file (time axis)
  end type field_type

end module lnd_comp_types
