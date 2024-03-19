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
  integer, parameter :: iMosaic = 1 ! mosaic grid, global
  integer, parameter :: iScrip  = 2 ! generic grid and mesh, defined with SCRIP grid definition file

  ! data type for domain related variables
  type domain_type
     type(ESMF_Grid)            :: grid      ! ESMF grid object, used only for mosaic case and converted to mesh
     type(ESMF_Mesh)            :: mesh      ! ESMF mesh object
     integer                    :: ntiles    ! number of tiles, just for mosaic grid
     integer                    :: layout(2) ! layout for domain decomposition, just for mosaic grid
     integer                    :: dims(2)   ! domain size in both direction, just for regional grid
     integer                    :: dtype     ! domain type: iMosaicG, iMosaicR, iMesh 
     integer                    :: begl      ! starting index of size in PE
     integer                    :: endl      ! ending index of size in PE
     integer                    :: im        ! number of grid points in each PE
     real(kind=r8), allocatable :: hgt  (:)  ! topography
     integer      , allocatable :: mask (:)  ! mesh land mask: 1 = land, 0 = ocean
     real(kind=r8), allocatable :: frac (:)  ! mesh fractional land
     real(kind=r8), allocatable :: latm (:)  ! mesh latitude
     real(kind=r8), allocatable :: lonm (:)  ! mesh longitude
     real(kind=r8), allocatable :: latg (:,:)! grid latitude
     real(kind=r8), allocatable :: long (:,:)! grid longitude
     real(kind=r8), allocatable :: garea(:)  ! mesh cell area
  end type domain_type

  ! data type for coupling
  type coupling_type
     real(kind=r8) :: dt ! coupling time step
  end type coupling_type

  ! data type for coupling level namelist options
  type namelist_type
     ! domain specific options
     character(len=cl) :: scrip_file  ! name of SCRIP grid definition file 
     character(len=cl) :: mosaic_file ! name of mosaic file
     character(len=cl) :: input_dir   ! input directory for tiled files
     logical           :: IsGlobal    ! global vs. regional, default is global

     ! output specific options
     character(len=cl) :: OutputMode  ! model output mode: all, low or mid
     integer           :: OutputFreq  ! model output interval in seconds
     logical           :: TransposeIO ! apply transpose to fields before write

     character(len=cl) :: CaseName   ! name of case
     integer           :: debug_level ! debug level
  end type namelist_type

  ! data type for forcing
  type forcing_type
    real(kind=kind_noahmp), allocatable :: SpecHumidityRefHeight  (:) ! specific humidity (water vapor/moist air) at reference height [kg/kg] 
    real(kind=kind_noahmp), allocatable :: TemperatureAirRefHeight(:) ! air temperature at reference height [K]
    real(kind=kind_noahmp), allocatable :: WindEastwardRefHeight  (:) ! wind speed in eastward direction at reference height [m/s]
    real(kind=kind_noahmp), allocatable :: WindNorthwardRefHeight (:) ! wind speed in northward direction at reference height [m/s]
    real(kind=kind_noahmp), allocatable :: TemperatureSoilBottom  (:) ! bottom boundary condition for soil temperature [K]
    real(kind=kind_noahmp), allocatable :: RadSwDownRefHeight     (:) ! downward shortwave radiation at reference height [W/m2]
    real(kind=kind_noahmp), allocatable :: RadLwDownRefHeight     (:) ! downward longwave radiation at reference height [W/m2]
    real(kind=kind_noahmp), allocatable :: PressureAirRefHeight   (:) ! air pressure at reference height [Pa]
    real(kind=kind_noahmp), allocatable :: PressureAirSurface     (:) ! air pressure at the surface-atmosphere interface level [Pa]
    real(kind=kind_noahmp), allocatable :: PrecipConvRefHeight    (:) ! convective precipitation rate at reference heighta [mm/s]
    real(kind=kind_noahmp), allocatable :: PrecipNonConvRefHeight (:) ! non-convective precipitation rate at reference height [mm/s]
    real(kind=kind_noahmp), allocatable :: PrecipShConvRefHeight  (:) ! shallow convective precipitation rate at reference height [mm/s]
    real(kind=kind_noahmp), allocatable :: PrecipSnowRefHeight    (:) ! snow rate at reference height [mm/s]
    real(kind=kind_noahmp), allocatable :: PrecipGraupelRefHeight (:) ! graupel rate at reference height [mm/s]
    real(kind=kind_noahmp), allocatable :: PrecipHailRefHeight    (:) ! hail rate at reference height [mm/s]
  end type forcing_type

  ! data type for I/O
  type io_type
    character(len=19) :: reference_date   ! reference date
  end type io_type 

  ! higher level data type for model 
  type model_type
     type(domain_type)   :: domain
     type(namelist_type) :: nmlist
     type(coupling_type) :: coupling
     type(forcing_type)  :: forcing
     type(noahmp_type)   :: noahmp
     type(io_type)       :: io

     contains

     procedure, public  :: AllocateInit
     procedure, private :: AllocateInitForcing
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

  ! field list type for import end export fields
  type fld_list_type
     character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
     logical :: connected = .false.
  end type fld_list_type

  ! import/export related variables
  integer, parameter     :: fldsMax = 100
  integer                :: fldsToLnd_num = 0
  integer                :: fldsFrLnd_num = 0
  type(fld_list_type)    :: fldsToLnd(fldsMax)
  type(fld_list_type)    :: fldsFrLnd(fldsMax)

  ! io related variables
  integer, parameter     :: fldsMaxIO = 200
  type(field_type)       :: histflds(fldsMaxIO)
  type(field_type)       :: restflds(fldsMaxIO)

!=============================================================================
contains
!=============================================================================

  subroutine AllocateInit(this)

    ! input/output variables
    class(model_type) :: this 
    !-------------------------------------------------------------------------

    call this%AllocateInitForcing(this%domain%begl, this%domain%endl)

  end subroutine AllocateInit

  !===========================================================================

  subroutine AllocateInitForcing(this, begl, endl)

    class(model_type) :: this
    integer :: begl, endl

    ! allocate
    if (.not. allocated(this%forcing%SpecHumidityRefHeight  )) allocate(this%forcing%SpecHumidityRefHeight  (begl:endl))
    if (.not. allocated(this%forcing%TemperatureAirRefHeight)) allocate(this%forcing%TemperatureAirRefHeight(begl:endl))
    if (.not. allocated(this%forcing%WindEastwardRefHeight  )) allocate(this%forcing%WindEastwardRefHeight  (begl:endl))
    if (.not. allocated(this%forcing%WindNorthwardRefHeight )) allocate(this%forcing%WindNorthwardRefHeight (begl:endl))
    if (.not. allocated(this%forcing%TemperatureSoilBottom  )) allocate(this%forcing%TemperatureSoilBottom  (begl:endl))
    if (.not. allocated(this%forcing%RadSwDownRefHeight     )) allocate(this%forcing%RadSwDownRefHeight     (begl:endl))
    if (.not. allocated(this%forcing%RadLwDownRefHeight     )) allocate(this%forcing%RadLwDownRefHeight     (begl:endl))
    if (.not. allocated(this%forcing%PressureAirRefHeight   )) allocate(this%forcing%PressureAirRefHeight   (begl:endl))
    if (.not. allocated(this%forcing%PressureAirSurface     )) allocate(this%forcing%PressureAirSurface     (begl:endl))
    if (.not. allocated(this%forcing%PrecipConvRefHeight    )) allocate(this%forcing%PrecipConvRefHeight    (begl:endl))
    if (.not. allocated(this%forcing%PrecipNonConvRefHeight )) allocate(this%forcing%PrecipNonConvRefHeight (begl:endl))
    if (.not. allocated(this%forcing%PrecipShConvRefHeight  )) allocate(this%forcing%PrecipShConvRefHeight  (begl:endl))
    if (.not. allocated(this%forcing%PrecipSnowRefHeight    )) allocate(this%forcing%PrecipSnowRefHeight    (begl:endl))
    if (.not. allocated(this%forcing%PrecipGraupelRefHeight )) allocate(this%forcing%PrecipGraupelRefHeight (begl:endl))
    if (.not. allocated(this%forcing%PrecipHailRefHeight    )) allocate(this%forcing%PrecipHailRefHeight    (begl:endl))

    ! init
    this%forcing%SpecHumidityRefHeight  (:) = undefined_real
    this%forcing%TemperatureAirRefHeight(:) = undefined_real
    this%forcing%WindEastwardRefHeight  (:) = undefined_real 
    this%forcing%WindNorthwardRefHeight (:) = undefined_real 
    this%forcing%TemperatureSoilBottom  (:) = undefined_real
    this%forcing%RadSwDownRefHeight     (:) = undefined_real 
    this%forcing%RadLwDownRefHeight     (:) = undefined_real 
    this%forcing%PressureAirRefHeight   (:) = undefined_real 
    this%forcing%PressureAirSurface     (:) = undefined_real 
    this%forcing%PrecipConvRefHeight    (:) = undefined_real 
    this%forcing%PrecipNonConvRefHeight (:) = undefined_real 
    this%forcing%PrecipShConvRefHeight  (:) = undefined_real 
    this%forcing%PrecipSnowRefHeight    (:) = undefined_real 
    this%forcing%PrecipGraupelRefHeight (:) = undefined_real 
    this%forcing%PrecipHailRefHeight    (:) = undefined_real 

  end subroutine AllocateInitForcing

end module lnd_comp_types

