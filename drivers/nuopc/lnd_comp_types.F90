module lnd_comp_types

  use ESMF         , only: ESMF_Grid, ESMF_Mesh
  use NoahmpVarType, only: noahmp_type
  use lnd_comp_kind, only: cl => shr_kind_cl

  !----------------------------------------------------------------------------
  ! Land NUOPC cap specific data types
  !----------------------------------------------------------------------------

  public

  ! data type for domain related variables
  type domain_type
     type(ESMF_Grid)      :: grid      ! ESMF grid object, used only for mosaic case and converted to mesh
     type(ESMF_Mesh)      :: mesh      ! ESMF mesh object
     integer, allocatable :: mask(:)   ! Mesh land-sea mask: 1 = land, 0 = ocean
     integer              :: ntiles    ! Number of tiles, just for mosaic grid
     integer              :: layout(2) ! Layout for domain decomposition, just for mosaic grid
  end type domain_type

  ! data type for coupling level namelist options
  type namelist_type
     character(len=cl) :: mesh_file   ! name of ESMF mesh file
     character(len=cl) :: domain_file ! name of domain file, incl. mask, area and also land fraction
     character(len=cl) :: mosaic_file ! name of mosaic file
     character(len=cl) :: input_dir   ! input directory for tiled files
     integer           :: debug_level ! debug level
  end type namelist_type

  ! higher level data type for model 
  type model_type
     type(domain_type)   :: domain
     type(namelist_type) :: nmlist
     type(noahmp_type)   :: noahmp
  end type model_type

end module lnd_comp_types
