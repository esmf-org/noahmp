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

  ! constants, noahmp table
  integer, parameter :: MVT         = 27     ! number of vegetation types
  integer, parameter :: MBAND       = 2      ! number of radiation bands
  integer, parameter :: MSC         = 8      ! number of soil texture
  integer, parameter :: MAX_SOILTYP = 30     ! max number of soil types
  integer, parameter :: NCROP       = 5      ! number of crop types
  integer, parameter :: NSTAGE      = 8      ! number of crop growth stages
  integer, parameter :: NUM_SLOPE   = 9      ! number of slope

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

  ! data type for noahmp table
  type table_type
    ! vegetation parameters
    character(len=256)         :: VEG_DATASET_DESCRIPTION_TABLE
    integer                    :: NVEG_TABLE                 ! number of vegetation types
    integer                    :: ISURBAN_TABLE              ! urban flag
    integer                    :: ISWATER_TABLE              ! water flag
    integer                    :: ISBARREN_TABLE             ! barren ground flag
    integer                    :: ISICE_TABLE                ! ice flag
    integer                    :: ISCROP_TABLE               ! cropland flag
    integer                    :: EBLFOREST_TABLE            ! evergreen broadleaf forest flag
    integer                    :: NATURAL_TABLE              ! natural vegetation type
    integer                    :: URBTYPE_BEG                ! urban type start number - 1
    integer                    :: LCZ_1_TABLE                ! urban LCZ 1
    integer                    :: LCZ_2_TABLE                ! urban LCZ 2
    integer                    :: LCZ_3_TABLE                ! urban LCZ 3
    integer                    :: LCZ_4_TABLE                ! urban LCZ 4
    integer                    :: LCZ_5_TABLE                ! urban LCZ 5
    integer                    :: LCZ_6_TABLE                ! urban LCZ 6
    integer                    :: LCZ_7_TABLE                ! urban LCZ 7
    integer                    :: LCZ_8_TABLE                ! urban LCZ 8
    integer                    :: LCZ_9_TABLE                ! urban LCZ 9
    integer                    :: LCZ_10_TABLE               ! urban LCZ 10
    integer                    :: LCZ_11_TABLE               ! urban LCZ 11
    real(kind=r8), allocatable :: CH2OP_TABLE    (:)         ! maximum intercepted h2o per unit lai+sai (mm)
    real(kind=r8), allocatable :: DLEAF_TABLE    (:)         ! characteristic leaf dimension (m)
    real(kind=r8), allocatable :: Z0MVT_TABLE    (:)         ! momentum roughness length (m)
    real(kind=r8), allocatable :: HVT_TABLE      (:)         ! top of canopy (m)
    real(kind=r8), allocatable :: HVB_TABLE      (:)         ! bottom of canopy (m)
    real(kind=r8), allocatable :: DEN_TABLE      (:)         ! tree density (no. of trunks per m2)
    real(kind=r8), allocatable :: RC_TABLE       (:)         ! tree crown radius (m)
    real(kind=r8), allocatable :: MFSNO_TABLE    (:)         ! snowmelt curve parameter
    real(kind=r8), allocatable :: SCFFAC_TABLE   (:)         ! snow cover factor (m) (replace original hard-coded 2.5*z0 in SCF formulation)
    real(kind=r8), allocatable :: CBIOM_TABLE    (:)         ! canopy biomass heat capacity parameter (m) 
    real(kind=r8), allocatable :: SAIM_TABLE     (:,:)       ! monthly stem area index, one-sided
    real(kind=r8), allocatable :: LAIM_TABLE     (:,:)       ! monthly leaf area index, one-sided
    real(kind=r8), allocatable :: SLA_TABLE      (:)         ! single-side leaf area per Kg [m2/kg]
    real(kind=r8), allocatable :: DILEFC_TABLE   (:)         ! coeficient for leaf stress death [1/s]
    real(kind=r8), allocatable :: DILEFW_TABLE   (:)         ! coeficient for leaf stress death [1/s]
    real(kind=r8), allocatable :: FRAGR_TABLE    (:)         ! fraction of growth respiration  !original was 0.3 
    real(kind=r8), allocatable :: LTOVRC_TABLE   (:)         ! leaf turnover [1/s]
    real(kind=r8), allocatable :: C3PSN_TABLE    (:)         ! photosynthetic pathway: 0. = c4, 1. = c3
    real(kind=r8), allocatable :: KC25_TABLE     (:)         ! co2 michaelis-menten constant at 25c (pa)
    real(kind=r8), allocatable :: AKC_TABLE      (:)         ! q10 for kc25
    real(kind=r8), allocatable :: KO25_TABLE     (:)         ! o2 michaelis-menten constant at 25c (pa)
    real(kind=r8), allocatable :: AKO_TABLE      (:)         ! q10 for ko25
    real(kind=r8), allocatable :: VCMX25_TABLE   (:)         ! maximum rate of carboxylation at 25c (umol co2/m2/s)
    real(kind=r8), allocatable :: AVCMX_TABLE    (:)         ! q10 for vcmx25
    real(kind=r8), allocatable :: BP_TABLE       (:)         ! minimum leaf conductance (umol/m2/s)
    real(kind=r8), allocatable :: MP_TABLE       (:)         ! slope of conductance-to-photosynthesis relationship
    real(kind=r8), allocatable :: QE25_TABLE     (:)         ! quantum efficiency at 25c (umol co2 / umol photon)
    real(kind=r8), allocatable :: AQE_TABLE      (:)         ! q10 for qe25
    real(kind=r8), allocatable :: RMF25_TABLE    (:)         ! leaf maintenance respiration at 25c (umol co2/m2/s)
    real(kind=r8), allocatable :: RMS25_TABLE    (:)         ! stem maintenance respiration at 25c (umol co2/kg bio/s)
    real(kind=r8), allocatable :: RMR25_TABLE    (:)         ! root maintenance respiration at 25c (umol co2/kg bio/s)
    real(kind=r8), allocatable :: ARM_TABLE      (:)         ! q10 for maintenance respiration
    real(kind=r8), allocatable :: FOLNMX_TABLE   (:)         ! foliage nitrogen concentration when f(n)=1 (%)
    real(kind=r8), allocatable :: TMIN_TABLE     (:)         ! minimum temperature for photosynthesis (k)
    real(kind=r8), allocatable :: XL_TABLE       (:)         ! leaf/stem orientation index
    real(kind=r8), allocatable :: RHOL_TABLE     (:,:)       ! leaf reflectance: 1=vis, 2=nir
    real(kind=r8), allocatable :: RHOS_TABLE     (:,:)       ! stem reflectance: 1=vis, 2=nir
    real(kind=r8), allocatable :: TAUL_TABLE     (:,:)       ! leaf transmittance: 1=vis, 2=nir
    real(kind=r8), allocatable :: TAUS_TABLE     (:,:)       ! stem transmittance: 1=vis, 2=nir
    real(kind=r8), allocatable :: MRP_TABLE      (:)         ! microbial respiration parameter (umol co2 /kg c/ s)
    real(kind=r8), allocatable :: CWPVT_TABLE    (:)         ! empirical canopy wind parameter
    real(kind=r8), allocatable :: WRRAT_TABLE    (:)         ! wood to non-wood ratio
    real(kind=r8), allocatable :: WDPOOL_TABLE   (:)         ! wood pool (switch 1 or 0) depending on woody or not [-]
    real(kind=r8), allocatable :: TDLEF_TABLE    (:)         ! characteristic T for leaf freezing [K]
    real(kind=r8), allocatable :: NROOT_TABLE    (:)         ! number of soil layers with root present
    real(kind=r8), allocatable :: RGL_TABLE      (:)         ! Parameter used in radiation stress function
    real(kind=r8), allocatable :: RS_TABLE       (:)         ! Minimum stomatal resistance [s m-1]
    real(kind=r8), allocatable :: HS_TABLE       (:)         ! Parameter used in vapor pressure deficit function
    real(kind=r8), allocatable :: TOPT_TABLE     (:)         ! Optimum transpiration air temperature [K]
    real(kind=r8), allocatable :: RSMAX_TABLE    (:)         ! Maximal stomatal resistance [s m-1]
    real(kind=r8), allocatable :: RTOVRC_TABLE   (:)         ! root turnover coefficient [1/s]
    real(kind=r8), allocatable :: RSWOODC_TABLE  (:)         ! wood respiration coeficient [1/s]
    real(kind=r8), allocatable :: BF_TABLE       (:)         ! parameter for present wood allocation [-]
    real(kind=r8), allocatable :: WSTRC_TABLE    (:)         ! water stress coeficient [-]
    real(kind=r8), allocatable :: LAIMIN_TABLE   (:)         ! minimum leaf area index [m2/m2]
    real(kind=r8), allocatable :: XSAMIN_TABLE   (:)         ! minimum stem area index [m2/m2]
    ! radiation parameters
    real(kind=r8), allocatable :: ALBSAT_TABLE   (:,:)       ! saturated soil albedos: 1=vis, 2=nir
    real(kind=r8), allocatable :: ALBDRY_TABLE   (:,:)       ! dry soil albedos: 1=vis, 2=nir
    real(kind=r8), allocatable :: ALBICE_TABLE   (:)         ! albedo land ice: 1=vis, 2=nir
    real(kind=r8), allocatable :: ALBLAK_TABLE   (:)         ! albedo frozen lakes: 1=vis, 2=nir
    real(kind=r8), allocatable :: OMEGAS_TABLE   (:)         ! two-stream parameter omega for snow
    real(kind=r8)              :: BETADS_TABLE               ! two-stream parameter betad for snow
    real(kind=r8)              :: BETAIS_TABLE               ! two-stream parameter betad for snow
    real(kind=r8), allocatable :: EG_TABLE       (:)         ! emissivity soil surface
    real(kind=r8)              :: EICE_TABLE
    ! global parameters
    real(kind=r8)              :: CO2_TABLE                  ! co2 partial pressure
    real(kind=r8)              :: O2_TABLE                   ! o2 partial pressure
    real(kind=r8)              :: TIMEAN_TABLE               ! gridcell mean topgraphic index (global mean)
    real(kind=r8)              :: FSATMX_TABLE               ! maximum surface saturated fraction (global mean)
    real(kind=r8)              :: Z0SNO_TABLE                ! snow surface roughness length (m) (0.002)
    real(kind=r8)              :: SSI_TABLE                  ! liquid water holding capacity for snowpack (m3/m3) (0.03)
    real(kind=r8)              :: SNOW_RET_FAC_TABLE         ! snowpack water release timescale factor (1/s)
    real(kind=r8)              :: SNOW_EMIS_TABLE            ! snow emissivity
    real(kind=r8)              :: SWEMX_TABLE                ! new snow mass to fully cover old snow (mm)
    real(kind=r8)              :: TAU0_TABLE                 ! tau0 from Yang97 eqn. 10a
    real(kind=r8)              :: GRAIN_GROWTH_TABLE         ! growth from vapor diffusion Yang97 eqn. 10b
    real(kind=r8)              :: EXTRA_GROWTH_TABLE         ! extra growth near freezing Yang97 eqn. 10c
    real(kind=r8)              :: DIRT_SOOT_TABLE            ! dirt and soot term Yang97 eqn. 10d
    real(kind=r8)              :: BATS_COSZ_TABLE            ! zenith angle snow albedo adjustment; b in Yang97 eqn. 15
    real(kind=r8)              :: BATS_VIS_NEW_TABLE         ! new snow visible albedo
    real(kind=r8)              :: BATS_NIR_NEW_TABLE         ! new snow NIR albedo
    real(kind=r8)              :: BATS_VIS_AGE_TABLE         ! age factor for diffuse visible snow albedo Yang97 eqn. 17
    real(kind=r8)              :: BATS_NIR_AGE_TABLE         ! age factor for diffuse NIR snow albedo Yang97 eqn. 18
    real(kind=r8)              :: BATS_VIS_DIR_TABLE         ! cosz factor for direct visible snow albedo Yang97 eqn. 15
    real(kind=r8)              :: BATS_NIR_DIR_TABLE         ! cosz factor for direct NIR snow albedo Yang97 eqn. 16
    real(kind=r8)              :: RSURF_SNOW_TABLE           ! surface resistance for snow(s/m)
    real(kind=r8)              :: RSURF_EXP_TABLE            ! exponent in the shape parameter for soil resistance option 1
    real(kind=r8)              :: C2_SNOWCOMPACT_TABLE       ! overburden snow compaction parameter (m3/kg)
    real(kind=r8)              :: C3_SNOWCOMPACT_TABLE       ! snow desctructive metamorphism compaction parameter1 [1/s]
    real(kind=r8)              :: C4_SNOWCOMPACT_TABLE       ! snow desctructive metamorphism compaction parameter2 [1/k]
    real(kind=r8)              :: C5_SNOWCOMPACT_TABLE       ! snow desctructive metamorphism compaction parameter3
    real(kind=r8)              :: DM_SNOWCOMPACT_TABLE       ! upper Limit on destructive metamorphism compaction [kg/m3]
    real(kind=r8)              :: ETA0_SNOWCOMPACT_TABLE     ! snow viscosity coefficient [kg-s/m2]
    real(kind=r8)              :: SNLIQMAXFRAC_TABLE         ! maximum liquid water fraction in snow
    real(kind=r8)              :: SWEMAXGLA_TABLE            ! Maximum SWE allowed at glaciers (mm)
    real(kind=r8)              :: WSLMAX_TABLE               ! maximum lake water storage (mm)
    real(kind=r8)              :: ROUS_TABLE                 ! specific yield [-] for Niu et al. 2007 groundwater scheme
    real(kind=r8)              :: CMIC_TABLE                 ! microprore content (0.0-1.0), 0.0: close to free drainage
    real(kind=r8)              :: SNOWDEN_MAX_TABLE          ! maximum fresh snowfall density (kg/m3)
    real(kind=r8)              :: CLASS_ALB_REF_TABLE        ! reference snow albedo in CLASS scheme
    real(kind=r8)              :: CLASS_SNO_AGE_TABLE        ! snow aging e-folding time (s) in CLASS albedo scheme
    real(kind=r8)              :: CLASS_ALB_NEW_TABLE        ! fresh snow albedo in CLASS scheme
    real(kind=r8)              :: PSIWLT_TABLE               ! soil metric potential for wilting point (m)
    real(kind=r8)              :: Z0SOIL_TABLE               ! Bare-soil roughness length (m) (i.e., under the canopy)
    real(kind=r8)              :: Z0LAKE_TABLE               ! Lake surface roughness length (m)
    ! irrigation parameters
    integer                    :: IRR_HAR_TABLE              ! number of days before harvest date to stop irrigation 
    real(kind=r8)              :: IRR_FRAC_TABLE             ! irrigation Fraction
    real(kind=r8)              :: IRR_LAI_TABLE              ! Minimum lai to trigger irrigation
    real(kind=r8)              :: IRR_MAD_TABLE              ! management allowable deficit (0-1)
    real(kind=r8)              :: FILOSS_TABLE               ! factor of flood irrigation loss
    real(kind=r8)              :: SPRIR_RATE_TABLE           ! mm/h, sprinkler irrigation rate
    real(kind=r8)              :: MICIR_RATE_TABLE           ! mm/h, micro irrigation rate
    real(kind=r8)              :: FIRTFAC_TABLE              ! flood application rate factor
    real(kind=r8)              :: IR_RAIN_TABLE              ! maximum precipitation to stop irrigation trigger
    ! tile drainage parameters
    integer                    :: DRAIN_LAYER_OPT_TABLE      ! tile drainage layer
    integer      , allocatable :: TD_DEPTH_TABLE (:)         ! tile drainage depth (layer number) from soil surface
    real(kind=r8), allocatable :: TDSMC_FAC_TABLE(:)         ! tile drainage soil moisture factor
    real(kind=r8), allocatable :: TD_DC_TABLE    (:)         ! tile drainage coefficient [mm/d]
    real(kind=r8), allocatable :: TD_DCOEF_TABLE (:)         ! tile drainage coefficient [mm/d]
    real(kind=r8), allocatable :: TD_D_TABLE     (:)         ! depth to impervious layer from drain water level [m]
    real(kind=r8), allocatable :: TD_ADEPTH_TABLE(:)         ! actual depth of impervious layer from land surface [m]
    real(kind=r8), allocatable :: TD_RADI_TABLE  (:)         ! effective radius of drain tubes [m]
    real(kind=r8), allocatable :: TD_SPAC_TABLE  (:)         ! distance between two drain tubes or tiles [m]
    real(kind=r8), allocatable :: TD_DDRAIN_TABLE(:)         ! tile drainage depth [m]
    real(kind=r8), allocatable :: KLAT_FAC_TABLE (:)         ! hydraulic conductivity mutiplification factor
    ! crop parameters
    integer                    :: DEFAULT_CROP_TABLE         ! Default crop index
    integer      , allocatable :: PLTDAY_TABLE   (:)         ! Planting date
    integer      , allocatable :: HSDAY_TABLE    (:)         ! Harvest date
    real(kind=r8), allocatable :: PLANTPOP_TABLE (:)         ! Plant density [per ha] - used?
    real(kind=r8), allocatable :: IRRI_TABLE     (:)         ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
    real(kind=r8), allocatable :: GDDTBASE_TABLE (:)         ! Base temperature for GDD accumulation [C]
    real(kind=r8), allocatable :: GDDTCUT_TABLE  (:)         ! Upper temperature for GDD accumulation [C]
    real(kind=r8), allocatable :: GDDS1_TABLE    (:)         ! GDD from seeding to emergence
    real(kind=r8), allocatable :: GDDS2_TABLE    (:)         ! GDD from seeding to initial vegetative 
    real(kind=r8), allocatable :: GDDS3_TABLE    (:)         ! GDD from seeding to post vegetative 
    real(kind=r8), allocatable :: GDDS4_TABLE    (:)         ! GDD from seeding to intial reproductive
    real(kind=r8), allocatable :: GDDS5_TABLE    (:)         ! GDD from seeding to pysical maturity 
    real(kind=r8), allocatable :: C3PSNI_TABLE   (:)         ! photosynthetic pathway: 0. = c4, 1. = c3 ! Zhe Zhang 2020-07-03
    real(kind=r8), allocatable :: KC25I_TABLE    (:)         ! co2 michaelis-menten constant at 25c (pa)
    real(kind=r8), allocatable :: AKCI_TABLE     (:)         ! q10 for kc25
    real(kind=r8), allocatable :: KO25I_TABLE    (:)         ! o2 michaelis-menten constant at 25c (pa)
    real(kind=r8), allocatable :: AKOI_TABLE     (:)         ! q10 for ko25
    real(kind=r8), allocatable :: VCMX25I_TABLE  (:)         ! maximum rate of carboxylation at 25c (umol co2/m2/s)
    real(kind=r8), allocatable :: AVCMXI_TABLE   (:)         ! q10 for vcmx25
    real(kind=r8), allocatable :: BPI_TABLE      (:)         ! minimum leaf conductance (umol/m2/s)
    real(kind=r8), allocatable :: MPI_TABLE      (:)         ! slope of conductance-to-photosynthesis relationship
    real(kind=r8), allocatable :: QE25I_TABLE    (:)         ! quantum efficiency at 25c (umol co2 / umol photon)
    real(kind=r8), allocatable :: FOLNMXI_TABLE  (:)         ! foliage nitrogen concentration when f(n)=1 (%)
    real(kind=r8), allocatable :: AREF_TABLE     (:)         ! reference maximum CO2 assimulation rate 
    real(kind=r8), allocatable :: PSNRF_TABLE    (:)         ! CO2 assimulation reduction factor(0-1) (caused by non-modeled part, pest,weeds)
    real(kind=r8), allocatable :: I2PAR_TABLE    (:)         ! Fraction of incoming solar radiation to photosynthetically active radiation
    real(kind=r8), allocatable :: TASSIM0_TABLE  (:)         ! Minimum temperature for CO2 assimulation [C]
    real(kind=r8), allocatable :: TASSIM1_TABLE  (:)         ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
    real(kind=r8), allocatable :: TASSIM2_TABLE  (:)         ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
    real(kind=r8), allocatable :: K_TABLE        (:)         ! light extinction coefficient
    real(kind=r8), allocatable :: EPSI_TABLE     (:)         ! initial light use efficiency
    real(kind=r8), allocatable :: Q10MR_TABLE    (:)         ! q10 for maintainance respiration
    real(kind=r8), allocatable :: LEFREEZ_TABLE  (:)         ! characteristic T for leaf freezing [K]
    real(kind=r8), allocatable :: DILE_FC_TABLE  (:,:)       ! coeficient for temperature leaf stress death [1/s]
    real(kind=r8), allocatable :: DILE_FW_TABLE  (:,:)       ! coeficient for water leaf stress death [1/s]
    real(kind=r8), allocatable :: FRA_GR_TABLE   (:)         ! fraction of growth respiration
    real(kind=r8), allocatable :: LF_OVRC_TABLE  (:,:)       ! fraction of leaf turnover  [1/s]
    real(kind=r8), allocatable :: ST_OVRC_TABLE  (:,:)       ! fraction of stem turnover  [1/s]
    real(kind=r8), allocatable :: RT_OVRC_TABLE  (:,:)       ! fraction of root tunrover  [1/s]
    real(kind=r8), allocatable :: LFMR25_TABLE   (:)         ! leaf maintenance respiration at 25C [umol CO2/m2/s]
    real(kind=r8), allocatable :: STMR25_TABLE   (:)         ! stem maintenance respiration at 25C [umol CO2/kg bio/s]
    real(kind=r8), allocatable :: RTMR25_TABLE   (:)         ! root maintenance respiration at 25C [umol CO2/kg bio/s]
    real(kind=r8), allocatable :: GRAINMR25_TABLE(:)         ! grain maintenance respiration at 25C [umol CO2/kg bio/s]
    real(kind=r8), allocatable :: LFPT_TABLE     (:,:)       ! fraction of carbohydrate flux to leaf
    real(kind=r8), allocatable :: STPT_TABLE     (:,:)       ! fraction of carbohydrate flux to stem
    real(kind=r8), allocatable :: RTPT_TABLE     (:,:)       ! fraction of carbohydrate flux to root
    real(kind=r8), allocatable :: GRAINPT_TABLE  (:,:)       ! fraction of carbohydrate flux to grain
    real(kind=r8), allocatable :: LFCT_TABLE     (:,:)       ! fraction of carbohydrate translocation from leaf to grain 
    real(kind=r8), allocatable :: STCT_TABLE     (:,:)       ! fraction of carbohydrate translocation from stem to grain
    real(kind=r8), allocatable :: RTCT_TABLE     (:,:)       ! fraction of carbohydrate translocation from root to grain
    real(kind=r8), allocatable :: BIO2LAI_TABLE  (:)         ! leaf area per living leaf biomass [m2/kg]
    ! soil parameters
    integer                    :: SLCATS_TABLE               ! number of soil categories
    real(kind=r8), allocatable :: BEXP_TABLE     (:)         ! soil B parameter
    real(kind=r8), allocatable :: SMCDRY_TABLE   (:)         ! dry soil moisture threshold
    real(kind=r8), allocatable :: SMCMAX_TABLE   (:)         ! porosity, saturated value of soil moisture (volumetric)
    real(kind=r8), allocatable :: SMCREF_TABLE   (:)         ! reference soil moisture (field capacity) (volumetric)
    real(kind=r8), allocatable :: PSISAT_TABLE   (:)         ! saturated soil matric potential
    real(kind=r8), allocatable :: DKSAT_TABLE    (:)         ! saturated soil hydraulic conductivity
    real(kind=r8), allocatable :: DWSAT_TABLE    (:)         ! saturated soil hydraulic diffusivity
    real(kind=r8), allocatable :: SMCWLT_TABLE   (:)         ! wilting point soil moisture (volumetric)
    real(kind=r8), allocatable :: QUARTZ_TABLE   (:)         ! soil quartz content
    real(kind=r8), allocatable :: BVIC_TABLE     (:)         ! VIC model infiltration parameter (-) for opt_run=6
    real(kind=r8), allocatable :: AXAJ_TABLE     (:)         ! Xinanjiang: Tension water distribution inflection parameter [-] for opt_run=7
    real(kind=r8), allocatable :: BXAJ_TABLE     (:)         ! Xinanjiang: Tension water distribution shape parameter [-] for opt_run=7
    real(kind=r8), allocatable :: XXAJ_TABLE     (:)         ! Xinanjiang: Free water distribution shape parameter [-] for opt_run=7
    real(kind=r8), allocatable :: BDVIC_TABLE    (:)         ! VIC model infiltration parameter (-)
    real(kind=r8), allocatable :: GDVIC_TABLE    (:)         ! mean capilary drive (m)
    real(kind=r8), allocatable :: BBVIC_TABLE    (:)         ! heterogeniety parameter for DVIC infiltration [-]
    ! general parameters
    real(kind=r8), allocatable :: SLOPE_TABLE    (:)         ! slope factor for soil drainage
    real(kind=r8)              :: CSOIL_TABLE                ! Soil heat capacity [J m-3 K-1]
    real(kind=r8)              :: REFDK_TABLE                ! Parameter in the surface runoff parameterization
    real(kind=r8)              :: REFKDT_TABLE               ! Parameter in the surface runoff parameterization
    real(kind=r8)              :: FRZK_TABLE                 ! Frozen ground parameter
    real(kind=r8)              :: ZBOT_TABLE                 ! Depth [m] of lower boundary soil temperature
    real(kind=r8)              :: CZIL_TABLE                 ! Parameter used in the calculation of the roughness length for heat
    ! optional parameters
    real(kind=r8)              :: SR2006_THETA_1500T_A_TABLE ! sand coefficient
    real(kind=r8)              :: SR2006_THETA_1500T_B_TABLE ! clay coefficient
    real(kind=r8)              :: SR2006_THETA_1500T_C_TABLE ! orgm coefficient
    real(kind=r8)              :: SR2006_THETA_1500T_D_TABLE ! sand*orgm coefficient
    real(kind=r8)              :: SR2006_THETA_1500T_E_TABLE ! clay*orgm coefficient
    real(kind=r8)              :: SR2006_THETA_1500T_F_TABLE ! sand*clay coefficient
    real(kind=r8)              :: SR2006_THETA_1500T_G_TABLE ! constant adjustment
    real(kind=r8)              :: SR2006_THETA_1500_A_TABLE  ! theta_1500t coefficient
    real(kind=r8)              :: SR2006_THETA_1500_B_TABLE  ! constant adjustment
    real(kind=r8)              :: SR2006_THETA_33T_A_TABLE   ! sand coefficient
    real(kind=r8)              :: SR2006_THETA_33T_B_TABLE   ! clay coefficient
    real(kind=r8)              :: SR2006_THETA_33T_C_TABLE   ! orgm coefficient
    real(kind=r8)              :: SR2006_THETA_33T_D_TABLE   ! sand*orgm coefficient
    real(kind=r8)              :: SR2006_THETA_33T_E_TABLE   ! clay*orgm coefficient
    real(kind=r8)              :: SR2006_THETA_33T_F_TABLE   ! sand*clay coefficient
    real(kind=r8)              :: SR2006_THETA_33T_G_TABLE   ! constant adjustment
    real(kind=r8)              :: SR2006_THETA_33_A_TABLE    ! theta_33t*theta_33t coefficient
    real(kind=r8)              :: SR2006_THETA_33_B_TABLE    ! theta_33t coefficient
    real(kind=r8)              :: SR2006_THETA_33_C_TABLE    ! constant adjustment
    real(kind=r8)              :: SR2006_THETA_S33T_A_TABLE  ! sand coefficient
    real(kind=r8)              :: SR2006_THETA_S33T_B_TABLE  ! clay coefficient
    real(kind=r8)              :: SR2006_THETA_S33T_C_TABLE  ! orgm coefficient
    real(kind=r8)              :: SR2006_THETA_S33T_D_TABLE  ! sand*orgm coefficient
    real(kind=r8)              :: SR2006_THETA_S33T_E_TABLE  ! clay*orgm coefficient
    real(kind=r8)              :: SR2006_THETA_S33T_F_TABLE  ! sand*clay coefficient
    real(kind=r8)              :: SR2006_THETA_S33T_G_TABLE  ! constant adjustment
    real(kind=r8)              :: SR2006_THETA_S33_A_TABLE   ! theta_s33t coefficient
    real(kind=r8)              :: SR2006_THETA_S33_B_TABLE   ! constant adjustment
    real(kind=r8)              :: SR2006_PSI_ET_A_TABLE      ! sand coefficient
    real(kind=r8)              :: SR2006_PSI_ET_B_TABLE      ! clay coefficient
    real(kind=r8)              :: SR2006_PSI_ET_C_TABLE      ! theta_s33 coefficient
    real(kind=r8)              :: SR2006_PSI_ET_D_TABLE      ! sand*theta_s33 coefficient
    real(kind=r8)              :: SR2006_PSI_ET_E_TABLE      ! clay*theta_s33 coefficient
    real(kind=r8)              :: SR2006_PSI_ET_F_TABLE      ! sand*clay coefficient
    real(kind=r8)              :: SR2006_PSI_ET_G_TABLE      ! constant adjustment
    real(kind=r8)              :: SR2006_PSI_E_A_TABLE       ! psi_et*psi_et coefficient
    real(kind=r8)              :: SR2006_PSI_E_B_TABLE       ! psi_et coefficient
    real(kind=r8)              :: SR2006_PSI_E_C_TABLE       ! constant adjustment
    real(kind=r8)              :: SR2006_SMCMAX_A_TABLE      ! sand adjustment
    real(kind=r8)              :: SR2006_SMCMAX_B_TABLE      ! constant adjustment
  end type table_type

  ! data type for coupling
  type coupling_type
     real(kind=r8)              :: MainTimeStep             ! coupling time step
     real(kind=r8), allocatable :: CosSolarZenithAngle(:)   ! cosine of solar zenith angle
     integer      , allocatable :: SoilType           (:)   ! soil type for each soil layer
     integer      , allocatable :: VegType            (:)   ! vegetation type 
     integer      , allocatable :: RunoffSlopeType    (:)   ! underground runoff slope term
     real(kind=r8), allocatable :: VegFracAnnMax      (:)   ! annual max vegetation fraction
     real(kind=r8), allocatable :: TemperatureSoilSnow(:,:) ! snow/soil temperature, K
     !real(kind=r8), allocatable :: SnowDepth          (:)   ! snow depth [mm]
     !real(kind=r8), allocatable :: SnowWaterEquiv     (:)   ! snow water equivalent (ice+liquid) [mm]
     integer      , allocatable :: NumSnowLayerNeg    (:)   ! actual number of snow layers (negative)

     ! water state variables
     real(kind=r8), allocatable :: CanopyLiqWater         (:) ! intercepted canopy liquid water [mm] 
     real(kind=r8), allocatable :: CanopyIce              (:) ! intercepted canopy ice [mm]
     real(kind=r8), allocatable :: CanopyWetFrac          (:) ! wetted or snowed fraction of the canopy
     real(kind=r8), allocatable :: SnowWaterEquiv         (:) ! snow water equivalent (ice+liquid) [mm]
     real(kind=r8), allocatable :: SnowWaterEquivPrev     (:) ! snow water equivalent at previous time step (mm)
     real(kind=r8), allocatable :: SnowDepth              (:) ! snow depth [m]
     real(kind=r8), allocatable :: IrrigationFracFlood    (:) ! fraction of grid under flood irrigation (0 to 1)
     real(kind=r8), allocatable :: IrrigationAmtFlood     (:) ! flood irrigation water amount [m]
     real(kind=r8), allocatable :: IrrigationFracMicro    (:) ! fraction of grid under micro irrigation (0 to 1)
     real(kind=r8), allocatable :: IrrigationAmtMicro     (:) ! micro irrigation water amount [m]
     real(kind=r8), allocatable :: IrrigationFracSprinkler(:) ! fraction of grid under sprinkler irrigation (0 to 1)
     real(kind=r8), allocatable :: IrrigationAmtSprinkler (:) ! sprinkler irrigation water amount [m]
     real(kind=r8), allocatable :: WaterTableDepth        (:) ! water table depth [m]
     real(kind=r8), allocatable :: SoilMoistureToWT       (:) ! soil moisture between bottom of the soil and the water table
     real(kind=r8), allocatable :: TileDrainFrac          (:) ! tile drainage fraction
     real(kind=r8), allocatable :: WaterStorageAquifer    (:) ! water storage in aquifer [mm]
     real(kind=r8), allocatable :: WaterStorageSoilAqf    (:) ! water storage in aquifer + saturated soil [mm]
     real(kind=r8), allocatable :: WaterStorageLake       (:) ! water storage in lake (can be negative) [mm] 
     real(kind=r8), allocatable :: IrrigationFracGrid     (:) ! total irrigation fraction from input for a grid
     integer      , allocatable :: IrrigationCntSprinkler (:) ! irrigation event number, Sprinkler
     integer      , allocatable :: IrrigationCntMicro     (:) ! irrigation event number, Micro
     integer      , allocatable :: IrrigationCntFlood     (:) ! irrigation event number, Flood 
     real(kind=r8), allocatable :: SnowIce                (:,:) ! snow layer ice [mm]
     real(kind=r8), allocatable :: SnowLiqWater           (:,:) ! snow layer liquid water [mm]
     real(kind=r8), allocatable :: SoilLiqWater           (:,:) ! soil liquid moisture [m3/m3]
     real(kind=r8), allocatable :: SoilMoisture           (:,:) ! total soil moisture [m3/m3]
     real(kind=r8), allocatable :: SoilMoistureEqui       (:,:) ! equilibrium soil water  content [m3/m3]
     real(kind=r8), allocatable :: RechargeGwDeepWT       (:)   ! groundwater recharge to or from the water table when deep [m]
     real(kind=r8), allocatable :: RechargeGwShallowWT    (:)   ! groundwater recharge to or from shallow water table [m]
     ! water flux variables
     real(kind=r8), allocatable :: EvapSoilSfcLiqAcc      (:)   ! accumulated soil surface water evaporation per soil timestep [m/s * dt_soil/dt_main]
     real(kind=r8), allocatable :: SoilSfcInflowAcc       (:)   ! accumulated water input on soil surface per soil timestep [m/s * dt_soil/dt_main]
     real(kind=r8), allocatable :: SfcWaterTotChgAcc      (:)   ! accumulated snow,soil,canopy water change per soil timestep [mm]
     real(kind=r8), allocatable :: PrecipTotAcc           (:)   ! accumulated precipitation per soil timestep [mm]
     real(kind=r8), allocatable :: EvapCanopyNetAcc       (:)   ! accumulated net evaporation of canopy intercepted water per soil timestep [mm]
     real(kind=r8), allocatable :: TranspirationAcc       (:)   ! accumulated transpiration per soil timestep [mm]
     real(kind=r8), allocatable :: EvapGroundNetAcc       (:)   ! accumulated net ground (soil/snow) evaporation per soil timestep [mm]
     real(kind=r8), allocatable :: TranspWatLossSoilAcc   (:,:) ! accumulated transpiration water loss from soil per soil timestep [m/s * dt_soil/dt_main]
     ! water parameter variables
     integer      , allocatable :: DrainSoilLayerInd      (:)   ! starting soil layer for drainage 
     real(kind=r8), allocatable :: CanopyLiqHoldCap       (:)   ! maximum canopy intercepted liquid water per unit veg area index [mm]
     real(kind=r8), allocatable :: SnowCompactBurdenFac   (:)   ! overburden snow compaction parameter [m3/kg]
     real(kind=r8), allocatable :: SnowCompactAgingFac1   (:)   ! snow desctructive metamorphism compaction parameter1 [1/s]
     real(kind=r8), allocatable :: SnowCompactAgingFac2   (:)   ! snow desctructive metamorphism compaction parameter2 [1/k]
     real(kind=r8), allocatable :: SnowCompactAgingFac3   (:)   ! snow desctructive metamorphism compaction parameter3
     real(kind=r8), allocatable :: SnowCompactAgingMax    (:)   ! upper Limit on destructive metamorphism compaction [kg/m3]
     real(kind=r8), allocatable :: SnowViscosityCoeff     (:)   ! snow viscosity coefficient [kg-s/m2], Anderson1979: 0.52e6~1.38e6
     real(kind=r8), allocatable :: SnowLiqFracMax         (:)   ! maximum liquid water fraction in snow
     real(kind=r8), allocatable :: SnowLiqHoldCap         (:)   ! liquid water holding capacity for snowpack [m3/m3]
     real(kind=r8), allocatable :: SnowLiqReleaseFac      (:)   ! snowpack water release timescale factor [1/s]
     real(kind=r8), allocatable :: IrriFloodRateFac       (:)   ! flood irrigation application rate factor
     real(kind=r8), allocatable :: IrriMicroRate          (:)   ! micro irrigation rate [mm/hr]
     real(kind=r8), allocatable :: SoilConductivityRef    (:)   ! reference Soil Conductivity parameter (used in runoff formulation)
     real(kind=r8), allocatable :: SoilInfilFacRef        (:)   ! reference Soil Infiltration Parameter (used in runoff formulation)
     real(kind=r8), allocatable :: GroundFrzCoeff         (:)   ! frozen ground parameter to compute frozen soil impervious fraction
     real(kind=r8), allocatable :: GridTopoIndex          (:)   ! gridcell mean topgraphic index (global mean)
     real(kind=r8), allocatable :: SoilSfcSatFracMax      (:)   ! maximum surface soil saturated fraction (global mean)
     real(kind=r8), allocatable :: SpecYieldGw            (:)   ! specific yield [-] for Niu et al. 2007 groundwater scheme
     real(kind=r8), allocatable :: MicroPoreContent       (:)   ! microprore content (0.0-1.0), 0.0: close to free drainage
     real(kind=r8), allocatable :: WaterStorageLakeMax    (:)   ! maximum lake water storage [mm]
     real(kind=r8), allocatable :: SnoWatEqvMaxGlacier    (:)   ! maximum SWE allowed at glaciers [mm]
     integer      , allocatable :: IrriStopDayBfHarvest   (:)   ! number of days before harvest date to stop irrigation
     real(kind=r8), allocatable :: IrriTriggerLaiMin      (:)   ! minimum lai to trigger irrigation
     real(kind=r8), allocatable :: SoilWatDeficitAllow    (:)   ! management allowable deficit (0-1)
     real(kind=r8), allocatable :: IrriFloodLossFrac      (:)   ! factor of flood irrigation loss
     real(kind=r8), allocatable :: IrriSprinklerRate      (:)   ! sprinkler irrigation rate [mm/h]
     real(kind=r8), allocatable :: IrriFracThreshold      (:)   ! irrigation Fraction threshold in a grid
     real(kind=r8), allocatable :: IrriStopPrecipThr      (:)   ! precipitation threshold [mm/hr] to stop irrigation trigger
     real(kind=r8), allocatable :: SnowfallDensityMax     (:)   ! maximum fresh snowfall density [kg/m3]
     real(kind=r8), allocatable :: SnowMassFullCoverOld   (:)   ! new snow mass to fully cover old snow [mm]
     real(kind=r8), allocatable :: SoilMatPotentialWilt   (:)   ! soil metric potential for wilting point [m]
     real(kind=r8), allocatable :: SnowMeltFac            (:)   ! snowmelt m parameter in snow cover fraction calculation
     real(kind=r8), allocatable :: SnowCoverFac           (:)   ! snow cover factor [m] (originally hard-coded 2.5*z0 in SCF formulation)
     real(kind=r8), allocatable :: InfilFacVic            (:)   ! VIC model infiltration parameter 
     real(kind=r8), allocatable :: TensionWatDistrInfl    (:)   ! tension water distribution inflection parameter
     real(kind=r8), allocatable :: TensionWatDistrShp     (:)   ! tension water distribution shape parameter
     real(kind=r8), allocatable :: FreeWatDistrShp        (:)   ! free water distribution shape parameter
     real(kind=r8), allocatable :: InfilHeteroDynVic      (:)   ! DVIC heterogeniety parameter for infiltration
     real(kind=r8), allocatable :: InfilCapillaryDynVic   (:)   ! DVIC Mean Capillary Drive (m) for infiltration models
     real(kind=r8), allocatable :: InfilFacDynVic         (:)   ! DVIC model infiltration parameter
     real(kind=r8), allocatable :: TileDrainCoeffSp       (:)   ! drainage coefficient [mm d^-1] for simple scheme
     integer      , allocatable :: TileDrainTubeDepth     (:)   ! depth [m] of drain tube from the soil surface for simple scheme
     real(kind=r8), allocatable :: DrainFacSoilWat        (:)   ! drainage factor for soil moisture
     real(kind=r8), allocatable :: TileDrainCoeff         (:)   ! drainage coefficent [m d^-1] for Hooghoudt scheme
     real(kind=r8), allocatable :: DrainDepthToImperv     (:)   ! actual depth of tile drainage to impermeable layer form surface
     real(kind=r8), allocatable :: LateralWatCondFac      (:)   ! multiplication factor to determine lateral hydraulic conductivity
     real(kind=r8), allocatable :: TileDrainDepth         (:)   ! depth of drain [m] for Hooghoudt scheme
     real(kind=r8), allocatable :: DrainTubeDist          (:)   ! distance between two drain tubes or tiles [m]
     real(kind=r8), allocatable :: DrainTubeRadius        (:)   ! effective radius of drain tubes [m]
     real(kind=r8), allocatable :: DrainWatDepToImperv    (:)   ! depth to impervious layer from drain water level [m]
     integer      , allocatable :: NumSoilLayerRoot       (:)   ! number of soil layers with root present
     real(kind=r8), allocatable :: SoilDrainSlope         (:)   ! slope index for soil drainage
     real(kind=r8), allocatable :: SoilMoistureSat        (:,:) ! saturated value of soil moisture [m3/m3]
     real(kind=r8), allocatable :: SoilMoistureWilt       (:,:) ! wilting point soil moisture [m3/m3]
     real(kind=r8), allocatable :: SoilMoistureFieldCap   (:,:) ! reference soil moisture (field capacity) [m3/m3]
     real(kind=r8), allocatable :: SoilMoistureDry        (:,:) ! dry soil moisture threshold [m3/m3]
     real(kind=r8), allocatable :: SoilWatDiffusivitySat  (:,:) ! saturated soil hydraulic diffusivity [m2/s]
     real(kind=r8), allocatable :: SoilWatConductivitySat (:,:) ! saturated soil hydraulic conductivity [m/s]
     real(kind=r8), allocatable :: SoilExpCoeffB          (:,:) ! soil exponent B parameter
     real(kind=r8), allocatable :: SoilMatPotentialSat    (:,:) ! saturated soil matric potential [m]
     real(kind=r8), allocatable :: SoilInfilMaxCoeff      (:)   ! parameter to calculate maximum soil infiltration rate
     real(kind=r8), allocatable :: SoilImpervFracCoeff    (:)   ! parameter to calculate frozen soil impermeable fraction
     real(kind=r8), allocatable :: SnowIceFracPrev        (:,:) ! ice fraction in snow layers at previous timestep
     ! table variables
     type(table_type)           :: table
  end type coupling_type

  ! data type for coupling level namelist options
  type namelist_type
     ! domain specific options
     character(len=cl) :: scrip_file   ! name of SCRIP grid definition file 
     character(len=cl) :: mosaic_file  ! name of mosaic file
     character(len=cl) :: input_dir    ! input directory for tiled files
     logical           :: IsGlobal     ! global vs. regional, default is global
     ! soil specific options
!     integer                    :: NumSoilLayer      ! number of soil layers
!     real(kind=r8), allocatable :: DepthSoilLayer(:) ! layer-bottom depth from soil surface
     ! input specific options
     character(len=cl) :: InputType          ! input type for static fields and initial condition: sfc, custom
     character(len=cl) :: InputIC            ! input file for initial condition or restart 
     character(len=cl) :: InputSoilType      ! input file for soil type data
     character(len=cl) :: InputVegType       ! input file for vegetation type data
     character(len=cl) :: InputSlopeType     ! input file for slope type 
     character(len=cl) :: InputBottomTemp    ! input file for bottom temperature
     character(len=cl) :: InputVegFracAnnMax ! input file annual max vegetation fraction
     ! output specific options
     character(len=cl) :: OutputMode  ! model output mode: all, low or mid
     integer           :: OutputFreq  ! model output interval in seconds
     ! generic options
     character(len=cl) :: CaseName    ! name of case
     integer           :: debug_level ! debug level
     character(len=cl) :: RestartType ! flag for restart run
     logical           :: IsRestart   ! flag if the run is restart or not
  end type namelist_type

  ! data type for forcing
  type forcing_type
    real(kind=r8), allocatable :: SpecHumidityRefHeight  (:) ! specific humidity (water vapor/moist air) at reference height [kg/kg] 
    real(kind=r8), allocatable :: TemperatureAirRefHeight(:) ! air temperature at reference height [K]
    real(kind=r8), allocatable :: WindEastwardRefHeight  (:) ! wind speed in eastward direction at reference height [m/s]
    real(kind=r8), allocatable :: WindNorthwardRefHeight (:) ! wind speed in northward direction at reference height [m/s]
    real(kind=r8), allocatable :: TemperatureSoilBottom  (:) ! bottom boundary condition for soil temperature [K]
    real(kind=r8), allocatable :: RadSwDownRefHeight     (:) ! downward shortwave radiation at reference height [W/m2]
    real(kind=r8), allocatable :: RadLwDownRefHeight     (:) ! downward longwave radiation at reference height [W/m2]
    real(kind=r8), allocatable :: PressureAirRefHeight   (:) ! air pressure at reference height [Pa]
    real(kind=r8), allocatable :: PressureAirSurface     (:) ! air pressure at the surface-atmosphere interface level [Pa]
    real(kind=r8), allocatable :: PrecipConvRefHeight    (:) ! convective precipitation rate at reference heighta [mm/s]
    real(kind=r8), allocatable :: PrecipNonConvRefHeight (:) ! non-convective precipitation rate at reference height [mm/s]
    real(kind=r8), allocatable :: PrecipShConvRefHeight  (:) ! shallow convective precipitation rate at reference height [mm/s]
    real(kind=r8), allocatable :: PrecipSnowRefHeight    (:) ! snow rate at reference height [mm/s]
    real(kind=r8), allocatable :: PrecipGraupelRefHeight (:) ! graupel rate at reference height [mm/s]
    real(kind=r8), allocatable :: PrecipHailRefHeight    (:) ! hail rate at reference height [mm/s]
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
     procedure, private :: AllocateInitCoupling
     procedure, private :: AllocateInitTable
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

    associate(begl  => this%domain%begl, &
              endl  => this%domain%endl, &
              nsoil => this%noahmp%config%domain%NumSoilLayer, &
              nsnow => this%noahmp%config%domain%NumSnowLayerMax)

    call this%AllocateInitForcing(begl, endl)
    call this%AllocateInitCoupling(begl, endl, nsoil, nsnow)

    end associate

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

  !===========================================================================

  subroutine AllocateInitCoupling(this, begl, endl, nsnow, nsoil)

    class(model_type) :: this
    integer :: begl, endl, nsnow, nsoil

    ! allocate
    if (.not. allocated(this%coupling%CosSolarZenithAngle)) allocate(this%coupling%CosSolarZenithAngle(begl:endl))
    if (.not. allocated(this%coupling%SoilType))            allocate(this%coupling%SoilType(begl:endl))
    if (.not. allocated(this%coupling%VegType))             allocate(this%coupling%VegType(begl:endl))
    if (.not. allocated(this%coupling%RunoffSlopeType))     allocate(this%coupling%RunoffSlopeType(begl:endl))
    if (.not. allocated(this%coupling%VegFracAnnMax))       allocate(this%coupling%VegFracAnnMax(begl:endl))
    if (.not. allocated(this%coupling%SoilMoisture))        allocate(this%coupling%SoilMoisture(begl:endl,nsoil))
    if (.not. allocated(this%coupling%TemperatureSoilSnow)) allocate(this%coupling%TemperatureSoilSnow(begl:endl,-nsnow+1:nsoil))
    if (.not. allocated(this%coupling%SnowDepth))           allocate(this%coupling%SnowDepth(begl:endl))
    if (.not. allocated(this%coupling%SnowWaterEquiv))      allocate(this%coupling%SnowWaterEquiv(begl:endl))
    if (.not. allocated(this%coupling%NumSnowLayerNeg))     allocate(this%coupling%NumSnowLayerNeg(begl:endl))

    ! init
    this%coupling%CosSolarZenithAngle(:)   = undefined_real
    this%coupling%SoilType(:)              = undefined_int
    this%coupling%VegType(:)               = undefined_int
    this%coupling%RunoffSlopeType(:)       = undefined_int
    this%coupling%VegFracAnnMax(:)         = undefined_real
    this%coupling%SoilMoisture(:,:)        = undefined_real
    this%coupling%TemperatureSoilSnow(:,:) = undefined_real
    this%coupling%SnowDepth(:)             = undefined_real
    this%coupling%SnowWaterEquiv(:)        = undefined_real
    this%coupling%NumSnowLayerNeg(:)       = undefined_int

    ! water state variables
    if (.not. allocated(this%coupling%CanopyLiqWater         )) allocate(this%coupling%CanopyLiqWater         (begl:endl)) 
    if (.not. allocated(this%coupling%CanopyIce              )) allocate(this%coupling%CanopyIce              (begl:endl))
    if (.not. allocated(this%coupling%CanopyWetFrac          )) allocate(this%coupling%CanopyWetFrac          (begl:endl))
    if (.not. allocated(this%coupling%SnowWaterEquiv         )) allocate(this%coupling%SnowWaterEquiv         (begl:endl))
    if (.not. allocated(this%coupling%SnowWaterEquivPrev     )) allocate(this%coupling%SnowWaterEquivPrev     (begl:endl))
    if (.not. allocated(this%coupling%SnowDepth              )) allocate(this%coupling%SnowDepth              (begl:endl))
    if (.not. allocated(this%coupling%IrrigationFracFlood    )) allocate(this%coupling%IrrigationFracFlood    (begl:endl))
    if (.not. allocated(this%coupling%IrrigationAmtFlood     )) allocate(this%coupling%IrrigationAmtFlood     (begl:endl))
    if (.not. allocated(this%coupling%IrrigationFracMicro    )) allocate(this%coupling%IrrigationFracMicro    (begl:endl))
    if (.not. allocated(this%coupling%IrrigationAmtMicro     )) allocate(this%coupling%IrrigationAmtMicro     (begl:endl))
    if (.not. allocated(this%coupling%IrrigationFracSprinkler)) allocate(this%coupling%IrrigationFracSprinkler(begl:endl))
    if (.not. allocated(this%coupling%IrrigationAmtSprinkler )) allocate(this%coupling%IrrigationAmtSprinkler (begl:endl))
    if (.not. allocated(this%coupling%WaterTableDepth        )) allocate(this%coupling%WaterTableDepth        (begl:endl))
    if (.not. allocated(this%coupling%SoilMoistureToWT       )) allocate(this%coupling%SoilMoistureToWT       (begl:endl))
    if (.not. allocated(this%coupling%TileDrainFrac          )) allocate(this%coupling%TileDrainFrac          (begl:endl))
    if (.not. allocated(this%coupling%WaterStorageAquifer    )) allocate(this%coupling%WaterStorageAquifer    (begl:endl))
    if (.not. allocated(this%coupling%WaterStorageSoilAqf    )) allocate(this%coupling%WaterStorageSoilAqf    (begl:endl))
    if (.not. allocated(this%coupling%WaterStorageLake       )) allocate(this%coupling%WaterStorageLake       (begl:endl))
    if (.not. allocated(this%coupling%IrrigationFracGrid     )) allocate(this%coupling%IrrigationFracGrid     (begl:endl))
    if (.not. allocated(this%coupling%IrrigationCntSprinkler )) allocate(this%coupling%IrrigationCntSprinkler (begl:endl))
    if (.not. allocated(this%coupling%IrrigationCntMicro     )) allocate(this%coupling%IrrigationCntMicro     (begl:endl))
    if (.not. allocated(this%coupling%IrrigationCntFlood     )) allocate(this%coupling%IrrigationCntFlood     (begl:endl))
    if (.not. allocated(this%coupling%SnowIce                )) allocate(this%coupling%SnowIce                (begl:endl,-nsnow+1:0))
    if (.not. allocated(this%coupling%SnowLiqWater           )) allocate(this%coupling%SnowLiqWater           (begl:endl,-nsnow+1:0))
    if (.not. allocated(this%coupling%SoilLiqWater           )) allocate(this%coupling%SoilLiqWater           (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilMoisture           )) allocate(this%coupling%SoilMoisture           (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilMoistureEqui       )) allocate(this%coupling%SoilMoistureEqui       (begl:endl,nsoil))
    if (.not. allocated(this%coupling%RechargeGwDeepWT       )) allocate(this%coupling%RechargeGwDeepWT       (begl:endl))
    if (.not. allocated(this%coupling%RechargeGwShallowWT    )) allocate(this%coupling%RechargeGwShallowWT    (begl:endl))
    if (.not. allocated(this%coupling%EvapSoilSfcLiqAcc      )) allocate(this%coupling%EvapSoilSfcLiqAcc      (begl:endl))
    if (.not. allocated(this%coupling%SoilSfcInflowAcc       )) allocate(this%coupling%SoilSfcInflowAcc       (begl:endl))
    if (.not. allocated(this%coupling%SfcWaterTotChgAcc      )) allocate(this%coupling%SfcWaterTotChgAcc      (begl:endl))
    if (.not. allocated(this%coupling%PrecipTotAcc           )) allocate(this%coupling%PrecipTotAcc           (begl:endl))
    if (.not. allocated(this%coupling%EvapCanopyNetAcc       )) allocate(this%coupling%EvapCanopyNetAcc       (begl:endl))
    if (.not. allocated(this%coupling%TranspirationAcc       )) allocate(this%coupling%TranspirationAcc       (begl:endl))
    if (.not. allocated(this%coupling%EvapGroundNetAcc       )) allocate(this%coupling%EvapGroundNetAcc       (begl:endl))
    if (.not. allocated(this%coupling%TranspWatLossSoilAcc   )) allocate(this%coupling%TranspWatLossSoilAcc   (begl:endl,nsoil))
    ! water parameter variables
    if (.not. allocated(this%coupling%DrainSoilLayerInd      )) allocate(this%coupling%DrainSoilLayerInd      (begl:endl))
    if (.not. allocated(this%coupling%CanopyLiqHoldCap       )) allocate(this%coupling%CanopyLiqHoldCap       (begl:endl))
    if (.not. allocated(this%coupling%SnowCompactBurdenFac   )) allocate(this%coupling%SnowCompactBurdenFac   (begl:endl))
    if (.not. allocated(this%coupling%SnowCompactAgingFac1   )) allocate(this%coupling%SnowCompactAgingFac1   (begl:endl))
    if (.not. allocated(this%coupling%SnowCompactAgingFac2   )) allocate(this%coupling%SnowCompactAgingFac2   (begl:endl))
    if (.not. allocated(this%coupling%SnowCompactAgingFac3   )) allocate(this%coupling%SnowCompactAgingFac3   (begl:endl))
    if (.not. allocated(this%coupling%SnowCompactAgingMax    )) allocate(this%coupling%SnowCompactAgingMax    (begl:endl))
    if (.not. allocated(this%coupling%SnowViscosityCoeff     )) allocate(this%coupling%SnowViscosityCoeff     (begl:endl))
    if (.not. allocated(this%coupling%SnowLiqFracMax         )) allocate(this%coupling%SnowLiqFracMax         (begl:endl))
    if (.not. allocated(this%coupling%SnowLiqHoldCap         )) allocate(this%coupling%SnowLiqHoldCap         (begl:endl))
    if (.not. allocated(this%coupling%SnowLiqReleaseFac      )) allocate(this%coupling%SnowLiqReleaseFac      (begl:endl))
    if (.not. allocated(this%coupling%IrriFloodRateFac       )) allocate(this%coupling%IrriFloodRateFac       (begl:endl))
    if (.not. allocated(this%coupling%IrriMicroRate          )) allocate(this%coupling%IrriMicroRate          (begl:endl))
    if (.not. allocated(this%coupling%SoilConductivityRef    )) allocate(this%coupling%SoilConductivityRef    (begl:endl))
    if (.not. allocated(this%coupling%SoilInfilFacRef        )) allocate(this%coupling%SoilInfilFacRef        (begl:endl))
    if (.not. allocated(this%coupling%GroundFrzCoeff         )) allocate(this%coupling%GroundFrzCoeff         (begl:endl))
    if (.not. allocated(this%coupling%GridTopoIndex          )) allocate(this%coupling%GridTopoIndex          (begl:endl))
    if (.not. allocated(this%coupling%SoilSfcSatFracMax      )) allocate(this%coupling%SoilSfcSatFracMax      (begl:endl))
    if (.not. allocated(this%coupling%SpecYieldGw            )) allocate(this%coupling%SpecYieldGw            (begl:endl))
    if (.not. allocated(this%coupling%MicroPoreContent       )) allocate(this%coupling%MicroPoreContent       (begl:endl))
    if (.not. allocated(this%coupling%WaterStorageLakeMax    )) allocate(this%coupling%WaterStorageLakeMax     (begl:endl))
    if (.not. allocated(this%coupling%SnoWatEqvMaxGlacier    )) allocate(this%coupling%SnoWatEqvMaxGlacier    (begl:endl))
    if (.not. allocated(this%coupling%IrriStopDayBfHarvest   )) allocate(this%coupling%IrriStopDayBfHarvest   (begl:endl))
    if (.not. allocated(this%coupling%IrriTriggerLaiMin      )) allocate(this%coupling%IrriTriggerLaiMin      (begl:endl))  
    if (.not. allocated(this%coupling%SoilWatDeficitAllow    )) allocate(this%coupling%SoilWatDeficitAllow    (begl:endl))
    if (.not. allocated(this%coupling%IrriFloodLossFrac      )) allocate(this%coupling%IrriFloodLossFrac      (begl:endl))
    if (.not. allocated(this%coupling%IrriSprinklerRate      )) allocate(this%coupling%IrriSprinklerRate      (begl:endl))
    if (.not. allocated(this%coupling%IrriFracThreshold      )) allocate(this%coupling%IrriFracThreshold      (begl:endl))
    if (.not. allocated(this%coupling%IrriStopPrecipThr      )) allocate(this%coupling%IrriStopPrecipThr      (begl:endl))
    if (.not. allocated(this%coupling%SnowfallDensityMax     )) allocate(this%coupling%SnowfallDensityMax     (begl:endl))
    if (.not. allocated(this%coupling%SnowMassFullCoverOld   )) allocate(this%coupling%SnowMassFullCoverOld   (begl:endl))
    if (.not. allocated(this%coupling%SoilMatPotentialWilt   )) allocate(this%coupling%SoilMatPotentialWilt   (begl:endl))
    if (.not. allocated(this%coupling%SnowMeltFac            )) allocate(this%coupling%SnowMeltFac            (begl:endl))
    if (.not. allocated(this%coupling%SnowCoverFac           )) allocate(this%coupling%SnowCoverFac           (begl:endl))
    if (.not. allocated(this%coupling%InfilFacVic            )) allocate(this%coupling%InfilFacVic            (begl:endl))
    if (.not. allocated(this%coupling%TensionWatDistrInfl    )) allocate(this%coupling%TensionWatDistrInfl    (begl:endl))
    if (.not. allocated(this%coupling%TensionWatDistrShp     )) allocate(this%coupling%TensionWatDistrShp     (begl:endl))
    if (.not. allocated(this%coupling%FreeWatDistrShp        )) allocate(this%coupling%FreeWatDistrShp        (begl:endl))
    if (.not. allocated(this%coupling%InfilHeteroDynVic      )) allocate(this%coupling%InfilHeteroDynVic      (begl:endl))
    if (.not. allocated(this%coupling%InfilCapillaryDynVic   )) allocate(this%coupling%InfilCapillaryDynVic   (begl:endl))
    if (.not. allocated(this%coupling%InfilFacDynVic         )) allocate(this%coupling%InfilFacDynVic         (begl:endl))
    if (.not. allocated(this%coupling%TileDrainCoeffSp       )) allocate(this%coupling%TileDrainCoeffSp       (begl:endl))
    if (.not. allocated(this%coupling%TileDrainTubeDepth     )) allocate(this%coupling%TileDrainTubeDepth     (begl:endl))
    if (.not. allocated(this%coupling%DrainFacSoilWat        )) allocate(this%coupling%DrainFacSoilWat        (begl:endl)) 
    if (.not. allocated(this%coupling%TileDrainCoeff         )) allocate(this%coupling%TileDrainCoeff         (begl:endl))
    if (.not. allocated(this%coupling%DrainDepthToImperv     )) allocate(this%coupling%DrainDepthToImperv     (begl:endl))
    if (.not. allocated(this%coupling%LateralWatCondFac      )) allocate(this%coupling%LateralWatCondFac      (begl:endl))
    if (.not. allocated(this%coupling%TileDrainDepth         )) allocate(this%coupling%TileDrainDepth         (begl:endl))
    if (.not. allocated(this%coupling%DrainTubeDist          )) allocate(this%coupling%DrainTubeDist          (begl:endl))
    if (.not. allocated(this%coupling%DrainTubeRadius        )) allocate(this%coupling%DrainTubeRadius        (begl:endl))
    if (.not. allocated(this%coupling%DrainWatDepToImperv    )) allocate(this%coupling%DrainWatDepToImperv    (begl:endl))
    if (.not. allocated(this%coupling%NumSoilLayerRoot       )) allocate(this%coupling%NumSoilLayerRoot       (begl:endl))
    if (.not. allocated(this%coupling%SoilDrainSlope         )) allocate(this%coupling%SoilDrainSlope         (begl:endl))
    if (.not. allocated(this%coupling%SoilMoistureSat        )) allocate(this%coupling%SoilMoistureSat        (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilMoistureWilt       )) allocate(this%coupling%SoilMoistureWilt       (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilMoistureFieldCap   )) allocate(this%coupling%SoilMoistureFieldCap   (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilMoistureDry        )) allocate(this%coupling%SoilMoistureDry        (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilWatDiffusivitySat  )) allocate(this%coupling%SoilWatDiffusivitySat  (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilWatConductivitySat )) allocate(this%coupling%SoilWatConductivitySat (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilExpCoeffB          )) allocate(this%coupling%SoilExpCoeffB          (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilMatPotentialSat    )) allocate(this%coupling%SoilMatPotentialSat    (begl:endl,nsoil))
    if (.not. allocated(this%coupling%SoilInfilMaxCoeff      )) allocate(this%coupling%SoilInfilMaxCoeff      (begl:endl))
    if (.not. allocated(this%coupling%SoilImpervFracCoeff    )) allocate(this%coupling%SoilImpervFracCoeff    (begl:endl))
    if (.not. allocated(this%coupling%SnowIceFracPrev        )) allocate(this%coupling%SnowIceFracPrev        (begl:endl,-nsnow+1:0))

    ! init water variables
    this%coupling%CanopyLiqWater          = undefined_real
    this%coupling%CanopyIce               = undefined_real
    this%coupling%CanopyWetFrac           = undefined_real
    this%coupling%SnowWaterEquiv          = undefined_real
    this%coupling%SnowWaterEquivPrev      = undefined_real
    this%coupling%SnowDepth               = undefined_real
    this%coupling%IrrigationFracFlood     = undefined_real
    this%coupling%IrrigationAmtFlood      = undefined_real
    this%coupling%IrrigationFracMicro     = undefined_real
    this%coupling%IrrigationAmtMicro      = undefined_real
    this%coupling%IrrigationFracSprinkler = undefined_real
    this%coupling%IrrigationAmtSprinkler  = undefined_real
    this%coupling%WaterTableDepth         = undefined_real
    this%coupling%SoilMoistureToWT        = undefined_real
    this%coupling%TileDrainFrac           = undefined_real
    this%coupling%WaterStorageAquifer     = undefined_real
    this%coupling%WaterStorageSoilAqf     = undefined_real
    this%coupling%WaterStorageLake        = undefined_real
    this%coupling%IrrigationFracGrid      = undefined_real
    this%coupling%IrrigationCntSprinkler  = undefined_int
    this%coupling%IrrigationCntMicro      = undefined_int
    this%coupling%IrrigationCntFlood      = undefined_int
    this%coupling%SnowIce                 = undefined_real 
    this%coupling%SnowLiqWater            = undefined_real
    this%coupling%SoilLiqWater            = undefined_real
    this%coupling%SoilMoisture            = undefined_real
    this%coupling%SoilMoistureEqui        = undefined_real
    this%coupling%RechargeGwDeepWT        = undefined_real 
    this%coupling%RechargeGwShallowWT     = undefined_real
    this%coupling%EvapSoilSfcLiqAcc       = undefined_real
    this%coupling%SoilSfcInflowAcc        = undefined_real
    this%coupling%SfcWaterTotChgAcc       = undefined_real
    this%coupling%PrecipTotAcc            = undefined_real
    this%coupling%EvapCanopyNetAcc        = undefined_real
    this%coupling%TranspirationAcc        = undefined_real
    this%coupling%EvapGroundNetAcc        = undefined_real
    this%coupling%TranspWatLossSoilAcc    = undefined_real
    this%coupling%DrainSoilLayerInd       = undefined_int
    this%coupling%CanopyLiqHoldCap        = undefined_real
    this%coupling%SnowCompactBurdenFac    = undefined_real
    this%coupling%SnowCompactAgingFac1    = undefined_real
    this%coupling%SnowCompactAgingFac2    = undefined_real
    this%coupling%SnowCompactAgingFac3    = undefined_real
    this%coupling%SnowCompactAgingMax     = undefined_real
    this%coupling%SnowViscosityCoeff      = undefined_real
    this%coupling%SnowLiqFracMax          = undefined_real
    this%coupling%SnowLiqHoldCap          = undefined_real
    this%coupling%SnowLiqReleaseFac       = undefined_real
    this%coupling%IrriFloodRateFac        = undefined_real
    this%coupling%IrriMicroRate           = undefined_real
    this%coupling%SoilConductivityRef     = undefined_real
    this%coupling%SoilInfilFacRef         = undefined_real
    this%coupling%GroundFrzCoeff          = undefined_real
    this%coupling%GridTopoIndex           = undefined_real
    this%coupling%SoilSfcSatFracMax       = undefined_real
    this%coupling%SpecYieldGw             = undefined_real
    this%coupling%MicroPoreContent        = undefined_real
    this%coupling%WaterStorageLakeMax     = undefined_real
    this%coupling%SnoWatEqvMaxGlacier     = undefined_real
    this%coupling%IrriStopDayBfHarvest    = undefined_int
    this%coupling%IrriTriggerLaiMin       = undefined_real 
    this%coupling%SoilWatDeficitAllow     = undefined_real
    this%coupling%IrriFloodLossFrac       = undefined_real
    this%coupling%IrriSprinklerRate       = undefined_real
    this%coupling%IrriFracThreshold       = undefined_real
    this%coupling%IrriStopPrecipThr       = undefined_real
    this%coupling%SnowfallDensityMax      = undefined_real
    this%coupling%SnowMassFullCoverOld    = undefined_real
    this%coupling%SoilMatPotentialWilt    = undefined_real
    this%coupling%SnowMeltFac             = undefined_real
    this%coupling%SnowCoverFac            = undefined_real
    this%coupling%InfilFacVic             = undefined_real
    this%coupling%TensionWatDistrInfl     = undefined_real
    this%coupling%TensionWatDistrShp      = undefined_real
    this%coupling%FreeWatDistrShp         = undefined_real
    this%coupling%InfilHeteroDynVic       = undefined_real
    this%coupling%InfilCapillaryDynVic    = undefined_real
    this%coupling%InfilFacDynVic          = undefined_real
    this%coupling%TileDrainCoeffSp        = undefined_real
    this%coupling%TileDrainTubeDepth      = undefined_int
    this%coupling%DrainFacSoilWat         = undefined_real 
    this%coupling%TileDrainCoeff          = undefined_real
    this%coupling%DrainDepthToImperv      = undefined_real
    this%coupling%LateralWatCondFac       = undefined_real
    this%coupling%TileDrainDepth          = undefined_real
    this%coupling%DrainTubeDist           = undefined_real
    this%coupling%DrainTubeRadius         = undefined_real
    this%coupling%DrainWatDepToImperv     = undefined_real
    this%coupling%NumSoilLayerRoot        = undefined_int
    this%coupling%SoilDrainSlope          = undefined_real

    ! init table variables
    call AllocateInitTable(this)

  end subroutine AllocateInitCoupling

  subroutine AllocateInitTable(this)

    class(model_type) :: this

    ! vegetation parameters
    if (.not. allocated(this%coupling%table%CH2OP_TABLE)  ) allocate(this%coupling%table%CH2OP_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%DLEAF_TABLE)  ) allocate(this%coupling%table%DLEAF_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%Z0MVT_TABLE)  ) allocate(this%coupling%table%Z0MVT_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%HVT_TABLE)    ) allocate(this%coupling%table%HVT_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%HVB_TABLE)    ) allocate(this%coupling%table%HVB_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%DEN_TABLE)    ) allocate(this%coupling%table%DEN_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%RC_TABLE)     ) allocate(this%coupling%table%RC_TABLE     (MVT))
    if (.not. allocated(this%coupling%table%MFSNO_TABLE)  ) allocate(this%coupling%table%MFSNO_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%SCFFAC_TABLE) ) allocate(this%coupling%table%SCFFAC_TABLE (MVT))
    if (.not. allocated(this%coupling%table%CBIOM_TABLE)  ) allocate(this%coupling%table%CBIOM_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%SAIM_TABLE)   ) allocate(this%coupling%table%SAIM_TABLE   (MVT,12))
    if (.not. allocated(this%coupling%table%LAIM_TABLE)   ) allocate(this%coupling%table%LAIM_TABLE   (MVT,12))
    if (.not. allocated(this%coupling%table%SLA_TABLE)    ) allocate(this%coupling%table%SLA_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%DILEFC_TABLE) ) allocate(this%coupling%table%DILEFC_TABLE (MVT))
    if (.not. allocated(this%coupling%table%DILEFW_TABLE) ) allocate(this%coupling%table%DILEFW_TABLE (MVT))
    if (.not. allocated(this%coupling%table%FRAGR_TABLE)  ) allocate(this%coupling%table%FRAGR_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%LTOVRC_TABLE) ) allocate(this%coupling%table%LTOVRC_TABLE (MVT))
    if (.not. allocated(this%coupling%table%C3PSN_TABLE)  ) allocate(this%coupling%table%C3PSN_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%KC25_TABLE)   ) allocate(this%coupling%table%KC25_TABLE   (MVT))
    if (.not. allocated(this%coupling%table%AKC_TABLE)    ) allocate(this%coupling%table%AKC_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%KO25_TABLE)   ) allocate(this%coupling%table%KO25_TABLE   (MVT))
    if (.not. allocated(this%coupling%table%AKO_TABLE)    ) allocate(this%coupling%table%AKO_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%VCMX25_TABLE) ) allocate(this%coupling%table%VCMX25_TABLE (MVT))
    if (.not. allocated(this%coupling%table%AVCMX_TABLE)  ) allocate(this%coupling%table%AVCMX_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%BP_TABLE)     ) allocate(this%coupling%table%BP_TABLE     (MVT))
    if (.not. allocated(this%coupling%table%MP_TABLE)     ) allocate(this%coupling%table%MP_TABLE     (MVT))
    if (.not. allocated(this%coupling%table%QE25_TABLE)   ) allocate(this%coupling%table%QE25_TABLE   (MVT))
    if (.not. allocated(this%coupling%table%AQE_TABLE)    ) allocate(this%coupling%table%AQE_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%RMF25_TABLE)  ) allocate(this%coupling%table%RMF25_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%RMS25_TABLE)  ) allocate(this%coupling%table%RMS25_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%RMR25_TABLE)  ) allocate(this%coupling%table%RMR25_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%ARM_TABLE)    ) allocate(this%coupling%table%ARM_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%FOLNMX_TABLE) ) allocate(this%coupling%table%FOLNMX_TABLE (MVT))
    if (.not. allocated(this%coupling%table%TMIN_TABLE)   ) allocate(this%coupling%table%TMIN_TABLE   (MVT))
    if (.not. allocated(this%coupling%table%XL_TABLE)     ) allocate(this%coupling%table%XL_TABLE     (MVT))
    if (.not. allocated(this%coupling%table%RHOL_TABLE)   ) allocate(this%coupling%table%RHOL_TABLE   (MVT,MBAND))
    if (.not. allocated(this%coupling%table%RHOS_TABLE)   ) allocate(this%coupling%table%RHOS_TABLE   (MVT,MBAND))
    if (.not. allocated(this%coupling%table%TAUL_TABLE)   ) allocate(this%coupling%table%TAUL_TABLE   (MVT,MBAND))
    if (.not. allocated(this%coupling%table%TAUS_TABLE)   ) allocate(this%coupling%table%TAUS_TABLE   (MVT,MBAND))
    if (.not. allocated(this%coupling%table%MRP_TABLE)    ) allocate(this%coupling%table%MRP_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%CWPVT_TABLE)  ) allocate(this%coupling%table%CWPVT_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%WRRAT_TABLE)  ) allocate(this%coupling%table%WRRAT_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%WDPOOL_TABLE) ) allocate(this%coupling%table%WDPOOL_TABLE (MVT))
    if (.not. allocated(this%coupling%table%TDLEF_TABLE)  ) allocate(this%coupling%table%TDLEF_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%NROOT_TABLE)  ) allocate(this%coupling%table%NROOT_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%RGL_TABLE)    ) allocate(this%coupling%table%RGL_TABLE    (MVT))
    if (.not. allocated(this%coupling%table%RS_TABLE)     ) allocate(this%coupling%table%RS_TABLE     (MVT))
    if (.not. allocated(this%coupling%table%HS_TABLE)     ) allocate(this%coupling%table%HS_TABLE     (MVT))
    if (.not. allocated(this%coupling%table%TOPT_TABLE)   ) allocate(this%coupling%table%TOPT_TABLE   (MVT))
    if (.not. allocated(this%coupling%table%RSMAX_TABLE)  ) allocate(this%coupling%table%RSMAX_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%RTOVRC_TABLE) ) allocate(this%coupling%table%RTOVRC_TABLE (MVT))
    if (.not. allocated(this%coupling%table%RSWOODC_TABLE)) allocate(this%coupling%table%RSWOODC_TABLE(MVT))
    if (.not. allocated(this%coupling%table%BF_TABLE)     ) allocate(this%coupling%table%BF_TABLE     (MVT))
    if (.not. allocated(this%coupling%table%WSTRC_TABLE)  ) allocate(this%coupling%table%WSTRC_TABLE  (MVT))
    if (.not. allocated(this%coupling%table%LAIMIN_TABLE) ) allocate(this%coupling%table%LAIMIN_TABLE (MVT))
    if (.not. allocated(this%coupling%table%XSAMIN_TABLE) ) allocate(this%coupling%table%XSAMIN_TABLE (MVT))

    ! soil parameters
    if (.not. allocated(this%coupling%table%BEXP_TABLE)   ) allocate(this%coupling%table%BEXP_TABLE  (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%SMCDRY_TABLE) ) allocate(this%coupling%table%SMCDRY_TABLE(MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%SMCMAX_TABLE) ) allocate(this%coupling%table%SMCMAX_TABLE(MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%SMCREF_TABLE) ) allocate(this%coupling%table%SMCREF_TABLE(MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%PSISAT_TABLE) ) allocate(this%coupling%table%PSISAT_TABLE(MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%DKSAT_TABLE)  ) allocate(this%coupling%table%DKSAT_TABLE (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%DWSAT_TABLE)  ) allocate(this%coupling%table%DWSAT_TABLE (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%SMCWLT_TABLE) ) allocate(this%coupling%table%SMCWLT_TABLE(MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%QUARTZ_TABLE) ) allocate(this%coupling%table%QUARTZ_TABLE(MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%BVIC_TABLE)   ) allocate(this%coupling%table%BVIC_TABLE  (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%AXAJ_TABLE)   ) allocate(this%coupling%table%AXAJ_TABLE  (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%BXAJ_TABLE)   ) allocate(this%coupling%table%BXAJ_TABLE  (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%XXAJ_TABLE)   ) allocate(this%coupling%table%XXAJ_TABLE  (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%BDVIC_TABLE)  ) allocate(this%coupling%table%BDVIC_TABLE (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%GDVIC_TABLE)  ) allocate(this%coupling%table%GDVIC_TABLE (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%BBVIC_TABLE)  ) allocate(this%coupling%table%BBVIC_TABLE (MAX_SOILTYP))

    ! general parameters
    if (.not. allocated(this%coupling%table%SLOPE_TABLE)  ) allocate(this%coupling%table%SLOPE_TABLE (NUM_SLOPE))

    ! radiation parameters
    if (.not. allocated(this%coupling%table%ALBSAT_TABLE) ) allocate(this%coupling%table%ALBSAT_TABLE(MSC,MBAND))
    if (.not. allocated(this%coupling%table%ALBDRY_TABLE) ) allocate(this%coupling%table%ALBDRY_TABLE(MSC,MBAND))
    if (.not. allocated(this%coupling%table%ALBICE_TABLE) ) allocate(this%coupling%table%ALBICE_TABLE(MBAND)    )
    if (.not. allocated(this%coupling%table%ALBLAK_TABLE) ) allocate(this%coupling%table%ALBLAK_TABLE(MBAND)    )
    if (.not. allocated(this%coupling%table%OMEGAS_TABLE) ) allocate(this%coupling%table%OMEGAS_TABLE(MBAND)    )
    if (.not. allocated(this%coupling%table%EG_TABLE)     ) allocate(this%coupling%table%EG_TABLE(2)            )

    ! tile drainage parameters
    if (.not. allocated(this%coupling%table%TDSMC_FAC_TABLE)) allocate(this%coupling%table%TDSMC_FAC_TABLE(MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%TD_DC_TABLE)    ) allocate(this%coupling%table%TD_DC_TABLE    (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%TD_DEPTH_TABLE) ) allocate(this%coupling%table%TD_DEPTH_TABLE (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%TD_DCOEF_TABLE) ) allocate(this%coupling%table%TD_DCOEF_TABLE (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%TD_D_TABLE)     ) allocate(this%coupling%table%TD_D_TABLE     (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%TD_ADEPTH_TABLE)) allocate(this%coupling%table%TD_ADEPTH_TABLE(MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%TD_RADI_TABLE)  ) allocate(this%coupling%table%TD_RADI_TABLE  (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%TD_SPAC_TABLE)  ) allocate(this%coupling%table%TD_SPAC_TABLE  (MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%TD_DDRAIN_TABLE)) allocate(this%coupling%table%TD_DDRAIN_TABLE(MAX_SOILTYP))
    if (.not. allocated(this%coupling%table%KLAT_FAC_TABLE) ) allocate(this%coupling%table%KLAT_FAC_TABLE (MAX_SOILTYP))

    ! crop parameters
    if (.not. allocated(this%coupling%table%PLTDAY_TABLE)   ) allocate(this%coupling%table%PLTDAY_TABLE   (NCROP))
    if (.not. allocated(this%coupling%table%HSDAY_TABLE)    ) allocate(this%coupling%table%HSDAY_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%PLANTPOP_TABLE) ) allocate(this%coupling%table%PLANTPOP_TABLE (NCROP))
    if (.not. allocated(this%coupling%table%IRRI_TABLE)     ) allocate(this%coupling%table%IRRI_TABLE     (NCROP))
    if (.not. allocated(this%coupling%table%GDDTBASE_TABLE) ) allocate(this%coupling%table%GDDTBASE_TABLE (NCROP))
    if (.not. allocated(this%coupling%table%GDDTCUT_TABLE)  ) allocate(this%coupling%table%GDDTCUT_TABLE  (NCROP))
    if (.not. allocated(this%coupling%table%GDDS1_TABLE)    ) allocate(this%coupling%table%GDDS1_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%GDDS2_TABLE)    ) allocate(this%coupling%table%GDDS2_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%GDDS3_TABLE)    ) allocate(this%coupling%table%GDDS3_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%GDDS4_TABLE)    ) allocate(this%coupling%table%GDDS4_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%GDDS5_TABLE)    ) allocate(this%coupling%table%GDDS5_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%C3PSNI_TABLE)   ) allocate(this%coupling%table%C3PSNI_TABLE   (NCROP))
    if (.not. allocated(this%coupling%table%KC25I_TABLE)    ) allocate(this%coupling%table%KC25I_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%AKCI_TABLE)     ) allocate(this%coupling%table%AKCI_TABLE     (NCROP))
    if (.not. allocated(this%coupling%table%KO25I_TABLE)    ) allocate(this%coupling%table%KO25I_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%AKOI_TABLE)     ) allocate(this%coupling%table%AKOI_TABLE     (NCROP))
    if (.not. allocated(this%coupling%table%VCMX25I_TABLE)  ) allocate(this%coupling%table%VCMX25I_TABLE  (NCROP))
    if (.not. allocated(this%coupling%table%AVCMXI_TABLE)   ) allocate(this%coupling%table%AVCMXI_TABLE   (NCROP))
    if (.not. allocated(this%coupling%table%BPI_TABLE)      ) allocate(this%coupling%table%BPI_TABLE      (NCROP))
    if (.not. allocated(this%coupling%table%MPI_TABLE)      ) allocate(this%coupling%table%MPI_TABLE      (NCROP))
    if (.not. allocated(this%coupling%table%QE25I_TABLE)    ) allocate(this%coupling%table%QE25I_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%FOLNMXI_TABLE)  ) allocate(this%coupling%table%FOLNMXI_TABLE  (NCROP))
    if (.not. allocated(this%coupling%table%AREF_TABLE)     ) allocate(this%coupling%table%AREF_TABLE     (NCROP))
    if (.not. allocated(this%coupling%table%PSNRF_TABLE)    ) allocate(this%coupling%table%PSNRF_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%I2PAR_TABLE)    ) allocate(this%coupling%table%I2PAR_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%TASSIM0_TABLE)  ) allocate(this%coupling%table%TASSIM0_TABLE  (NCROP))
    if (.not. allocated(this%coupling%table%TASSIM1_TABLE)  ) allocate(this%coupling%table%TASSIM1_TABLE  (NCROP))
    if (.not. allocated(this%coupling%table%TASSIM2_TABLE)  ) allocate(this%coupling%table%TASSIM2_TABLE  (NCROP))
    if (.not. allocated(this%coupling%table%K_TABLE)        ) allocate(this%coupling%table%K_TABLE        (NCROP))
    if (.not. allocated(this%coupling%table%EPSI_TABLE)     ) allocate(this%coupling%table%EPSI_TABLE     (NCROP))
    if (.not. allocated(this%coupling%table%Q10MR_TABLE)    ) allocate(this%coupling%table%Q10MR_TABLE    (NCROP))
    if (.not. allocated(this%coupling%table%LEFREEZ_TABLE)  ) allocate(this%coupling%table%LEFREEZ_TABLE  (NCROP))
    if (.not. allocated(this%coupling%table%DILE_FC_TABLE)  ) allocate(this%coupling%table%DILE_FC_TABLE  (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%DILE_FW_TABLE)  ) allocate(this%coupling%table%DILE_FW_TABLE  (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%FRA_GR_TABLE)   ) allocate(this%coupling%table%FRA_GR_TABLE   (NCROP))
    if (.not. allocated(this%coupling%table%LF_OVRC_TABLE)  ) allocate(this%coupling%table%LF_OVRC_TABLE  (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%ST_OVRC_TABLE)  ) allocate(this%coupling%table%ST_OVRC_TABLE  (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%RT_OVRC_TABLE)  ) allocate(this%coupling%table%RT_OVRC_TABLE  (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%LFMR25_TABLE)   ) allocate(this%coupling%table%LFMR25_TABLE   (NCROP))
    if (.not. allocated(this%coupling%table%STMR25_TABLE)   ) allocate(this%coupling%table%STMR25_TABLE   (NCROP))
    if (.not. allocated(this%coupling%table%RTMR25_TABLE)   ) allocate(this%coupling%table%RTMR25_TABLE   (NCROP))
    if (.not. allocated(this%coupling%table%GRAINMR25_TABLE)) allocate(this%coupling%table%GRAINMR25_TABLE(NCROP))
    if (.not. allocated(this%coupling%table%LFPT_TABLE)     ) allocate(this%coupling%table%LFPT_TABLE     (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%STPT_TABLE)     ) allocate(this%coupling%table%STPT_TABLE     (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%RTPT_TABLE)     ) allocate(this%coupling%table%RTPT_TABLE     (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%GRAINPT_TABLE)  ) allocate(this%coupling%table%GRAINPT_TABLE  (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%LFCT_TABLE)     ) allocate(this%coupling%table%LFCT_TABLE     (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%STCT_TABLE)     ) allocate(this%coupling%table%STCT_TABLE     (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%RTCT_TABLE)     ) allocate(this%coupling%table%RTCT_TABLE     (NCROP,NSTAGE))
    if (.not. allocated(this%coupling%table%BIO2LAI_TABLE)  ) allocate(this%coupling%table%BIO2LAI_TABLE  (NCROP))

    ! vegetation parameters
    this%coupling%table%ISURBAN_TABLE      = undefined_int
    this%coupling%table%ISWATER_TABLE      = undefined_int
    this%coupling%table%ISBARREN_TABLE     = undefined_int
    this%coupling%table%ISICE_TABLE        = undefined_int
    this%coupling%table%ISCROP_TABLE       = undefined_int
    this%coupling%table%EBLFOREST_TABLE    = undefined_int
    this%coupling%table%NATURAL_TABLE      = undefined_int
    this%coupling%table%URBTYPE_BEG        = undefined_int
    this%coupling%table%LCZ_1_TABLE        = undefined_int
    this%coupling%table%LCZ_2_TABLE        = undefined_int
    this%coupling%table%LCZ_3_TABLE        = undefined_int
    this%coupling%table%LCZ_4_TABLE        = undefined_int
    this%coupling%table%LCZ_5_TABLE        = undefined_int
    this%coupling%table%LCZ_6_TABLE        = undefined_int
    this%coupling%table%LCZ_7_TABLE        = undefined_int
    this%coupling%table%LCZ_8_TABLE        = undefined_int
    this%coupling%table%LCZ_9_TABLE        = undefined_int
    this%coupling%table%LCZ_10_TABLE       = undefined_int
    this%coupling%table%LCZ_11_TABLE       = undefined_int
    this%coupling%table%CH2OP_TABLE        = undefined_real
    this%coupling%table%DLEAF_TABLE        = undefined_real
    this%coupling%table%Z0MVT_TABLE        = undefined_real
    this%coupling%table%HVT_TABLE          = undefined_real
    this%coupling%table%HVB_TABLE          = undefined_real
    this%coupling%table%DEN_TABLE          = undefined_real
    this%coupling%table%RC_TABLE           = undefined_real
    this%coupling%table%MFSNO_TABLE        = undefined_real
    this%coupling%table%SCFFAC_TABLE       = undefined_real
    this%coupling%table%CBIOM_TABLE        = undefined_real
    this%coupling%table%RHOL_TABLE         = undefined_real
    this%coupling%table%RHOS_TABLE         = undefined_real
    this%coupling%table%TAUL_TABLE         = undefined_real
    this%coupling%table%TAUS_TABLE         = undefined_real
    this%coupling%table%XL_TABLE           = undefined_real
    this%coupling%table%CWPVT_TABLE        = undefined_real
    this%coupling%table%C3PSN_TABLE        = undefined_real
    this%coupling%table%KC25_TABLE         = undefined_real
    this%coupling%table%AKC_TABLE          = undefined_real
    this%coupling%table%KO25_TABLE         = undefined_real
    this%coupling%table%AKO_TABLE          = undefined_real
    this%coupling%table%AVCMX_TABLE        = undefined_real
    this%coupling%table%AQE_TABLE          = undefined_real
    this%coupling%table%LTOVRC_TABLE       = undefined_real
    this%coupling%table%DILEFC_TABLE       = undefined_real
    this%coupling%table%DILEFW_TABLE       = undefined_real
    this%coupling%table%RMF25_TABLE        = undefined_real
    this%coupling%table%SLA_TABLE          = undefined_real
    this%coupling%table%FRAGR_TABLE        = undefined_real
    this%coupling%table%TMIN_TABLE         = undefined_real
    this%coupling%table%VCMX25_TABLE       = undefined_real
    this%coupling%table%TDLEF_TABLE        = undefined_real
    this%coupling%table%BP_TABLE           = undefined_real
    this%coupling%table%MP_TABLE           = undefined_real
    this%coupling%table%QE25_TABLE         = undefined_real
    this%coupling%table%RMS25_TABLE        = undefined_real
    this%coupling%table%RMR25_TABLE        = undefined_real
    this%coupling%table%ARM_TABLE          = undefined_real
    this%coupling%table%FOLNMX_TABLE       = undefined_real
    this%coupling%table%WDPOOL_TABLE       = undefined_real
    this%coupling%table%WRRAT_TABLE        = undefined_real
    this%coupling%table%MRP_TABLE          = undefined_real
    this%coupling%table%SAIM_TABLE         = undefined_real
    this%coupling%table%LAIM_TABLE         = undefined_real
    this%coupling%table%NROOT_TABLE        = undefined_real
    this%coupling%table%RGL_TABLE          = undefined_real
    this%coupling%table%RS_TABLE           = undefined_real
    this%coupling%table%HS_TABLE           = undefined_real
    this%coupling%table%TOPT_TABLE         = undefined_real
    this%coupling%table%RSMAX_TABLE        = undefined_real
    this%coupling%table%RTOVRC_TABLE       = undefined_real
    this%coupling%table%RSWOODC_TABLE      = undefined_real
    this%coupling%table%BF_TABLE           = undefined_real
    this%coupling%table%WSTRC_TABLE        = undefined_real
    this%coupling%table%LAIMIN_TABLE       = undefined_real
    this%coupling%table%XSAMIN_TABLE       = undefined_real

    ! soil parameters
    this%coupling%table%SLCATS_TABLE       = undefined_int
    this%coupling%table%BEXP_TABLE         = undefined_real
    this%coupling%table%SMCDRY_TABLE       = undefined_real
    this%coupling%table%SMCMAX_TABLE       = undefined_real
    this%coupling%table%SMCREF_TABLE       = undefined_real
    this%coupling%table%PSISAT_TABLE       = undefined_real
    this%coupling%table%DKSAT_TABLE        = undefined_real
    this%coupling%table%DWSAT_TABLE        = undefined_real
    this%coupling%table%SMCWLT_TABLE       = undefined_real
    this%coupling%table%QUARTZ_TABLE       = undefined_real
    this%coupling%table%BVIC_TABLE         = undefined_real
    this%coupling%table%AXAJ_TABLE         = undefined_real
    this%coupling%table%BXAJ_TABLE         = undefined_real
    this%coupling%table%XXAJ_TABLE         = undefined_real
    this%coupling%table%BDVIC_TABLE        = undefined_real
    this%coupling%table%GDVIC_TABLE        = undefined_real
    this%coupling%table%BBVIC_TABLE        = undefined_real

    ! general parameters
    this%coupling%table%SLOPE_TABLE        = undefined_real
    this%coupling%table%CSOIL_TABLE        = undefined_real
    this%coupling%table%REFDK_TABLE        = undefined_real
    this%coupling%table%REFKDT_TABLE       = undefined_real
    this%coupling%table%FRZK_TABLE         = undefined_real
    this%coupling%table%ZBOT_TABLE         = undefined_real
    this%coupling%table%CZIL_TABLE         = undefined_real

    ! radiation parameters
    this%coupling%table%ALBSAT_TABLE       = undefined_real
    this%coupling%table%ALBDRY_TABLE       = undefined_real
    this%coupling%table%ALBICE_TABLE       = undefined_real
    this%coupling%table%ALBLAK_TABLE       = undefined_real
    this%coupling%table%OMEGAS_TABLE       = undefined_real
    this%coupling%table%BETADS_TABLE       = undefined_real
    this%coupling%table%BETAIS_TABLE       = undefined_real
    this%coupling%table%EG_TABLE           = undefined_real
    this%coupling%table%EICE_TABLE         = undefined_real

    ! global parameters
    this%coupling%table%CO2_TABLE              = undefined_real
    this%coupling%table%O2_TABLE               = undefined_real
    this%coupling%table%TIMEAN_TABLE           = undefined_real
    this%coupling%table%FSATMX_TABLE           = undefined_real
    this%coupling%table%Z0SNO_TABLE            = undefined_real
    this%coupling%table%SSI_TABLE              = undefined_real
    this%coupling%table%SNOW_RET_FAC_TABLE     = undefined_real
    this%coupling%table%SNOW_EMIS_TABLE        = undefined_real
    this%coupling%table%SWEMX_TABLE            = undefined_real
    this%coupling%table%TAU0_TABLE             = undefined_real
    this%coupling%table%GRAIN_GROWTH_TABLE     = undefined_real
    this%coupling%table%EXTRA_GROWTH_TABLE     = undefined_real
    this%coupling%table%DIRT_SOOT_TABLE        = undefined_real
    this%coupling%table%BATS_COSZ_TABLE        = undefined_real
    this%coupling%table%BATS_VIS_NEW_TABLE     = undefined_real
    this%coupling%table%BATS_NIR_NEW_TABLE     = undefined_real
    this%coupling%table%BATS_VIS_AGE_TABLE     = undefined_real
    this%coupling%table%BATS_NIR_AGE_TABLE     = undefined_real
    this%coupling%table%BATS_VIS_DIR_TABLE     = undefined_real
    this%coupling%table%BATS_NIR_DIR_TABLE     = undefined_real
    this%coupling%table%RSURF_SNOW_TABLE       = undefined_real
    this%coupling%table%RSURF_EXP_TABLE        = undefined_real
    this%coupling%table%C2_SNOWCOMPACT_TABLE   = undefined_real
    this%coupling%table%C3_SNOWCOMPACT_TABLE   = undefined_real
    this%coupling%table%C4_SNOWCOMPACT_TABLE   = undefined_real
    this%coupling%table%C5_SNOWCOMPACT_TABLE   = undefined_real
    this%coupling%table%DM_SNOWCOMPACT_TABLE   = undefined_real
    this%coupling%table%ETA0_SNOWCOMPACT_TABLE = undefined_real
    this%coupling%table%SNLIQMAXFRAC_TABLE     = undefined_real
    this%coupling%table%SWEMAXGLA_TABLE        = undefined_real
    this%coupling%table%WSLMAX_TABLE           = undefined_real
    this%coupling%table%ROUS_TABLE             = undefined_real
    this%coupling%table%CMIC_TABLE             = undefined_real
    this%coupling%table%SNOWDEN_MAX_TABLE      = undefined_real
    this%coupling%table%CLASS_ALB_REF_TABLE    = undefined_real
    this%coupling%table%CLASS_SNO_AGE_TABLE    = undefined_real
    this%coupling%table%CLASS_ALB_NEW_TABLE    = undefined_real
    this%coupling%table%PSIWLT_TABLE           = undefined_real
    this%coupling%table%Z0SOIL_TABLE           = undefined_real
    this%coupling%table%Z0LAKE_TABLE           = undefined_real

    ! irrigation parameters
    this%coupling%table%IRR_HAR_TABLE          = undefined_int
    this%coupling%table%IRR_FRAC_TABLE         = undefined_real
    this%coupling%table%IRR_LAI_TABLE          = undefined_real
    this%coupling%table%IRR_MAD_TABLE          = undefined_real
    this%coupling%table%FILOSS_TABLE           = undefined_real
    this%coupling%table%SPRIR_RATE_TABLE       = undefined_real
    this%coupling%table%MICIR_RATE_TABLE       = undefined_real
    this%coupling%table%FIRTFAC_TABLE          = undefined_real
    this%coupling%table%IR_RAIN_TABLE          = undefined_real

    ! crop parameters
    this%coupling%table%DEFAULT_CROP_TABLE     = undefined_int
    this%coupling%table%PLTDAY_TABLE           = undefined_int
    this%coupling%table%HSDAY_TABLE            = undefined_int
    this%coupling%table%PLANTPOP_TABLE         = undefined_real
    this%coupling%table%IRRI_TABLE             = undefined_real
    this%coupling%table%GDDTBASE_TABLE         = undefined_real
    this%coupling%table%GDDTCUT_TABLE          = undefined_real
    this%coupling%table%GDDS1_TABLE            = undefined_real
    this%coupling%table%GDDS2_TABLE            = undefined_real
    this%coupling%table%GDDS3_TABLE            = undefined_real
    this%coupling%table%GDDS4_TABLE            = undefined_real
    this%coupling%table%GDDS5_TABLE            = undefined_real
    this%coupling%table%C3PSNI_TABLE           = undefined_real
    this%coupling%table%KC25I_TABLE            = undefined_real
    this%coupling%table%AKCI_TABLE             = undefined_real
    this%coupling%table%KO25I_TABLE            = undefined_real
    this%coupling%table%AKOI_TABLE             = undefined_real
    this%coupling%table%AVCMXI_TABLE           = undefined_real
    this%coupling%table%VCMX25I_TABLE          = undefined_real
    this%coupling%table%BPI_TABLE              = undefined_real
    this%coupling%table%MPI_TABLE              = undefined_real
    this%coupling%table%FOLNMXI_TABLE          = undefined_real
    this%coupling%table%QE25I_TABLE            = undefined_real
    this%coupling%table%AREF_TABLE             = undefined_real
    this%coupling%table%PSNRF_TABLE            = undefined_real
    this%coupling%table%I2PAR_TABLE            = undefined_real
    this%coupling%table%TASSIM0_TABLE          = undefined_real
    this%coupling%table%TASSIM1_TABLE          = undefined_real
    this%coupling%table%TASSIM2_TABLE          = undefined_real
    this%coupling%table%K_TABLE                = undefined_real
    this%coupling%table%EPSI_TABLE             = undefined_real
    this%coupling%table%Q10MR_TABLE            = undefined_real
    this%coupling%table%LEFREEZ_TABLE          = undefined_real
    this%coupling%table%DILE_FC_TABLE          = undefined_real
    this%coupling%table%DILE_FW_TABLE          = undefined_real
    this%coupling%table%FRA_GR_TABLE           = undefined_real
    this%coupling%table%LF_OVRC_TABLE          = undefined_real
    this%coupling%table%ST_OVRC_TABLE          = undefined_real
    this%coupling%table%RT_OVRC_TABLE          = undefined_real
    this%coupling%table%LFMR25_TABLE           = undefined_real
    this%coupling%table%STMR25_TABLE           = undefined_real
    this%coupling%table%RTMR25_TABLE           = undefined_real
    this%coupling%table%GRAINMR25_TABLE        = undefined_real
    this%coupling%table%LFPT_TABLE             = undefined_real
    this%coupling%table%STPT_TABLE             = undefined_real
    this%coupling%table%RTPT_TABLE             = undefined_real
    this%coupling%table%GRAINPT_TABLE          = undefined_real
    this%coupling%table%LFCT_TABLE             = undefined_real
    this%coupling%table%STCT_TABLE             = undefined_real
    this%coupling%table%RTCT_TABLE             = undefined_real
    this%coupling%table%BIO2LAI_TABLE          = undefined_real

    ! tile drainage parameters
    this%coupling%table%DRAIN_LAYER_OPT_TABLE  = undefined_int
    this%coupling%table%TD_DEPTH_TABLE         = undefined_int
    this%coupling%table%TDSMC_FAC_TABLE        = undefined_real 
    this%coupling%table%TD_DC_TABLE            = undefined_real
    this%coupling%table%TD_DCOEF_TABLE         = undefined_real
    this%coupling%table%TD_D_TABLE             = undefined_real
    this%coupling%table%TD_ADEPTH_TABLE        = undefined_real
    this%coupling%table%TD_RADI_TABLE          = undefined_real
    this%coupling%table%TD_SPAC_TABLE          = undefined_real
    this%coupling%table%TD_DDRAIN_TABLE        = undefined_real
    this%coupling%table%KLAT_FAC_TABLE         = undefined_real

    ! optional parameters
    this%coupling%table%SR2006_THETA_1500T_A_TABLE = undefined_real
    this%coupling%table%SR2006_THETA_1500T_B_TABLE = undefined_real
    this%coupling%table%SR2006_THETA_1500T_C_TABLE = undefined_real
    this%coupling%table%SR2006_THETA_1500T_D_TABLE = undefined_real
    this%coupling%table%SR2006_THETA_1500T_E_TABLE = undefined_real
    this%coupling%table%SR2006_THETA_1500T_F_TABLE = undefined_real
    this%coupling%table%SR2006_THETA_1500T_G_TABLE = undefined_real
    this%coupling%table%SR2006_THETA_1500_A_TABLE  = undefined_real
    this%coupling%table%SR2006_THETA_1500_B_TABLE  = undefined_real
    this%coupling%table%SR2006_THETA_33T_A_TABLE   = undefined_real
    this%coupling%table%SR2006_THETA_33T_B_TABLE   = undefined_real
    this%coupling%table%SR2006_THETA_33T_C_TABLE   = undefined_real
    this%coupling%table%SR2006_THETA_33T_D_TABLE   = undefined_real
    this%coupling%table%SR2006_THETA_33T_E_TABLE   = undefined_real
    this%coupling%table%SR2006_THETA_33T_F_TABLE   = undefined_real
    this%coupling%table%SR2006_THETA_33T_G_TABLE   = undefined_real
    this%coupling%table%SR2006_THETA_33_A_TABLE    = undefined_real
    this%coupling%table%SR2006_THETA_33_B_TABLE    = undefined_real
    this%coupling%table%SR2006_THETA_33_C_TABLE    = undefined_real
    this%coupling%table%SR2006_THETA_S33T_A_TABLE  = undefined_real
    this%coupling%table%SR2006_THETA_S33T_B_TABLE  = undefined_real
    this%coupling%table%SR2006_THETA_S33T_C_TABLE  = undefined_real
    this%coupling%table%SR2006_THETA_S33T_D_TABLE  = undefined_real
    this%coupling%table%SR2006_THETA_S33T_E_TABLE  = undefined_real
    this%coupling%table%SR2006_THETA_S33T_F_TABLE  = undefined_real
    this%coupling%table%SR2006_THETA_S33T_G_TABLE  = undefined_real
    this%coupling%table%SR2006_THETA_S33_A_TABLE   = undefined_real
    this%coupling%table%SR2006_THETA_S33_B_TABLE   = undefined_real
    this%coupling%table%SR2006_PSI_ET_A_TABLE      = undefined_real
    this%coupling%table%SR2006_PSI_ET_B_TABLE      = undefined_real
    this%coupling%table%SR2006_PSI_ET_C_TABLE      = undefined_real
    this%coupling%table%SR2006_PSI_ET_D_TABLE      = undefined_real
    this%coupling%table%SR2006_PSI_ET_E_TABLE      = undefined_real
    this%coupling%table%SR2006_PSI_ET_F_TABLE      = undefined_real
    this%coupling%table%SR2006_PSI_ET_G_TABLE      = undefined_real
    this%coupling%table%SR2006_PSI_E_A_TABLE       = undefined_real
    this%coupling%table%SR2006_PSI_E_B_TABLE       = undefined_real
    this%coupling%table%SR2006_PSI_E_C_TABLE       = undefined_real
    this%coupling%table%SR2006_SMCMAX_A_TABLE      = undefined_real
    this%coupling%table%SR2006_SMCMAX_B_TABLE      = undefined_real

  end subroutine AllocateInitTable

end module lnd_comp_types

