module lnd_comp_driver

  ! This file contains the NoahMP land surface model driver

  use ESMF, only: operator(+), operator(-), operator(/)
  use ESMF, only: ESMF_GridComp, ESMF_GridCompGet
  use ESMF, only: ESMF_VM, ESMF_VMGet, ESMF_UtilIOUnitGet
  use ESMF, only: ESMF_Clock, ESMF_Time, ESMF_TimeSet
  use ESMF, only: ESMF_TimeInterval, ESMF_TimeGet, ESMF_ClockGet
  use ESMF, only: ESMF_TimeIntervalGet, ESMF_KIND_R8
  use ESMF, only: ESMF_TimeIsLeapYear, ESMF_LOGMSG_INFO
  use ESMF, only: ESMF_LOGMSG_ERROR, ESMF_FAILURE
  use ESMF, only: ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO

  use lnd_comp_types, only: model_type
  use lnd_comp_shr  , only: ChkErr
  use lnd_comp_kind , only: cl => shr_kind_cl
  use lnd_comp_kind , only: r8 => shr_kind_r8
  use lnd_comp_types, only: histflds, restflds

  use lnd_comp_io, only: SetupWriteFields
  use lnd_comp_io, only: WriteFile
  use lnd_comp_io, only: ReadStatic
  use lnd_comp_io, only: ReadIC

  use ConfigVarInitMod , only: ConfigVarInitDefault
  use ForcingVarInitMod, only: ForcingVarInitDefault
  use EnergyVarInitMod , only: EnergyVarInitDefault
  use WaterVarInitMod  , only: WaterVarInitDefault
  use BiochemVarInitMod, only: BiochemVarInitDefault
  use NoahmpMainGlacierMod, only: NoahmpMainGlacier
  use NoahmpMainMod, only: NoahmpMain

  implicit none
  private

  public :: NoahmpDriverInit
  public :: NoahmpDriverMain
  public :: NoahmpDriverFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(ESMF_Time)        :: epoc
  type(ESMF_Time)        :: dummTime, prevTime, nextTime
  character(*),parameter :: modName =  "(lnd_comp_driver)"

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=============================================================================
contains
!=============================================================================

  subroutine NoahmpDriverInit(gcomp, model, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    type(model_type)   , intent(inout) :: model
    integer            , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)           :: vm
    integer                 :: localPet
    character(len=*), parameter :: subname=trim(modName)//':(NoahmpDriverInit) '
    !-------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ----------------------
    ! Query component
    ! ----------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Set epoc as 1970-01-01 00:00:00
    !----------------------

    model%io%reference_date = "1970-01-01 00:00:00"
    call ESMF_TimeSet(epoc, yy=1970, mm=1, dd=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------
    ! Allocate and initialize model data
    ! ----------------------

    call model%AllocateInit()

    ! ----------------------
    ! Read in Noah-MP Parameters from MPTABLE
    ! ----------------------

    call NoahmpReadTable(model, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------
    ! Read static information
    ! ----------------------

    call ReadStatic(model, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Read initial condition / restart file
    ! ---------------------

    if (model%nmlist%IsRestart) then
    !   call ReadRestart()
    else
       call ReadIC(model, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ---------------------
    ! Initialize initial state
    ! ---------------------

    call InitializeStates(model, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------
    ! Setup data structures to write model output
    ! ----------------------

    ! history
    call SetupWriteFields(model, 'hist', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! restart
    call SetupWriteFields(model, 'rest', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine NoahmpDriverInit

  !===============================================================================

  subroutine NoahmpDriverMain(gcomp, model, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    type(model_type)   , intent(inout) :: model
    integer            , intent(out)   :: rc

    ! local variables
    integer                 :: i
    integer                 :: year, month, day
    integer                 :: hour, minute, second
    logical, save           :: first_time = .true.
    real(ESMF_KIND_R8)      :: now_time, dofy
    real(ESMF_KIND_R8)      :: tol = 1.0d-20
    character(len=cl)       :: filename
    type(ESMF_VM)           :: vm
    type(ESMF_Clock)        :: clock
    type(ESMF_Time)         :: startTime
    type(ESMF_Time)         :: currTime
    type(ESMF_TimeInterval) :: timeStep
    logical, save           :: firstTime = .true.
    character(len=cl) :: msg
    character(len=*), parameter :: subname=trim(modName)//':(NoahmpDriverMain) '
    !-------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! Query component clock and extract required information for model run
    !----------------------

    call ESMF_GridCompGet(gcomp, vm=vm, clock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, &
      h=hour, m=minute, s=second, dayOfYear_r8=dofy, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------
    ! Set model specific parameters
    ! ----------------------

    ! get coupling time step
    call ESMF_TimeIntervalGet(timeStep, s_r8=model%coupling%MainTimeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    model%noahmp%config%domain%MainTimeStep = model%coupling%MainTimeStep

    ! check model time step
    if (model%noahmp%config%domain%SoilTimeStep < 0.0) then
       ! use coupling time step and internal time step of model
       model%noahmp%config%domain%SoilTimeStep = model%noahmp%config%domain%MainTimeStep
       model%noahmp%config%domain%NumSoilTimeStep = 1
    else
       ! check specified time step divides evenly with coupling time step
       if (abs(mod(model%coupling%MainTimeStep, model%noahmp%config%domain%SoilTimeStep)) > tol) then
          call ESMF_LogWrite(trim(subname)//': model time step does not divide evenly with coupling time step. Exiting ... ', ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if

       ! calculate number of soil time step based on coupling and soil time steps
       model%noahmp%config%domain%NumSoilTimeStep = int(model%noahmp%config%domain%MainTimeStep/model%noahmp%config%domain%SoilTimeStep)
       write(msg, fmt='(I5)') model%noahmp%config%domain%NumSoilTimeStep
       call ESMF_LogWrite(trim(subname)//': NumSoilTimeStep = '//trim(msg), ESMF_LOGMSG_INFO)
    end if

    ! number of days in the particular year
    model%noahmp%config%domain%NumDayInYear = 365
    if (ESMF_TimeIsLeapYear(currTime, rc=rc)) then
       model%noahmp%config%domain%NumDayInYear = 366
    end if

    ! Julian day of year
    model%noahmp%config%domain%DayJulianInYear = dofy

    ! cosine solar zenith angle (0-1)
    call CalcCosZ(currTime, model%domain%latm, model%domain%lonm, model%coupling%CosSolarZenithAngle, dofy, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! main NoahMP model time step
    model%noahmp%config%domain%MainTimeStep = model%coupling%MainTimeStep

    ! set surface type (1 - soil and 2 - lake)
    model%noahmp%config%domain%SurfaceType = 1

    ! number of shortwave radiation bands
    model%noahmp%config%domain%NumSwRadBand = 2

    ! soil color type for albedo
    ! TODO: This is read from file in UFS
    model%noahmp%config%domain%SoilColor = 4

    ! number of crop growth stages
    model%noahmp%config%domain%NumCropGrowStage = 8

    ! model grid spacing size
    ! TODO: maybe model could be calculated from grid area?
    model%noahmp%config%domain%GridSize = -9999.0

    ! thickness of lowest atmosphere layer
    ! TODO: namelist option? or retrieved from active atmosphere
    model%noahmp%config%domain%ThicknessAtmosBotLayer = -9999.0

    ! ----------------------
    ! Loop over grid points and run NoahMP model
    ! ----------------------

    do i = model%domain%begl, model%domain%endl
       ! ----------------------
       ! Skip any open water points
       ! ----------------------

       if (model%domain%mask(i) < 1) cycle

       ! ----------------------
       ! Set config related variables 
       ! NOTE: soil type same in all layers
       ! ----------------------

       model%noahmp%config%domain%GridIndexI            = i
       model%noahmp%config%domain%GridIndexJ            = -9999
       model%noahmp%config%domain%Latitude              = model%domain%latm(i)
       model%noahmp%config%domain%CosSolarZenithAngle   = model%coupling%CosSolarZenithAngle(i)
       model%noahmp%config%domain%SoilType(:)           = model%coupling%SoilType(i)
       model%noahmp%config%domain%VegType               = model%coupling%VegType(i)
       model%noahmp%config%domain%RunoffSlopeType       = model%coupling%RunoffSlopeType(i)
       model%noahmp%config%domain%NumSnowLayerNeg       = model%coupling%NumSnowLayerNeg(i)

       ! ----------------------
       ! Transfer all the inputs from couplign layer to the model 
       ! ----------------------

       call ConfigVarInTransfer (model)
       call ForcingVarInTransfer(model)
       call EnergyVarInTransfer (model)
       call WaterVarInTransfer  (model)
       call BiochemVarInTransfer(model)

       ! ----------------------
       ! Set energy related variables 
       ! ----------------------

       !model%noahmp%energy%param%VegFracAnnMax          = model%coupling%VegFracAnnMax(i)
       !model%noahmp%energy%state%TemperatureSoilSnow(:) = model%coupling%TemperatureSoilSnow(i,:)

       ! ----------------------
       ! Set water related variables 
       ! ----------------------
       
       !model%noahmp%water%state%SoilMoisture(:)         = model%coupling%SoilMoisture(i,:)
       !model%noahmp%water%state%SnowWaterEquiv          = model%coupling%SnowWaterEquiv(i)
       !model%noahmp%water%state%SnowDepth               = model%coupling%SnowDepth(i)

    !   ! ----------------------
    !   ! Skip any open water points
    !   ! ----------------------

    !   if (model%domain%mask(i) < 1) cycle

    !   ! ----------------------
    !   ! Transfer model variables to 1-D Noah-MP column variables
    !   ! ----------------------

    !   !call ConfigVarInTransfer(model)
    !   call ForcingVarInTransfer(model, i)

    !   ! ----------------------
    !   ! Call 1D Noah-MP LSM
    !   ! ----------------------

    !   if (model%noahmp%config%domain%VegType == model%noahmp%config%domain%IndexIcePoint ) then
    !      ! land ice point
    !      model%noahmp%config%domain%IndicatorIceSfc = -1

    !      ! set deep glaicer temp to >= -10C
    !      model%noahmp%forcing%TemperatureSoilBottom = min(model%noahmp%forcing%TemperatureSoilBottom, 263.15)

    !      ! call glacier land
    !      call NoahmpMainGlacier(model%noahmp)
    !   else
    !      ! land soil point
    !      model%noahmp%config%domain%IndicatorIceSfc = 0

    !      ! call non-glacier land
    !      call NoahmpMain(model%noahmp)
    !   endif

    !   ! ----------------------
    !   ! Transfer 1-D Noah-MP column variables to model variables
    !   ! ----------------------

    end do

    ! ----------------------
    ! Write model output and restart files
    ! ----------------------

    ! Query date to create file name
    call ESMF_TimeGet(currTime+timeStep, yy=year, mm=month, dd=day, &
      h=hour, m=minute, s=second, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get date as second
    call ESMF_TimeIntervalGet(currTime-epoc+timeStep, s_r8=now_time, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Check the output frequency before calling write method
    if (mod(int(now_time), model%nmlist%OutputFreq) == 0) then
       write(filename, fmt='(a,i4,a1,i2.2,a1,i2.2,a1,i5.5)') &
          trim(model%nmlist%CaseName)//'.lnd.out.', &
          year, '-', month, '-', day, '-', hour*60*60+minute*60+second
       call WriteFile(model, filename, histflds, now_time, vm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine NoahmpDriverMain

  !===============================================================================

  subroutine NoahmpDriverFinalize(gcomp, model, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)   :: gcomp
    type(model_type)  , intent(inout) :: model
    integer            , intent(out)  :: rc

    ! local variables
    integer            :: localPet
    type(ESMF_VM)      :: vm
    character(len=*), parameter :: subname=trim(modName)//':(NoahmpDriverFinalize) '
    !-------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine NoahmpDriverFinalize

  !===============================================================================

  subroutine InitializeStates(model, rc)

    ! input/output variables
    type(model_type), intent(inout) :: model
    integer         , intent(out)   :: rc

    ! local variables
    integer :: i
    real(ESMF_KIND_R8) :: depthm
    character(len=*), parameter :: subname=trim(modName)//':(InitializeStates) '
    !-------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! loop over all grid points
    do i = model%domain%begl, model%domain%endl

       ! snow depth in meters
       depthm = model%coupling%SnowDepth(i)/1000.0_r8

       ! set number of snow layer based on depth
       !dzsno = 0.0
       if (depthm < 0.025 ) then
          model%coupling%NumSnowLayerNeg(i) = 0
          !dzsno(-2:0) = 0.0
       elseif (depthm >= 0.025 .and. depthm <= 0.05 ) then
          model%coupling%NumSnowLayerNeg(i) = -1
          !dzsno(0) = depthm
       elseif (depthm > 0.05 .and. depthm <= 0.10 ) then
          model%coupling%NumSnowLayerNeg(i) = -2
          !dzsno(-1) = 0.5*depthm
          !dzsno(0) = 0.5*depthm
       elseif (depthm > 0.10 .and. depthm <= 0.25 ) then
          model%coupling%NumSnowLayerNeg(i) = -2
          !dzsno(-1) = 0.05
          !dzsno(0) = depthm - 0.05
       elseif (depthm > 0.25 .and. depthm <= 0.45 ) then
          model%coupling%NumSnowLayerNeg(i) = -3
          !dzsno(-2) = 0.05
          !dzsno(-1) = 0.5*(depthm-0.05)
          !dzsno(0) = 0.5*(depthm-0.05)
       elseif (depthm > 0.45) then 
          model%coupling%NumSnowLayerNeg(i) = -3
          !dzsno(-2) = 0.05
          !dzsno(-1) = 0.20
          !dzsno(0) = depthm - 0.05 - 0.20
       endif

    end do

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeStates

  !===============================================================================

  subroutine NoahmpReadTable(model, rcode)

    ! Adapted from HRLDAS driver

    use lnd_comp_types, only: MVT, MBAND, MSC, MAX_SOILTYP, NCROP, NSTAGE, NUM_SLOPE

    ! input/output variables
    type(model_type), target, intent(inout) :: model
    integer         , intent(out)   :: rcode

    ! local variables
    ! vegetation parameters
    character(len=256)                  :: VEG_DATASET_DESCRIPTION
    integer                             :: IERR, IK, IM
    integer                             :: NVEG, ISURBAN, ISWATER, ISBARREN, ISICE, ISCROP, EBLFOREST, NATURAL, URBTYPE_beg
    integer                             :: LCZ_1, LCZ_2, LCZ_3, LCZ_4, LCZ_5, LCZ_6, LCZ_7, LCZ_8, LCZ_9, LCZ_10, LCZ_11    
    real(kind=r8), dimension(MVT)       :: SAI_JAN, SAI_FEB, SAI_MAR, SAI_APR, SAI_MAY, SAI_JUN, SAI_JUL, SAI_AUG,        &
                                           SAI_SEP, SAI_OCT, SAI_NOV, SAI_DEC, LAI_JAN, LAI_FEB, LAI_MAR, LAI_APR,        &
                                           LAI_MAY, LAI_JUN, LAI_JUL, LAI_AUG, LAI_SEP, LAI_OCT, LAI_NOV, LAI_DEC,        &
                                           RHOL_VIS, RHOL_NIR, RHOS_VIS, RHOS_NIR, TAUL_VIS, TAUL_NIR, TAUS_VIS, TAUS_NIR,&
                                           CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, MFSNO, SCFFAC, XL, CWPVT, C3PSN, KC25, &
                                           AKC, KO25, AKO, AVCMX, AQE, LTOVRC, DILEFC, DILEFW, RMF25, SLA, FRAGR, TMIN,   &
                                           VCMX25, TDLEF, BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,    &
                                           NROOT, RGL, RS, HS, TOPT, RSMAX, RTOVRC, RSWOODC, BF, WSTRC, LAIMIN, CBIOM,    &
                                           XSAMIN    
    namelist / noahmp_usgs_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    namelist / noahmp_usgs_parameters     / ISURBAN, ISWATER, ISBARREN, ISICE, ISCROP, EBLFOREST, NATURAL, URBTYPE_beg,   &
                                            LCZ_1, LCZ_2, LCZ_3, LCZ_4, LCZ_5, LCZ_6, LCZ_7, LCZ_8, LCZ_9, LCZ_10, LCZ_11,&
                                            CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, MFSNO, SCFFAC, XL, CWPVT, C3PSN, KC25,&
                                            AKC, KO25, AKO, AVCMX, AQE, LTOVRC, DILEFC, DILEFW, RMF25, SLA, FRAGR, TMIN,  &
                                            VCMX25, TDLEF, BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
                                            NROOT, RGL, RS, HS, TOPT, RSMAX, RTOVRC, RSWOODC, BF, WSTRC, LAIMIN, CBIOM,   &
                                            XSAMIN, SAI_JAN, SAI_FEB, SAI_MAR, SAI_APR, SAI_MAY,                          &
                                            SAI_JUN, SAI_JUL, SAI_AUG, SAI_SEP, SAI_OCT, SAI_NOV, SAI_DEC, LAI_JAN,       &
                                            LAI_FEB, LAI_MAR, LAI_APR, LAI_MAY, LAI_JUN, LAI_JUL, LAI_AUG, LAI_SEP,       &
                                            LAI_OCT, LAI_NOV, LAI_DEC, RHOL_VIS, RHOL_NIR, RHOS_VIS, RHOS_NIR, TAUL_VIS,  &
                                            TAUL_NIR, TAUS_VIS, TAUS_NIR
    namelist / noahmp_modis_veg_categories /VEG_DATASET_DESCRIPTION, NVEG
    namelist / noahmp_modis_parameters     /ISURBAN, ISWATER, ISBARREN, ISICE, ISCROP, EBLFOREST, NATURAL, URBTYPE_beg,   &
                                            LCZ_1, LCZ_2, LCZ_3, LCZ_4, LCZ_5, LCZ_6, LCZ_7, LCZ_8, LCZ_9, LCZ_10, LCZ_11,&
                                            CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, MFSNO, SCFFAC, XL, CWPVT, C3PSN, KC25,&
                                            AKC, KO25, AKO, AVCMX, AQE, LTOVRC, DILEFC, DILEFW, RMF25, SLA, FRAGR, TMIN,  &
                                            VCMX25, TDLEF, BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
                                            NROOT, RGL, RS, HS, TOPT, RSMAX, RTOVRC, RSWOODC, BF, WSTRC, LAIMIN, CBIOM,   &
                                            XSAMIN, SAI_JAN, SAI_FEB, SAI_MAR, SAI_APR, SAI_MAY,                          &
                                            SAI_JUN, SAI_JUL, SAI_AUG, SAI_SEP, SAI_OCT, SAI_NOV, SAI_DEC, LAI_JAN,       &
                                            LAI_FEB, LAI_MAR, LAI_APR, LAI_MAY, LAI_JUN, LAI_JUL, LAI_AUG, LAI_SEP,       &
                                            LAI_OCT, LAI_NOV, LAI_DEC, RHOL_VIS, RHOL_NIR, RHOS_VIS, RHOS_NIR, TAUL_VIS,  &
                                            TAUL_NIR, TAUS_VIS, TAUS_NIR
    ! soil parameters
    character(len=10)                    :: SLTYPE
    integer                              :: SLCATS
    real(kind=r8), dimension(MAX_SOILTYP):: BB, DRYSMC, MAXSMC, REFSMC, SATPSI, SATDK, SATDW, WLTSMC, QTZ,    &
                                            BVIC, AXAJ, BXAJ, XXAJ, BDVIC, BBVIC, GDVIC, HC
    namelist / noahmp_stas_soil_categories /SLTYPE, SLCATS
    namelist / noahmp_soil_stas_parameters /BB, DRYSMC, MAXSMC, REFSMC, SATPSI, SATDK, SATDW, WLTSMC, QTZ,    &
                                            BVIC, AXAJ, BXAJ, XXAJ, BDVIC, BBVIC, GDVIC
    namelist / noahmp_soil_stas_ruc_parameters / BB, DRYSMC, HC, MAXSMC, REFSMC, SATPSI, SATDK, SATDW, WLTSMC, QTZ,    &
                                            BVIC, AXAJ, BXAJ, XXAJ, BDVIC, BBVIC, GDVIC
    ! general parameters
    real(kind=r8)                        :: CSOIL_DATA, REFDK_DATA, REFKDT_DATA, FRZK_DATA, ZBOT_DATA, CZIL_DATA
    real(kind=r8), dimension(NUM_SLOPE)  :: SLOPE_DATA
    namelist / noahmp_general_parameters /  SLOPE_DATA, CSOIL_DATA, REFDK_DATA, REFKDT_DATA, FRZK_DATA, ZBOT_DATA,   &
                                            CZIL_DATA
    ! radiation parameters
    real(kind=r8)                        :: BETADS, BETAIS, EICE
    real(kind=r8), dimension(MBAND)      :: ALBICE, ALBLAK, OMEGAS 
    real(kind=r8), dimension(2)          :: EG
    real(kind=r8), dimension(MSC)        :: ALBSAT_VIS, ALBSAT_NIR, ALBDRY_VIS, ALBDRY_NIR
    namelist / noahmp_rad_parameters /      ALBSAT_VIS, ALBSAT_NIR, ALBDRY_VIS, ALBDRY_NIR, ALBICE, ALBLAK, OMEGAS,      &
                                            BETADS, BETAIS, EG, EICE
    ! global parameters
    real(kind=r8)                        :: CO2, O2, TIMEAN, FSATMX, Z0SNO, SSI, SNOW_RET_FAC ,SNOW_EMIS, SWEMX, TAU0,   &
                                            GRAIN_GROWTH, EXTRA_GROWTH, DIRT_SOOT, BATS_COSZ, BATS_VIS_NEW,              &
                                            BATS_NIR_NEW, BATS_VIS_AGE, BATS_NIR_AGE, BATS_VIS_DIR, BATS_NIR_DIR,        &
                                            RSURF_SNOW, RSURF_EXP, C2_SNOWCOMPACT, C3_SNOWCOMPACT, C4_SNOWCOMPACT,       &
                                            C5_SNOWCOMPACT, DM_SNOWCOMPACT, ETA0_SNOWCOMPACT, SNLIQMAXFRAC, SWEMAXGLA,   &
                                            WSLMAX, ROUS, CMIC, SNOWDEN_MAX, CLASS_ALB_REF, CLASS_SNO_AGE, CLASS_ALB_NEW,&
                                            PSIWLT, Z0SOIL, Z0LAKE
    namelist / noahmp_global_parameters /   CO2, O2, TIMEAN, FSATMX, Z0SNO, SSI, SNOW_RET_FAC ,SNOW_EMIS, SWEMX, TAU0,   &
                                            GRAIN_GROWTH, EXTRA_GROWTH, DIRT_SOOT, BATS_COSZ, BATS_VIS_NEW,              &
                                            BATS_NIR_NEW, BATS_VIS_AGE, BATS_NIR_AGE, BATS_VIS_DIR, BATS_NIR_DIR,        &
                                            RSURF_SNOW, RSURF_EXP, C2_SNOWCOMPACT, C3_SNOWCOMPACT, C4_SNOWCOMPACT,       &
                                            C5_SNOWCOMPACT, DM_SNOWCOMPACT, ETA0_SNOWCOMPACT, SNLIQMAXFRAC, SWEMAXGLA,   &
                                            WSLMAX, ROUS, CMIC, SNOWDEN_MAX, CLASS_ALB_REF, CLASS_SNO_AGE, CLASS_ALB_NEW,&
                                            PSIWLT, Z0SOIL, Z0LAKE
    ! irrigation parameters
    integer                              :: IRR_HAR
    real(kind=r8)                        :: IRR_FRAC, IRR_LAI, IRR_MAD, FILOSS, SPRIR_RATE, MICIR_RATE, FIRTFAC, IR_RAIN
    namelist / noahmp_irrigation_parameters / IRR_FRAC, IRR_HAR, IRR_LAI, IRR_MAD, FILOSS, SPRIR_RATE, MICIR_RATE,       &
                                              FIRTFAC, IR_RAIN
    ! crop parameters
    integer                              :: DEFAULT_CROP
    integer      , dimension(NCROP)      :: PLTDAY, HSDAY
    real(kind=r8), dimension(NCROP)      :: PLANTPOP, IRRI, GDDTBASE, GDDTCUT, GDDS1, GDDS2, GDDS3, GDDS4, GDDS5, C3PSNI,&
                                            KC25I, AKCI, KO25I, AKOI, AVCMXI, VCMX25I, BPI, MPII, FOLNMXI, QE25I, AREF,  &
                                            PSNRF, I2PAR, TASSIM0, TASSIM1, TASSIM2, K, EPSI, Q10MR, LEFREEZ,            &
                                            DILE_FC_S1, DILE_FC_S2, DILE_FC_S3, DILE_FC_S4, DILE_FC_S5, DILE_FC_S6,      &
                                            DILE_FC_S7, DILE_FC_S8, DILE_FW_S1, DILE_FW_S2, DILE_FW_S3, DILE_FW_S4,      &
                                            DILE_FW_S5, DILE_FW_S6, DILE_FW_S7, DILE_FW_S8, FRA_GR, LF_OVRC_S1,          &
                                            LF_OVRC_S2, LF_OVRC_S3, LF_OVRC_S4, LF_OVRC_S5, LF_OVRC_S6, LF_OVRC_S7,      &
                                            LF_OVRC_S8, ST_OVRC_S1, ST_OVRC_S2, ST_OVRC_S3, ST_OVRC_S4, ST_OVRC_S5,      &
                                            ST_OVRC_S6, ST_OVRC_S7, ST_OVRC_S8, RT_OVRC_S1, RT_OVRC_S2, RT_OVRC_S3,      &
                                            RT_OVRC_S4, RT_OVRC_S5, RT_OVRC_S6, RT_OVRC_S7, RT_OVRC_S8, LFMR25, STMR25,  &
                                            RTMR25, GRAINMR25, LFPT_S1, LFPT_S2, LFPT_S3, LFPT_S4, LFPT_S5, LFPT_S6,     &
                                            LFPT_S7, LFPT_S8, STPT_S1, STPT_S2, STPT_S3, STPT_S4, STPT_S5, STPT_S6,      &
                                            STPT_S7, STPT_S8, RTPT_S1, RTPT_S2, RTPT_S3, RTPT_S4, RTPT_S5, RTPT_S6,      &
                                            RTPT_S7, RTPT_S8, GRAINPT_S1, GRAINPT_S2, GRAINPT_S3, GRAINPT_S4, GRAINPT_S5,&
                                            GRAINPT_S6, GRAINPT_S7, GRAINPT_S8, LFCT_S1, LFCT_S2, LFCT_S3, LFCT_S4,      &
                                            LFCT_S5, LFCT_S6, LFCT_S7, LFCT_S8, STCT_S1, STCT_S2, STCT_S3, STCT_S4,      &
                                            STCT_S5, STCT_S6, STCT_S7, STCT_S8, RTCT_S1, RTCT_S2, RTCT_S3, RTCT_S4,      &
                                            RTCT_S5, RTCT_S6, RTCT_S7, RTCT_S8, BIO2LAI
    namelist / noahmp_crop_parameters /     DEFAULT_CROP, PLTDAY, HSDAY, PLANTPOP, IRRI, GDDTBASE, GDDTCUT, GDDS1, GDDS2,&
                                            GDDS3, GDDS4, GDDS5, C3PSNI, KC25I, AKCI, KO25I, AKOI, AVCMXI, VCMX25I, BPI, &
                                            MPII, FOLNMXI, QE25I, AREF, PSNRF, I2PAR, TASSIM0, TASSIM1, TASSIM2, K,      &
                                            EPSI,Q10MR, LEFREEZ, DILE_FC_S1, DILE_FC_S2, DILE_FC_S3, DILE_FC_S4,         &
                                            DILE_FC_S5, DILE_FC_S6, DILE_FC_S7, DILE_FC_S8, DILE_FW_S1, DILE_FW_S2,      &
                                            DILE_FW_S3, DILE_FW_S4, DILE_FW_S5, DILE_FW_S6, DILE_FW_S7, DILE_FW_S8,      &
                                            FRA_GR, LF_OVRC_S1, LF_OVRC_S2, LF_OVRC_S3, LF_OVRC_S4, LF_OVRC_S5,          &
                                            LF_OVRC_S6, LF_OVRC_S7, LF_OVRC_S8, ST_OVRC_S1, ST_OVRC_S2, ST_OVRC_S3,      &
                                            ST_OVRC_S4, ST_OVRC_S5, ST_OVRC_S6, ST_OVRC_S7, ST_OVRC_S8, RT_OVRC_S1,      &
                                            RT_OVRC_S2, RT_OVRC_S3, RT_OVRC_S4, RT_OVRC_S5, RT_OVRC_S6, RT_OVRC_S7,      &
                                            RT_OVRC_S8, LFMR25, STMR25, RTMR25, GRAINMR25, LFPT_S1, LFPT_S2, LFPT_S3,    &
                                            LFPT_S4, LFPT_S5, LFPT_S6, LFPT_S7, LFPT_S8, STPT_S1, STPT_S2, STPT_S3,      &
                                            STPT_S4, STPT_S5, STPT_S6, STPT_S7, STPT_S8, RTPT_S1, RTPT_S2, RTPT_S3,      &
                                            RTPT_S4, RTPT_S5, RTPT_S6, RTPT_S7, RTPT_S8, GRAINPT_S1, GRAINPT_S2,         &
                                            GRAINPT_S3, GRAINPT_S4, GRAINPT_S5, GRAINPT_S6, GRAINPT_S7, GRAINPT_S8,      &
                                            LFCT_S1, LFCT_S2, LFCT_S3, LFCT_S4, LFCT_S5, LFCT_S6, LFCT_S7, LFCT_S8,      &
                                            STCT_S1, STCT_S2, STCT_S3, STCT_S4, STCT_S5, STCT_S6, STCT_S7, STCT_S8,      &
                                            RTCT_S1, RTCT_S2, RTCT_S3, RTCT_S4, RTCT_S5, RTCT_S6, RTCT_S7, RTCT_S8,      &
                                            BIO2LAI
    ! tile drainage parameters
    integer                               :: NSOILTYPE, DRAIN_LAYER_OPT
    integer      , dimension(MAX_SOILTYP) :: TD_DEPTH
    real(kind=r8), dimension(MAX_SOILTYP) :: TDSMC_FAC, TD_DC, TD_DCOEF, TD_D, TD_ADEPTH, TD_RADI, TD_SPAC,         &
                                             TD_DDRAIN, KLAT_FAC
    namelist / noahmp_tiledrain_parameters / NSOILTYPE, DRAIN_LAYER_OPT, TDSMC_FAC, TD_DEPTH, TD_DC, TD_DCOEF,      &
                                             TD_D, TD_ADEPTH, TD_RADI, TD_SPAC, TD_DDRAIN, KLAT_FAC
    ! optional parameters
    real(kind=r8)                         :: SR2006_THETA_1500T_A, SR2006_THETA_1500T_B, SR2006_THETA_1500T_C,      &
                                             SR2006_THETA_1500T_D, SR2006_THETA_1500T_E, SR2006_THETA_1500T_F,      &
                                             SR2006_THETA_1500T_G, SR2006_THETA_1500_A , SR2006_THETA_1500_B,       &
                                             SR2006_THETA_33T_A, SR2006_THETA_33T_B, SR2006_THETA_33T_C,            &
                                             SR2006_THETA_33T_D, SR2006_THETA_33T_E, SR2006_THETA_33T_F,            &
                                             SR2006_THETA_33T_G, SR2006_THETA_33_A, SR2006_THETA_33_B,              &
                                             SR2006_THETA_33_C, SR2006_THETA_S33T_A, SR2006_THETA_S33T_B,           &
                                             SR2006_THETA_S33T_C, SR2006_THETA_S33T_D, SR2006_THETA_S33T_E,         &
                                             SR2006_THETA_S33T_F, SR2006_THETA_S33T_G, SR2006_THETA_S33_A,          &
                                             SR2006_THETA_S33_B, SR2006_PSI_ET_A, SR2006_PSI_ET_B, SR2006_PSI_ET_C, &
                                             SR2006_PSI_ET_D, SR2006_PSI_ET_E, SR2006_PSI_ET_F, SR2006_PSI_ET_G,    &
                                             SR2006_PSI_E_A, SR2006_PSI_E_B, SR2006_PSI_E_C, SR2006_SMCMAX_A,       &
                                             SR2006_SMCMAX_B
    namelist / noahmp_optional_parameters /  SR2006_THETA_1500T_A, SR2006_THETA_1500T_B, SR2006_THETA_1500T_C,      &
                                             SR2006_THETA_1500T_D, SR2006_THETA_1500T_E, SR2006_THETA_1500T_F,      &
                                             SR2006_THETA_1500T_G, SR2006_THETA_1500_A, SR2006_THETA_1500_B,        &
                                             SR2006_THETA_33T_A, SR2006_THETA_33T_B, SR2006_THETA_33T_C,            &
                                             SR2006_THETA_33T_D, SR2006_THETA_33T_E, SR2006_THETA_33T_F,            &
                                             SR2006_THETA_33T_G, SR2006_THETA_33_A, SR2006_THETA_33_B,              &
                                             SR2006_THETA_33_C, SR2006_THETA_S33T_A, SR2006_THETA_S33T_B,           &
                                             SR2006_THETA_S33T_C, SR2006_THETA_S33T_D, SR2006_THETA_S33T_E,         &
                                             SR2006_THETA_S33T_F, SR2006_THETA_S33T_G, SR2006_THETA_S33_A,          &
                                             SR2006_THETA_S33_B, SR2006_PSI_ET_A, SR2006_PSI_ET_B, SR2006_PSI_ET_C, &
                                             SR2006_PSI_ET_D, SR2006_PSI_ET_E, SR2006_PSI_ET_F, SR2006_PSI_ET_G,    &
                                             SR2006_PSI_E_A, SR2006_PSI_E_B, SR2006_PSI_E_C, SR2006_SMCMAX_A,       &
                                             SR2006_SMCMAX_B
    logical :: file_exist
    integer :: noahmp_table_unit
    character(len=*), parameter :: subname=trim(modName)//':(NoahmpReadTable) '
    !-------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ----------------------
    ! Open NoahmpTable.TBL
    ! Since we are not opening/closing file in each read, the order of read is important
    ! ----------------------

    inquire(file='NoahmpTable.TBL', exist=file_exist)

    if (file_exist) then
       call ESMF_UtilIOUnitGet(unit=noahmp_table_unit, rc=rcode)
       open(noahmp_table_unit, file="NoahmpTable.TBL", status='old', form='formatted', action='read', iostat=ierr)
    else
       call ESMF_LogWrite(trim(subname)//': Cannot find file NoahmpTable.TBL.', ESMF_LOGMSG_ERROR) 
       rcode = ESMF_FAILURE
       return
    end if

    ! ----------------------
    ! Read vegetation parameters 
    ! ----------------------

    if (trim(model%noahmp%config%domain%LandUseDataName) == "USGS") then
       read(noahmp_table_unit, noahmp_usgs_veg_categories)
       read(noahmp_table_unit, noahmp_usgs_parameters)
    elseif (trim(model%noahmp%config%domain%LandUseDataName) == "MODIFIED_IGBP_MODIS_NOAH") then
       read(noahmp_table_unit, noahmp_modis_veg_categories)
       read(noahmp_table_unit, noahmp_modis_parameters)
    else
       call ESMF_LogWrite(trim(subname)//': Unrecognized LandUseDataName!', ESMF_LOGMSG_ERROR)
       rcode = ESMF_FAILURE
       close(noahmp_table_unit)
       return
    endif

    ! assign values
    model%coupling%table%ISURBAN_TABLE         = ISURBAN
    model%coupling%table%ISWATER_TABLE         = ISWATER
    model%coupling%table%ISBARREN_TABLE        = ISBARREN
    model%coupling%table%ISICE_TABLE           = ISICE
    model%coupling%table%ISCROP_TABLE          = ISCROP
    model%coupling%table%EBLFOREST_TABLE       = EBLFOREST
    model%coupling%table%NATURAL_TABLE         = NATURAL
    model%coupling%table%URBTYPE_beg           = URBTYPE_beg
    model%coupling%table%LCZ_1_TABLE           = LCZ_1
    model%coupling%table%LCZ_2_TABLE           = LCZ_2
    model%coupling%table%LCZ_3_TABLE           = LCZ_3
    model%coupling%table%LCZ_4_TABLE           = LCZ_4
    model%coupling%table%LCZ_5_TABLE           = LCZ_5
    model%coupling%table%LCZ_6_TABLE           = LCZ_6
    model%coupling%table%LCZ_7_TABLE           = LCZ_7
    model%coupling%table%LCZ_8_TABLE           = LCZ_8
    model%coupling%table%LCZ_9_TABLE           = LCZ_9
    model%coupling%table%LCZ_10_TABLE          = LCZ_10
    model%coupling%table%LCZ_11_TABLE          = LCZ_11
    model%coupling%table%CH2OP_TABLE  (1:NVEG) = CH2OP  (1:NVEG)
    model%coupling%table%DLEAF_TABLE  (1:NVEG) = DLEAF  (1:NVEG)
    model%coupling%table%Z0MVT_TABLE  (1:NVEG) = Z0MVT  (1:NVEG)
    model%coupling%table%HVT_TABLE    (1:NVEG) = HVT    (1:NVEG)
    model%coupling%table%HVB_TABLE    (1:NVEG) = HVB    (1:NVEG)
    model%coupling%table%DEN_TABLE    (1:NVEG) = DEN    (1:NVEG)
    model%coupling%table%RC_TABLE     (1:NVEG) = RC     (1:NVEG)
    model%coupling%table%MFSNO_TABLE  (1:NVEG) = MFSNO  (1:NVEG)
    model%coupling%table%SCFFAC_TABLE (1:NVEG) = SCFFAC (1:NVEG)
    model%coupling%table%CBIOM_TABLE  (1:NVEG) = CBIOM  (1:NVEG)
    model%coupling%table%XL_TABLE     (1:NVEG) = XL     (1:NVEG)
    model%coupling%table%CWPVT_TABLE  (1:NVEG) = CWPVT  (1:NVEG)
    model%coupling%table%C3PSN_TABLE  (1:NVEG) = C3PSN  (1:NVEG)
    model%coupling%table%KC25_TABLE   (1:NVEG) = KC25   (1:NVEG)
    model%coupling%table%AKC_TABLE    (1:NVEG) = AKC    (1:NVEG)
    model%coupling%table%KO25_TABLE   (1:NVEG) = KO25   (1:NVEG)
    model%coupling%table%AKO_TABLE    (1:NVEG) = AKO    (1:NVEG)
    model%coupling%table%AVCMX_TABLE  (1:NVEG) = AVCMX  (1:NVEG)
    model%coupling%table%AQE_TABLE    (1:NVEG) = AQE    (1:NVEG)
    model%coupling%table%LTOVRC_TABLE (1:NVEG) = LTOVRC (1:NVEG)
    model%coupling%table%DILEFC_TABLE (1:NVEG) = DILEFC (1:NVEG)
    model%coupling%table%DILEFW_TABLE (1:NVEG) = DILEFW (1:NVEG)
    model%coupling%table%RMF25_TABLE  (1:NVEG) = RMF25  (1:NVEG)
    model%coupling%table%SLA_TABLE    (1:NVEG) = SLA    (1:NVEG)
    model%coupling%table%FRAGR_TABLE  (1:NVEG) = FRAGR  (1:NVEG)
    model%coupling%table%TMIN_TABLE   (1:NVEG) = TMIN   (1:NVEG)
    model%coupling%table%VCMX25_TABLE (1:NVEG) = VCMX25 (1:NVEG)
    model%coupling%table%TDLEF_TABLE  (1:NVEG) = TDLEF  (1:NVEG)
    model%coupling%table%BP_TABLE     (1:NVEG) = BP     (1:NVEG)
    model%coupling%table%MP_TABLE     (1:NVEG) = MP     (1:NVEG)
    model%coupling%table%QE25_TABLE   (1:NVEG) = QE25   (1:NVEG)
    model%coupling%table%RMS25_TABLE  (1:NVEG) = RMS25  (1:NVEG)
    model%coupling%table%RMR25_TABLE  (1:NVEG) = RMR25  (1:NVEG)
    model%coupling%table%ARM_TABLE    (1:NVEG) = ARM    (1:NVEG)
    model%coupling%table%FOLNMX_TABLE (1:NVEG) = FOLNMX (1:NVEG)
    model%coupling%table%WDPOOL_TABLE (1:NVEG) = WDPOOL (1:NVEG)
    model%coupling%table%WRRAT_TABLE  (1:NVEG) = WRRAT  (1:NVEG)
    model%coupling%table%MRP_TABLE    (1:NVEG) = MRP    (1:NVEG)
    model%coupling%table%NROOT_TABLE  (1:NVEG) = NROOT  (1:NVEG)
    model%coupling%table%RGL_TABLE    (1:NVEG) = RGL    (1:NVEG)
    model%coupling%table%RS_TABLE     (1:NVEG) = RS     (1:NVEG)
    model%coupling%table%HS_TABLE     (1:NVEG) = HS     (1:NVEG)
    model%coupling%table%TOPT_TABLE   (1:NVEG) = TOPT   (1:NVEG)
    model%coupling%table%RSMAX_TABLE  (1:NVEG) = RSMAX  (1:NVEG)
    model%coupling%table%RTOVRC_TABLE (1:NVEG) = RTOVRC (1:NVEG)
    model%coupling%table%RSWOODC_TABLE(1:NVEG) = RSWOODC(1:NVEG)
    model%coupling%table%BF_TABLE     (1:NVEG) = BF     (1:NVEG)
    model%coupling%table%WSTRC_TABLE  (1:NVEG) = WSTRC  (1:NVEG)
    model%coupling%table%LAIMIN_TABLE (1:NVEG) = LAIMIN (1:NVEG)
    model%coupling%table%XSAMIN_TABLE (1:NVEG) = XSAMIN (1:NVEG)

    model%coupling%table%SAIM_TABLE(1:NVEG, 1) = SAI_JAN(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG, 2) = SAI_FEB(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG, 3) = SAI_MAR(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG, 4) = SAI_APR(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG, 5) = SAI_MAY(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG, 6) = SAI_JUN(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG, 7) = SAI_JUL(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG, 8) = SAI_AUG(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG, 9) = SAI_SEP(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG,10) = SAI_OCT(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG,11) = SAI_NOV(1:NVEG)
    model%coupling%table%SAIM_TABLE(1:NVEG,12) = SAI_DEC(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG, 1) = LAI_JAN(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG, 2) = LAI_FEB(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG, 3) = LAI_MAR(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG, 4) = LAI_APR(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG, 5) = LAI_MAY(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG, 6) = LAI_JUN(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG, 7) = LAI_JUL(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG, 8) = LAI_AUG(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG, 9) = LAI_SEP(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG,10) = LAI_OCT(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG,11) = LAI_NOV(1:NVEG)
    model%coupling%table%LAIM_TABLE(1:NVEG,12) = LAI_DEC(1:NVEG)
    model%coupling%table%RHOL_TABLE(1:NVEG,1)  = RHOL_VIS(1:NVEG) !leaf reflectance: 1=vis, 2=nir
    model%coupling%table%RHOL_TABLE(1:NVEG,2)  = RHOL_NIR(1:NVEG) !leaf reflectance: 1=vis, 2=nir
    model%coupling%table%RHOS_TABLE(1:NVEG,1)  = RHOS_VIS(1:NVEG) !stem reflectance: 1=vis, 2=nir
    model%coupling%table%RHOS_TABLE(1:NVEG,2)  = RHOS_NIR(1:NVEG) !stem reflectance: 1=vis, 2=nir
    model%coupling%table%TAUL_TABLE(1:NVEG,1)  = TAUL_VIS(1:NVEG) !leaf transmittance: 1=vis, 2=nir
    model%coupling%table%TAUL_TABLE(1:NVEG,2)  = TAUL_NIR(1:NVEG) !leaf transmittance: 1=vis, 2=nir
    model%coupling%table%TAUS_TABLE(1:NVEG,1)  = TAUS_VIS(1:NVEG) !stem transmittance: 1=vis, 2=nir
    model%coupling%table%TAUS_TABLE(1:NVEG,2)  = TAUS_NIR(1:NVEG) !stem transmittance: 1=vis, 2=nir

    ! ----------------------
    ! Read radiation parameters 
    ! ----------------------

    read(noahmp_table_unit, noahmp_rad_parameters)

    ! assign values
    model%coupling%table%ALBSAT_TABLE(:,1) = ALBSAT_VIS ! saturated soil albedos: 1=vis, 2=nir
    model%coupling%table%ALBSAT_TABLE(:,2) = ALBSAT_NIR ! saturated soil albedos: 1=vis, 2=nir
    model%coupling%table%ALBDRY_TABLE(:,1) = ALBDRY_VIS ! dry soil albedos: 1=vis, 2=nir
    model%coupling%table%ALBDRY_TABLE(:,2) = ALBDRY_NIR ! dry soil albedos: 1=vis, 2=nir
    model%coupling%table%ALBICE_TABLE      = ALBICE
    model%coupling%table%ALBLAK_TABLE      = ALBLAK
    model%coupling%table%OMEGAS_TABLE      = OMEGAS
    model%coupling%table%BETADS_TABLE      = BETADS
    model%coupling%table%BETAIS_TABLE      = BETAIS
    model%coupling%table%EG_TABLE          = EG
    model%coupling%table%EICE_TABLE        = EICE

    ! ----------------------
    ! Read global parameters 
    ! ----------------------

    read(noahmp_table_unit, noahmp_global_parameters)

    ! assign values
    model%coupling%table%CO2_TABLE              = CO2
    model%coupling%table%O2_TABLE               = O2
    model%coupling%table%TIMEAN_TABLE           = TIMEAN
    model%coupling%table%FSATMX_TABLE           = FSATMX
    model%coupling%table%Z0SNO_TABLE            = Z0SNO
    model%coupling%table%SSI_TABLE              = SSI
    model%coupling%table%SNOW_RET_FAC_TABLE     = SNOW_RET_FAC
    model%coupling%table%SNOW_EMIS_TABLE        = SNOW_EMIS
    model%coupling%table%SWEMX_TABLE            = SWEMX
    model%coupling%table%TAU0_TABLE             = TAU0
    model%coupling%table%GRAIN_GROWTH_TABLE     = GRAIN_GROWTH
    model%coupling%table%EXTRA_GROWTH_TABLE     = EXTRA_GROWTH
    model%coupling%table%DIRT_SOOT_TABLE        = DIRT_SOOT
    model%coupling%table%BATS_COSZ_TABLE        = BATS_COSZ
    model%coupling%table%BATS_VIS_NEW_TABLE     = BATS_VIS_NEW
    model%coupling%table%BATS_NIR_NEW_TABLE     = BATS_NIR_NEW
    model%coupling%table%BATS_VIS_AGE_TABLE     = BATS_VIS_AGE
    model%coupling%table%BATS_NIR_AGE_TABLE     = BATS_NIR_AGE
    model%coupling%table%BATS_VIS_DIR_TABLE     = BATS_VIS_DIR
    model%coupling%table%BATS_NIR_DIR_TABLE     = BATS_NIR_DIR
    model%coupling%table%RSURF_SNOW_TABLE       = RSURF_SNOW
    model%coupling%table%RSURF_EXP_TABLE        = RSURF_EXP
    model%coupling%table%C2_SNOWCOMPACT_TABLE   = C2_SNOWCOMPACT
    model%coupling%table%C3_SNOWCOMPACT_TABLE   = C3_SNOWCOMPACT
    model%coupling%table%C4_SNOWCOMPACT_TABLE   = C4_SNOWCOMPACT
    model%coupling%table%C5_SNOWCOMPACT_TABLE   = C5_SNOWCOMPACT
    model%coupling%table%DM_SNOWCOMPACT_TABLE   = DM_SNOWCOMPACT
    model%coupling%table%ETA0_SNOWCOMPACT_TABLE = ETA0_SNOWCOMPACT
    model%coupling%table%SNLIQMAXFRAC_TABLE     = SNLIQMAXFRAC
    model%coupling%table%SWEMAXGLA_TABLE        = SWEMAXGLA
    model%coupling%table%WSLMAX_TABLE           = WSLMAX
    model%coupling%table%ROUS_TABLE             = ROUS
    model%coupling%table%CMIC_TABLE             = CMIC
    model%coupling%table%SNOWDEN_MAX_TABLE      = SNOWDEN_MAX
    model%coupling%table%CLASS_ALB_REF_TABLE    = CLASS_ALB_REF
    model%coupling%table%CLASS_SNO_AGE_TABLE    = CLASS_SNO_AGE
    model%coupling%table%CLASS_ALB_NEW_TABLE    = CLASS_ALB_NEW
    model%coupling%table%PSIWLT_TABLE           = PSIWLT
    model%coupling%table%Z0SOIL_TABLE           = Z0SOIL
    model%coupling%table%Z0LAKE_TABLE           = Z0LAKE

    ! ----------------------
    ! Read irrigation parameters
    ! ----------------------

    read(noahmp_table_unit, noahmp_irrigation_parameters)
    if ((FILOSS < 0.0) .or. (FILOSS > 0.99)) then
       write(*,'("WARNING: FILOSS should be >=0.0 and <=0.99")')
       stop "STOP in NoahMP_irrigation_parameters"
       call ESMF_LogWrite(trim(subname)//': FILOSS should be >=0.0 and <=0.99!', ESMF_LOGMSG_ERROR)
       call ESMF_LogWrite(trim(subname)//': STOP while reading NoahMP irrigation parameters', ESMF_LOGMSG_ERROR)
       rcode = ESMF_FAILURE
       close(noahmp_table_unit)
       return
    endif

    ! assign values
    model%coupling%table%IRR_FRAC_TABLE   = IRR_FRAC
    model%coupling%table%IRR_HAR_TABLE    = IRR_HAR
    model%coupling%table%IRR_LAI_TABLE    = IRR_LAI
    model%coupling%table%IRR_MAD_TABLE    = IRR_MAD
    model%coupling%table%FILOSS_TABLE     = FILOSS  
    model%coupling%table%SPRIR_RATE_TABLE = SPRIR_RATE
    model%coupling%table%MICIR_RATE_TABLE = MICIR_RATE
    model%coupling%table%FIRTFAC_TABLE    = FIRTFAC
    model%coupling%table%IR_RAIN_TABLE    = IR_RAIN 

    ! ----------------------
    ! Read crop parameters
    ! ----------------------

    read(noahmp_table_unit, noahmp_crop_parameters)

    ! assign values
    model%coupling%table%DEFAULT_CROP_TABLE     = DEFAULT_CROP
    model%coupling%table%PLTDAY_TABLE           = PLTDAY
    model%coupling%table%HSDAY_TABLE            = HSDAY
    model%coupling%table%PLANTPOP_TABLE         = PLANTPOP
    model%coupling%table%IRRI_TABLE             = IRRI
    model%coupling%table%GDDTBASE_TABLE         = GDDTBASE
    model%coupling%table%GDDTCUT_TABLE          = GDDTCUT
    model%coupling%table%GDDS1_TABLE            = GDDS1
    model%coupling%table%GDDS2_TABLE            = GDDS2
    model%coupling%table%GDDS3_TABLE            = GDDS3
    model%coupling%table%GDDS4_TABLE            = GDDS4
    model%coupling%table%GDDS5_TABLE            = GDDS5
    model%coupling%table%C3PSNI_TABLE (1:5)     = C3PSNI (1:5)
    model%coupling%table%KC25I_TABLE  (1:5)     = KC25I  (1:5)
    model%coupling%table%AKCI_TABLE   (1:5)     = AKCI   (1:5)
    model%coupling%table%KO25I_TABLE  (1:5)     = KO25I  (1:5)
    model%coupling%table%AKOI_TABLE   (1:5)     = AKOI   (1:5)
    model%coupling%table%AVCMXI_TABLE (1:5)     = AVCMXI (1:5)
    model%coupling%table%VCMX25I_TABLE(1:5)     = VCMX25I(1:5)
    model%coupling%table%BPI_TABLE    (1:5)     = BPI    (1:5)
    model%coupling%table%MPI_TABLE    (1:5)     = MPII   (1:5)
    model%coupling%table%FOLNMXI_TABLE(1:5)     = FOLNMXI(1:5)
    model%coupling%table%QE25I_TABLE  (1:5)     = QE25I  (1:5)
    model%coupling%table%AREF_TABLE             = AREF
    model%coupling%table%PSNRF_TABLE            = PSNRF
    model%coupling%table%I2PAR_TABLE            = I2PAR
    model%coupling%table%TASSIM0_TABLE          = TASSIM0
    model%coupling%table%TASSIM1_TABLE          = TASSIM1
    model%coupling%table%TASSIM2_TABLE          = TASSIM2
    model%coupling%table%K_TABLE                = K
    model%coupling%table%EPSI_TABLE             = EPSI
    model%coupling%table%Q10MR_TABLE            = Q10MR
    model%coupling%table%LEFREEZ_TABLE          = LEFREEZ
    model%coupling%table%FRA_GR_TABLE           = FRA_GR
    model%coupling%table%LFMR25_TABLE           = LFMR25
    model%coupling%table%STMR25_TABLE           = STMR25
    model%coupling%table%RTMR25_TABLE           = RTMR25
    model%coupling%table%GRAINMR25_TABLE        = GRAINMR25
    model%coupling%table%BIO2LAI_TABLE          = BIO2LAI
    model%coupling%table%DILE_FC_TABLE(:,1)     = DILE_FC_S1
    model%coupling%table%DILE_FC_TABLE(:,2)     = DILE_FC_S2
    model%coupling%table%DILE_FC_TABLE(:,3)     = DILE_FC_S3
    model%coupling%table%DILE_FC_TABLE(:,4)     = DILE_FC_S4
    model%coupling%table%DILE_FC_TABLE(:,5)     = DILE_FC_S5
    model%coupling%table%DILE_FC_TABLE(:,6)     = DILE_FC_S6
    model%coupling%table%DILE_FC_TABLE(:,7)     = DILE_FC_S7
    model%coupling%table%DILE_FC_TABLE(:,8)     = DILE_FC_S8
    model%coupling%table%DILE_FW_TABLE(:,1)     = DILE_FW_S1
    model%coupling%table%DILE_FW_TABLE(:,2)     = DILE_FW_S2
    model%coupling%table%DILE_FW_TABLE(:,3)     = DILE_FW_S3
    model%coupling%table%DILE_FW_TABLE(:,4)     = DILE_FW_S4
    model%coupling%table%DILE_FW_TABLE(:,5)     = DILE_FW_S5
    model%coupling%table%DILE_FW_TABLE(:,6)     = DILE_FW_S6
    model%coupling%table%DILE_FW_TABLE(:,7)     = DILE_FW_S7
    model%coupling%table%DILE_FW_TABLE(:,8)     = DILE_FW_S8
    model%coupling%table%LF_OVRC_TABLE(:,1)     = LF_OVRC_S1
    model%coupling%table%LF_OVRC_TABLE(:,2)     = LF_OVRC_S2
    model%coupling%table%LF_OVRC_TABLE(:,3)     = LF_OVRC_S3
    model%coupling%table%LF_OVRC_TABLE(:,4)     = LF_OVRC_S4
    model%coupling%table%LF_OVRC_TABLE(:,5)     = LF_OVRC_S5
    model%coupling%table%LF_OVRC_TABLE(:,6)     = LF_OVRC_S6
    model%coupling%table%LF_OVRC_TABLE(:,7)     = LF_OVRC_S7
    model%coupling%table%LF_OVRC_TABLE(:,8)     = LF_OVRC_S8
    model%coupling%table%ST_OVRC_TABLE(:,1)     = ST_OVRC_S1
    model%coupling%table%ST_OVRC_TABLE(:,2)     = ST_OVRC_S2
    model%coupling%table%ST_OVRC_TABLE(:,3)     = ST_OVRC_S3
    model%coupling%table%ST_OVRC_TABLE(:,4)     = ST_OVRC_S4
    model%coupling%table%ST_OVRC_TABLE(:,5)     = ST_OVRC_S5
    model%coupling%table%ST_OVRC_TABLE(:,6)     = ST_OVRC_S6
    model%coupling%table%ST_OVRC_TABLE(:,7)     = ST_OVRC_S7
    model%coupling%table%ST_OVRC_TABLE(:,8)     = ST_OVRC_S8
    model%coupling%table%RT_OVRC_TABLE(:,1)     = RT_OVRC_S1
    model%coupling%table%RT_OVRC_TABLE(:,2)     = RT_OVRC_S2
    model%coupling%table%RT_OVRC_TABLE(:,3)     = RT_OVRC_S3
    model%coupling%table%RT_OVRC_TABLE(:,4)     = RT_OVRC_S4
    model%coupling%table%RT_OVRC_TABLE(:,5)     = RT_OVRC_S5
    model%coupling%table%RT_OVRC_TABLE(:,6)     = RT_OVRC_S6
    model%coupling%table%RT_OVRC_TABLE(:,7)     = RT_OVRC_S7
    model%coupling%table%RT_OVRC_TABLE(:,8)     = RT_OVRC_S8
    model%coupling%table%LFPT_TABLE   (:,1)     = LFPT_S1
    model%coupling%table%LFPT_TABLE   (:,2)     = LFPT_S2
    model%coupling%table%LFPT_TABLE   (:,3)     = LFPT_S3
    model%coupling%table%LFPT_TABLE   (:,4)     = LFPT_S4
    model%coupling%table%LFPT_TABLE   (:,5)     = LFPT_S5
    model%coupling%table%LFPT_TABLE   (:,6)     = LFPT_S6
    model%coupling%table%LFPT_TABLE   (:,7)     = LFPT_S7
    model%coupling%table%LFPT_TABLE   (:,8)     = LFPT_S8
    model%coupling%table%STPT_TABLE   (:,1)     = STPT_S1
    model%coupling%table%STPT_TABLE   (:,2)     = STPT_S2
    model%coupling%table%STPT_TABLE   (:,3)     = STPT_S3
    model%coupling%table%STPT_TABLE   (:,4)     = STPT_S4
    model%coupling%table%STPT_TABLE   (:,5)     = STPT_S5
    model%coupling%table%STPT_TABLE   (:,6)     = STPT_S6
    model%coupling%table%STPT_TABLE   (:,7)     = STPT_S7
    model%coupling%table%STPT_TABLE   (:,8)     = STPT_S8
    model%coupling%table%RTPT_TABLE   (:,1)     = RTPT_S1
    model%coupling%table%RTPT_TABLE   (:,2)     = RTPT_S2
    model%coupling%table%RTPT_TABLE   (:,3)     = RTPT_S3
    model%coupling%table%RTPT_TABLE   (:,4)     = RTPT_S4
    model%coupling%table%RTPT_TABLE   (:,5)     = RTPT_S5
    model%coupling%table%RTPT_TABLE   (:,6)     = RTPT_S6
    model%coupling%table%RTPT_TABLE   (:,7)     = RTPT_S7
    model%coupling%table%RTPT_TABLE   (:,8)     = RTPT_S8
    model%coupling%table%GRAINPT_TABLE(:,1)     = GRAINPT_S1
    model%coupling%table%GRAINPT_TABLE(:,2)     = GRAINPT_S2
    model%coupling%table%GRAINPT_TABLE(:,3)     = GRAINPT_S3
    model%coupling%table%GRAINPT_TABLE(:,4)     = GRAINPT_S4
    model%coupling%table%GRAINPT_TABLE(:,5)     = GRAINPT_S5
    model%coupling%table%GRAINPT_TABLE(:,6)     = GRAINPT_S6
    model%coupling%table%GRAINPT_TABLE(:,7)     = GRAINPT_S7
    model%coupling%table%GRAINPT_TABLE(:,8)     = GRAINPT_S8
    model%coupling%table%LFCT_TABLE   (:,1)     = LFCT_S1
    model%coupling%table%LFCT_TABLE   (:,2)     = LFCT_S2
    model%coupling%table%LFCT_TABLE   (:,3)     = LFCT_S3
    model%coupling%table%LFCT_TABLE   (:,4)     = LFCT_S4
    model%coupling%table%LFCT_TABLE   (:,5)     = LFCT_S5
    model%coupling%table%LFCT_TABLE   (:,6)     = LFCT_S6
    model%coupling%table%LFCT_TABLE   (:,7)     = LFCT_S7
    model%coupling%table%LFCT_TABLE   (:,8)     = LFCT_S8
    model%coupling%table%STCT_TABLE   (:,1)     = STCT_S1
    model%coupling%table%STCT_TABLE   (:,2)     = STCT_S2
    model%coupling%table%STCT_TABLE   (:,3)     = STCT_S3
    model%coupling%table%STCT_TABLE   (:,4)     = STCT_S4
    model%coupling%table%STCT_TABLE   (:,5)     = STCT_S5
    model%coupling%table%STCT_TABLE   (:,6)     = STCT_S6
    model%coupling%table%STCT_TABLE   (:,7)     = STCT_S7
    model%coupling%table%STCT_TABLE   (:,8)     = STCT_S8
    model%coupling%table%RTCT_TABLE   (:,1)     = RTCT_S1
    model%coupling%table%RTCT_TABLE   (:,2)     = RTCT_S2
    model%coupling%table%RTCT_TABLE   (:,3)     = RTCT_S3
    model%coupling%table%RTCT_TABLE   (:,4)     = RTCT_S4
    model%coupling%table%RTCT_TABLE   (:,5)     = RTCT_S5
    model%coupling%table%RTCT_TABLE   (:,6)     = RTCT_S6
    model%coupling%table%RTCT_TABLE   (:,7)     = RTCT_S7
    model%coupling%table%RTCT_TABLE   (:,8)     = RTCT_S8

    ! ----------------------
    ! Read drainage parameters
    ! ----------------------

    read(noahmp_table_unit, noahmp_tiledrain_parameters)

    ! assign values
    model%coupling%table%DRAIN_LAYER_OPT_TABLE        = DRAIN_LAYER_OPT
    model%coupling%table%TDSMC_FAC_TABLE(1:NSOILTYPE) = TDSMC_FAC(1:NSOILTYPE)
    model%coupling%table%TD_DEPTH_TABLE (1:NSOILTYPE) = TD_DEPTH (1:NSOILTYPE)
    model%coupling%table%TD_DC_TABLE    (1:NSOILTYPE) = TD_DC    (1:NSOILTYPE)
    model%coupling%table%TD_DCOEF_TABLE (1:NSOILTYPE) = TD_DCOEF (1:NSOILTYPE)
    model%coupling%table%TD_D_TABLE     (1:NSOILTYPE) = TD_D     (1:NSOILTYPE)
    model%coupling%table%TD_ADEPTH_TABLE(1:NSOILTYPE) = TD_ADEPTH(1:NSOILTYPE)
    model%coupling%table%TD_RADI_TABLE  (1:NSOILTYPE) = TD_RADI  (1:NSOILTYPE)
    model%coupling%table%TD_SPAC_TABLE  (1:NSOILTYPE) = TD_SPAC  (1:NSOILTYPE)
    model%coupling%table%TD_DDRAIN_TABLE(1:NSOILTYPE) = TD_DDRAIN(1:NSOILTYPE)
    model%coupling%table%KLAT_FAC_TABLE (1:NSOILTYPE) = KLAT_FAC (1:NSOILTYPE)

    ! ----------------------
    ! Read optional parameters
    ! ----------------------

    read(noahmp_table_unit, noahmp_optional_parameters)

    ! assign values
    MODEL%COUPLING%TABLE%SR2006_THETA_1500T_A_TABLE = SR2006_THETA_1500T_A
    MODEL%COUPLING%TABLE%SR2006_THETA_1500T_B_TABLE = SR2006_THETA_1500T_B
    MODEL%COUPLING%TABLE%SR2006_THETA_1500T_C_TABLE = SR2006_THETA_1500T_C
    MODEL%COUPLING%TABLE%SR2006_THETA_1500T_D_TABLE = SR2006_THETA_1500T_D
    MODEL%COUPLING%TABLE%SR2006_THETA_1500T_E_TABLE = SR2006_THETA_1500T_E
    MODEL%COUPLING%TABLE%SR2006_THETA_1500T_F_TABLE = SR2006_THETA_1500T_F
    MODEL%COUPLING%TABLE%SR2006_THETA_1500T_G_TABLE = SR2006_THETA_1500T_G
    MODEL%COUPLING%TABLE%SR2006_THETA_1500_A_TABLE  = SR2006_THETA_1500_A
    MODEL%COUPLING%TABLE%SR2006_THETA_1500_B_TABLE  = SR2006_THETA_1500_B
    MODEL%COUPLING%TABLE%SR2006_THETA_33T_A_TABLE   = SR2006_THETA_33T_A
    MODEL%COUPLING%TABLE%SR2006_THETA_33T_B_TABLE   = SR2006_THETA_33T_B
    MODEL%COUPLING%TABLE%SR2006_THETA_33T_C_TABLE   = SR2006_THETA_33T_C
    MODEL%COUPLING%TABLE%SR2006_THETA_33T_D_TABLE   = SR2006_THETA_33T_D
    MODEL%COUPLING%TABLE%SR2006_THETA_33T_E_TABLE   = SR2006_THETA_33T_E
    MODEL%COUPLING%TABLE%SR2006_THETA_33T_F_TABLE   = SR2006_THETA_33T_F
    MODEL%COUPLING%TABLE%SR2006_THETA_33T_G_TABLE   = SR2006_THETA_33T_G
    MODEL%COUPLING%TABLE%SR2006_THETA_33_A_TABLE    = SR2006_THETA_33_A
    MODEL%COUPLING%TABLE%SR2006_THETA_33_B_TABLE    = SR2006_THETA_33_B
    MODEL%COUPLING%TABLE%SR2006_THETA_33_C_TABLE    = SR2006_THETA_33_C
    MODEL%COUPLING%TABLE%SR2006_THETA_S33T_A_TABLE  = SR2006_THETA_S33T_A
    MODEL%COUPLING%TABLE%SR2006_THETA_S33T_B_TABLE  = SR2006_THETA_S33T_B
    MODEL%COUPLING%TABLE%SR2006_THETA_S33T_C_TABLE  = SR2006_THETA_S33T_C
    MODEL%COUPLING%TABLE%SR2006_THETA_S33T_D_TABLE  = SR2006_THETA_S33T_D
    MODEL%COUPLING%TABLE%SR2006_THETA_S33T_E_TABLE  = SR2006_THETA_S33T_E
    MODEL%COUPLING%TABLE%SR2006_THETA_S33T_F_TABLE  = SR2006_THETA_S33T_F
    MODEL%COUPLING%TABLE%SR2006_THETA_S33T_G_TABLE  = SR2006_THETA_S33T_G
    MODEL%COUPLING%TABLE%SR2006_THETA_S33_A_TABLE   = SR2006_THETA_S33_A
    MODEL%COUPLING%TABLE%SR2006_THETA_S33_B_TABLE   = SR2006_THETA_S33_B
    MODEL%COUPLING%TABLE%SR2006_PSI_ET_A_TABLE      = SR2006_PSI_ET_A
    MODEL%COUPLING%TABLE%SR2006_PSI_ET_B_TABLE      = SR2006_PSI_ET_B
    MODEL%COUPLING%TABLE%SR2006_PSI_ET_C_TABLE      = SR2006_PSI_ET_C
    MODEL%COUPLING%TABLE%SR2006_PSI_ET_D_TABLE      = SR2006_PSI_ET_D
    MODEL%COUPLING%TABLE%SR2006_PSI_ET_E_TABLE      = SR2006_PSI_ET_E
    MODEL%COUPLING%TABLE%SR2006_PSI_ET_F_TABLE      = SR2006_PSI_ET_F
    MODEL%COUPLING%TABLE%SR2006_PSI_ET_G_TABLE      = SR2006_PSI_ET_G
    MODEL%COUPLING%TABLE%SR2006_PSI_E_A_TABLE       = SR2006_PSI_E_A
    MODEL%COUPLING%TABLE%SR2006_PSI_E_B_TABLE       = SR2006_PSI_E_B
    MODEL%COUPLING%TABLE%SR2006_PSI_E_C_TABLE       = SR2006_PSI_E_C
    MODEL%COUPLING%TABLE%SR2006_SMCMAX_A_TABLE      = SR2006_SMCMAX_A
    MODEL%COUPLING%TABLE%SR2006_SMCMAX_B_TABLE      = SR2006_SMCMAX_B

    ! ----------------------
    ! Read general parameters 
    ! ----------------------

    read(noahmp_table_unit, noahmp_general_parameters)

    ! assign values
    model%coupling%table%SLOPE_TABLE(1:NUM_SLOPE) = SLOPE_DATA(1:NUM_SLOPE)
    model%coupling%table%CSOIL_TABLE              = CSOIL_DATA
    model%coupling%table%REFDK_TABLE              = REFDK_DATA
    model%coupling%table%REFKDT_TABLE             = REFKDT_DATA
    model%coupling%table%FRZK_TABLE               = FRZK_DATA
    model%coupling%table%ZBOT_TABLE               = ZBOT_DATA
    model%coupling%table%CZIL_TABLE               = CZIL_DATA

    ! ----------------------
    ! Read soil parameters 
    ! ----------------------

    read(noahmp_table_unit, noahmp_stas_soil_categories)
    if (trim(SLTYPE) == "STAS") then
       read(noahmp_table_unit, noahmp_soil_stas_parameters)
    elseif (trim(SLTYPE) == "STAS_RUC") then
       read(noahmp_table_unit, noahmp_soil_stas_ruc_parameters)
    else
       call ESMF_LogWrite(trim(subname)//': Unrecognized SOILTYPE! '//' It is set to '//trim(SLTYPE), ESMF_LOGMSG_ERROR)
       rcode = ESMF_FAILURE
       close(noahmp_table_unit)
       return
    endif

    ! assign values
    model%coupling%table%SLCATS_TABLE           = SLCATS
    model%coupling%table%BEXP_TABLE  (1:SLCATS) = BB    (1:SLCATS)
    model%coupling%table%SMCDRY_TABLE(1:SLCATS) = DRYSMC(1:SLCATS)
    model%coupling%table%SMCMAX_TABLE(1:SLCATS) = MAXSMC(1:SLCATS)
    model%coupling%table%SMCREF_TABLE(1:SLCATS) = REFSMC(1:SLCATS)
    model%coupling%table%PSISAT_TABLE(1:SLCATS) = SATPSI(1:SLCATS)
    model%coupling%table%DKSAT_TABLE (1:SLCATS) = SATDK (1:SLCATS)
    model%coupling%table%DWSAT_TABLE (1:SLCATS) = SATDW (1:SLCATS)
    model%coupling%table%SMCWLT_TABLE(1:SLCATS) = WLTSMC(1:SLCATS)
    model%coupling%table%QUARTZ_TABLE(1:SLCATS) = QTZ   (1:SLCATS)
    model%coupling%table%BVIC_TABLE  (1:SLCATS) = BVIC  (1:SLCATS)
    model%coupling%table%AXAJ_TABLE  (1:SLCATS) = AXAJ  (1:SLCATS)
    model%coupling%table%BXAJ_TABLE  (1:SLCATS) = BXAJ  (1:SLCATS)
    model%coupling%table%XXAJ_TABLE  (1:SLCATS) = XXAJ  (1:SLCATS)
    model%coupling%table%BDVIC_TABLE (1:SLCATS) = BDVIC (1:SLCATS)
    model%coupling%table%GDVIC_TABLE (1:SLCATS) = GDVIC (1:SLCATS)
    model%coupling%table%BBVIC_TABLE (1:SLCATS) = BBVIC (1:SLCATS)

    ! ----------------------
    ! Close NoahmpTable.TBL
    ! ----------------------

    close(noahmp_table_unit)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine NoahmpReadTable

  !===============================================================================

  subroutine BiochemVarInTransfer(model)

    ! input/output variables
    type(model_type), target, intent(inout) :: model

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//': (BiochemVarInTransfer) '
    !-------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    associate(                                                           &
              I             => model%noahmp%config%domain%GridIndexI   , &
              VegType       => model%noahmp%config%domain%VegType      , &
              CropType      => model%noahmp%config%domain%CropType     , &
              OptCropModel  => model%noahmp%config%nmlist%OptCropModel , &
              OptIrrigation => model%noahmp%config%nmlist%OptIrrigation  &
             )

    !----------------------
    ! State variables
    !----------------------

    !model%noahmp%biochem%state%PlantGrowStage        = model%coupling%table%PGSXY   (I,J)   
    !model%noahmp%biochem%state%LeafMass              = model%coupling%table%LFMASSXY(I,J)
    !model%noahmp%biochem%state%RootMass              = model%coupling%table%RTMASSXY(I,J)
    !model%noahmp%biochem%state%StemMass              = model%coupling%table%STMASSXY(I,J) 
    !model%noahmp%biochem%state%WoodMass              = model%coupling%table%WOODXY  (I,J) 
    !model%noahmp%biochem%state%CarbonMassDeepSoil    = model%coupling%table%STBLCPXY(I,J) 
    !model%noahmp%biochem%state%CarbonMassShallowSoil = model%coupling%table%FASTCPXY(I,J)
    !model%noahmp%biochem%state%GrainMass             = model%coupling%table%GRAINXY (I,J)  
    !model%noahmp%biochem%state%GrowDegreeDay         = model%coupling%table%GDDXY   (I,J)  
    !model%noahmp%biochem%state%NitrogenConcFoliage   = 1.0 

    !----------------------
    ! Transfer parameters
    !----------------------

    model%noahmp%biochem%param%NitrogenConcFoliageMax = model%coupling%table%FOLNMX_TABLE (VegType)
    model%noahmp%biochem%param%QuantumEfficiency25C   = model%coupling%table%QE25_TABLE   (VegType)
    model%noahmp%biochem%param%CarboxylRateMax25C     = model%coupling%table%VCMX25_TABLE (VegType)
    model%noahmp%biochem%param%CarboxylRateMaxQ10     = model%coupling%table%AVCMX_TABLE  (VegType)
    model%noahmp%biochem%param%PhotosynPathC3         = model%coupling%table%C3PSN_TABLE  (VegType)
    model%noahmp%biochem%param%SlopeConductToPhotosyn = model%coupling%table%MP_TABLE     (VegType)
    model%noahmp%biochem%param%RespMaintQ10           = model%coupling%table%ARM_TABLE    (VegType)
    model%noahmp%biochem%param%RespMaintLeaf25C       = model%coupling%table%RMF25_TABLE  (VegType)
    model%noahmp%biochem%param%RespMaintStem25C       = model%coupling%table%RMS25_TABLE  (VegType)
    model%noahmp%biochem%param%RespMaintRoot25C       = model%coupling%table%RMR25_TABLE  (VegType)
    model%noahmp%biochem%param%WoodToRootRatio        = model%coupling%table%WRRAT_TABLE  (VegType)
    model%noahmp%biochem%param%WoodPoolIndex          = model%coupling%table%WDPOOL_TABLE (VegType)
    model%noahmp%biochem%param%TurnoverCoeffLeafVeg   = model%coupling%table%LTOVRC_TABLE (VegType)
    model%noahmp%biochem%param%TemperaureLeafFreeze   = model%coupling%table%TDLEF_TABLE  (VegType)
    model%noahmp%biochem%param%LeafDeathWaterCoeffVeg = model%coupling%table%DILEFW_TABLE (VegType)
    model%noahmp%biochem%param%LeafDeathTempCoeffVeg  = model%coupling%table%DILEFC_TABLE (VegType)
    model%noahmp%biochem%param%GrowthRespFrac         = model%coupling%table%FRAGR_TABLE  (VegType)
    model%noahmp%biochem%param%MicroRespCoeff         = model%coupling%table%MRP_TABLE    (VegType)
    model%noahmp%biochem%param%TemperatureMinPhotosyn = model%coupling%table%TMIN_TABLE   (VegType)
    model%noahmp%biochem%param%LeafAreaPerMass1side   = model%coupling%table%SLA_TABLE    (VegType)
    model%noahmp%biochem%param%StemAreaIndexMin       = model%coupling%table%XSAMIN_TABLE (VegType)
    model%noahmp%biochem%param%WoodAllocFac           = model%coupling%table%BF_TABLE     (VegType)
    model%noahmp%biochem%param%WaterStressCoeff       = model%coupling%table%WSTRC_TABLE  (VegType)
    model%noahmp%biochem%param%LeafAreaIndexMin       = model%coupling%table%LAIMIN_TABLE (VegType)
    model%noahmp%biochem%param%TurnoverCoeffRootVeg   = model%coupling%table%RTOVRC_TABLE (VegType)
    model%noahmp%biochem%param%WoodRespCoeff          = model%coupling%table%RSWOODC_TABLE(VegType)

    !----------------------
    ! Crop model specific parameters
    !----------------------

    if ((OptCropModel > 0) .and. (CropType > 0)) then
       model%noahmp%biochem%param%DatePlanting            = model%coupling%table%PLTDAY_TABLE   (CropType)
       model%noahmp%biochem%param%DateHarvest             = model%coupling%table%HSDAY_TABLE    (CropType)
       model%noahmp%biochem%param%NitrogenConcFoliageMax  = model%coupling%table%FOLNMXI_TABLE  (CropType)
       model%noahmp%biochem%param%QuantumEfficiency25C    = model%coupling%table%QE25I_TABLE    (CropType)
       model%noahmp%biochem%param%CarboxylRateMax25C      = model%coupling%table%VCMX25I_TABLE  (CropType)
       model%noahmp%biochem%param%CarboxylRateMaxQ10      = model%coupling%table%AVCMXI_TABLE   (CropType)
       model%noahmp%biochem%param%PhotosynPathC3          = model%coupling%table%C3PSNI_TABLE   (CropType)
       model%noahmp%biochem%param%SlopeConductToPhotosyn  = model%coupling%table%MPI_TABLE      (CropType)
       model%noahmp%biochem%param%RespMaintQ10            = model%coupling%table%Q10MR_TABLE    (CropType)
       model%noahmp%biochem%param%RespMaintLeaf25C        = model%coupling%table%LFMR25_TABLE   (CropType)
       model%noahmp%biochem%param%RespMaintStem25C        = model%coupling%table%STMR25_TABLE   (CropType)
       model%noahmp%biochem%param%RespMaintRoot25C        = model%coupling%table%RTMR25_TABLE   (CropType)
       model%noahmp%biochem%param%GrowthRespFrac          = model%coupling%table%FRA_GR_TABLE   (CropType)
       model%noahmp%biochem%param%TemperaureLeafFreeze    = model%coupling%table%LEFREEZ_TABLE  (CropType)
       model%noahmp%biochem%param%LeafAreaPerBiomass      = model%coupling%table%BIO2LAI_TABLE  (CropType)
       model%noahmp%biochem%param%TempBaseGrowDegDay      = model%coupling%table%GDDTBASE_TABLE (CropType)
       model%noahmp%biochem%param%TempMaxGrowDegDay       = model%coupling%table%GDDTCUT_TABLE  (CropType)
       model%noahmp%biochem%param%GrowDegDayEmerg         = model%coupling%table%GDDS1_TABLE    (CropType)
       model%noahmp%biochem%param%GrowDegDayInitVeg       = model%coupling%table%GDDS2_TABLE    (CropType)
       model%noahmp%biochem%param%GrowDegDayPostVeg       = model%coupling%table%GDDS3_TABLE    (CropType)
       model%noahmp%biochem%param%GrowDegDayInitReprod    = model%coupling%table%GDDS4_TABLE    (CropType)
       model%noahmp%biochem%param%GrowDegDayMature        = model%coupling%table%GDDS5_TABLE    (CropType)
       model%noahmp%biochem%param%PhotosynRadFrac         = model%coupling%table%I2PAR_TABLE    (CropType)
       model%noahmp%biochem%param%TempMinCarbonAssim      = model%coupling%table%TASSIM0_TABLE  (CropType)
       model%noahmp%biochem%param%TempMaxCarbonAssim      = model%coupling%table%TASSIM1_TABLE  (CropType)
       model%noahmp%biochem%param%TempMaxCarbonAssimMax   = model%coupling%table%TASSIM2_TABLE  (CropType)
       model%noahmp%biochem%param%CarbonAssimRefMax       = model%coupling%table%AREF_TABLE     (CropType)
       model%noahmp%biochem%param%LightExtCoeff           = model%coupling%table%K_TABLE        (CropType)
       model%noahmp%biochem%param%LightUseEfficiency      = model%coupling%table%EPSI_TABLE     (CropType)
       model%noahmp%biochem%param%CarbonAssimReducFac     = model%coupling%table%PSNRF_TABLE    (CropType)
       model%noahmp%biochem%param%RespMaintGrain25C       = model%coupling%table%GRAINMR25_TABLE(CropType)
       model%noahmp%biochem%param%LeafDeathTempCoeffCrop  = model%coupling%table%DILE_FC_TABLE  (CropType,:)
       model%noahmp%biochem%param%LeafDeathWaterCoeffCrop = model%coupling%table%DILE_FW_TABLE  (CropType,:)
       model%noahmp%biochem%param%CarbohydrLeafToGrain    = model%coupling%table%LFCT_TABLE     (CropType,:)
       model%noahmp%biochem%param%CarbohydrStemToGrain    = model%coupling%table%STCT_TABLE     (CropType,:)
       model%noahmp%biochem%param%CarbohydrRootToGrain    = model%coupling%table%RTCT_TABLE     (CropType,:)
       model%noahmp%biochem%param%CarbohydrFracToLeaf     = model%coupling%table%LFPT_TABLE     (CropType,:)
       model%noahmp%biochem%param%CarbohydrFracToStem     = model%coupling%table%STPT_TABLE     (CropType,:)
       model%noahmp%biochem%param%CarbohydrFracToRoot     = model%coupling%table%RTPT_TABLE     (CropType,:)
       model%noahmp%biochem%param%CarbohydrFracToGrain    = model%coupling%table%GRAINPT_TABLE  (CropType,:)
       model%noahmp%biochem%param%TurnoverCoeffLeafCrop   = model%coupling%table%LF_OVRC_TABLE  (CropType,:)
       model%noahmp%biochem%param%TurnoverCoeffStemCrop   = model%coupling%table%ST_OVRC_TABLE  (CropType,:)
       model%noahmp%biochem%param%TurnoverCoeffRootCrop   = model%coupling%table%RT_OVRC_TABLE  (CropType,:)

       !if (OptCropModel == 1) then
       !   model%noahmp%biochem%param%DatePlanting         = model%coupling%table%PLANTING(I,J)
       !   model%noahmp%biochem%param%DateHarvest          = model%coupling%table%HARVEST(I,J)
       !   model%noahmp%biochem%param%GrowDegDayEmerg      = model%coupling%table%SEASON_GDD(I,J) / 1770.0 * &
       !                                                     model%noahmp%biochem%param%GrowDegDayEmerg
       !   model%noahmp%biochem%param%GrowDegDayInitVeg    = model%coupling%table%SEASON_GDD(I,J) / 1770.0 * &
       !                                                     model%noahmp%biochem%param%GrowDegDayInitVeg
       !   model%noahmp%biochem%param%GrowDegDayPostVeg    = model%coupling%table%SEASON_GDD(I,J) / 1770.0 * &
       !                                                     model%noahmp%biochem%param%GrowDegDayPostVeg
       !   model%noahmp%biochem%param%GrowDegDayInitReprod = model%coupling%table%SEASON_GDD(I,J) / 1770.0 * &
       !                                                     model%noahmp%biochem%param%GrowDegDayInitReprod
       !   model%noahmp%biochem%param%GrowDegDayMature     = model%coupling%table%SEASON_GDD(I,J) / 1770.0 * &
       !                                                     model%noahmp%biochem%param%GrowDegDayMature
       ! endif
    end if

    !if (OptIrrigation == 2) then
    !   model%noahmp%biochem%param%DatePlanting = NoahmpIO%PLANTING(I,J)
    !   model%noahmp%biochem%param%DateHarvest  = NoahmpIO%HARVEST (I,J)
    !endif

    end associate

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine BiochemVarInTransfer

  !===============================================================================

  subroutine ConfigVarInTransfer(model)

    ! input/output variables
    type(model_type), target, intent(inout) :: model

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//': (ConfigVarInTransfer) '
    !-------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ConfigVarInTransfer

  !===============================================================================

  subroutine EnergyVarInTransfer(model)

    ! input/output variables
    type(model_type), target, intent(inout) :: model

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//': (EnergyVarInTransfer) '
    !-------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine EnergyVarInTransfer

  !===============================================================================

  subroutine ForcingVarInTransfer(model)

    ! input/output variables
    type(model_type), target, intent(inout) :: model

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//': (ForcingVarInTransfer) '
    !-------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    associate(I => model%noahmp%config%domain%GridIndexI)

    !----------------------
    ! Transfer variables
    !---------------------- 

    ! TODO: check conversion from mixing ratio to specific humidity
    model%noahmp%forcing%SpecHumidityRefHeight   = model%forcing%SpecHumidityRefHeight  (I) 

    model%noahmp%forcing%TemperatureAirRefHeight = model%forcing%TemperatureAirRefHeight(I) 
    model%noahmp%forcing%WindEastwardRefHeight   = model%forcing%WindEastwardRefHeight  (I) 
    model%noahmp%forcing%WindNorthwardRefHeight  = model%forcing%WindNorthwardRefHeight (I) 
    model%noahmp%forcing%TemperatureSoilBottom   = model%forcing%TemperatureSoilBottom  (I)
    model%noahmp%forcing%RadSwDownRefHeight      = model%forcing%RadSwDownRefHeight     (I) 
    model%noahmp%forcing%RadLwDownRefHeight      = model%forcing%RadLwDownRefHeight     (I) 
    model%noahmp%forcing%PressureAirRefHeight    = model%forcing%PressureAirRefHeight   (I) 
    model%noahmp%forcing%PressureAirSurface      = model%forcing%PressureAirSurface     (I) 

    model%noahmp%forcing%PrecipConvRefHeight     = model%forcing%PrecipConvRefHeight    (I) / model%coupling%MainTimeStep 
    model%noahmp%forcing%PrecipNonConvRefHeight  = model%forcing%PrecipNonConvRefHeight (I) / model%coupling%MainTimeStep
    model%noahmp%forcing%PrecipShConvRefHeight   = model%forcing%PrecipShConvRefHeight  (I) / model%coupling%MainTimeStep
    model%noahmp%forcing%PrecipSnowRefHeight     = model%forcing%PrecipSnowRefHeight    (I) / model%coupling%MainTimeStep 
    model%noahmp%forcing%PrecipGraupelRefHeight  = model%forcing%PrecipGraupelRefHeight (I) / model%coupling%MainTimeStep
    model%noahmp%forcing%PrecipHailRefHeight     = model%forcing%PrecipHailRefHeight    (I) / model%coupling%MainTimeStep

    ! TODO: treat other precipitation (e.g. fog) contained in total precipitation

    end associate

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ForcingVarInTransfer

  !===============================================================================

  subroutine WaterVarInTransfer(model)

    ! input/output variables
    type(model_type), target, intent(inout) :: model

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//': (WaterVarInTransfer) '
    !-------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    associate(                                                          &
              I             => model%noahmp%config%domain%GridIndexI  , &
              nsoil         => model%noahmp%config%domain%NumSoilLayer, &
              nsnow         => model%noahmp%config%domain%NumSnowLayerMax
             )

    ! water state variables
    model%noahmp%water%state%CanopyLiqWater          = model%coupling%CanopyLiqWater          (I) 
    model%noahmp%water%state%CanopyIce               = model%coupling%CanopyIce               (I) 
    model%noahmp%water%state%CanopyWetFrac           = model%coupling%CanopyWetFrac           (I) 
    model%noahmp%water%state%SnowWaterEquiv          = model%coupling%SnowWaterEquiv          (I) 
    model%noahmp%water%state%SnowWaterEquivPrev      = model%coupling%SnowWaterEquivPrev      (I) 
    model%noahmp%water%state%SnowDepth               = model%coupling%SnowDepth               (I) 
    model%noahmp%water%state%IrrigationFracFlood     = model%coupling%IrrigationFracFlood     (I) 
    model%noahmp%water%state%IrrigationAmtFlood      = model%coupling%IrrigationAmtFlood      (I) 
    model%noahmp%water%state%IrrigationFracMicro     = model%coupling%IrrigationFracMicro     (I) 
    model%noahmp%water%state%IrrigationAmtMicro      = model%coupling%IrrigationAmtMicro      (I) 
    model%noahmp%water%state%IrrigationFracSprinkler = model%coupling%IrrigationFracSprinkler (I) 
    model%noahmp%water%state%IrrigationAmtSprinkler  = model%coupling%IrrigationAmtSprinkler  (I) 
    model%noahmp%water%state%WaterTableDepth         = model%coupling%WaterTableDepth         (I) 
    model%noahmp%water%state%SoilMoistureToWT        = model%coupling%SoilMoistureToWT        (I) 
    model%noahmp%water%state%TileDrainFrac           = model%coupling%TileDrainFrac           (I) 
    model%noahmp%water%state%WaterStorageAquifer     = model%coupling%WaterStorageAquifer     (I) 
    model%noahmp%water%state%WaterStorageSoilAqf     = model%coupling%WaterStorageSoilAqf     (I) 
    model%noahmp%water%state%WaterStorageLake        = model%coupling%WaterStorageLake        (I) 
    model%noahmp%water%state%IrrigationFracGrid      = model%coupling%IrrigationFracGrid      (I) 
    model%noahmp%water%state%IrrigationCntSprinkler  = model%coupling%IrrigationCntSprinkler  (I) 
    model%noahmp%water%state%IrrigationCntMicro      = model%coupling%IrrigationCntMicro      (I) 
    model%noahmp%water%state%IrrigationCntFlood      = model%coupling%IrrigationCntFlood      (I) 
    model%noahmp%water%state%SnowIce(-nsnow+1:0)     = model%coupling%SnowIce                 (I,-nsnow+1:0) 
    model% model%coupling%SnowLiqWater            (I) 
    model% model%coupling%SoilLiqWater            (I) 
    model% model%coupling%SoilMoisture            (I) 
    model% model%coupling%SoilMoistureEqui        (I) 
    model% model%coupling%RechargeGwDeepWT        (I) 
    model% model%coupling%RechargeGwShallowWT     (I) 
    model% model%coupling%EvapSoilSfcLiqAcc       (I) 
    model% model%coupling%SoilSfcInflowAcc        (I) 
    model% model%coupling%SfcWaterTotChgAcc       (I) 
    model% model%coupling%PrecipTotAcc            (I) 
    model% model%coupling%EvapCanopyNetAcc        (I) 
    model% model%coupling%TranspirationAcc        (I) 
    model% model%coupling%EvapGroundNetAcc        (I) 
    model% model%coupling%TranspWatLossSoilAcc    (I) 
    model% model%coupling%DrainSoilLayerInd       (I) 
    model% model%coupling%CanopyLiqHoldCap        (I) 
    model% model%coupling%SnowCompactBurdenFac    (I) 
    model% model%coupling%SnowCompactAgingFac1    (I) 
    model% model%coupling%SnowCompactAgingFac2    (I) 
    model% model%coupling%SnowCompactAgingFac3    (I) 
    model% model%coupling%SnowCompactAgingMax     (I) 
    model% model%coupling%SnowViscosityCoeff      (I) 
    model% model%coupling%SnowLiqFracMax          (I) 
    model% model%coupling%SnowLiqHoldCap          (I) 
    model% model%coupling%SnowLiqReleaseFac       (I) 
    model% model%coupling%IrriFloodRateFac        (I) 
    model% model%coupling%IrriMicroRate           (I) 
    model% model%coupling%SoilConductivityRef     (I) 
    model% model%coupling%SoilInfilFacRef         (I) 
    model% model%coupling%GroundFrzCoeff          (I) 
    model% model%coupling%GridTopoIndex           (I) 
    model% model%coupling%SoilSfcSatFracMax       (I) 
    model% model%coupling%SpecYieldGw             (I) 
    model% model%coupling%MicroPoreContent        (I) 
    model% model%coupling%WaterStorageLakeMax     (I) 
    model% model%coupling%SnoWatEqvMaxGlacier     (I) 
    model% model%coupling%IrriStopDayBfHarvest    (I) 
    model% model%coupling%IrriTriggerLaiMin       (I) 
    model% model%coupling%SoilWatDeficitAllow     (I) 
    model% model%coupling%IrriFloodLossFrac       (I) 
    model% model%coupling%IrriSprinklerRate       (I) 
    model% model%coupling%IrriFracThreshold       (I) 
    model% model%coupling%IrriStopPrecipThr       (I) 
    model% model%coupling%SnowfallDensityMax      (I) 
    model% model%coupling%SnowMassFullCoverOld    (I) 
    model% model%coupling%SoilMatPotentialWilt    (I) 
    model% model%coupling%SnowMeltFac             (I) 
    model% model%coupling%SnowCoverFac            (I) 
    model% model%coupling%InfilFacVic             (I) 
    model% model%coupling%TensionWatDistrInfl     (I) 
    model% model%coupling%TensionWatDistrShp      (I) 
    model% model%coupling%FreeWatDistrShp         (I) 
    model% model%coupling%InfilHeteroDynVic       (I) 
    model% model%coupling%InfilCapillaryDynVic    (I) 
    model% model%coupling%InfilFacDynVic          (I) 
    model% model%coupling%TileDrainCoeffSp        (I) 
    model% model%coupling%TileDrainTubeDepth      (I) 
    model% model%coupling%DrainFacSoilWat         (I) 
    model% model%coupling%TileDrainCoeff          (I) 
    model% model%coupling%DrainDepthToImperv      (I) 
    model% model%coupling%LateralWatCondFac       (I) 
    model% model%coupling%TileDrainDepth          (I) 
    model% model%coupling%DrainTubeDist           (I) 
    model% model%coupling%DrainTubeRadius         (I) 
    model% model%coupling%DrainWatDepToImperv     (I) 
    model% model%coupling%NumSoilLayerRoot        (I) 
    model% model%coupling%SoilDrainSlope          (I) 

    end associate

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine WaterVarInTransfer

  !===============================================================================

  subroutine CalcCosZ(currTime, lat, long, cosz, julian, rc)
    implicit none

    ! input/output variables
    type(ESMF_Time), intent(in) :: currTime
    real(r8), intent(in)        :: lat(:)
    real(r8), intent(in)        :: long(:)
    real(r8), intent(inout)     :: cosz(:)
    real(r8), intent(inout)     :: julian
    integer, intent(inout)      :: rc

    ! local variables
    real(r8)                    :: sec_since
    real(r8)                    :: obecl, sinob, sxlong, arg, tloctim, hrang, declin
    integer                     :: iloc, iyear, imonth, iday, ihour, iminute, isecond
    real(r8), parameter         :: degrad = 3.14159265/180.0
    real(r8), parameter         :: dpd    = 360.0/365.0
    character(len=*), parameter :: subname=trim(modName)//':(CalcCosZ) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------
    ! query current time
    !---------------------- 

    call ESMF_TimeGet(currTime, yy=iyear, mm=imonth, dd=iday, h=ihour, m=iminute, s=isecond, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeSet(dummTime, yy=iyear, mm=1, dd=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet(currTime-dummTime, s_r8=sec_since, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! obliquity = 23.5 degree.
    !----------------------

    obecl = 23.5*degrad
    sinob = sin(obecl)

    !----------------------
    ! calculate longitude of the sun from vernal equinox
    !----------------------

    julian = sec_since/86400.0
    if (julian >= 80.0) sxlong = dpd*(julian-80.0)*degrad
    if (julian <  80.0) sxlong = dpd*(julian+285.0)*degrad
    arg = sinob*sin(sxlong)
    declin = asin(arg)

    do iloc = 1, size(lat)
       tloctim = real(ihour)+real(iminute)/60.0+real(isecond)/3600.0+ long(iloc)/15.0 ! local time in hours
       tloctim = mod(tloctim+24.0, 24.0)
       hrang = 15.*(tloctim-12.0)*degrad
       cosz(iloc) = sin(lat(iloc)*degrad)*sin(declin)+cos(lat(iloc)*degrad)*cos(declin)*cos(hrang)
    end do

  end subroutine CalcCosZ

end module lnd_comp_driver
