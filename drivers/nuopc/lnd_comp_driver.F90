module lnd_comp_driver

  ! This file contains the NoahMP land surface model driver

  use ESMF, only: operator(+), operator(-), operator(/)
  use ESMF, only: ESMF_GridComp, ESMF_GridCompGet
  use ESMF, only: ESMF_VM, ESMF_VMGet
  use ESMF, only: ESMF_Clock, ESMF_Time, ESMF_TimeSet
  use ESMF, only: ESMF_TimeInterval, ESMF_TimeGet, ESMF_ClockGet
  use ESMF, only: ESMF_TimeIntervalGet, ESMF_KIND_R8
  use ESMF, only: ESMF_TimeIsLeapYear
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

    !call NoahmpReadTable(model)

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
    character(len=cl)       :: filename
    type(ESMF_VM)           :: vm
    type(ESMF_Clock)        :: clock
    type(ESMF_Time)         :: startTime
    type(ESMF_Time)         :: currTime
    type(ESMF_TimeInterval) :: timeStep
    logical, save           :: firstTime = .true.
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

    ! use coupling time step and internal time step of model
    call ESMF_TimeIntervalGet(timeStep, s_r8=model%coupling%MainTimeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

    ! number of crop growth stages
    model%noahmp%config%domain%NumCropGrowStage = 8

    ! model grid spacing size
    ! TODO: maybe this could be calculated from grid area?
    model%noahmp%config%domain%GridSize = -9999.0

    ! thickness of lowest atmosphere layer
    ! TODO: namelist option? or retrieved from active atmosphere
    model%noahmp%config%domain%ThicknessAtmosBotLayer = -9999.0

!    if (firstTime) then
!    ! history
!    call SetupWriteFields(model, 'hist', rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    firstTime = .false.
!    end if

    ! restart
    !call SetupWriteFields(model, 'rest', rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------
    ! Loop over grid points and run NoahMP model
    ! ----------------------

    do i = model%domain%begl, model%domain%endl
       ! ----------------------
       ! Set domain specific model parameters
       ! ----------------------

       model%noahmp%config%domain%GridIndexI            = i
       model%noahmp%config%domain%GridIndexJ            = -9999
       model%noahmp%config%domain%Latitude              = model%domain%latm(i)
       model%noahmp%config%domain%CosSolarZenithAngle   = model%coupling%CosSolarZenithAngle(i)
       model%noahmp%config%domain%SoilType(:)           = model%coupling%SoilType(i)
       model%noahmp%config%domain%VegType               = model%coupling%VegType(i)
       model%noahmp%config%domain%RunoffSlopeType       = model%coupling%RunoffSlopeType(i)

       model%noahmp%energy%param%VegFracAnnMax          = model%coupling%VegFracAnnMax(i)
       model%noahmp%energy%state%TemperatureSoilSnow(:) = model%coupling%TemperatureSoilSnow(i,:)
       
       model%noahmp%water%state%SoilMoisture(:)         = model%coupling%SoilMoisture(i,:)
       model%noahmp%water%state%SnowWaterEquiv          = model%coupling%SnowWaterEquiv(i)
       model%noahmp%water%state%SnowDepth               = model%coupling%SnowDepth(i)

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

  !subroutine NoahmpReadTable(model)

  !  ! input/output variables
  !  type(model_type), target, intent(inout) :: model

  !  ! local variables
  !  integer, parameter :: NumVegType = 27 ! number of vegetation types
  !  character(len=*), parameter :: subname=trim(modName)//':(NoahmpReadTable) '
  !  !-------------------------------------------------------------------------

  !  call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

  !  ! ----------------------
  !  ! Associate variables
  !  ! ----------------------

  !  ! vegetation variables
  !  associate(model%noahmp%water%param%CanopyLiqHoldCap       => CH2OP_TABLE, &
  !            model%noahmp%energy%param%LeafDimLength         => DLEAF_TABLE, &
  !            model%noahmp%energy%param%RoughLenMomVeg        => Z0MVT_TABLE, &
  !            model%noahmp%energy%param%HeightCanopyTop       => HVT_TABLE, &
  !            model%noahmp%energy%param%HeightCanopyBot       => HVB_TABLE, &
  !            model%noahmp%energy%param%TreeDensity           => DEN_TABLE, &
  !            model%noahmp%energy%param%TreeCrownRadius       => RC_TABLE, &
  !            model%noahmp%water%param%SnowMeltFac            => MFSNO_TABLE, &
  !            model%noahmp%water%param%SnowCoverFac           => SCFFAC_TABLE, &
  !            model%noahmp%energy%param%HeatCapacCanFac       => CBIOM_TABLE, &
  !            model%noahmp%energy%param%StemAreaIndexMon      => SAIM_TABLE, &
  !            model%noahmp%energy%param%LeafAreaIndexMon      => LAIM_TABLE, &
  !            model%noahmp%biochem%param%LeafAreaPerMass1side => SLA_TABLE, &
  !            model% => DILEFC_TABLE, &
  !            model% => DILEFW_TABLE, &
  !            model% => FRAGR_TABLE, &
  !            model% => LTOVRC_TABLE, &
  !            model% => C3PSN_TABLE, &
  !            model% => KC25_TABLE, &
  !            model% => AKC_TABLE, &
  !            model% => KO25_TABLE, &
  !            model% => AKO_TABLE, &
  !            model% => VCMX25_TABLE, &
  !            model% => AVCMX_TABLE, &
  !            model% => BP_TABLE, &
  !            model% => MP_TABLE, &
  !            model% => QE25_TABLE, &
  !            model% => AQE_TABLE, &
  !            model% => RMF25_TABLE, &
  !            model% => RMS25_TABLE, &
  !            model% => RMR25_TABLE, &
  !            model% => ARM_TABLE, &
  !            model% => FOLNMX_TABLE, &
  !            model% => TMIN_TABLE, &
  !            model% => XL_TABLE, &
  !            model% => RHOL_TABLE, &
  !            model% => RHOS_TABLE, &
  !            model% => TAUL_TABLE, &
  !            model% => TAUS_TABLE, &
  !            model% => MRP_TABLE, &
  !            model% => CWPVT_TABLE, &
  !            model% => WRRAT_TABLE, &
  !            model% => WDPOOL_TABLE, &
  !            model% => TDLEF_TABLE, &
  !            model% => NROOT_TABLE, &
  !            model% => RGL_TABLE, &
  !            model% => RS_TABLE, &
  !            model% => HS_TABLE, &
  !            model% => TOPT_TABLE, &
  !            model% => RSMAX_TABLE, &
  !            model% => RTOVRC_TABLE, &
  !            model% => RSWOODC_TABLE, &
  !            model% => BF_TABLE, &
  !            model% => WSTRC_TABLE, &
  !            model% => LAIMIN_TABLE, &
  !            model% => XSAMIN_TABLE)

  !  ! ----------------------
  !  ! Allocate table variables
  !  ! ----------------------

  !  if (.not. allocated(CH2OP_TABLE )) allocate(CH2OP_TABLE(NumVegType))
  !  if (.not. allocated(DLEAF_TABLE )) allocate(DLEAF_TABLE(NumVegType))
  !  if (.not. allocated(Z0MVT_TABLE )) allocate(Z0MVT_TABLE
  !  if (.not. allocated(HVT_TABLE   )) allocate(HVT_TABLE
  !  if (.not. allocated(HVB_TABLE   )) allocate(HVB_TABLE
  !  if (.not. allocated(DEN_TABLE   )) allocate(DEN_TABLE
  !  if (.not. allocated(RC_TABLE    )) allocate(RC_TABLE
  !  if (.not. allocated(MFSNO_TABLE )) allocate(MFSNO_TABLE
  !  if (.not. allocated(SCFFAC_TABLE)) allocate(SCFFAC_TABLE
  !  if (.not. allocated(CBIOM_TABLE )) allocate(CBIOM_TABLE
  !  if (.not. allocated(SAIM_TABLE  )) allocate(SAIM_TABLE
  !  if (.not. allocated(LAIM_TABLE  )) allocate(LAIM_TABLE
  !  if (.not. allocated(SLA_TABLE   )) allocate(SLA_TABLE




  !  call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  !end subroutine NoahmpReadTable

  !===============================================================================

  subroutine ForcingVarInTransfer(model, indx)

    ! input/output variables
    type(model_type), target, intent(inout) :: model
    integer :: indx

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//': (ForcingVarInTransfer) '
    !-------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! Transfer variables
    !---------------------- 

    ! TODO: check conversion from mixing ratio to specific humidity
    model%noahmp%forcing%SpecHumidityRefHeight   = model%forcing%SpecHumidityRefHeight  (indx) 

    model%noahmp%forcing%TemperatureAirRefHeight = model%forcing%TemperatureAirRefHeight(indx) 
    model%noahmp%forcing%WindEastwardRefHeight   = model%forcing%WindEastwardRefHeight  (indx) 
    model%noahmp%forcing%WindNorthwardRefHeight  = model%forcing%WindNorthwardRefHeight (indx) 
    model%noahmp%forcing%TemperatureSoilBottom   = model%forcing%TemperatureSoilBottom  (indx)
    model%noahmp%forcing%RadSwDownRefHeight      = model%forcing%RadSwDownRefHeight     (indx) 
    model%noahmp%forcing%RadLwDownRefHeight      = model%forcing%RadLwDownRefHeight     (indx) 
    model%noahmp%forcing%PressureAirRefHeight    = model%forcing%PressureAirRefHeight   (indx) 
    model%noahmp%forcing%PressureAirSurface      = model%forcing%PressureAirSurface     (indx) 

    model%noahmp%forcing%PrecipConvRefHeight     = model%forcing%PrecipConvRefHeight    (indx) / model%coupling%MainTimeStep 
    model%noahmp%forcing%PrecipNonConvRefHeight  = model%forcing%PrecipNonConvRefHeight (indx) / model%coupling%MainTimeStep
    model%noahmp%forcing%PrecipShConvRefHeight   = model%forcing%PrecipShConvRefHeight  (indx) / model%coupling%MainTimeStep
    model%noahmp%forcing%PrecipSnowRefHeight     = model%forcing%PrecipSnowRefHeight    (indx) / model%coupling%MainTimeStep 
    model%noahmp%forcing%PrecipGraupelRefHeight  = model%forcing%PrecipGraupelRefHeight (indx) / model%coupling%MainTimeStep
    model%noahmp%forcing%PrecipHailRefHeight     = model%forcing%PrecipHailRefHeight    (indx) / model%coupling%MainTimeStep

    ! TODO: treat other precipitation (e.g. fog) contained in total precipitation

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ForcingVarInTransfer

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
