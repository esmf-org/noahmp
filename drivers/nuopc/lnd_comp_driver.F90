module lnd_comp_driver

  ! This file contains the NoahMP land surface model driver

  use ESMF, only: operator(+), operator(-), operator(/)
  use ESMF, only: ESMF_GridComp, ESMF_GridCompGet
  use ESMF, only: ESMF_VM, ESMF_VMGet
  use ESMF, only: ESMF_Clock, ESMF_Time, ESMF_TimeSet
  use ESMF, only: ESMF_TimeInterval, ESMF_TimeGet, ESMF_ClockGet
  use ESMF, only: ESMF_TimeIntervalGet, ESMF_KIND_R8
  use ESMF, only: ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO

  use lnd_comp_types, only: model_type
  use lnd_comp_shr  , only: ChkErr
  use lnd_comp_kind , only: cl => shr_kind_cl
  use lnd_comp_types, only: histflds, restflds

  use lnd_comp_io, only: SetupWriteFields
  use lnd_comp_io, only: WriteFile

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
    type(ESMF_GridComp), intent(in)   :: gcomp
    type(model_type)  , intent(inout) :: model
    integer            , intent(out)  :: rc

    ! local variables
    integer            :: localPet
    type(ESMF_VM)      :: vm
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
    ! Allocate and initialize noahmp model data
    ! ----------------------

    !call InitNoahMP(model, model%domain%begl, model%domain%endl)


    ! get coupling time step 
    !call ESMF_TimeIntervalGet(timeStep, s_r8=noahmp%static%delt, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! history
    !call SetupWriteFields(model, 'hist', rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! restart
    !call SetupWriteFields(model, 'rest', rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return


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
    real(ESMF_KIND_R8)      :: now_time
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
    ! Query component
    !----------------------

    call ESMF_GridCompGet(gcomp, vm=vm, clock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (firstTime) then
    ! history
    call SetupWriteFields(model, 'hist', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    firstTime = .false.
    end if

    ! restart
    !call SetupWriteFields(model, 'rest', rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------
    ! Loop over grid points and run NoahMP model
    ! ----------------------

    !do i = model%domain%begl, model%domain%endl

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

    !end do

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

    model%noahmp%forcing%PrecipConvRefHeight     = model%forcing%PrecipConvRefHeight    (indx) / model%coupling%dt 
    model%noahmp%forcing%PrecipNonConvRefHeight  = model%forcing%PrecipNonConvRefHeight (indx) / model%coupling%dt 
    model%noahmp%forcing%PrecipShConvRefHeight   = model%forcing%PrecipShConvRefHeight  (indx) / model%coupling%dt 
    model%noahmp%forcing%PrecipSnowRefHeight     = model%forcing%PrecipSnowRefHeight    (indx) / model%coupling%dt 
    model%noahmp%forcing%PrecipGraupelRefHeight  = model%forcing%PrecipGraupelRefHeight (indx) / model%coupling%dt 
    model%noahmp%forcing%PrecipHailRefHeight     = model%forcing%PrecipHailRefHeight    (indx) / model%coupling%dt 

    ! TODO: treat other precipitation (e.g. fog) contained in total precipitation

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ForcingVarInTransfer

end module lnd_comp_driver
