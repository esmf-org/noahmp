module lnd_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for NoahMP
  !----------------------------------------------------------------------------

  use ESMF             , only: operator(+), operator(-), operator(==)
  use ESMF             , only: ESMF_GridCompSetEntryPoint, ESMF_GridComp
  use ESMF             , only: ESMF_GridCompGet, ESMF_VMGet, ESMF_Mesh
  use ESMF             , only: ESMF_MethodRemove, ESMF_LogWrite, ESMF_LOGMSG_INFO
  use ESMF             , only: ESMF_METHOD_INITIALIZE, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF             , only: ESMF_State, ESMF_Clock, ESMF_Time, ESMF_VM
  use ESMF             , only: ESMF_Array, ESMF_ArrayRead, ESMF_ArrayGet, ESMF_ArrayDestroy
  use ESMF             , only: ESMF_TimeInterval, ESMF_Alarm, ESMF_ClockGet
  use ESMF             , only: ESMF_ClockGetAlarmList, ESMF_Clock, ESMF_Time
  use ESMF             , only: ESMF_ClockSet, ESMF_TimeInterval, ESMF_ALARMLIST_ALL 
  use ESMF             , only: ESMF_AlarmSet, ESMF_ClockAdvance
  use ESMF             , only: ESMF_TimeGet, ESMF_TimeInterval
  use ESMF             , only: ESMF_GEOMTYPE_GRID, ESMF_GEOMTYPE_MESH
  use ESMF             , only: ESMF_MeshCreate, ESMF_Grid, ESMF_GeomType_Flag
  use ESMF             , only: ESMF_VMGet, ESMF_VMGetCurrent
  use NUOPC            , only: NUOPC_CompDerive, NUOPC_CompAttributeGet
  use NUOPC            , only: NUOPC_CompFilterPhaseMap, NUOPC_CompSetEntryPoint 
  use NUOPC            , only: NUOPC_CompSpecialize
  use NUOPC_Model      , only: NUOPC_ModelGet
  use NUOPC_Model      , only: model_routine_SS => SetServices
  use NUOPC_Model      , only: SetVM
  use NUOPC_Model      , only: model_label_Advance => label_Advance
  use NUOPC_Model      , only: model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model      , only: model_label_Finalize => label_Finalize

  use lnd_comp_domain  , only: SetDomain 
  use lnd_comp_shr     , only: ChkErr
  use lnd_comp_shr     , only: ReadNamelist
  use lnd_comp_types   , only: model_type
  use lnd_comp_import_export, only: AdvertiseFields, RealizeFields, ImportFields

  use lnd_comp_driver  , only: NoahmpDriverInit
  use lnd_comp_driver  , only: NoahmpDriverMain
  use lnd_comp_driver  , only: NoahmpDriverFinalize

  implicit none
  private ! except

  ! Module public routines
  public  :: SetServices, SetVM

  ! Module private routines
  private :: InitializeP0        ! Phase zero of initialization
  private :: InitializeAdvertise ! Advertise the fields that can be passed
  private :: InitializeRealize   ! Realize the list of fields that will be exchanged
  private :: ModelAdvance        ! Advance the model
  private :: ModelSetRunClock    ! Set the run clock

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(model_type)       :: model
  character(*),parameter :: modName =  "(lnd_comp_nuopc)"

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    ! Setup the pointers to the function calls for the different models phases (initialize, run, finalize)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! Phase zero initialization
    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! Advertise the fields that can be exchanged
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Advertise fields
    ! ---------------------

    call AdvertiseFields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise

  !===============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! Realize the list of fields that will be exchanged
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Read coupling specific namelist options
    ! ---------------------

    call ReadNamelist(gcomp, model, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Create grid and convert it to mesh 
    ! ---------------------

    call SetDomain(gcomp, model, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Initialize NoahMP
    !----------------------

    call NoahmpDriverInit(gcomp, model, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return     

    ! ---------------------
    ! Realize the actively coupled fields
    ! ---------------------

    call RealizeFields(importState, exportState, model%domain%mesh, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize

  !===============================================================================
  subroutine ModelAdvance(gcomp, rc)

    !------------------------
    ! Run NoahMP
    !------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !-----------------------
    ! Query the Component for its clock, importState and exportState
    !-----------------------

    !call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------
    ! import state
    !-----------------------

    call ImportFields(gcomp, model, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! diagnostics
    !----------------------

    !if (dbug > 1) then
    !   call state_diagnose(importState, subname//': ImportState ',rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !endif

    !----------------------
    ! Run NoahMP
    !----------------------

    call NoahmpDriverMain(gcomp, model, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! export state
    !----------------------

    !if (noahmp%nmlist%has_export) then
    !   call export_fields(gcomp, noahmp, rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !end if

    !----------------------
    ! diagnostics
    !----------------------

    !if (dbug > 1) then
    !   call state_diagnose(exportState, subname//': ExportState ',rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !endif

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelAdvance

  !===============================================================================
  subroutine ModelSetRunClock(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: mclock, dclock
    type(ESMF_Time)         :: mcurrtime, dcurrtime
    type(ESMF_Time)         :: mstoptime
    type(ESMF_TimeInterval) :: mtimestep, dtimestep
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !--------------------------------
    ! query the Component for its clocks
    !--------------------------------

    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! call model finalize routine
    !call e(gcomp, noahmp, rc)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

end module lnd_comp_nuopc
