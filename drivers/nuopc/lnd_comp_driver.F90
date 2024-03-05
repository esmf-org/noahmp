module lnd_comp_driver

  ! This file contains the NoahMP land surface model driver

  use ESMF, only: ESMF_GridComp, ESMF_GridCompGet
  use ESMF, only: ESMF_VM, ESMF_VMGet
  use ESMF, only: ESMF_Time
  use ESMF, only: ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO

  use lnd_comp_types, only: model_type
  use lnd_comp_shr  , only: ChkErr

  implicit none
  private

  public :: drv_init!, drv_run, drv_finalize

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

  subroutine drv_init(gcomp, model, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)   :: gcomp
    type(model_type)  , intent(inout) :: model
    integer            , intent(out)  :: rc

    ! local variables
    integer            :: localPet
    type(ESMF_VM)      :: vm
    character(len=*), parameter :: subname=trim(modName)//':(drv_init) '
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

    ! ----------------------
    ! Allocate and initialize noahmp model data
    ! ----------------------

    call InitNoahMP(model, model%domain%begl, model%domain%endl)



  end subroutine drv_init

  !===============================================================================

  subroutine InitNoahMP(model, begl, endl)

    use Machine

    ! input/output variables
    type(model_type), target, intent(inout) :: model
    integer, intent(in) :: begl
    integer, intent(in) :: endl

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(InitNoahMP) '
    !-------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ----------------------
    ! Allocate forcing data
    ! ----------------------

    allocate(model%forcing%q1(begl:endl))
    allocate(model%forcing%t1(begl:endl))
    allocate(model%forcing%u1(begl:endl))
    allocate(model%forcing%v1(begl:endl))
    allocate(model%forcing%dswsfc(begl:endl))
    allocate(model%forcing%dlwflx(begl:endl))
    allocate(model%forcing%pbot(begl:endl))
    allocate(model%forcing%ps(begl:endl))
    allocate(model%forcing%tprcpc(begl:endl))
    allocate(model%forcing%tprcpl(begl:endl))
    allocate(model%forcing%snowc(begl:endl))
    allocate(model%forcing%snowl(begl:endl))
    allocate(model%forcing%graupel(begl:endl))
    allocate(model%forcing%hail(begl:endl))
    !?allocate(model%forcing%TemperatureSoilBottom  (begl:endl))

    ! ----------------------
    ! Initialize forcing data
    ! ----------------------

    model%forcing%q1 = undefined_real 
    model%forcing%t1 = undefined_real
    model%forcing%u1 = undefined_real
    model%forcing%v1 = undefined_real
    model%forcing%dswsfc = undefined_real
    model%forcing%dlwflx = undefined_real
    model%forcing%pbot = undefined_real
    model%forcing%ps = undefined_real
    model%forcing%tprcpc = undefined_real
    model%forcing%tprcpl = undefined_real
    model%forcing%snowc = undefined_real
    model%forcing%snowl = undefined_real
    model%forcing%graupel = undefined_real
    model%forcing%hail = undefined_real
    !? = undefined_real

    ! ----------------------
    ! Allocate domain data
    ! ----------------------

    !associate(NumSoilLayer => model%noahmp%config%domain%NumSoilLayer)
    !associate(Domain => model%noahmp%config%domain)

    !allocate(Domain%DepthSoilLayer(1:NumSoilLayer))
    !allocate(Domain%ThicknessSoilLayer(1:NumSoilLayer))






    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitNoahMP 


end module lnd_comp_driver
