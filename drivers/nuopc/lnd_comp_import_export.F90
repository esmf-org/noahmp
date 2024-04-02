module lnd_comp_import_export

  !----------------------------------------------------------------------------
  ! NoahMP import and export fields exchanged with the coupler
  !----------------------------------------------------------------------------

  use ESMF, only: operator(==)
  use ESMF, only: ESMF_GridComp, ESMF_Mesh
  use ESMF, only: ESMF_State, ESMF_StateItem_Flag, ESMF_StateGet
  use ESMF, only: ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF, only: ESMF_LOGMSG_ERROR, ESMF_END_ABORT, ESMF_Finalize
  use ESMF, only: ESMF_Field, ESMF_FieldStatus_Flag
  use ESMF, only: ESMF_FieldGet, ESMF_MeshGet
  use ESMF, only: ESMF_STATEITEM_FIELD, ESMF_FIELDSTATUS_COMPLETE
  use ESMF, only: ESMF_FAILURE

  use NUOPC, only: NUOPC_Advertise
  use NUOPC_Model, only: NUOPC_ModelGet

  use lnd_comp_shr  , only: ChkErr
  use lnd_comp_kind , only: r8 => shr_kind_r8
  use lnd_comp_types, only: fld_list_type, fldsMax
  use lnd_comp_types, only: fldsToLnd, fldsToLnd_num
  use lnd_comp_types, only: fldsFrLnd, fldsFrLnd_num
  use lnd_comp_types, only: model_type

  implicit none
  private ! except

  public :: AdvertiseFields
  public :: ImportFields
  public :: RealizeFields

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(*),parameter :: modName =  "(lnd_comp_import_export)"
  character(*),parameter :: u_FILE_u = __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine AdvertiseFields(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    integer,             intent(out) :: rc

    ! local variables
    integer           :: n
    type(ESMF_State)  :: importState
    type(ESMF_State)  :: exportState
    character(len=*), parameter :: subname=trim(modName)//':(AdvertiseFields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !--------------------------------
    ! Query import and export states
    !--------------------------------

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    ! export to med
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Sl_lfrin')

    ! export to atm
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Sl_sfrac')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Fall_lat')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Fall_sen')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Fall_evap')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Sl_tref')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Sl_qref')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Sl_q')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Fall_gflx')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Fall_roff')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Fall_soff')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Sl_cmm')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Sl_chh')
    call fldListAdd(fldsFrLnd_num, fldsFrlnd, 'Sl_zvfun')

    ! Now advertise above export fields
    do n = 1,fldsFrLnd_num
       call NUOPC_Advertise(exportState, standardName=fldsFrLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    ! import from atm
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_z')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_tbot')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_ta')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_tskn')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_pslv')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_prsl')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_pbot')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_shum')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_qa')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_u')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_v')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_exner')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_ustar')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Faxa_swdn')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Faxa_lwdn')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Faxa_swnet')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Faxa_rainc')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Faxa_rainl')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Faxa_rain')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Faxa_snow')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Faxa_snowc')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Faxa_snowl')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_vfrac')
    call fldListAdd(fldsToLnd_num, fldsToLnd, 'Sa_zorl')

    ! Now advertise import fields
    do n = 1,fldsToLnd_num
       call NUOPC_Advertise(importState, standardName=fldsToLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine AdvertiseFields

  !===============================================================================

  subroutine RealizeFields(importState, exportState, Emesh, rc)

    ! input/output variables
    type(ESMF_State) , intent(inout) :: importState
    type(ESMF_State) , intent(inout) :: exportState
    type(ESMF_Mesh)  , intent(in)    :: Emesh
    integer          , intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(RealizeFields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! export state
    call fldListRealize(exportState, fldsFrLnd, fldsFrLnd_num, trim(subname)//':NoahMP_Export', mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! import state
    call fldListRealize(importState, fldsToLnd, fldsToLnd_num, trim(subname)//':NoahMP_Import', mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine RealizeFields

  !===============================================================================

  subroutine fldListAdd(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    integer,          optional, intent(in)    :: ungridded_lbound
    integer,          optional, intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(fldListAdd)'
    !-------------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Set up a list of field information
    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine fldListAdd

  !===============================================================================

  subroutine fldListRealize(state, fldList, numflds, tag, mesh, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(inout) :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname=trim(modName)//':(fldListRealize)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    do n = 1, numflds
       stdname = trim(fldList(n)%stdname)
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          ! Create the field
          if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                  ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                  ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                  gridToFieldMap=(/2/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
             ESMF_LOGMSG_INFO)

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Set flag for connected fields
          fldList(n)%connected = .true.
       else
          call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
             ESMF_LOGMSG_INFO)
          call ESMF_StateRemove(state, (/stdname/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine fldListRealize

  !===============================================================================

  subroutine ImportFields(gcomp, model, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    type(model_type)   , intent(inout) :: model 
    integer            , intent(out)   :: rc

    ! local variables
    type(ESMF_State)            :: importState
    character(len=*), parameter :: subname=trim(modName)//': (ImportFields)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get export state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! atm input fields
    ! -----------------------

    !call StateGetImport(importState, 'Sa_z'      , noahmp%forc%hgt, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Sa_ta'     , model%forcing%TemperatureAirRefHeight, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Sa_tbot'   , model%forcing%TemperatureAirRefHeight, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Sa_tskn'   , model%forcing%TemperatureSoilBottom, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Sa_prsl'   , model%forc%t, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Sa_pbot'   , model%forcing%PressureAirSurface, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Sa_pslv'   , model%forcing%PressureAirRefHeight, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Sa_shum'   , model%forcing%SpecHumidityRefHeight, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Sa_qa'     , model%forcing%SpecHumidityRefHeight, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Faxa_swdn' , model%forcing%RadSwDownRefHeight, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Faxa_lwdn' , model%forcing%RadLwDownRefHeight, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Faxa_swnet', noahmp%forc%snet, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Sa_u'      , model%forcing%WindEastwardRefHeight, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call StateGetImport(importState, 'Sa_v'      , model%forcing%WindNorthwardRefHeight, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Sa_exner'  , noahmp%forc%prslk1, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Sa_ustar'  , noahmp%forc%ustar1, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Faxa_rain' , noahmp%forc%tprcp, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Faxa_rainc', noahmp%forc%tprcpc, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Faxa_rainl', noahmp%forc%tprcpl, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Faxa_snow' , noahmp%forc%snow, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Faxa_snowc', noahmp%forc%snowc, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Faxa_snowl', noahmp%forc%snowl, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Sa_vfrac'  , noahmp%forc%vegfrac, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call StateGetImport(importState, 'Sa_zorl'   , noahmp%forc%zorl, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ImportFields

  !===============================================================================

  subroutine StateGetImport(state, fldname, arr1d, rc)

    ! fill in noahmp import data for 1d field

    ! input/output variabes
    type(ESMF_State) , intent(in)    :: state
    character(len=*) , intent(in)    :: fldname
    real(r8)         , intent(inout) :: arr1d(:)
    integer          , intent(out)   :: rc

    ! local variables
    real(r8), pointer           :: fldPtr1d(:)
    type(ESMF_StateItem_Flag)   :: itemType
    character(len=*), parameter :: subname=trim(modName)//':(StateGetImport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=trim(fldname), itemType=itemType, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (itemType == ESMF_STATEITEM_FIELD) then
       call StateGetFldPtr(State, trim(fldname), fldptr1d=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       arr1d(:) = fldptr1d(:)
       call CheckForNaNs(arr1d, trim(fldname), 1, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(subname//' '//trim(fldname)//' is not in the state!', ESMF_LOGMSG_INFO)
    end if

  end subroutine StateGetImport

  !===============================================================================

  subroutine StateGetFldPtr(state, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State),            intent(in)    :: State
    character(len=*),            intent(in)    :: fldname
    real(R8), pointer, optional, intent(out)   :: fldptr1d(:)
    real(R8), pointer, optional, intent(out)   :: fldptr2d(:,:)
    integer,                     intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    character(len=*), parameter :: subname=trim(modName)//':(StateGetFldPtr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (present(fldptr1d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (present(fldptr2d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)//": either fldptr1d or fldptr2d must be an input argument", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

  end subroutine StateGetFldPtr

  !=============================================================================

  subroutine CheckForNaNs(array, fname, begg, rc)

    ! input/output variables
    real(r8)        , intent(in)  :: array(:)
    character(len=*), intent(in)  :: fname
    integer         , intent(in)  :: begg
    integer         , intent(out) :: rc

    ! local variables
    integer :: i
    character(len=*), parameter :: subname='(CheckForNaNs)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Check if any input from mediator or output to mediator is NaN
    if (any(isnan(array))) then
       write(*,*) '# of NaNs = ', count(isnan(array))
       write(*,*) 'Which are NaNs = ', isnan(array)
       do i = 1, size(array)
          if (isnan(array(i))) then
             write(*,*) "NaN found in field ", trim(fname), ' at gridcell index ',begg+i-1
          end if
       end do
       call ESMF_LogWrite(trim(subname)//": one or more of the output from NoahMP to the coupler are NaN", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine CheckForNaNs

end module lnd_comp_import_export
