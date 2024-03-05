module lnd_comp_import_export

  !----------------------------------------------------------------------------
  ! NoahMP import and export fields exchanged with the coupler
  !----------------------------------------------------------------------------

  use ESMF, only: ESMF_GridComp
  use ESMF, only: ESMF_State, ESMF_Mesh
  use ESMF, only: ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS

  use NUOPC_Model, only: NUOPC_ModelGet

  use lnd_comp_shr, only: ChkErr

  implicit none
  private ! except

  public :: advertise_fields
  public :: realize_fields

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(*),parameter :: modName =  "(lnd_comp_import_export)"
  character(*),parameter :: u_FILE_u = __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    integer,             intent(out) :: rc

    ! local variables
    integer           :: n
    type(ESMF_State)  :: importState
    type(ESMF_State)  :: exportState
    character(len=*), parameter :: subname=trim(modName)//':(advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine advertise_fields

  !===============================================================================

  subroutine realize_fields(importState, exportState, Emesh, rc)

    ! input/output variables
    type(ESMF_State) , intent(inout) :: importState
    type(ESMF_State) , intent(inout) :: exportState
    type(ESMF_Mesh)  , intent(in)    :: Emesh
    integer          , intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(realize_fields)'
    !---------------------------------------------------------------------------

  end subroutine realize_fields

end module lnd_comp_import_export
