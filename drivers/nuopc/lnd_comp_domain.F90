module lnd_comp_domain

  ! This file contains domain related coupling routines 

  use ESMF,  only: ESMF_GridComp, ESMF_GridCreateMosaic
  use ESMF,  only: ESMF_MeshCreate, ESMF_MeshGet, ESMF_MeshSet, ESMF_FILEFORMAT_ESMFMESH
  use ESMF,  only: ESMF_Field, ESMF_FieldGet, ESMF_FieldRead 
  use ESMF,  only: ESMF_LogWrite, ESMF_INDEX_GLOBAL, ESMF_Decomp_Flag
  use ESMF,  only: ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_DECOMP_SYMMEDGEMAX
  use ESMF,  only: ESMF_STAGGERLOC_CORNER, ESMF_STAGGERLOC_CENTER

  use lnd_comp_shr  , only: ChkErr
  use lnd_comp_shr  , only: ChkErrNc
  use lnd_comp_shr  , only: FieldRead
  use lnd_comp_types, only: model_type 

  implicit none
  private

  public :: SetDomain

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(*), parameter :: modName = "(lnd_comp_domain)"

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetDomain(gcomp, model, rc)

    use netcdf

    ! input/output variables
    type(ESMF_GridComp), intent(in) :: gcomp
    type(model_type), intent(inout) :: model 
    integer, intent(out)            :: rc

    ! local variables
    type(ESMF_Field)                    :: field
    integer                             :: ncid, dimid, ncerr
    integer                             :: n, lsize
    integer, allocatable                :: decomptile(:,:)
    integer                             :: maxIndex(2)
    type(ESMF_Decomp_Flag), allocatable :: decompflagPTile(:,:)
    integer         , allocatable       :: tmp1d(:)
    character(len=*), parameter         :: subname = trim(modName)//':(SetDomain) '

    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Create ESMF mesh
    ! ---------------------

    if (trim(model%nmlist%mesh_file) /= '') then
       ! read mesh
       model%domain%mesh = ESMF_MeshCreate(trim(model%nmlist%mesh_file), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc) 
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! read mask from domain file
       call FieldRead(model%domain%mesh, trim(model%nmlist%domain_file), 'mask', tmp1d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
       ! set mask in the mesh
       call ESMF_MeshSet(model%domain%mesh, elementMask=tmp1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (trim(model%nmlist%mosaic_file) /= '') then
       ! open file
       ncerr = nf90_open(trim(model%nmlist%mosaic_file), NF90_NOWRITE, ncid=ncid)
       if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
       
       ! query dimension size of ntiles 
       ncerr = nf90_inq_dimid(ncid, 'ntiles', dimid)
       if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
       ncerr = nf90_inquire_dimension(ncid, dimid, len=model%domain%ntiles)
       if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

       ! close file
       ncerr = nf90_close(ncid)
       if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

       ! allocate required data structures and fill them
       if (.not. allocated(decomptile)) allocate(decomptile(2,model%domain%ntiles))
       if (.not. allocated(decompflagPTile)) allocate(decompflagPTile(2,model%domain%ntiles))

       do n = 1, model%domain%ntiles
          decomptile(1,n) = model%domain%layout(1)
          decomptile(2,n) = model%domain%layout(2)
          decompflagPTile(:,n) = (/ ESMF_DECOMP_SYMMEDGEMAX, ESMF_DECOMP_SYMMEDGEMAX /)
       end do

       ! create grid
       model%domain%grid = ESMF_GridCreateMosaic(filename=trim(model%nmlist%mosaic_file), &
          regDecompPTile=decomptile, tileFilePath=trim(model%nmlist%input_dir), &
          decompflagPTile=decompflagPTile, &
          staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
          indexflag=ESMF_INDEX_GLOBAL, &
          name='lnd_grid', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return



       print*, model%domain%ntiles

    end if 

    ! ---------------------
    ! Get fraction from orography file
    ! ---------------------




 

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetDomain

end module lnd_comp_domain
