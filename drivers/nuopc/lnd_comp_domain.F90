module lnd_comp_domain

  ! This file contains domain related coupling routines 

  use ESMF, only: operator(/=)
  use ESMF, only: ESMF_GridComp, ESMF_GridCreateMosaic
  use ESMF, only: ESMF_MeshCreate, ESMF_MeshGet, ESMF_MeshSet, ESMF_FILEFORMAT_ESMFMESH
  use ESMF, only: ESMF_Field, ESMF_FieldGet, ESMF_FieldRead 
  use ESMF, only: ESMF_LogWrite, ESMF_INDEX_GLOBAL, ESMF_Decomp_Flag
  use ESMF, only: ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_DECOMP_SYMMEDGEMAX
  use ESMF, only: ESMF_STAGGERLOC_CORNER, ESMF_STAGGERLOC_CENTER
  use ESMF, only: ESMF_MeshGetFieldBounds, ESMF_MESHLOC_ELEMENT
  use ESMF, only: ESMF_FieldCreate, ESMF_FieldRegridGetArea
  use ESMF, only: ESMF_CoordSys_Flag, ESMF_COORDSYS_CART
  use ESMF, only: ESMF_TYPEKIND_R8, ESMF_KIND_R8, ESMF_MeshWriteVTK
  use ESMF, only: ESMF_GridAddItem, ESMF_GRIDITEM_MASK, ESMF_STAGGERLOC_CENTER
  use ESMF, only: ESMF_GridCreate, ESMF_GridGetCoord, ESMF_INDEX_DELOCAL
  use ESMF, only: ESMF_GridCreateNoPeriDim, ESMF_COORDSYS_SPH_RAD
  use ESMF, only: ESMF_FAILURE, ESMF_FILEFORMAT_SCRIP

  use lnd_comp_shr  , only: ChkErr
  use lnd_comp_shr  , only: ChkErrNc
  use lnd_comp_shr  , only: FieldRead
  use lnd_comp_types, only: model_type
  use lnd_comp_types, only: field_type
  use lnd_comp_types, only: iMosaic, iScrip
  use lnd_comp_kind , only: cl => shr_kind_cl
  use lnd_comp_kind , only: r4 => shr_kind_r4
  use lnd_comp_kind , only: r8 => shr_kind_r8
  use lnd_comp_io   , only: ReadFile 
  use lnd_comp_io   , only: GetNumTiles

  implicit none
  private

  public :: SetDomain

  !-----------------------------------------------------------------------------
  ! Private module data
  !-----------------------------------------------------------------------------

  real(r8), parameter     :: rearth = 6.3712e+6_r8
  character(*), parameter :: modName = "(lnd_comp_domain)"

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetDomain(gcomp, model, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in) :: gcomp
    type(model_type), intent(inout) :: model
    integer, intent(out)            :: rc

    character(len=*), parameter         :: subname = trim(modName)//':(SetDomain) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! create domain
    if (model%domain%dtype == iMosaic) then
       call SetDomainMosaicGlobal(gcomp, model, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (model%domain%dtype == iScrip) then
       call SetDomainScrip(gcomp, model, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! check mesh for debugging purposes
    if (model%nmlist%debug_level > 2) then
       call ESMF_MeshWriteVTK(model%domain%mesh, "lnd_mesh", rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetDomain

  !=============================================================================

  subroutine SetDomainMosaicGlobal(gcomp, model, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in) :: gcomp
    type(model_type), intent(inout) :: model 
    integer, intent(out)            :: rc

    ! local variables
    type(ESMF_Field)                    :: field
    type(ESMF_Field)                    :: farea
    type(ESMF_CoordSys_Flag)            :: coordSys
    type(field_type)                    :: flds(1)
    integer                             :: n, lsize
    integer                             :: tlb(1), tub(1), tc(1)
    integer                             :: spatialDim, numOwnedElements
    integer, allocatable                :: decomptile(:,:)
    type(ESMF_Decomp_Flag), allocatable :: decompflagPTile(:,:)
    real(r4), target, allocatable       :: tmpr4(:)
    real(ESMF_KIND_R8), pointer         :: ptr1d(:)
    real(r8), allocatable               :: ownedElemCoords(:)
    character(len=cl)                   :: filename, msg
    character(len=*), parameter         :: subname = trim(modName)//':(SetDomainMosaicGlobal) '
    !----------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Create ESMF grid
    ! ---------------------

    ! allocate required data structures and fill them
    if (.not. allocated(decomptile)) allocate(decomptile(2,model%domain%ntiles))
    if (.not. allocated(decompflagPTile)) allocate(decompflagPTile(2,model%domain%ntiles))

    do n = 1, model%domain%ntiles
       decomptile(1,n) = model%domain%layout(1)
       decomptile(2,n) = model%domain%layout(2)
       decompflagPTile(:,n) = (/ ESMF_DECOMP_SYMMEDGEMAX, ESMF_DECOMP_SYMMEDGEMAX /)
    end do

    ! create mosaic grid
    model%domain%grid = ESMF_GridCreateMosaic(filename=trim(model%nmlist%mosaic_file), &
       regDecompPTile=decomptile, tileFilePath=trim(model%nmlist%input_dir), &
       decompflagPTile=decompflagPTile, &
       staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
       indexflag=ESMF_INDEX_GLOBAL, &
       name='lnd_grid', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Add mask to the grid
    ! ---------------------

    call ESMF_GridAddItem(model%domain%grid, itemflag=ESMF_GRIDITEM_MASK, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Convert ESMF grid to mesh 
    ! ---------------------

    model%domain%mesh = ESMF_MeshCreate(model%domain%grid, 'lnd_mesh_from_grid', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Query sizes from mesh
    ! ---------------------

    call ESMF_MeshGetFieldBounds(model%domain%mesh, meshloc=ESMF_MESHLOC_ELEMENT, &
      totalLBound=tlb, totalUBound=tub, totalCount=tc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    model%domain%begl = tlb(1)
    model%domain%endl = tub(1)
    model%domain%im = tc(1)
    write(msg, fmt='(A,3I5)') trim(subname)//' : begl, endl, im = ', model%domain%begl, model%domain%endl, model%domain%im
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    !----------------------
    ! Allocate temporary data structures
    !----------------------

    if (.not. allocated(tmpr4)) then
       allocate(tmpr4(model%domain%begl:model%domain%endl))
       tmpr4(:) = 0.0
    end if

    ! ---------------------
    ! Get land fraction information from orography file
    ! ---------------------

    if (.not. allocated(model%domain%frac)) then
       allocate(model%domain%frac(model%domain%begl:model%domain%endl))
    end if

    filename = trim(model%nmlist%input_dir)//'oro_data.tile*.nc'
    flds(1)%short_name = 'land_frac'
    flds(1)%ptr1r4 => tmpr4
    call ReadFile(model, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    model%domain%frac = dble(tmpr4)

    ! ---------------------
    ! Calculate mask from land-sea fraction
    ! ---------------------

    if (.not. allocated(model%domain%mask)) then
       allocate(model%domain%mask(model%domain%begl:model%domain%endl))
    end if

    where (model%domain%frac(:) > 0.0_r8)
       model%domain%mask(:) = 1
    elsewhere
       model%domain%mask(:) = 0
    end where

    ! ---------------------
    ! Set mask in mesh 
    ! ---------------------

    call ESMF_MeshSet(model%domain%mesh, elementMask=model%domain%mask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Get height from orography file
    ! ---------------------

    ! read field
    filename = trim(model%nmlist%input_dir)//'oro_data.tile*.nc'
    flds(1)%short_name = 'orog_raw'
    flds(1)%ptr1r4 => tmpr4
    call ReadFile(model, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! allocate data
    if (.not. allocated(model%domain%hgt)) then
       allocate(model%domain%hgt(model%domain%begl:model%domain%endl))
    end if
    model%domain%hgt = dble(tmpr4)

    ! ---------------------
    ! Query cell area 
    ! ---------------------

    ! create field in R8 type
    farea = ESMF_FieldCreate(model%domain%mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get cell area to the field
    call ESMF_FieldRegridGetArea(farea, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (.not. allocated(model%domain%garea)) then
       allocate(model%domain%garea(model%domain%begl:model%domain%endl))
    end if

    ! retrieve pointer and fill area array
    call ESMF_FieldGet(farea, farrayPtr=ptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    model%domain%garea(:) = ptr1d(:)

    ! make unit conversion from square radians to square meters if it is
    ! required
    call ESMF_MeshGet(model%domain%mesh, coordSys=coordSys, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (coordSys /= ESMF_COORDSYS_CART) then
       model%domain%garea(:) = model%domain%garea(:)*(rearth**2)
    end if

    ! ---------------------
    ! Query coordiates from ESMF mesh
    ! ---------------------

    ! determine dimensions in mesh
    call ESMF_MeshGet(model%domain%mesh, spatialDim=spatialDim, &
      numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! allocate data
    if (.not. allocated(model%domain%lats)) then
       allocate(model%domain%lats(numOwnedElements))
    end if
    if (.not. allocated(model%domain%lons)) then
       allocate(model%domain%lons(numOwnedElements))
    end if

    ! allocate temporary array
    if (.not. allocated(ownedElemCoords)) then
       allocate(ownedElemCoords(spatialDim*numOwnedElements))
    end if

    call ESMF_MeshGet(model%domain%mesh, ownedElemCoords=ownedElemCoords)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1,numOwnedElements
       model%domain%lons(n) = ownedElemCoords(2*n-1)
       model%domain%lats(n) = ownedElemCoords(2*n)
    end do

    ! clean memory
    deallocate(ownedElemCoords)
    deallocate(tmpr4)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetDomainMosaicGlobal

  !=============================================================================

  subroutine SetDomainScrip(gcomp, model, rc)

    use netcdf

    ! input/output variables
    type(ESMF_GridComp), intent(in) :: gcomp
    type(model_type), intent(inout) :: model
    integer, intent(out)            :: rc

    ! local variables
    integer                       :: n
    integer                       :: tlb(1), tub(1), tc(1)
    integer                       :: spatialDim, numOwnedElements
    character(len=cl)             :: filename, msg
    logical                       :: isFound
    type(ESMF_Field)              :: farea
    type(ESMF_CoordSys_Flag)      :: coordSys
    type(field_type)              :: flds(1)
    real(ESMF_KIND_R8), pointer   :: ptr1d(:)
    real(r4), target, allocatable :: tmpr4(:)
    real(r8), allocatable         :: ownedElemCoords(:)
    character(len=*), parameter   :: subname = trim(modName)//':(SetDomainScrip) '
    !----------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Read in SCRIP file and create grid
    ! ---------------------

    model%domain%grid = ESMF_GridCreate(filename=trim(model%nmlist%scrip_file), fileformat=ESMF_FILEFORMAT_SCRIP, &
       isSphere=.false., addCornerStagger=.true., indexflag=ESMF_INDEX_GLOBAL, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Convert ESMF grid to mesh 
    ! ---------------------

    model%domain%mesh = ESMF_MeshCreate(model%domain%grid, 'lnd_mesh_from_grid', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Query sizes from mesh
    ! ---------------------

    call ESMF_MeshGetFieldBounds(model%domain%mesh, meshloc=ESMF_MESHLOC_ELEMENT, &
      totalLBound=tlb, totalUBound=tub, totalCount=tc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    model%domain%begl = tlb(1)
    model%domain%endl = tub(1)
    model%domain%im = tc(1)
    write(msg, fmt='(A,3I5)') trim(subname)//' : begl, endl, im = ', model%domain%begl, model%domain%endl, model%domain%im
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    !----------------------
    ! Allocate temporary data structures
    !----------------------

    if (.not. allocated(tmpr4)) then
       allocate(tmpr4(model%domain%begl:model%domain%endl))
       tmpr4(:) = 0.0
    end if

    ! ---------------------
    ! Get height from orography file
    ! ---------------------

    ! read field
    filename = trim(model%nmlist%input_dir)//'oro_data.nc'
    flds(1)%short_name = 'orog_raw'
    flds(1)%ptr1r4 => tmpr4
    call ReadFile(model, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! allocate data
    if (.not. allocated(model%domain%hgt)) then
       allocate(model%domain%hgt(model%domain%begl:model%domain%endl))
    end if
    model%domain%hgt = dble(tmpr4)

    ! ---------------------
    ! Query cell area 
    ! ---------------------

    ! create field in R8 type
    farea = ESMF_FieldCreate(model%domain%mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get cell area to the field
    call ESMF_FieldRegridGetArea(farea, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (.not. allocated(model%domain%garea)) then
       allocate(model%domain%garea(model%domain%begl:model%domain%endl))
    end if

    ! retrieve pointer and fill area array
    call ESMF_FieldGet(farea, farrayPtr=ptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    model%domain%garea(:) = ptr1d(:)

    ! make unit conversion from square radians to square meters if it is
    ! required
    call ESMF_MeshGet(model%domain%mesh, coordSys=coordSys, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (coordSys /= ESMF_COORDSYS_CART) then
       model%domain%garea(:) = model%domain%garea(:)*(rearth**2)
    end if

    ! ---------------------
    ! Query coordiates from ESMF mesh
    ! ---------------------

    ! determine dimensions in mesh
    call ESMF_MeshGet(model%domain%mesh, spatialDim=spatialDim, &
      numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! allocate data
    if (.not. allocated(model%domain%lats)) then
       allocate(model%domain%lats(numOwnedElements))
    end if
    if (.not. allocated(model%domain%lons)) then
       allocate(model%domain%lons(numOwnedElements))
    end if

    ! allocate temporary array
    if (.not. allocated(ownedElemCoords)) then
       allocate(ownedElemCoords(spatialDim*numOwnedElements))
    end if

    call ESMF_MeshGet(model%domain%mesh, ownedElemCoords=ownedElemCoords)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1,numOwnedElements
       model%domain%lons(n) = ownedElemCoords(2*n-1)
       model%domain%lats(n) = ownedElemCoords(2*n)
    end do

    ! clean memory
    deallocate(ownedElemCoords)
    deallocate(tmpr4)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetDomainScrip

end module lnd_comp_domain
