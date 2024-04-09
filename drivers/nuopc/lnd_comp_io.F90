module lnd_comp_io

  ! This file contains I/O routines for the NoahMP land surface model

  use ESMF          , only: operator(==), operator(/=)
  use ESMF          , only: ESMF_FieldBundle, ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF          , only: ESMF_FieldBundleGet, ESMF_FieldBundleRead, ESMF_FieldBundleWrite
  use ESMF          , only: ESMF_FieldBundleRemove, ESMF_FieldBundleDestroy
  use ESMF          , only: ESMF_FieldBundleRedistStore, ESMF_FieldBundleRedist
  use ESMF          , only: ESMF_RouteHandleDestroy, ESMF_RouteHandle, ESMF_FieldWrite
  use ESMF          , only: ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet, ESMF_FieldWriteVTK
  use ESMF          , only: ESMF_FieldDestroy, ESMF_ArraySpec, ESMF_ArraySpecSet
  use ESMF          , only: ESMF_LogWrite, ESMF_FieldStatus_Flag, ESMF_GeomType_Flag
  use ESMF          , only: ESMF_Mesh, ESMF_MeshGet, ESMF_FieldBundleAddReplace
  use ESMF          , only: ESMF_KIND_R8, ESMF_TYPEKIND_R8
  use ESMF          , only: ESMF_KIND_R4, ESMF_TYPEKIND_R4
  use ESMF          , only: ESMF_KIND_I4, ESMF_TYPEKIND_I4
  use ESMF          , only: ESMF_STAGGERLOC_CENTER, ESMF_MESHLOC_ELEMENT
  use ESMF          , only: ESMF_INDEX_DELOCAL, ESMF_INDEX_GLOBAL
  use ESMF          , only: ESMF_MAXSTR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF          , only: ESMF_LOGMSG_ERROR, ESMF_LOGMSG_INFO
  use ESMF          , only: ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_COMPLETE
  use ESMF          , only: ESMF_AttributeAdd, ESMF_AttributeSet, ESMF_Grid
  use ESMF          , only: ESMF_LocStream, ESMF_LocStreamCreate, ESMF_LocStreamDestroy
  use ESMF          , only: ESMF_GridCreate, ESMF_FILESTATUS_OLD, ESMF_TypeKind_Flag
  use ESMF          , only: ESMF_DistGrid, ESMF_DistGridCreate
  use ESMF          , only: ESMF_VMGet, ESMF_VMBarrier, ESMF_VM

  use lnd_comp_types, only: model_type
  use lnd_comp_types, only: field_type
  use lnd_comp_types, only: iMosaic, iScrip
  use lnd_comp_types, only: fldsMaxIO, histflds, restflds
  use lnd_comp_kind , only: cl => shr_kind_cl
  use lnd_comp_kind , only: r4 => shr_kind_r4
  use lnd_comp_kind , only: r8 => shr_kind_r8
  use lnd_comp_kind , only: i4 => shr_kind_i4
  use lnd_comp_shr  , only: chkerr
  use lnd_comp_shr  , only: chkerrnc

  implicit none
  private

  public :: GetNumTiles
  public :: ReadFile
  public :: ReadStatic
  public :: ReadIC
  public :: SetupWriteFields
  public :: WriteFile

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer(ESMF_KIND_I4)        :: missing_i4 = -999
  real(ESMF_KIND_R4)           :: missing_r4 = 1.0e20
  real(ESMF_KIND_R8)           :: missing_r8 = 1.0d20

  character(*), parameter      :: modName = "(lnd_comp_io)"
  character(len=*) , parameter :: u_FILE_u = __FILE__

!===============================================================================
contains
!===============================================================================

  integer function GetNumTiles(filename)

    use netcdf

    ! input/output variables
    character(len=*), intent(in) :: filename

    ! local variables
    integer :: ncid, dimid, ncerr
    !-----------------------------------------------------------------------------

    ! open file
    ncerr = nf90_open(trim(filename), NF90_NOWRITE, ncid=ncid)
    if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

    ! query dimension size of ntiles 
    ncerr = nf90_inq_dimid(ncid, 'ntiles', dimid)
    if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
    ncerr = nf90_inquire_dimension(ncid, dimid, len=GetNumTiles)
    if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

    ! close file
    ncerr = nf90_close(ncid)
    if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

  end function GetNumTiles

  !===============================================================================

  subroutine ReadStatic(model, rc)

    ! input/output variables
    type(model_type), target, intent(inout) :: model 
    integer         ,         intent(inout) :: rc

    ! local variables
    type(field_type), allocatable :: flds(:)
    real(r4), target, allocatable :: tmpr4(:)
    real(r8), target, allocatable :: tmpr8(:)
    real(r4), target, allocatable :: tmp2r4(:,:)
    real(r8), target, allocatable :: tmp2r8(:,:)
    character(len=CL)             :: filename
    character(len=*), parameter   :: subname=trim(modName)//':(ReadStatic) '
    !-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! allocate temporary data structures
    !----------------------

    if (.not. allocated(tmpr4)) then
       allocate(tmpr4(model%domain%begl:model%domain%endl))
       tmpr4(:) = 0.0
    end if

    if (.not. allocated(tmpr8)) then
       allocate(tmpr8(model%domain%begl:model%domain%endl))
       tmpr8(:) = 0.0_r8
    end if

    if (.not. allocated(tmp2r4)) then
       allocate(tmp2r4(model%domain%begl:model%domain%endl,12))
       tmp2r4(:,:) = 0.0
    end if

    if (.not. allocated(tmp2r8)) then
       allocate(tmp2r8(model%domain%begl:model%domain%endl,1))
       tmp2r8(:,:) = 0.0_r8
    end if

    !----------------------
    ! Set file name for sfc 
    !----------------------

    if (trim(model%nmlist%InputType) == 'sfc') then
       if (model%domain%ntiles == 1) then
          filename = trim(model%nmlist%input_dir)//'sfc_data.nc'
       else
          filename = trim(model%nmlist%input_dir)//'sfc_data.tile*.nc'
       end if
    end if

    !----------------------
    ! Read soil type
    !----------------------

    if (trim(model%nmlist%InputType) == 'sfc') then
       allocate(flds(1))
       flds(1)%short_name = 'stype'
       flds(1)%ptr1r8 => tmpr8
       call ReadFile(model, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       model%coupling%SoilType = int(tmpr8)
      deallocate(flds)
    else
       write(filename, fmt="(A)") trim(model%nmlist%input_dir)//trim(model%nmlist%InputSoilType)
       allocate(flds(1))
       flds(1)%short_name = 'soil_type'
       flds(1)%ptr1r4 => tmpr4
       call ReadFile(model, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       model%coupling%SoilType = int(tmpr4)
       deallocate(flds)
    end if

    !----------------------
    ! Read vegetation type
    !----------------------

    if (trim(model%nmlist%InputType) == 'sfc') then
       allocate(flds(1))
       flds(1)%short_name = 'vtype'
       flds(1)%ptr1r8 => tmpr8
       call ReadFile(model, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       model%coupling%VegType = int(tmpr8)
       deallocate(flds)
    else
       write(filename, fmt="(A)") trim(model%nmlist%input_dir)//trim(model%nmlist%InputVegType)
       allocate(flds(1))
       flds(1)%short_name = 'vegetation_type'
       flds(1)%ptr1r4 => tmpr4
       call ReadFile(model, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       model%coupling%VegType = int(tmpr4)
       deallocate(flds)
    end if

    !----------------------
    ! Read slope type
    !----------------------

    if (trim(model%nmlist%InputType) == 'sfc') then
       allocate(flds(1))
       flds(1)%short_name = 'slope'
       flds(1)%ptr1r8 => tmpr8
       call ReadFile(model, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       model%coupling%RunoffSlopeType = int(tmpr8)
       deallocate(flds)
    else
       write(filename, fmt="(A)") trim(model%nmlist%input_dir)//trim(model%nmlist%InputSlopeType)
       allocate(flds(1))
       flds(1)%short_name = 'slope_type'
       flds(1)%ptr1r4 => tmpr4
       call ReadFile(model, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       model%coupling%RunoffSlopeType = int(tmpr4)
       deallocate(flds)
    end if

    !----------------------
    ! Read deep soil temperature
    !----------------------

    if (trim(model%nmlist%InputType) == 'sfc') then
       allocate(flds(1))
       flds(1)%short_name = 'tg3'
       flds(1)%nrec = 1; flds(1)%ptr2r8 => tmp2r8
       call ReadFile(model, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       model%forcing%TemperatureSoilBottom = tmp2r8(:,1)
       deallocate(flds)
    else
       write(filename, fmt="(A)") trim(model%nmlist%input_dir)//trim(model%nmlist%InputBottomTemp)
       allocate(flds(1))
       flds(1)%short_name = 'substrate_temperature'
       flds(1)%ptr1r4 => tmpr4
       call ReadFile(model, filename, flds, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       model%forcing%TemperatureSoilBottom = dble(tmpr4)
       deallocate(flds)
    end if

    !----------------------
    ! Read upper bound on max albedo over deep snow
    !----------------------

    !if (trim(noahmp%nmlist%ic_type) == 'sfc') then
    !   allocate(flds(1))
    !   flds(1)%short_name = 'snoalb'
    !   flds(1)%nrec = 1; flds(1)%ptr2r8 => tmp2r8
    !   call read_tiled_file(noahmp, filename, flds, rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !   model%coupling%AlbedoSnowPrev = tmp2r8(:,1)
    !   deallocate(flds)
    !else
    !   allocate(flds(1))
    !   flds(1)%short_name = 'maximum_snow_albedo'
    !   flds(1)%ptr1r4 => tmpr4
    !   call read_tiled_file(noahmp, filename, flds, rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !   noahmp%model%snoalb = dble(tmpr4)
    !   deallocate(flds)
    !end if

    !----------------------
    ! Read vegetation greenness, monthly average 
    !----------------------

    write(filename, fmt="(A)") trim(model%nmlist%input_dir)//trim(model%nmlist%InputVegFracAnnMax)
    allocate(flds(1))
    flds(1)%short_name = 'vegetation_greenness'
    flds(1)%nrec = 12; flds(1)%ptr2r4 => tmp2r4
    call ReadFile(model, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    model%coupling%VegFracAnnMax(:) = maxval(dble(tmp2r4(:,:)), dim=2)
    deallocate(flds)

    !----------------------
    ! Soil color
    !----------------------

    !allocate(flds(1))
    !write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', noahmp%domain%ni, '.soil_color.tile*.nc'
    !flds(1)%short_name = 'soil_color'
    !flds(1)%ptr1r4 => tmpr4
    !call read_tiled_file(noahmp, filename, flds, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !noahmp%model%soilcol = int(tmpr4)
    !deallocate(flds)

    call ESMF_LogWrite(subname//' done for '//trim(filename), ESMF_LOGMSG_INFO)

  end subroutine ReadStatic

  !===============================================================================

  subroutine ReadIC(model, rc)

    ! input/output variables
    type(model_type), target, intent(inout) :: model
    integer         ,         intent(inout) :: rc

    ! local variables
    type(field_type), allocatable :: flds(:)
    character(len=CL)             :: filename
    character(len=*), parameter   :: subname=trim(modName)//':(ReadIC) '
    !-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! Read information from input file
    !----------------------

    write(filename, fmt="(A)") trim(model%nmlist%input_dir)//trim(model%nmlist%InputIC)

    if (trim(model%nmlist%InputType) == 'sfc') then

    else
       !----------------------
       ! Allocate field list
       !----------------------

       allocate(flds(3))

       !----------------------
       ! Snow albedo at last time step (CLASS type), CanopyTotalWater
       !----------------------

       !----------------------
       ! Surface skin temperature
       !----------------------

       !----------------------
       ! Water equivalent accumulated snow depth
       !----------------------

       flds(1)%short_name = 'snow_water_equivalent'
       flds(1)%ptr1r8 => model%coupling%SnowWaterEquiv

       !----------------------
       ! Snow depth
       !----------------------

       flds(2)%short_name = 'snow_depth'
       flds(2)%ptr1r8 => model%coupling%SnowDepth

       !----------------------
       ! Surface soil temperature
       !----------------------

       !flds(3)%short_name = 'soil_temperature'
       !flds(3)%nrec = model%noahmp%config%domain%NumSoilLayer
       !flds(3)%ptr2r8 => model%coupling%TemperatureSoilSnow

       !----------------------
       ! Surface soil moisture
       !----------------------

       flds(3)%short_name = 'soil_moisture'
       flds(3)%nrec = model%noahmp%config%domain%NumSoilLayer
       flds(3)%ptr2r8 => model%coupling%SoilMoisture
    end if

    !----------------------
    ! Read data 
    !----------------------

    call ReadFile(model, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Clean memory
    !----------------------

    if (allocated(flds)) deallocate(flds)

    call ESMF_LogWrite(subname//' done for '//trim(filename), ESMF_LOGMSG_INFO)

  end subroutine ReadIC

  !===============================================================================

  subroutine ReadFile(model, filename, flds, maskflag, notransfer, rh, rc)

    ! input/output variables
    type(model_type) , intent(inout) :: model
    character(len=*) , intent(in)    :: filename
    type(field_type) , intent(in)    :: flds(:)
    logical, optional, intent(in)    :: maskflag
    logical, optional, intent(in)    :: notransfer
    type(ESMF_RouteHandle), optional, intent(in) :: rh
    integer, optional, intent(inout) :: rc

    ! local variables
    logical                     :: amask, transfer_flag
    integer                     :: i, j, k, rank, fieldCount
    integer , pointer           :: ptr1i4(:)
    real(r4), pointer           :: ptr1r4(:)
    real(r8), pointer           :: ptr1r8(:)
    integer , pointer           :: ptr2i4(:,:)
    real(r4), pointer           :: ptr2r4(:,:)
    real(r8), pointer           :: ptr2r8(:,:)
    type(ESMF_RouteHandle)      :: rh_local
    type(ESMF_FieldBundle)      :: FBgrid, FBmesh
    type(ESMF_ArraySpec)        :: arraySpec
    type(ESMF_Field)            :: fgrid, fmesh, ftmp
    type(ESMF_TypeKind_Flag)    :: typekind
    character(len=cl)           :: fname
    character(len=cl), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname = trim(modName)//': (ReadFile) '
    !-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called for '//trim(filename), ESMF_LOGMSG_INFO)

    !----------------------
    ! Check mask flag
    !----------------------

    if (.not. present(maskflag)) then
       amask = .false.
    else
       amask = maskflag
    end if

    if (.not. allocated(model%domain%mask) .and. amask) then
       call ESMF_LogWrite(trim(subname)//' maskflag = .true. but model%domain%mask is not allocated yet! Skip applying mask ...', ESMF_LOGMSG_INFO)
    end if

    !----------------------
    ! Check grid to mesh transfer flag 
    !----------------------

    transfer_flag = .true.
    if (present(notransfer)) then
       if (notransfer) transfer_flag = .false.
    end if

    !----------------------
    ! Create field bundles
    !----------------------

    ! create empty field bundle on grid
    FBgrid = ESMF_FieldBundleCreate(name="fields_on_grid", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create empty field bundle on mesh
    if (transfer_flag) then
       FBmesh = ESMF_FieldBundleCreate(name="fields_on_mesh", rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! Loop over fields and add them to the field bundles
    !----------------------

    do i = 1, size(flds)
       ! 2d/r8 field (x,y)
       if (associated(flds(i)%ptr1r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          if (transfer_flag) then
             fmesh = ESMF_FieldCreate(model%domain%mesh, flds(i)%ptr1r8, meshloc=ESMF_MESHLOC_ELEMENT, &
                name=trim(flds(i)%short_name), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

       ! 2d/r4 field (x,y)
       else if (associated(flds(i)%ptr1r4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R4, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          if (transfer_flag) then
             fmesh = ESMF_FieldCreate(model%domain%mesh, flds(i)%ptr1r4, meshloc=ESMF_MESHLOC_ELEMENT, &
                name=trim(flds(i)%short_name), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

       ! 2d/i4 field (x,y)
       else if (associated(flds(i)%ptr1i4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_I4, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          if (transfer_flag) then
             fmesh = ESMF_FieldCreate(model%domain%mesh, flds(i)%ptr1i4, meshloc=ESMF_MESHLOC_ELEMENT, &
                name=trim(flds(i)%short_name), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

       ! 3d/r8 field (x,y,rec)
       else if (associated(flds(i)%ptr2r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/flds(i)%nrec/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          if (transfer_flag) then
             fmesh = ESMF_FieldCreate(model%domain%mesh, flds(i)%ptr2r8, meshloc=ESMF_MESHLOC_ELEMENT, &
                name=trim(flds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

       ! 3d/r4 field (x,y,rec)
       else if (associated(flds(i)%ptr2r4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R4, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/flds(i)%nrec/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          if (transfer_flag) then
             fmesh = ESMF_FieldCreate(model%domain%mesh, flds(i)%ptr2r4, meshloc=ESMF_MESHLOC_ELEMENT, &
                name=trim(flds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

       ! 3d/i4 field (x,y,rec)
       else if (associated(flds(i)%ptr2i4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_I4, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/flds(i)%nrec/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          if (transfer_flag) then
             fmesh = ESMF_FieldCreate(model%domain%mesh, flds(i)%ptr2i4, meshloc=ESMF_MESHLOC_ELEMENT, &
                name=trim(flds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if

       ! debug print
       call ESMF_LogWrite(trim(subname)//' adding '//trim(flds(i)%short_name)//' to FB', ESMF_LOGMSG_INFO)

       ! add it to the field bundle on grid
       call ESMF_FieldBundleAdd(FBgrid, [fgrid], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! add it to the field bundle on mesh
       if (transfer_flag) then
          call ESMF_FieldBundleAdd(FBmesh, [fmesh], rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    !----------------------
    ! Read data
    !----------------------

    call ESMF_FieldBundleRead(FBgrid, fileName=trim(filename), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !do i = 1, size(flds)
    !   ! get field from FB
    !   call ESMF_FieldBundleGet(FBgrid, fieldName=trim(flds(i)%short_name), field=fgrid, rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !   call ESMF_FieldWriteVTK(fgrid, trim(flds(i)%short_name), rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !end do

    !----------------------
    ! Create routehandle if it is not provided to transfer data from grid to mesh and move data
    !----------------------

    if (transfer_flag) then
       ! use provided routehandle
       if (present(rh)) then
          rh_local = rh
       ! create routehandle
       else
          call ESMF_FieldBundleRedistStore(FBgrid, FBmesh, routehandle=rh_local, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! move data from ESMF grid to mesh
       call ESMF_FieldBundleRedist(FBgrid, FBmesh, rh_local, checkflag=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! Apply masking
    !----------------------

    do i = 1, size(flds)
       ! get field from FB
       call ESMF_FieldBundleGet(FBmesh, fieldName=trim(flds(i)%short_name), field=fmesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! check its rank
       call ESMF_FieldGet(fmesh, rank=rank, typekind=typekind, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! query pointer and apply mask 
       if (rank .eq. 1) then
          if (amask .and. allocated(model%domain%mask)) then
             if (typekind == ESMF_TYPEKIND_R4) then
                call ESMF_FieldGet(fmesh, farrayPtr=ptr1r4, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                where (model%domain%mask < 1) ptr1r4 = 0.0_r4
                nullify(ptr1r4)
             else if (typekind == ESMF_TYPEKIND_R8) then
                call ESMF_FieldGet(fmesh, farrayPtr=ptr1r8, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                where (model%domain%mask < 1) ptr1r8 = 0.0_r8
                nullify(ptr1r8)
             else if (typekind == ESMF_TYPEKIND_I4) then
                call ESMF_FieldGet(fmesh, farrayPtr=ptr1i4, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                where (model%domain%mask < 1) ptr1i4 = 0
                nullify(ptr1i4)
             end if
          end if

          ! write field to VTK file
          if (model%nmlist%debug_level > 10) then
             call ESMF_FieldWriteVTK(fmesh, trim(flds(i)%short_name), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       else
          if (typekind == ESMF_TYPEKIND_R4) then
             call ESMF_FieldGet(fmesh, farrayPtr=ptr2r4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (amask .and. allocated(model%domain%mask)) then
                do j = 1, flds(i)%nrec
                   where (model%domain%mask < 1) ptr2r4(:,j) = 0.0_r4
                end do
             end if
          else if (typekind == ESMF_TYPEKIND_R8) then
             call ESMF_FieldGet(fmesh, farrayPtr=ptr2r8, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (amask .and. allocated(model%domain%mask)) then
                do j = 1, flds(i)%nrec
                   where (model%domain%mask < 1) ptr2r8(:,j) = 0.0_r8
                end do
             end if
          else if (typekind == ESMF_TYPEKIND_I4) then
             call ESMF_FieldGet(fmesh, farrayPtr=ptr2i4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (amask .and. allocated(model%domain%mask)) then
                do j = 1, flds(i)%nrec
                   where (model%domain%mask < 1) ptr2i4(:,j) = 0
                end do
             end if
          end if
 
          ! write field to VTK file, each record seperate file
          if (model%nmlist%debug_level > 10) then
             do j = 1, flds(i)%nrec
                ! file name
                write(fname, fmt='(A,I2.2)') trim(flds(i)%short_name)//'_rec_', j

                ! 2d/r4 field (element,layer)
                if (typekind == ESMF_TYPEKIND_R4) then
                   ! create temporary field and write it
                   ftmp = ESMF_FieldCreate(model%domain%mesh, typekind=ESMF_TYPEKIND_R4, &
                     name=trim(flds(i)%short_name), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return

                   ! get pointer and fill it
                   call ESMF_FieldGet(ftmp, localDe=0, farrayPtr=ptr1r4, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   ptr1r4(:) = ptr2r4(:,j)
                   nullify(ptr1r4)

                ! 2d/r8 field (element,layer)
                else if (typekind == ESMF_TYPEKIND_R8) then
                   ! create temporary field and write it
                   ftmp = ESMF_FieldCreate(model%domain%mesh, typekind=ESMF_TYPEKIND_R8, &
                     name=trim(flds(i)%short_name), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return

                   ! get pointer and fill it
                   call ESMF_FieldGet(ftmp, localDe=0, farrayPtr=ptr1r8, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   ptr1r8(:) = ptr2r8(:,j)
                   nullify(ptr1r8)

                ! 2d/i4 field (element,layer)
                else if (typekind == ESMF_TYPEKIND_I4) then
                   ! create temporary field and write it
                   ftmp = ESMF_FieldCreate(model%domain%mesh, typekind=ESMF_TYPEKIND_I4, &
                     name=trim(flds(i)%short_name), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return

                   ! get pointer and fill it
                   call ESMF_FieldGet(ftmp, localDe=0, farrayPtr=ptr1i4, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   ptr1i4(:) = ptr2i4(:,j)
                   nullify(ptr1i4)

                end if

                ! write
                call ESMF_FieldWriteVTK(ftmp, trim(fname), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return

                ! delete temporary field
                call ESMF_FieldDestroy(ftmp, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end do
          end if

          ! nullify pointers
          if (typekind == ESMF_TYPEKIND_R4) nullify(ptr2r4)
          if (typekind == ESMF_TYPEKIND_R8) nullify(ptr2r8)
          if (typekind == ESMF_TYPEKIND_I4) nullify(ptr2i4)
       end if
    end do

    !----------------------
    ! Empty FBs and destroy them 
    !----------------------

    ! FB grid
    call ESMF_FieldBundleGet(FBgrid, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(fieldNameList(fieldCount))
    call ESMF_FieldBundleGet(FBgrid, fieldNameList=fieldNameList, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do i = 1, fieldCount
       ! pull field from FB
       call ESMF_FieldBundleGet(FBgrid, fieldName=trim(fieldNameList(i)), field=ftmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! destroy field
       call ESMF_FieldDestroy(ftmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! remove field from FB
       call ESMF_FieldBundleRemove(FBgrid, fieldNameList=[trim(fieldNameList(i))], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do
    deallocate(fieldNameList)

    ! destroy grid FB
    call ESMF_FieldBundleDestroy(FBgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! FB mesh 
    call ESMF_FieldBundleGet(FBmesh, fieldCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(fieldNameList(fieldCount))
    call ESMF_FieldBundleGet(FBmesh, fieldNameList=fieldNameList, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do i = 1, fieldCount
       ! pull field from FB
       call ESMF_FieldBundleGet(FBmesh, fieldName=trim(fieldNameList(i)), field=ftmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! destroy field
       call ESMF_FieldDestroy(ftmp, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! remove field from FB
       call ESMF_FieldBundleRemove(FBmesh, fieldNameList=[trim(fieldNameList(i))], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do
    deallocate(fieldNameList)

    ! destroy grid FB
    call ESMF_FieldBundleDestroy(FBmesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Destroy route handle if it is created locally 
    !----------------------

    if (.not. present(rh)) then
       call ESMF_RouteHandleDestroy(rh_local, rc=rc)
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ReadFile

  !===============================================================================

  subroutine WriteFile(model, prefix, outflds, now_time, vm, rc)
    ! use statement
    use netcdf

    ! input/output variables
    type(model_type)  , target, intent(inout) :: model
    character(len=*)  , intent(in)    :: prefix 
    type(field_type)  , intent(in)    :: outflds(:)
    real(ESMF_KIND_R8), intent(in)    :: now_time
    type(ESMF_VM)     , intent(in)    :: vm
    integer, optional , intent(inout) :: rc

    ! local variables
    integer                     :: i, j, k
    integer                     :: localPet, rank, nlev, nfld, sub_str_indx
    integer                     :: ncerr, ncid, varid, dimid, max_indx
    real(r4), pointer           :: ptr2r4(:,:)
    real(r8), pointer           :: ptr2r8(:,:)
    integer , pointer           :: ptr2i4(:,:)
    integer , pointer           :: ptrMask(:,:)
    real(r4), pointer           :: ptr3r4(:,:,:)
    real(r8), pointer           :: ptr3r8(:,:,:)
    integer , pointer           :: ptr3i4(:,:,:)
    character(cl)               :: zaxis_name
    character(cl), allocatable  :: fieldNameList(:)
    character(cl)               :: filename, filename_tile
    logical                     :: isPresent, file_exists
    logical                     :: flag_soil_levels
    logical                     :: flag_snow_levels
    logical                     :: flag_snso_levels
    type(ESMF_RouteHandle)      :: rh_local
    type(ESMF_FieldBundle)      :: FBgrid, FBmesh
    type(ESMF_ArraySpec)        :: arraySpecI4, arraySpecR4, arraySpecR8
    type(ESMF_Field)            :: fgrid, fmesh
    type(ESMF_TypeKind_Flag)    :: typekind
    type(ESMF_LocStream)        :: locs
    type(ESMF_Grid)             :: grid
    type(ESMF_DistGrid)         :: distgrid
    character(len=*), parameter :: subname = trim(modName)//': (WriteFile) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------
    ! Query vm 
    !----------------------

    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Add file extension
    !----------------------

    if (model%domain%dtype == iMosaic) then 
       filename = trim(prefix)//'.tile*.nc'
    else
       filename = trim(prefix)//'.nc'
    end if
    call ESMF_LogWrite(trim(subname)//' called for '//trim(filename), ESMF_LOGMSG_INFO)

    !----------------------
    ! Remove existing file
    !----------------------

    if (localPet == 0) then
       if (model%domain%ntiles > 1) then
          ! loop over tiles
          do i = 1, model%domain%ntiles
             ! file name for tile
             sub_str_indx = index(trim(filename), "*", .true.)
             write(filename_tile, fmt='(a,i1,a)') trim(filename(:sub_str_indx-1)), i , trim(filename(sub_str_indx+1:))

             ! check file and delete
             inquire(file=trim(filename_tile), exist=file_exists)
             if (file_exists) then
                call ESMF_LogWrite(trim(subname)//' deleting '//trim(filename_tile), ESMF_LOGMSG_INFO)
                call system('rm -f '//trim(filename_tile))
             end if
          end do
       else
          ! check file and delete
          inquire(file=trim(filename), exist=file_exists)
          if (file_exists) then
             call ESMF_LogWrite(trim(subname)//'deleting '//trim(filename), ESMF_LOGMSG_INFO)
             call system('rm -f '//trim(filename))
          end if
       end if
    end if

    ! wait for sync
    call ESMF_VMBarrier(vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Create field bundles
    !----------------------

    ! create empty field bundle on grid
    FBgrid = ESMF_FieldBundleCreate(name="fields_on_grid", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create empty field bundle on mesh
    FBmesh = ESMF_FieldBundleCreate(name="fields_on_mesh", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Add metadata to grid
    !----------------------

    ! add coordinate dimensions: grid_xt and grid_yt
    call ESMF_AttributeAdd(model%domain%grid, convention="NetCDF", purpose="NOAHMP", attrList=(/"ESMF:gridded_dim_labels"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(model%domain%grid, convention="NetCDF", purpose="NOAHMP", name="ESMF:gridded_dim_labels", valueList=(/"grid_xt", "grid_yt"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Add coordinate variables (horizontal) 
    !----------------------

    ! create field for x coordinate on grid
    fgrid = ESMF_FieldCreate(model%domain%grid, farray=model%domain%long, indexflag=ESMF_INDEX_GLOBAL, &
       staggerloc=ESMF_STAGGERLOC_CENTER, name="grid_xt", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add coordinate attributes to the field
    call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/"cartesian_axis"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name="cartesian_axis", value="X", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/"long_name"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name="long_name", value="T-cell longitude", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/"units"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name="units", value="degrees_E", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! write to the file
    call ESMF_FieldWrite(fgrid, fileName=trim(filename), convention="NetCDF", purpose="NOAHMP", overwrite=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! destroy field
    call ESMF_FieldDestroy(fgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create field for y coordinate on grid
    fgrid = ESMF_FieldCreate(model%domain%grid, farray=model%domain%latg, indexflag=ESMF_INDEX_GLOBAL, &
       staggerloc=ESMF_STAGGERLOC_CENTER, name="grid_yt", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add coordinate attributes to the field
    call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/"cartesian_axis"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name="cartesian_axis", value="Y", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/"long_name"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name="long_name", value="T-cell latitude", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/"units"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name="units", value="degrees_N", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! write to the file
    call ESMF_FieldWrite(fgrid, fileName=trim(filename), convention="NetCDF", purpose="NOAHMP", &
       status=ESMF_FILESTATUS_OLD, overwrite=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! destroy field
    call ESMF_FieldDestroy(fgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Loop over fields and add them to the field bundles
    !----------------------

    ! init flags
    flag_soil_levels = .false.
    flag_snow_levels = .false.
    flag_snso_levels = .false.

    ! get number of fields
    max_indx = count(outflds(:)%id /= -999, 1)

    do i = 1, max_indx
       ! debug information
       if (model%nmlist%debug_level > 0) then
          call ESMF_LogWrite(trim(subname)//' adding '//trim(outflds(i)%short_name), ESMF_LOGMSG_INFO)
       end if

       ! set size and name of z-axis
       nlev = 0
       if (trim(outflds(i)%zaxis) == "z") then
          nlev = model%noahmp%config%domain%NumSoilLayer
          zaxis_name = "soil_levels"
          flag_soil_levels = .true.
       !else if (trim(outflds(i)%zaxis) == "z1") then
       !   nlev = abs(noahmp%static%lsnowl)+1
       !   zaxis_name = "snow_levels"
       !   flag_snow_levels = .true.
       !else if (trim(outflds(i)%zaxis) == "z2") then
       !   nlev = size(noahmp%nmlist%soil_level_nodes)+abs(noahmp%static%lsnowl)+1
       !   zaxis_name = "snso_levels"
       !   flag_snso_levels = .true.
       end if

       ! 2d/r8 field (x,y)
       if (associated(outflds(i)%ptr1r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpecR8, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpecR8, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=missing_r8, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(model%domain%mesh, outflds(i)%ptr1r8, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 2d/r4 field (x,y)
       else if (associated(outflds(i)%ptr1r4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpecR4, typekind=ESMF_TYPEKIND_R4, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpecR4, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=missing_r4, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(model%domain%mesh, outflds(i)%ptr1r4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 2d/i4 field (x,y)
       else if (associated(outflds(i)%ptr1i4)) then
          call ESMF_LogWrite(trim(subname)//' hoho', ESMF_LOGMSG_INFO)
          ! set field type
          call ESMF_ArraySpecSet(arraySpecI4, typekind=ESMF_TYPEKIND_I4, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpecI4, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=missing_i4, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(model%domain%mesh, outflds(i)%ptr1i4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 3d/r8 field (x,y,z)
       else if (associated(outflds(i)%ptr2r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpecR8, typekind=ESMF_TYPEKIND_R8, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpecR8, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/nlev/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=missing_r8, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(model%domain%mesh, outflds(i)%ptr2r8, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 3d/r4 field (x,y,z)
       else if (associated(outflds(i)%ptr2r4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpecR4, typekind=ESMF_TYPEKIND_R4, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpecR4, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/nlev/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=missing_r4, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(model%domain%mesh, outflds(i)%ptr2r4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 3d/i4 field (x,y,z)
       else if (associated(outflds(i)%ptr2i4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpecI4, typekind=ESMF_TYPEKIND_I4, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
          ! create field on grid
          fgrid = ESMF_FieldCreate(model%domain%grid, arraySpecI4, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/nlev/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=missing_i4, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(model%domain%mesh, outflds(i)%ptr2i4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_LogWrite(trim(subname)//" ERROR!!! The pointer is not associated with actual data.", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       ! write data on mesh for debugging
       !if (model%nmlist%debug_level > 10) then
       !   call ESMF_FieldWriteVTK(fmesh, trim(outflds(i)%short_name), rc=rc)
       !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
       !end if

       ! add long_name and units attributes to the field on grid
       call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'long_name'/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='long_name', value=trim(outflds(i)%long_name), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'units'/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='units', value=trim(outflds(i)%units), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! add vertical dimension name to the field on grid if it has ungridded dimension
       if (nlev > 0) then
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/"ESMF:ungridded_dim_labels"/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name="ESMF:ungridded_dim_labels", valueList=(/trim(zaxis_name)/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! add it to the field bundle on grid
       call ESMF_FieldBundleAdd(FBgrid, [fgrid], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! add it to the field bundle on mesh
       call ESMF_FieldBundleAdd(FBmesh, [fmesh], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !----------------------
    ! Add metadata to FB, global attributes
    !----------------------

    !call ESMF_AttributeAdd(FBgrid, convention="NetCDF", purpose="NOAHMP", &
    !  attrList=(/ "delt     ", &
    !              "idveg    ", & 
    !              "iopt_crs ", & 
    !              "iopt_btr ", & 
    !              "iopt_run ", &
    !              "iopt_sfc ", &
    !              "iopt_frz ", &
    !              "iopt_inf ", &
    !              "iopt_rad ", &
    !              "iopt_alb ", &
    !              "iopt_snf ", &
    !              "iopt_tbot", &
    !              "iopt_stc ", &
    !              "iopt_trs " /), rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="delt"     , value=noahmp%static%delt, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="idveg"    , value=noahmp%static%idveg, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_crs" , value=noahmp%static%iopt_crs, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_btr" , value=noahmp%static%iopt_btr, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_run" , value=noahmp%static%iopt_run, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_sfc" , value=noahmp%static%iopt_sfc, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_frz" , value=noahmp%static%iopt_frz, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_inf" , value=noahmp%static%iopt_inf, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_rad" , value=noahmp%static%iopt_rad, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_alb" , value=noahmp%static%iopt_alb, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_snf" , value=noahmp%static%iopt_snf, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_tbot", value=noahmp%static%iopt_tbot, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_stc" , value=noahmp%static%iopt_stc, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="iopt_trs" , value=noahmp%static%iopt_trs, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Create routehandle if it is not provided to transfer data from mesh to grid
    !----------------------

    !if (present(rh)) then
    !   rh_local = rh
    !else
       call ESMF_FieldBundleRedistStore(FBmesh, FBgrid, routehandle=rh_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !end if

    !----------------------
    ! Move data from ESMF grid to mesh
    !----------------------

    call ESMF_FieldBundleRedist(FBmesh, FBgrid, rh_local, checkflag=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Loop over fields on grid and apply mask and transpose if it is requested
    !----------------------

    ! check mask is in FB or not
    call ESMF_FieldBundleGet(FBgrid, fieldName="mask", isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (isPresent) then
       ! query mask information
       call ESMF_FieldBundleGet(FBgrid, fieldName="mask", field=fgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! get pointer from field
       call ESMF_FieldGet(fgrid, farrayPtr=ptrMask, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! loop over fields 
    do i = 1, max_indx
       ! get field from FB
       call ESMF_FieldBundleGet(FBgrid, fieldName=trim(outflds(i)%short_name), field=fgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! query type of the field
       call ESMF_FieldGet(fgrid, rank=rank, typekind=typekind, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! apply mask
       if (rank .eq. 2) then
          if (typekind == ESMF_TYPEKIND_R4) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr2r4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (associated(ptrMask)) then
                where(ptrMask < 1) ptr2r4 = missing_r4
             end if
             nullify(ptr2r4)
          else if (typekind == ESMF_TYPEKIND_R8) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr2r8, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (associated(ptrMask)) then
                where(ptrMask < 1) ptr2r8 = missing_r8
             end if
             nullify(ptr2r8)
          else if (typekind == ESMF_TYPEKIND_I4) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr2i4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (associated(ptrMask)) then
                where(ptrMask < 1) ptr2i4 = missing_i4
             end if
             nullify(ptr2i4)
          end if
       else if (rank .eq. 3) then
          if (typekind == ESMF_TYPEKIND_R4) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr3r4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (associated(ptrMask)) then
                do k = 1, ubound(ptr3r4, dim=3)
                   where(ptrMask < 1) ptr3r4(:,:,k) = missing_r4
                end do
             end if
             nullify(ptr3r4)
          else if (typekind == ESMF_TYPEKIND_R8) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr3r8, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (associated(ptrMask)) then
                do k = 1, ubound(ptr3r8, dim=3)
                   where(ptrMask < 1) ptr3r8(:,:,k) = missing_r8
                end do
             end if
             nullify(ptr3r8)
          else if (typekind == ESMF_TYPEKIND_I4) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr3i4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (associated(ptrMask)) then
                do k = 1, ubound(ptr3i4, dim=3)
                   where(ptrMask < 1) ptr3i4(:,:,k) = missing_i4
                end do
             end if
             nullify(ptr3i4)
          end if
       end if
    end do

    if (associated(ptrMask)) nullify(ptrMask)

    !----------------------
    ! Append fields to file
    !----------------------

    call ESMF_FieldBundleWrite(FBgrid, fileName=trim(filename), convention="NetCDF", purpose="NOAHMP", timeslice=1, overwrite=.true., status=ESMF_FILESTATUS_OLD, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMBarrier(vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Append coordinate variables (vertical and time)
    ! TODO: It uses serial netcdf library, once ESMF I/O layer extended to support
    ! adding coordinate variables and their attributes, this part will be removed.
    !----------------------

    ! only on the root pet
    if (localPet == 0) then
       ! loop over tiles
       do i = 1, model%domain%ntiles
          ! file name for tile
          sub_str_indx = index(trim(filename), "*", .true.)
          write(filename_tile, fmt='(a,i1,a)') trim(filename(:sub_str_indx-1)), i , trim(filename(sub_str_indx+1:))
       
          ! open file
          ncerr = nf90_open(trim(filename_tile), NF90_WRITE, ncid=ncid)
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

          !----------------------

          ! enter define mode
          ncerr = nf90_redef(ncid=ncid)
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

          ! check time dimension
          ncerr = nf90_inq_dimid(ncid, "time", dimid=dimid)
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

          ! if the time dimension does not exist, add it
          if (ncerr /= NF90_NOERR) then
             ncerr = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dimid=dimid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
          end if

          ! define variable
          ncerr = nf90_def_var(ncid, "time", NF90_DOUBLE, dimids=(/dimid/), varid=varid)
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

          ! add attributes
          ncerr = nf90_put_att(ncid, varid, "long_name", "valid time")
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
          ncerr = nf90_put_att(ncid, varid, "units", "seconds since "//trim(model%io%reference_date))
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
          ncerr = nf90_put_att(ncid, varid, "calendar", "gregorian")
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
          ncerr = nf90_put_att(ncid, varid, "cartesian_axis", "T")
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

          ! exit from define mode
          ncerr = nf90_enddef(ncid=ncid)
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

          ! add value to time
          ncerr = nf90_put_var(ncid, varid, values=[now_time])
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

          !----------------------

          if (flag_soil_levels) then
             ! enter define mode
             ncerr = nf90_redef(ncid=ncid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! check soil_levels dimension
             ncerr = nf90_inq_dimid(ncid, "soil_levels", dimid=dimid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! define variable
             ncerr = nf90_def_var(ncid, "soil_levels", NF90_DOUBLE, dimids=(/dimid/), varid=varid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! add attributes
             ncerr = nf90_put_att(ncid, varid, "long_name", "soil levels")
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
             ncerr = nf90_put_att(ncid, varid, "units", "meters")
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! exit from define mode
             ncerr = nf90_enddef(ncid=ncid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! add value to soil levels
             ncerr = nf90_put_var(ncid, varid, values=model%noahmp%config%domain%DepthSoilLayer)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
          end if

          !----------------------

          if (flag_snow_levels) then
             ! enter define mode
             ncerr = nf90_redef(ncid=ncid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! check snow_levels dimension
             ncerr = nf90_inq_dimid(ncid, "snow_levels", dimid=dimid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! define variable
             ncerr = nf90_def_var(ncid, "snow_levels", NF90_DOUBLE, dimids=(/dimid/), varid=varid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! add attributes
             ncerr = nf90_put_att(ncid, varid, "long_name", "snow levels")
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
             ncerr = nf90_put_att(ncid, varid, "units", "unitless")
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! exit from define mode
             ncerr = nf90_enddef(ncid=ncid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! add value to snow levels
             !ncerr = nf90_put_var(ncid, varid, values=(/(j*1.0d0,j=noahmp%static%lsnowl,0)/))
             !if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
          end if

          !----------------------

          if (flag_snso_levels) then
             ! enter define mode
             ncerr = nf90_redef(ncid=ncid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! check snso_levels dimension
             ncerr = nf90_inq_dimid(ncid, "snso_levels", dimid=dimid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
             
             ! define variable
             ncerr = nf90_def_var(ncid, "snso_levels", NF90_DOUBLE, dimids=(/dimid/), varid=varid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! add attributes
             ncerr = nf90_put_att(ncid, varid, "long_name", "snow soil levels")
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
             ncerr = nf90_put_att(ncid, varid, "units", "unitless")
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! exit from define mode
             ncerr = nf90_enddef(ncid=ncid)
             if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

             ! add value to snow levels
             !ncerr = nf90_put_var(ncid, varid, values=(/(j*1.0d0,j=noahmp%static%lsnowl,noahmp%nmlist%num_soil_levels)/))
             !if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
          end if

          ! close file
          ncerr = nf90_close(ncid=ncid)
          if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return
       end do
    end if

    !----------------------
    ! Empty FBs and destroy them 
    !----------------------

    ! loop over FB and remove fields
    call ESMF_FieldBundleGet(FBgrid, fieldCount=nfld, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! allocate field list
    allocate(fieldNameList(nfld))

    ! get field names
    call ESMF_FieldBundleGet(FBgrid, fieldNameList=fieldNameList, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do i = 1, nfld
       ! debug information
       if (model%nmlist%debug_level > 0) then
          call ESMF_LogWrite(trim(subname)//' removing '//trim(fieldNameList(i))//' from FBgrid and FBmesh', ESMF_LOGMSG_INFO)
       end if

       ! get field on grid
       call ESMF_FieldBundleGet(FBgrid, fieldName=trim(fieldNameList(i)), field=fgrid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! destroy it
       call ESMF_FieldDestroy(fgrid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! get field on mesh
       call ESMF_FieldBundleGet(FBmesh, fieldName=trim(fieldNameList(i)), field=fmesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! destroy it
       call ESMF_FieldDestroy(fmesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! deallocate temporrary field name array
    deallocate(fieldNameList)

    ! destroy field bundles
    call ESMF_FieldBundleDestroy(FBgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleDestroy(FBmesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! destroy routehandle if it is created locally
    !if (.not. present(rh)) then
       call ESMF_RouteHandleDestroy(rh_local, rc=rc)
    !end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine WriteFile

  !=============================================================================
  
  subroutine SetupWriteFields(model, outtype, rc)

    ! input/output variables
    type(model_type), target, intent(inout)   :: model 
    character(len=*), intent(in)              :: outtype
    integer         , optional, intent(inout) :: rc

    ! local variables
    character(len=*), parameter :: subname = trim(modName)//': (SetupWriteFields) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//' called ', ESMF_LOGMSG_INFO)

    !----------------------
    ! Prepare data structure to write the data 
    !----------------------

    ! history
    if (trim(outtype) == 'hist') then
       call FieldAdd('mask', 'Land-sea mask', '1', histflds, ptr1i4=model%domain%mask)

       ! mode = all
       if (trim(model%nmlist%OutputMode) == 'all') then
          ! forcing fields
          call FieldAdd('SpecHumidityRefHeight'  , 'specific humidity at reference height'                    , 'kg/kg'   , histflds, ptr1r8=model%forcing%SpecHumidityRefHeight)
          call FieldAdd('TemperatureAirRefHeight', 'air temperature at reference height'                      , 'K'       , histflds, ptr1r8=model%forcing%TemperatureAirRefHeight)
          call FieldAdd('WindEastwardRefHeight'  , 'wind speed in eastward direction at reference height'     , 'm s-1'   , histflds, ptr1r8=model%forcing%WindEastwardRefHeight)
          call FieldAdd('WindNorthwardRefHeight' , 'wind speed in northward direction at reference height'    , 'm s-1'   , histflds, ptr1r8=model%forcing%WindNorthwardRefHeight)
          call FieldAdd('TemperatureSoilBottom'  , 'bottom boundary condition for soil temperature'           , 'K'       , histflds, ptr1r8=model%forcing%TemperatureSoilBottom)
          call FieldAdd('RadSwDownRefHeight'     , 'downward shortwave radiation at reference height'         , 'W m-2'   , histflds, ptr1r8=model%forcing%RadSwDownRefHeight)
          call FieldAdd('RadLwDownRefHeight'     , 'downward longwave radiation at reference height'          , 'W m-2'   , histflds, ptr1r8=model%forcing%RadLwDownRefHeight)
          call FieldAdd('PressureAirRefHeight'   , 'air pressure at reference height'                         , 'Pa'      , histflds, ptr1r8=model%forcing%PressureAirRefHeight)
          call FieldAdd('PressureAirSurface'     , 'air pressure at the surface-atmosphere interface level'   , 'Pa'      , histflds, ptr1r8=model%forcing%PressureAirSurface)
          call FieldAdd('PrecipConvRefHeight'    , 'convective precipitation rate at reference height'        , 'mm s-1'  , histflds, ptr1r8=model%forcing%PrecipConvRefHeight)
          call FieldAdd('PrecipNonConvRefHeight' , 'non-convective precipitation rate at reference height'    , 'mm s-1'  , histflds, ptr1r8=model%forcing%PrecipNonConvRefHeight)
          call FieldAdd('PrecipShConvRefHeight'  , 'shallow convective precipitation rate at reference height', 'mm s-1'  , histflds, ptr1r8=model%forcing%PrecipShConvRefHeight)
          call FieldAdd('PrecipSnowRefHeight'    , 'snow rate at reference height'                            , 'mm s-1'  , histflds, ptr1r8=model%forcing%PrecipSnowRefHeight)
          call FieldAdd('PrecipGraupelRefHeight' , 'graupel rate at reference height'                         , 'mm s-1'  , histflds, ptr1r8=model%forcing%PrecipGraupelRefHeight)
          call FieldAdd('PrecipHailRefHeight'    , 'hail rate at reference height'                            , 'mm s-1'  , histflds, ptr1r8=model%forcing%PrecipHailRefHeight)

          ! static fields
          call FieldAdd('SoilType'               , 'soil type for each soil layer'                            , 'unitless', histflds, ptr1i4=model%coupling%SoilType)
          call FieldAdd('VegType'                , 'vegetation type'                                          , 'unitless', histflds, ptr1i4=model%coupling%VegType)
          call FieldAdd('RunoffSlopeType'        , 'underground runoff slope term'                            , 'unitless', histflds, ptr1i4=model%coupling%RunoffSlopeType)
          call FieldAdd('VegFracAnnMax'          , 'annual max vegetation fraction'                           , 'unitless', histflds, ptr1r8=model%coupling%VegFracAnnMax)
          call FieldAdd('SnowWaterEquiv'         , 'snow water equivalent'                                    , 'mm'      , histflds, ptr1r8=model%coupling%SnowWaterEquiv)
          call FieldAdd('SnowDepth'              , 'snow depth'                                               , 'mm'      , histflds, ptr1r8=model%coupling%SnowDepth)

          ! other fields
          call FieldAdd('CosSolarZenithAngle'    , 'cosine solar zenith angle'                                , '1'       , histflds, ptr1r8=model%coupling%CosSolarZenithAngle)

          ! 3d fields
          call FieldAdd('SoilMoisture'           , 'soil moisture'                                            , 'm3/m3'   , histflds, ptr2r8=model%coupling%SoilMoisture, zaxis="z")

       ! mode = mid
       else if (trim(model%nmlist%OutputMode) == 'mid') then

       ! mode = low
       else if (trim(model%nmlist%OutputMode) == 'low') then

       end if

    ! restart
    else if (trim(outtype) == 'rest') then

    endif

    call ESMF_LogWrite(trim(subname)//' done ', ESMF_LOGMSG_INFO)

  end subroutine SetupWriteFields 

  !===============================================================================

  subroutine FieldAdd(varName, varLName, varUnit, outflds, ptr1r4, ptr1r8, ptr1i4, ptr2r4, ptr2r8, ptr2i4, zAxis)

    ! input/output variables
    character(len=*)  , intent(in)           :: varName
    character(len=*)  , intent(in)           :: varUnit
    character(len=*)  , intent(in)           :: varLName
    type(field_type)  , intent(inout)        :: outflds(:) 
    real(r4), optional, pointer , intent(in) :: ptr1r4(:)
    real(r8), optional, pointer , intent(in) :: ptr1r8(:)
    integer , optional, pointer , intent(in) :: ptr1i4(:)
    real(r4), optional, pointer , intent(in) :: ptr2r4(:,:)
    real(r8), optional, pointer , intent(in) :: ptr2r8(:,:)
    integer , optional, pointer , intent(in) :: ptr2i4(:,:)
    character(len=*)  , optional, intent(in) :: zAxis

    ! local variables
    integer                     :: i, max_indx, indx
    logical                     :: found, restart
    character(len=*), parameter :: subname=trim(modName)//': (FieldAdd) '
    !-------------------------------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called for "//trim(varName), ESMF_LOGMSG_INFO)

    ! find out indices
    indx = 0
    found = .false.
    do i = 1, fldsMaxIO
       ! do not add to the list if it is found
       if (trim(outflds(i)%short_name) == trim(varName)) then
          indx = i
          found = .true.
          exit
       end if
    end do

    ! if it is a new entry, increment max_indx
    if (.not. found) then
       max_indx = count(outflds(:)%id /= -999, 1)
       indx = max_indx+1
       max_indx = max_indx+1
       if (max_indx > fldsMaxIO) then
          call ESMF_LogWrite(trim(subname)//": max_indx > fldsMaxIO could not add more variable! increase fldsMaxIO!", ESMF_LOGMSG_INFO)
          return
       end if
    end if

    ! add field metadata 
    outflds(indx)%id = indx
    outflds(indx)%short_name = trim(varName)
    outflds(indx)%units = trim(varUnit)
    outflds(indx)%long_name = trim(varLName)

    ! assign pointers
    if (present(ptr1r4)) then
       outflds(indx)%ptr1r4 => ptr1r4
    else if (present(ptr1r8)) then
       outflds(indx)%ptr1r8 => ptr1r8
    else if (present(ptr1i4)) then
       outflds(indx)%ptr1i4 => ptr1i4
    else if (present(ptr2r4)) then
       outflds(indx)%ptr2r4 => ptr2r4
    else if (present(ptr2r8)) then
       outflds(indx)%ptr2r8 => ptr2r8
    else if (present(ptr2i4)) then
       outflds(indx)%ptr2i4 => ptr2i4
    end if

    ! add extra metadata for the fields with z-axis
    if (present(zAxis)) then
       if (present(ptr2r4)) then
          outflds(indx)%zaxis = trim(zAxis)
          outflds(indx)%nlev = size(ptr2r4, dim=2)
       else if (present(ptr2r8)) then
          outflds(indx)%zaxis = trim(zAxis)
          outflds(indx)%nlev = size(ptr2r8, dim=2)
       else if (present(ptr2i4)) then
          outflds(indx)%zaxis = trim(zAxis)
          outflds(indx)%nlev = size(ptr2i4, dim=2)
      end if
    end if

    call ESMF_LogWrite(trim(subname)//' done for '//trim(varName), ESMF_LOGMSG_INFO)

  end subroutine FieldAdd

end module lnd_comp_io
