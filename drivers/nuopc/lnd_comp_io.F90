module lnd_comp_io

  ! This file contains I/O routines for the NoahMP land surface model

  use ESMF          , only : operator(==), operator(/=)
  use ESMF          , only : ESMF_FieldBundle, ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF          , only : ESMF_FieldBundleGet, ESMF_FieldBundleRead, ESMF_FieldBundleWrite
  use ESMF          , only : ESMF_FieldBundleRemove, ESMF_FieldBundleDestroy
  use ESMF          , only : ESMF_FieldBundleRedistStore, ESMF_FieldBundleRedist
  use ESMF          , only : ESMF_RouteHandleDestroy, ESMF_RouteHandle, ESMF_FieldWrite
  use ESMF          , only : ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet, ESMF_FieldWriteVTK
  use ESMF          , only : ESMF_FieldDestroy, ESMF_ArraySpec, ESMF_ArraySpecSet
  use ESMF          , only : ESMF_LogWrite, ESMF_FieldStatus_Flag, ESMF_GeomType_Flag
  use ESMF          , only : ESMF_Mesh, ESMF_MeshGet, ESMF_FieldBundleAddReplace
  use ESMF          , only : ESMF_KIND_R8, ESMF_TYPEKIND_R8
  use ESMF          , only : ESMF_KIND_R4, ESMF_TYPEKIND_R4
  use ESMF          , only : ESMF_KIND_I4, ESMF_TYPEKIND_I4
  use ESMF          , only : ESMF_STAGGERLOC_CENTER, ESMF_MESHLOC_ELEMENT
  use ESMF          , only : ESMF_INDEX_DELOCAL, ESMF_INDEX_GLOBAL
  use ESMF          , only : ESMF_MAXSTR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF          , only : ESMF_LOGMSG_ERROR, ESMF_LOGMSG_INFO
  use ESMF          , only : ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_COMPLETE
  use ESMF          , only : ESMF_AttributeAdd, ESMF_AttributeSet
  use ESMF          , only : ESMF_TypeKind_Flag

  use lnd_comp_types, only : noahmp_type, field_type
  use lnd_comp_kind , only : cl => shr_kind_cl
  use lnd_comp_kind , only : r4 => shr_kind_r4
  use lnd_comp_kind , only : r8 => shr_kind_r8
  use lnd_comp_kind , only : i4 => shr_kind_i4
  use lnd_comp_shr  , only : chkerr

  implicit none
  private

  public :: read_initial
  public :: read_restart
  public :: read_static
  public :: read_tiled_file
  public :: write_mosaic_output

  !-----------------------------------------------------------------------------
  ! Private module data
  !-----------------------------------------------------------------------------

  integer, parameter           :: dbug = 0
  integer, parameter           :: iswater = 17
  integer, parameter           :: max_num_variables = 200
  integer                      :: max_indx = 0
  type(field_type)             :: outflds(max_num_variables)
  character(len=1024)          :: msgString
  character(*), parameter      :: modName = "(lnd_comp_io)"
  character(len=*) , parameter :: u_FILE_u = __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine read_initial(noahmp, rc)

    ! input/output variables
    type(noahmp_type), target, intent(inout) :: noahmp
    integer,                   intent(inout) :: rc

    ! local variables
    type(field_type), allocatable    :: flds(:)
    character(len=cl)                :: filename
    character(len=*), parameter      :: subname=trim(modName)//':(read_initial) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------
    ! Prepare data structure to read the data 
    !----------------------

    if (noahmp%nmlist%ic_type == 'sfc') then
       ! input file name
       filename = trim(noahmp%nmlist%input_dir)//'sfc_data.tile#.nc'

       ! create field list
       allocate(flds(9))
       flds(1)%short_name = 'sheleg'; flds(1)%ptr1r8 => noahmp%init%snow_water_equivalent
       flds(2)%short_name = 'snwdph'; flds(2)%ptr1r8 => noahmp%init%snow_depth
       flds(3)%short_name = 'canopy'; flds(3)%ptr1r8 => noahmp%init%canopy_water
       flds(4)%short_name = 'tsea'  ; flds(4)%ptr1r8 => noahmp%init%skin_temperature
       flds(5)%short_name = 'stc'   ; flds(5)%nrec = noahmp%nmlist%num_soil_levels; flds(5)%ptr2r8 => noahmp%init%soil_temperature
       flds(6)%short_name = 'smc'   ; flds(6)%nrec = noahmp%nmlist%num_soil_levels; flds(6)%ptr2r8 => noahmp%init%soil_moisture
       flds(7)%short_name = 'slc'   ; flds(7)%nrec = noahmp%nmlist%num_soil_levels; flds(7)%ptr2r8 => noahmp%init%soil_liquid
       flds(8)%short_name = 'zorl'  ; flds(8)%ptr1r8 => noahmp%init%surface_roughness
       flds(9)%short_name = 'uustar'; flds(9)%ptr1r8 => noahmp%init%friction_velocity
    else 
       ! input file name
       write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.initial.tile#.nc'

       ! create field list
       allocate(flds(7))

       flds(1)%short_name = 'snow_water_equivalent'
       !flds(1)%input_type = ESMF_TYPEKIND_R4 
       flds(1)%ptr1r8 => noahmp%init%snow_water_equivalent

       flds(2)%short_name = 'snow_depth'
       !flds(2)%input_type = ESMF_TYPEKIND_R4
       flds(2)%ptr1r8 => noahmp%init%snow_depth

       flds(3)%short_name = 'canopy_water'
       !flds(3)%input_type = ESMF_TYPEKIND_R4
       flds(3)%ptr1r8 => noahmp%init%canopy_water

       flds(4)%short_name = 'skin_temperature'
       !flds(4)%input_type = ESMF_TYPEKIND_R4
       flds(4)%ptr1r8 => noahmp%init%skin_temperature

       flds(5)%short_name = 'soil_temperature'
       !flds(5)%input_type = ESMF_TYPEKIND_R4
       flds(5)%nrec = noahmp%nmlist%num_soil_levels
       flds(5)%ptr2r8 => noahmp%init%soil_temperature

       flds(6)%short_name = 'soil_moisture'
       !flds(6)%input_type = ESMF_TYPEKIND_R4
       flds(6)%nrec = noahmp%nmlist%num_soil_levels
       flds(6)%ptr2r8 => noahmp%init%soil_moisture

       flds(7)%short_name = 'soil_liquid'
       !flds(7)%input_type = ESMF_TYPEKIND_R4
       flds(7)%nrec = noahmp%nmlist%num_soil_levels
       flds(7)%ptr2r8 => noahmp%init%soil_liquid
    end if

    call ESMF_LogWrite(subname//' called for '//trim(filename), ESMF_LOGMSG_INFO)

    !----------------------
    ! Read file
    !----------------------

    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Clean memory
    !----------------------

    if (allocated(flds)) deallocate(flds)

    call ESMF_LogWrite(subname//' done for '//trim(filename), ESMF_LOGMSG_INFO)

  end subroutine read_initial

  !=============================================================================

  subroutine read_restart(noahmp, rc)

    ! input/output variables
    type(noahmp_type), target, intent(inout) :: noahmp
    integer,                   intent(inout) :: rc

    ! local variables
    type(field_type), allocatable    :: flds(:)
    character(len=cl)                :: filename
    character(len=*), parameter      :: subname=trim(modName)//':(read_restart) '
    !---------------------------------------------------------------------------

    ! input file name
    filename = trim(noahmp%nmlist%restart_dir)//trim(noahmp%nmlist%restart_file)
    call ESMF_LogWrite(subname//' called for '//trim(filename), ESMF_LOGMSG_INFO)

    ! create field list
    allocate(flds(126))


  end subroutine read_restart

  !=============================================================================

  subroutine read_static(noahmp, rc)

    ! input/output variables
    type(noahmp_type), target, intent(inout) :: noahmp
    integer,                   intent(inout) :: rc

    ! local variables
    type(field_type), allocatable :: flds(:)
    real(r4), target, allocatable :: tmpr4(:)
    real(r4), target, allocatable :: tmp2r4(:,:)
    character(len=cl)             :: filename
    real(r8)        , parameter   :: pi_8 = 3.14159265358979323846_r8
    character(len=*), parameter   :: subname=trim(modName)//':(read_static) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! allocate teemporary data structures
    !----------------------

    if (.not. allocated(tmpr4)) then
       allocate(tmpr4(noahmp%domain%begl:noahmp%domain%endl))
       tmpr4(:) = 0.0
    end if

    if (.not. allocated(tmp2r4)) then
       allocate(tmp2r4(noahmp%domain%begl:noahmp%domain%endl,12))
       tmp2r4(:,:) = 0.0
    end if

    !----------------------
    ! Set data sources
    !----------------------

    noahmp%static%isot = 1
    noahmp%static%ivegsrc = 1

    !----------------------
    ! Read latitude
    !----------------------

    allocate(flds(1))
    filename = trim(noahmp%nmlist%input_dir)//'oro_data.tile#.nc'
    flds(1)%short_name = 'geolat'; flds(1)%ptr1r8 => noahmp%model%xlatin
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    deallocate(flds)

    ! convert it to radian
    noahmp%model%xlatin(:) = noahmp%model%xlatin(:)*pi_8/180.0_r8

    !----------------------
    ! Read soil type
    !----------------------

    allocate(flds(1))
    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.soil_type.tile#.nc'
    flds(1)%short_name = 'soil_type'; flds(1)%ptr1r4 => tmpr4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%soiltyp = int(tmpr4)
    deallocate(flds)

    !----------------------
    ! Read vegetation type
    !----------------------

    allocate(flds(1))
    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.vegetation_type.tile#.nc'
    flds(1)%short_name = 'vegetation_type'; flds(1)%ptr1r4 => tmpr4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%vegtype = int(tmpr4)
    deallocate(flds)

    !----------------------
    ! Read slope type
    !----------------------

    allocate(flds(1))
    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.slope_type.tile#.nc'
    flds(1)%short_name = 'slope_type'; flds(1)%ptr1r4 => tmpr4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%slopetyp = int(tmpr4)
    deallocate(flds)

    !----------------------
    ! Read deep soil temperature
    !----------------------

    allocate(flds(1))
    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.substrate_temperature.tile#.nc'
    flds(1)%short_name = 'substrate_temperature'; flds(1)%ptr1r4 => tmpr4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%tg3 = dble(tmpr4)
    deallocate(flds)

    !----------------------
    ! Set emissivity
    !----------------------

    ! TODO: this needs to be option in nems.configure
    noahmp%model%emiss(:) = 0.95

    !----------------------
    ! Set albedo
    !----------------------

    ! TODO: this needs to be option in nems.configure
    noahmp%model%alb_monthly(:,:) = 0.25

    !----------------------
    ! Read maximum snow albedo
    !----------------------

    allocate(flds(1))
    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.maximum_snow_albedo.tile#.nc'
    flds(1)%short_name = 'maximum_snow_albedo'; flds(1)%ptr1r4 => tmpr4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%snoalb = dble(tmpr4)
    deallocate(flds)

    !----------------------
    ! Read vegetation greenness, monthly average, 12 months
    !----------------------

    allocate(flds(1))
    write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.vegetation_greenness.tile#.nc'
    flds(1)%short_name = 'vegetation_greenness'; flds(1)%nrec = 12; flds(1)%ptr2r4 => tmp2r4
    call read_tiled_file(noahmp, filename, flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%model%gvf_monthly(:,:) = dble(tmp2r4)
    noahmp%model%shdmin(:) = minval(noahmp%model%gvf_monthly(:,:), dim=2)
    noahmp%model%shdmax(:) = maxval(noahmp%model%gvf_monthly(:,:), dim=2)
    deallocate(flds)

    !----------------------
    ! Set dry
    !----------------------

    noahmp%model%dry(:) = .false.
    where(noahmp%model%vegtype(:) /= iswater) noahmp%model%dry(:) = .true.

    call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

  end subroutine read_static

  !=============================================================================

  subroutine read_tiled_file(noahmp, filename, flds, rh, rc)

    ! input/output variables
    type(noahmp_type), intent(inout) :: noahmp
    character(len=*),  intent(in)    :: filename
    type(field_type),  intent(in)    :: flds(:)
    type(ESMF_RouteHandle), optional, intent(in) :: rh
    integer, optional, intent(inout) :: rc

    ! local variables
    integer                     :: i, j, k, rank
    integer, pointer            :: ptr_i4(:)
    real(r4), pointer           :: ptr_r4(:)
    real(r8), pointer           :: ptr_r8(:)
    type(ESMF_RouteHandle)      :: rh_local
    type(ESMF_FieldBundle)      :: FBgrid, FBmesh
    type(ESMF_ArraySpec)        :: arraySpec
    type(ESMF_Field)            :: fgrid, fmesh, ftmp
    character(len=cl)           :: fname
    character(len=cl), allocatable :: fields(:)
    character(len=*), parameter :: subname = trim(modName)//': (read_tiled_file) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//' called for '//trim(filename), ESMF_LOGMSG_INFO)

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
    ! Loop over fields and add them to the field bundles
    !----------------------

    do i = 1, size(flds)
       ! 2d/r8 field (x,y)
       if (associated(flds(i)%ptr1r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(noahmp%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(noahmp%domain%mesh, flds(i)%ptr1r8, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 2d/r4 field (x,y)
       else if (associated(flds(i)%ptr1r4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R4, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(noahmp%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(noahmp%domain%mesh, flds(i)%ptr1r4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 3d/r8 field (x,y,rec)
       else if (associated(flds(i)%ptr2r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(noahmp%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/flds(i)%nrec/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(noahmp%domain%mesh, flds(i)%ptr2r8, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 3d/r4 field (x,y,rec)
       else if (associated(flds(i)%ptr2r4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R4, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(noahmp%domain%grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/flds(i)%nrec/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(noahmp%domain%mesh, flds(i)%ptr2r4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(flds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! debug print
       call ESMF_LogWrite(trim(subname)//' adding '//trim(flds(i)%short_name)//' to FB', ESMF_LOGMSG_INFO)

       ! add it to the field bundle on grid
       call ESMF_FieldBundleAdd(FBgrid, [fgrid], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! add it to the field bundle on mesh
       call ESMF_FieldBundleAdd(FBmesh, [fmesh], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !----------------------
    ! Read data
    !----------------------

    call ESMF_FieldBundleRead(FBgrid, fileName=trim(filename), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return     

    !----------------------
    ! Create routehandle if it is not provided to transfer data from grid to mesh
    !----------------------

    if (present(rh)) then
       rh_local = rh
    else
       call ESMF_FieldBundleRedistStore(FBgrid, FBmesh, routehandle=rh_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! Move data from ESMF grid to mesh
    !----------------------

    call ESMF_FieldBundleRedist(FBgrid, FBmesh, rh_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !call FB_diagnose(FBmesh, trim(subname), rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Clean memory
    !----------------------

    if (dbug > 0) then
       do i = 1, size(flds)
          ! get field from FB
          call ESMF_FieldBundleGet(FBmesh, fieldName=trim(flds(i)%short_name), field=fmesh, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! check its rank
          call ESMF_FieldGet(fmesh, rank=rank, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! TODO: ESMF_FieldWriteVTK() call does not support ungridded dimension
          ! The workaround is implemented in here but it would be nice to extend
          ! ESMF_FieldWriteVTK() call to handle it.  
          if (rank > 1) then
             ! create temporary field
             if (associated(flds(i)%ptr2r4)) then
                ftmp = ESMF_FieldCreate(noahmp%domain%mesh, typekind=ESMF_TYPEKIND_R4, &
                  name=trim(flds(i)%short_name), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return

                call ESMF_FieldGet(ftmp, localDe=0, farrayPtr=ptr_r4, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else if (associated(flds(i)%ptr2r8)) then
                ftmp = ESMF_FieldCreate(noahmp%domain%mesh, typekind=ESMF_TYPEKIND_R8, &
                  name=trim(flds(i)%short_name), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return

                call ESMF_FieldGet(ftmp, localDe=0, farrayPtr=ptr_r8, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else if (associated(flds(i)%ptr2i4)) then
                ftmp = ESMF_FieldCreate(noahmp%domain%mesh, typekind=ESMF_TYPEKIND_I4, &
                  name=trim(flds(i)%short_name), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return

                call ESMF_FieldGet(ftmp, localDe=0, farrayPtr=ptr_i4, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             ! write all record to seperate VTK file
             do j = 1, flds(i)%nrec
                if (associated(flds(i)%ptr2i4)) ptr_i4(:) = flds(i)%ptr2i4(:,j)
                if (associated(flds(i)%ptr2r4)) ptr_r4(:) = flds(i)%ptr2r4(:,j)
                if (associated(flds(i)%ptr2r8)) ptr_r8(:) = flds(i)%ptr2r8(:,j)
                write(fname, fmt='(A,I2.2)') trim(flds(i)%short_name)//'_rec', j
                call ESMF_FieldWriteVTK(ftmp, trim(fname), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end do

             ! delete temporary field
             call ESMF_FieldDestroy(ftmp, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             ! write field to VTK file
             call ESMF_FieldWriteVTK(fmesh, trim(flds(i)%short_name), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
    end if

    !----------------------
    ! Empty FBs and destroy them 
    !----------------------

    do i = 1, size(flds)
       ! remove field from both FBs
       call ESMF_FieldBundleRemove(FBgrid, fieldNameList=[trim(flds(i)%short_name)], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldBundleRemove(FBmesh, fieldNameList=[trim(flds(i)%short_name)], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    call ESMF_FieldBundleDestroy(FBgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleDestroy(FBmesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (.not. present(rh)) then
       call ESMF_RouteHandleDestroy(rh_local, rc=rc)
    end if

    call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

  end subroutine read_tiled_file

  !=============================================================================

  subroutine write_mosaic_output(filename, noahmp, now_time, rc)

    ! input/output variables
    character(len=*)  , intent(in)    :: filename
    type(noahmp_type) , target, intent(inout) :: noahmp
    real(ESMF_KIND_R8), intent(in)    :: now_time
    integer           , intent(inout) :: rc

    ! local variables
    integer                     :: i, k, rank, nlev
    real(r4), pointer           :: ptr2r4(:,:)
    real(r8), pointer           :: ptr2r8(:,:)
    integer , pointer           :: ptr2i4(:,:)
    integer , pointer           :: ptrMask(:,:)
    real(r4), pointer           :: ptr3r4(:,:,:)
    real(r8), pointer           :: ptr3r8(:,:,:)
    integer , pointer           :: ptr3i4(:,:,:)
    character(cl)               :: zaxis_name   
    type(ESMF_RouteHandle)      :: rh_local
    type(ESMF_FieldBundle)      :: FBgrid, FBmesh
    type(ESMF_ArraySpec)        :: arraySpecI4, arraySpecR4, arraySpecR8
    type(ESMF_Field)            :: fgrid, fmesh
    type(ESMF_TypeKind_Flag)    :: typekind
    character(len=*), parameter :: subname = trim(modName)//': (write_mosaic_output) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//' called for '//trim(filename), ESMF_LOGMSG_INFO)

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

    call ESMF_AttributeAdd(noahmp%domain%grid, convention="NetCDF", purpose="NOAHMP", attrList=(/"ESMF:gridded_dim_labels"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(noahmp%domain%grid, convention="NetCDF", purpose="NOAHMP", name="ESMF:gridded_dim_labels", valueList=(/"grid_xt", "grid_yt"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Prepare data structure to write the data 
    !----------------------

    call fld_add("ps"        , "surface pressure"                                                  , "Pa"     , ptr1r8=noahmp%forc%ps)
    call fld_add("u1"        , "u-component of wind"                                               , "m/s"    , ptr1r8=noahmp%model%u1)
    call fld_add("v1"        , "v-component of wind"                                               , "m/s"    , ptr1r8=noahmp%model%v1)
    call fld_add("t1"        , "forcing air temperature"                                           , "K"      , ptr1r8=noahmp%forc%t1)
    call fld_add("q1"        , "forcing specific humidity"                                         , "kg/kg"  , ptr1r8=noahmp%forc%q1)
    call fld_add("soiltyp"   , "soil type"                                                         , "1"      , ptr1i4=noahmp%model%soiltyp)
    call fld_add("vegtype"   , "vegetation type"                                                   , "1"      , ptr1i4=noahmp%model%vegtype)
    call fld_add("sigmaf"    , "green vegetation fraction"                                         , "1"      , ptr1r8=noahmp%model%sigmaf) 
    call fld_add("dlwflx"    , "forcing longwave downward flux"                                    , "W/m2"   , ptr1r8=noahmp%forc%dlwflx)
    call fld_add("dswsfc"    , "forcing shortwave downward flux"                                   , "W/m2"   , ptr1r8=noahmp%forc%dswsfc)
    call fld_add("snet"      , "forcing net shortwave flux"                                        , "W/m2"   , ptr1r8=noahmp%model%snet) 
    call fld_add("tg3"       , "deep soil temperature"                                             , "K"      , ptr1r8=noahmp%model%tg3)
    call fld_add("cm"        , "surface exchange coeff for momentum"                               , "m/s"    , ptr1r8=noahmp%model%cm)
    call fld_add("ch"        , "surface exchange coeff for heat and moisture"                      , "m/s"    , ptr1r8=noahmp%model%ch) 
    call fld_add("prsl1"     , "mean pressure at lowest model layer"                               , "Pa"     , ptr1r8=noahmp%model%prsl1)
    call fld_add("prslk1"    , "dimensionless Exner function at the lowest model layer"            , "1"      , ptr1r8=noahmp%model%prslk1)
    call fld_add("prslki"    , "Exner function ratio bt midlayer and interface at 1st layer"       , "1"      , ptr1r8=noahmp%model%prslki)
    call fld_add("prsik1"    , "dimensionless Exner function at the ground surface"                , "1"      , ptr1r8=noahmp%model%prsik1)
    call fld_add("zf"        , "height of bottom layer"                                            , "m"      , ptr1r8=noahmp%model%zf)
    call fld_add("mask"      , "land-sea mask"                                                     , "1"      , ptr1i4=noahmp%domain%mask)
    call fld_add("wind"      , "wind speed"                                                        , "m/s"    , ptr1r8=noahmp%forc%wind)
    call fld_add("slopetyp"  , "class of sfc slope"                                                , "1"      , ptr1i4=noahmp%model%slopetyp)
    call fld_add("shdmin"    , "min fractional coverage of green veg"                              , "1"      , ptr1r8=noahmp%model%shdmin)
    call fld_add("shdmax"    , "max fractional coverage of green veg"                              , "1"      , ptr1r8=noahmp%model%shdmax)
    call fld_add("snoalb"    , "upper bound on max albedo over deep snow"                          , "1"      , ptr1r8=noahmp%model%snoalb)
    call fld_add("sfalb"     , "mean sfc diffuse sw albedo"                                        , "1"      , ptr1r8=noahmp%model%sfalb)
    call fld_add("xlatin"    , "latitude"                                                          , "radian" , ptr1r8=noahmp%model%xlatin)
    call fld_add("xcoszin"   , "cosine of zenith angle"                                            , "degree" , ptr1r8=noahmp%model%xcoszin)
    call fld_add("garea"     , "area of the grid cell"                                             , "m2"     , ptr1r8=noahmp%domain%garea)
    call fld_add("rainn_mp"  , "microphysics non-convective precipitation"                         , "mm"     , ptr1r8=noahmp%model%rainn_mp)
    call fld_add("rainc_mp"  , "microphysics convective precipitation"                             , "mm"     , ptr1r8=noahmp%model%rainc_mp)
    call fld_add("snow_mp"   , "microphysics snow"                                                 , "mm"     , ptr1r8=noahmp%model%snow_mp)
    call fld_add("graupel_mp", "microphysics graupel"                                              , "mm"     , ptr1r8=noahmp%model%graupel_mp)
    call fld_add("ice_mp"    , "microphysics ice/hail"                                             , "mm"     , ptr1r8=noahmp%model%ice_mp)
    call fld_add("weasd"     , "water equivalent accumulated snow depth"                           , "mm"     , ptr1r8=noahmp%model%weasd)
    call fld_add("snwdph"    , "snow depth (water equiv) over land"                                , "m"      , ptr1r8=noahmp%model%snwdph)
    call fld_add("tskin"     , "ground surface skin temperature"                                   , "K"      , ptr1r8=noahmp%model%tskin)
    call fld_add("tprcp"     , "total precipitation"                                               , "mm"     , ptr1r8=noahmp%forc%tprcp)
    call fld_add("srflag"    , "snow/rain flag for precipitation"                                  , "1"      , ptr1r8=noahmp%model%srflag)
    call fld_add("smc"       ,"total soil moisture content"                                        ,"m3/m3"   , ptr2r8=noahmp%model%smc, zaxis="z")
    !!call fld_add("stc"       ,"soil temperature"                                                  ,"K"      , v2r8=noahmp%model%stc, zaxis="z")
    !!call fld_add("slc"       ,"liquid soil moisture"                                              ,"m3/m3"  , v2r8=noahmp%model%slc, zaxis="z")
    call fld_add("canopy"    , "canopy moisture content"                                           , "m"      , ptr1r8=noahmp%model%canopy)
    call fld_add("trans"     , "total plant transpiration"                                         , "m/2"    , ptr1r8=noahmp%model%trans)
    call fld_add("tsurf"     , "surface skin temperature (after iteration)"                        , "K"      , ptr1r8=noahmp%model%tsurf)
    call fld_add("zorl"      , "surface roughness"                                                 , "m"      , ptr1r8=noahmp%model%zorl)
    call fld_add("rb1"       , "bulk Richardson number at the surface over land"                   , "1"      , ptr1r8=noahmp%model%rb1)
    call fld_add("fm1"       , "Monin-Obukhov similarity function for momentum over land"          , "1"      , ptr1r8=noahmp%model%fm1)
    call fld_add("fh1"       , "Monin-Obukhov similarity function for heat over land"              , "1"      , ptr1r8=noahmp%model%fh1)
    call fld_add("ustar1"    , "surface friction velocity over land"                               , "m/s"    , ptr1r8=noahmp%model%ustar1)
    call fld_add("stress1"   , "surface wind stress over land"                                     , "m2/s2"  , ptr1r8=noahmp%model%stress1)
    call fld_add("fm101"     , "Monin-Obukhov similarity parameter for momentum at 10m over land"  , "1"      , ptr1r8=noahmp%model%fm101)
    call fld_add("fh21"      , "Monin-Obukhov similarity parameter for heat at 2m over land"       , "1"      , ptr1r8=noahmp%model%fh21)
    call fld_add("snowxy"    , "actual no. of snow layers"                                         , "1"      , ptr1r8=noahmp%model%snowxy)
    call fld_add("tvxy"      , "vegetation leaf temperature"                                       , "K"      , ptr1r8=noahmp%model%tvxy)
    call fld_add("tgxy"      , "bulk ground surface temperature"                                   , "K"      , ptr1r8=noahmp%model%tgxy)
    call fld_add("canicexy"  , "canopy-intercepted ice"                                            , "mm"     , ptr1r8=noahmp%model%canicexy)
    call fld_add("canliqxy"  , "canopy-intercepted liquid water"                                   , "mm"     , ptr1r8=noahmp%model%canliqxy)
    call fld_add("eahxy"     , "canopy air vapor pressure"                                         , "Pa"     , ptr1r8=noahmp%model%eahxy)
    call fld_add("tahxy"     , "canopy air temperature"                                            , "K"      , ptr1r8=noahmp%model%tahxy)
    call fld_add("cmxy"      , "bulk momentum drag coefficient"                                    , "m/s"    , ptr1r8=noahmp%model%cmxy)
    call fld_add("chxy"      , "bulk sensible heat exchange coefficient"                           , "m/s"    , ptr1r8=noahmp%model%chxy)
    call fld_add("fwetxy"    , "wetted or snowed fraction of the canopy"                           , "1"      , ptr1r8=noahmp%model%fwetxy)
    call fld_add("sneqvoxy"  , "snow mass at last time step"                                       , "mm"     , ptr1r8=noahmp%model%sneqvoxy)
    call fld_add("alboldxy"  , "snow albedo at last time step"                                     , "1"      , ptr1r8=noahmp%model%alboldxy)
    call fld_add("qsnowxy"   , "snowfall on the ground"                                            , "mm/s"   , ptr1r8=noahmp%model%qsnowxy)
    call fld_add("wslakexy"  , "lake water storage"                                                , "mm"     , ptr1r8=noahmp%model%wslakexy)
    call fld_add("zwtxy"     , "water table depth"                                                 , "m"      , ptr1r8=noahmp%model%zwtxy)
    call fld_add("waxy"      , "water in the aquifer"                                              , "mm"     , ptr1r8=noahmp%model%waxy)
    call fld_add("wtxy"      , "groundwater storage"                                               , "mm"     , ptr1r8=noahmp%model%wtxy)
    !!call fld_add("tsnoxy"    ,"temperature in surface snow"                                       ,"K"      , v2r8=noahmp%model%tsnoxy , zaxis="z1")
    !!call fld_add("zsnsoxy"   ,"depth from the top of the snow surface at the bottom of the layer" ,"m"      , v2r8=noahmp%model%zsnsoxy, zaxis="z2")
    !!call fld_add("snicexy"   ,"lwe thickness of ice in surface snow"                              ,"mm"     , v2r8=noahmp%model%snicexy, zaxis="z1")
    !!call fld_add("snliqxy"   ,"snow layer liquid water"                                           ,"mm"     , v2r8=noahmp%model%snliqxy, zaxis="z1")
    call fld_add("lfmassxy"  , "leaf mass"                                                         , "g/m2"   , ptr1r8=noahmp%model%lfmassxy)
    call fld_add("rtmassxy"  , "mass of fine roots"                                                , "g/m2"   , ptr1r8=noahmp%model%rtmassxy)
    call fld_add("stmassxy"  , "stem mas"                                                          , "g/m2"   , ptr1r8=noahmp%model%stmassxy)
    call fld_add("woodxy"    , "mass of wood incl woody roots"                                     , "g/m2"   , ptr1r8=noahmp%model%woodxy)
    call fld_add("stblcpxy"  , "stable carbon in deep soil"                                        , "g/m2"   , ptr1r8=noahmp%model%stblcpxy)
    call fld_add("fastcpxy"  , "short-lived carbon, shallow soil"                                  , "g/m2"   , ptr1r8=noahmp%model%fastcpxy)
    call fld_add("xlaixy"    , "leaf area index"                                                   , "1"      , ptr1r8=noahmp%model%xlaixy)
    call fld_add("xsaixy"    , "stem area index"                                                   , "1"      , ptr1r8=noahmp%model%xsaixy)
    call fld_add("taussxy"   , "snow age factor"                                                   , "1"      , ptr1r8=noahmp%model%taussxy)
    !!call fld_add("smoiseq"   ,"equilibrium soil water content"                                    ,"m3/m3"  , v2r8=noahmp%model%smoiseq, zaxis="z")
    call fld_add("smcwtdxy"  , "soil moisture content in the layer to the water table when deep"   , "mm"     , ptr1r8=noahmp%model%smcwtdxy)
    call fld_add("deeprechxy", "recharge to the water table when deep"                             , "1"      , ptr1r8=noahmp%model%deeprechxy)
    call fld_add("rechxy"    , "recharge to the water table (diagnostic)"                          , "1"      , ptr1r8=noahmp%model%rechxy)
    call fld_add("albdvis"   , "albedo - direct visible"                                           , "1"      , ptr1r8=noahmp%model%albdvis)
    call fld_add("albdnir"   , "albedo - direct NIR"                                               , "1"      , ptr1r8=noahmp%model%albdnir)
    call fld_add("albivis"   , "albedo - diffuse visible"                                          , "1"      , ptr1r8=noahmp%model%albivis)
    call fld_add("albinir"   , "albedo - diffuse NIR"                                              , "1"      , ptr1r8=noahmp%model%albinir)
    call fld_add("emiss"     , "surface emissivity"                                                , "1"      , ptr1r8=noahmp%model%emiss)
    call fld_add("sncovr1"   , "snow cover over land"                                              , "1"      , ptr1r8=noahmp%model%sncovr1)
    call fld_add("qsurf"     , "specific humidity at sfc"                                          , "kg/kg"  , ptr1r8=noahmp%model%qsurf)
    call fld_add("gflux"     , "soil heat flux"                                                    , "W/m2"   , ptr1r8=noahmp%model%gflux)
    call fld_add("drain"     , "subsurface runoff"                                                 , "mm/s"   , ptr1r8=noahmp%model%drain)
    call fld_add("evap"      , "evaporation from latent heat flux"                                 , "mm/s"   , ptr1r8=noahmp%model%evap)
    call fld_add("hflx"      , "sensible heat flux"                                                , "W/m2"   , ptr1r8=noahmp%model%hflx)
    call fld_add("ep"        , "potential evaporation"                                             , "W/m2"   , ptr1r8=noahmp%model%ep)
    call fld_add("runoff"    , "surface runoff"                                                    , "m/s"    , ptr1r8=noahmp%model%runoff)
    call fld_add("cmm"       , "cm * rho"                                                          , "m/s"    , ptr1r8=noahmp%model%cmm)
    call fld_add("chh"       , "ch * rho"                                                          , "kg/m2/s", ptr1r8=noahmp%model%chh)
    call fld_add("evbs"      , "direct soil evaporation"                                           , "m/s"    , ptr1r8=noahmp%model%evbs)
    call fld_add("evcw"      , "canopy water evaporation"                                          , "m/s"    , ptr1r8=noahmp%model%evcw)
    call fld_add("sbsno"     , "sublimation/deposit from snopack"                                  , "m/s"    , ptr1r8=noahmp%model%sbsno)
    call fld_add("pah"       , "precipitation advected heat - total"                               , "W/m2"   , ptr1r8=noahmp%model%pah)
    call fld_add("ecan"      , "evaporation of intercepted water"                                  , "kg/m2/s", ptr1r8=noahmp%model%ecan)
    call fld_add("etran"     , "transpiration rate"                                                , "kg/m2/s", ptr1r8=noahmp%model%etran)
    call fld_add("edir"      , "soil surface evaporation rate"                                     , "kg/m2/s", ptr1r8=noahmp%model%edir)
    call fld_add("snowc"     , "fractional snow cover"                                             , "1"      , ptr1r8=noahmp%model%snowc)
    call fld_add("stm"       , "total soil column moisture content"                                , "m"      , ptr1r8=noahmp%model%stm)
    call fld_add("snohf"     , "snow/freezing-rain latent heat flux"                               , "W/m2"   , ptr1r8=noahmp%model%snohf)
    call fld_add("smcwlt2"   , "dry soil moisture threshold"                                       , "m3/m3"  , ptr1r8=noahmp%model%smcwlt2)
    call fld_add("smcref2"   , "soil moisture threshold"                                           , "m3/m3"  , ptr1r8=noahmp%model%smcref2)
    call fld_add("wet1"      , "normalized soil wetness"                                           , "1"      , ptr1r8=noahmp%model%wet1)
    call fld_add("t2mmp"     , "combined T2m from tiles"                                           , "K"      , ptr1r8=noahmp%model%t2mmp)
    call fld_add("q2mp"      , "combined q2m from tiles"                                           , "kg/kg"  , ptr1r8=noahmp%model%q2mp)
    call fld_add("zvfun"     , "function of surface roughness length and green vegetation fraction", "1"      , ptr1r8=noahmp%model%zvfun)
    !!call fld_add("rho"       ,"density"                                                           ,"kg/m3"  , v1r8=noahmp%model%rho)
    !!call fld_add("hgt"       ,"forcing height"                                                    ,"m"      , v1r8=noahmp%forc%hgt)
    !!call fld_add("pblh"      ,"height of pbl"                                                     ,"m"      , v1r8=noahmp%model%pblh)

    !----------------------
    ! Loop over fields and add them to the field bundles
    !----------------------

    do i = 1, max_indx
       ! query z axis
       nlev = 0
       if (trim(outflds(i)%zaxis) == "z") then
          nlev = size(noahmp%nmlist%soil_level_nodes)
          zaxis_name = "soil_levels"
       end if

       ! 2d/r8 field (x,y)
       if (associated(outflds(i)%ptr1r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpecR8, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(noahmp%domain%grid, arraySpecR8, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=1.0d20, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(noahmp%domain%mesh, outflds(i)%ptr1r8, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 2d/r4 field (x,y)
       else if (associated(outflds(i)%ptr1r4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpecR4, typekind=ESMF_TYPEKIND_R4, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(noahmp%domain%grid, arraySpecR4, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=1.0e20, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(noahmp%domain%mesh, outflds(i)%ptr1r4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 2d/i4 field (x,y)
       else if (associated(outflds(i)%ptr1i4)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpecI4, typekind=ESMF_TYPEKIND_I4, rank=2, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(noahmp%domain%grid, arraySpecI4, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=-999, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(noahmp%domain%mesh, outflds(i)%ptr1i4, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! 3d/r8 field (x,y,z)
       else if (associated(outflds(i)%ptr2r8)) then
          ! set field type
          call ESMF_ArraySpecSet(arraySpecR8, typekind=ESMF_TYPEKIND_R8, rank=3, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on grid
          fgrid = ESMF_FieldCreate(noahmp%domain%grid, arraySpecR8, staggerloc=ESMF_STAGGERLOC_CENTER, &
             indexflag=ESMF_INDEX_GLOBAL, name=trim(outflds(i)%short_name), ungriddedLbound=(/1/), &
             ungriddedUbound=(/nlev/), gridToFieldMap=(/1,2/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! add missing value attribute to the field
          call ESMF_AttributeAdd(fgrid, convention="NetCDF", purpose="NOAHMP", attrList=(/'missing_value'/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_AttributeSet(fgrid, convention="NetCDF", purpose="NOAHMP", name='missing_value', value=1.0d20, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! create field on mesh
          fmesh = ESMF_FieldCreate(noahmp%domain%mesh, outflds(i)%ptr2r8, meshloc=ESMF_MESHLOC_ELEMENT, &
             name=trim(outflds(i)%short_name), gridToFieldMap=(/1/), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

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
    !    attrList=(/ "time               ", &
    !                "time:long_name     ", &
    !                "time:units         ", &
    !                "time:cartesian_axis", &
    !                "time:calendar_type ", &
    !                "time:calendar      " /), rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !call ESMF_AttributeSet(FBgrid, convention="NetCDF", purpose="NOAHMP", name="time", value=0, rc=rc)
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

    call ESMF_FieldBundleRedist(FBmesh, FBgrid, rh_local, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Loop over fields on grid and apply mask
    !----------------------

    ! get mask information
    call ESMF_FieldBundleGet(FBgrid, fieldName="mask", field=fgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(fgrid, farrayPtr=ptrMask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! apply mask
    do i = 1, max_indx
       call ESMF_FieldBundleGet(FBgrid, fieldName=trim(outflds(i)%short_name), field=fgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldGet(fgrid, rank=rank, typekind=typekind, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (rank .eq. 2) then
          if (typekind == ESMF_TYPEKIND_R4) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr2r4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             where(ptrMask < 1) ptr2r4 = 1.0e20
             nullify(ptr2r4)
          else if (typekind == ESMF_TYPEKIND_R8) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr2r8, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             where(ptrMask < 1) ptr2r8 = 1.0d20
             nullify(ptr2r8)
          else if (typekind == ESMF_TYPEKIND_I4) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr2i4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             where(ptrMask < 1) ptr2i4 = -999
             nullify(ptr2i4)
          end if
       else if (rank .eq. 3) then
          if (typekind == ESMF_TYPEKIND_R4) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr3r4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do k = 1, ubound(ptr3r4, dim=3)
                where(ptrMask < 1) ptr3r4(:,:,k) = 1.0e20
             end do
             nullify(ptr3r4)
          else if (typekind == ESMF_TYPEKIND_R8) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr3r8, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do k = 1, ubound(ptr3r8, dim=3)
                where(ptrMask < 1) ptr3r8(:,:,k) = 1.0d20
             end do
             nullify(ptr3r8)
          else if (typekind == ESMF_TYPEKIND_I4) then
             call ESMF_FieldGet(fgrid, farrayPtr=ptr3i4, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do k = 1, ubound(ptr3i4, dim=3)
                where(ptrMask < 1) ptr3i4(:,:,k) = -999
             end do
             nullify(ptr3i4)
          end if
       end if
    end do

    !----------------------
    ! Write data
    !----------------------

    call ESMF_FieldBundleWrite(FBgrid, fileName=trim(filename), convention="NetCDF", purpose="NOAHMP", overwrite=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Empty FBs and destroy them 
    !----------------------

    do i = 1, max_indx
       ! remove field from both FBs
       call ESMF_FieldBundleRemove(FBgrid, fieldNameList=[trim(outflds(i)%short_name)], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldBundleRemove(FBmesh, fieldNameList=[trim(outflds(i)%short_name)], rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    call ESMF_FieldBundleDestroy(FBgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleDestroy(FBmesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !if (.not. present(rh)) then
       call ESMF_RouteHandleDestroy(rh_local, rc=rc)
    !end if

    call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

  end subroutine write_mosaic_output

  !===============================================================================

  subroutine fld_add(varName, varLName, varUnit, ptr1r4, ptr1r8, ptr1i4, ptr2r4, ptr2r8, ptr2i4, zAxis)

    ! input/output variables
    character(len=*)  , intent(in)           :: varName
    character(len=*)  , intent(in)           :: varUnit
    character(len=*)  , intent(in)           :: varLName
    real(r4), optional, pointer , intent(in) :: ptr1r4(:)
    real(r8), optional, pointer , intent(in) :: ptr1r8(:)
    integer , optional, pointer , intent(in) :: ptr1i4(:)
    real(r4), optional, pointer , intent(in) :: ptr2r4(:,:)
    real(r8), optional, pointer , intent(in) :: ptr2r8(:,:)
    integer , optional, pointer , intent(in) :: ptr2i4(:,:)
    character(len=*)  , optional, intent(in) :: zAxis

    ! local variables
    integer                     :: i, indx
    logical                     :: found, restart
    character(len=*), parameter :: subname=trim(modName)//': (fld_add) '
    !-------------------------------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called for "//trim(varName), ESMF_LOGMSG_INFO)

    ! find out indices
    indx = 0
    found = .false.
    do i = 1, max_num_variables
       ! do not add to the list if it is found
       if (trim(outflds(i)%short_name) == trim(varName)) then
          indx = i
          found = .true.
          exit
       end if
    end do

    ! if it is a new entry, increment max_indx
    if (.not. found) then
       indx = max_indx+1
       max_indx = max_indx+1
       if (max_indx > max_num_variables) then
          print*, "max_indx > max_num_variables could not add more variable! increase max_num_variables ..."
          return
       end if
    end if

    ! add field metadata 
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
 
    call ESMF_LogWrite(trim(subname)//' done for '//trim(varName), ESMF_LOGMSG_INFO)

  end subroutine fld_add  

  !=============================================================================

  subroutine FB_diagnose(FB, string, rc)

    ! input/output variables
    type(ESMF_FieldBundle) , intent(inout)        :: FB
    character(len=*)       , intent(in), optional :: string
    integer                , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR), pointer :: lfieldnamelist(:)
    character(len=CL)               :: lstring
    real(R8), pointer               :: dataPtr1d(:)
    real(R8), pointer               :: dataPtr2d(:,:)
    type(ESMF_Field)                :: lfield
    character(len=*), parameter     :: subname='(FB_diagnose)'
    !---------------------------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) then
       lstring = trim(string) // ' '
    endif

    ! Determine number of fields in field bundle and allocate memory for lfieldnamelist
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    ! Get the fields in the field bundle
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! For each field in the bundle, get its memory location and print out the field
    do n = 1, fieldCount
       call ESMF_FieldBundleGet(FB, fieldName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call Field_GetFldPtr(lfield, fldptr1=dataptr1d, fldptr2=dataptr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data

       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), " no data"
          endif

       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif

       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    ! Deallocate memory
    deallocate(lfieldnamelist)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine FB_diagnose

  !=============================================================================

  subroutine Field_GetFldPtr(field, fldptr1, fldptr2, rank, abort, rc)

    ! ----------------------------------------------
    ! for a field, determine rank and return fldptr1 or fldptr2
    ! abort is true by default and will abort if fldptr is not yet allocated in field
    ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_Field)  , intent(in)              :: field
    real(R8), pointer , intent(inout), optional :: fldptr1(:)
    real(R8), pointer , intent(inout), optional :: fldptr2(:,:)
    integer           , intent(out)  , optional :: rank
    logical           , intent(in)   , optional :: abort
    integer           , intent(out)  , optional :: rc

    ! local variables
    type(ESMF_Mesh)             :: lmesh
    integer                     :: lrank, nnodes, nelements
    logical                     :: labort
    type(ESMF_GeomType_Flag)    :: geomtype
    type(ESMF_FieldStatus_Flag) :: status
    character(len=*), parameter :: subname='(Field_GetFldPtr)'
    !---------------------------------------------------------------------------

    if (dbug > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    if (.not.present(rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR rc not present ", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    rc = ESMF_SUCCESS

    labort = .true.
    if (present(abort)) then
      labort = abort
    endif
    lrank = -99

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
      lrank = 0
      if (labort) then
        call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO)
        rc = ESMF_FAILURE
        return
      else
        call ESMF_LogWrite(trim(subname)//": WARNING data not allocated ", ESMF_LOGMSG_INFO)
      endif
    else

      call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      if (geomtype == ESMF_GEOMTYPE_GRID) then
        call ESMF_FieldGet(field, rank=lrank, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      elseif (geomtype == ESMF_GEOMTYPE_MESH) then
        call ESMF_FieldGet(field, rank=lrank, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        if (nnodes == 0 .and. nelements == 0) lrank = 0
      else
         call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", ESMF_LOGMSG_INFO)
        rc = ESMF_FAILURE
        return
      endif ! geomtype

      if (lrank == 0) then
         call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", &
              ESMF_LOGMSG_INFO)

      elseif (lrank == 1) then
        if (.not.present(fldptr1)) then
           call ESMF_LogWrite(trim(subname)//": ERROR missing rank=1 array ", &
                ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
        endif
        call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

      elseif (lrank == 2) then
        if (.not.present(fldptr2)) then
           call ESMF_LogWrite(trim(subname)//": ERROR missing rank=2 array ", &
                ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
        endif
        call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

      else
         call ESMF_LogWrite(trim(subname)//": ERROR in rank ", &
              ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
        rc = ESMF_FAILURE
        return
      endif

    endif  ! status

    if (present(rank)) then
      rank = lrank
    endif

    if (dbug > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine Field_GetFldPtr

end module lnd_comp_io
