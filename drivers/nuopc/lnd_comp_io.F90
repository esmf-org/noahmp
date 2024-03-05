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
  use ESMF          , only : ESMF_AttributeAdd, ESMF_AttributeSet, ESMF_Grid
  use ESMF          , only : ESMF_LocStream, ESMF_LocStreamCreate, ESMF_LocStreamDestroy
  use ESMF          , only : ESMF_GridCreate, ESMF_FILESTATUS_OLD, ESMF_TypeKind_Flag
  use ESMF          , only : ESMF_DistGrid, ESMF_DistGridCreate
  use ESMF          , only : ESMF_VMBarrier, ESMF_VM

  use lnd_comp_types, only : model_type
  use lnd_comp_types, only : field_type
  use lnd_comp_kind , only : cl => shr_kind_cl
  use lnd_comp_kind , only : r4 => shr_kind_r4
  use lnd_comp_kind , only : r8 => shr_kind_r8
  use lnd_comp_kind , only : i4 => shr_kind_i4
  use lnd_comp_shr  , only : chkerr
  use lnd_comp_shr  , only : chkerrnc

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(*), parameter      :: modName = "(lnd_comp_io)"
  character(len=*) , parameter :: u_FILE_u = __FILE__

!===============================================================================
contains
!===============================================================================

  integer function get_num_tiles(filename)

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
    ncerr = nf90_inquire_dimension(ncid, dimid, len=get_num_tiles)
    if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

    ! close file
    ncerr = nf90_close(ncid)
    if (ChkErrNc(ncerr,__LINE__,u_FILE_u)) return

  end function get_num_tiles

  !===============================================================================

  subroutine read_tiled_file(model, filename, flds, maskflag, notransfer, rh, rc)

    ! input/output variables
    type(model_type), intent(inout) :: model
    character(len=*),  intent(in)    :: filename
    type(field_type),  intent(in)    :: flds(:)
    logical, optional, intent(in)    :: maskflag
    logical, optional, intent(in)    :: notransfer
    type(ESMF_RouteHandle), optional, intent(in) :: rh
    integer, optional, intent(inout) :: rc

    ! local variables
    logical                     :: amask, transfer_flag
    integer                     :: i, j, k, rank, fieldCount
    integer, pointer            :: ptr1i4(:)
    real(r4), pointer           :: ptr1r4(:)
    real(r8), pointer           :: ptr1r8(:)
    integer, pointer            :: ptr2i4(:,:)
    real(r4), pointer           :: ptr2r4(:,:)
    real(r8), pointer           :: ptr2r8(:,:)
    type(ESMF_RouteHandle)      :: rh_local
    type(ESMF_FieldBundle)      :: FBgrid, FBmesh
    type(ESMF_ArraySpec)        :: arraySpec
    type(ESMF_Field)            :: fgrid, fmesh, ftmp
    type(ESMF_TypeKind_Flag)    :: typekind
    character(len=cl)           :: fname
    character(len=cl), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname = trim(modName)//': (read_tiled_file) '
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
       call ESMF_FieldBundleRedist(FBgrid, FBmesh, rh_local, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! Apply masking
    !----------------------

    return

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
          if (dbug > 0) then
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
          if (dbug > 0) then
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

  end subroutine read_tiled_file

end module lnd_comp_io
