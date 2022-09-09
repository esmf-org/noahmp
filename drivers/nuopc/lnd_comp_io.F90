module lnd_comp_io

  ! This file contains I/O routines for the NoahMP land surface model

  use ESMF          , only : operator(==), operator(/=)
  use ESMF          , only : ESMF_FieldBundle, ESMF_FieldBundleCreate, ESMF_FieldBundleAdd
  use ESMF          , only : ESMF_FieldBundleGet, ESMF_FieldBundleRead, ESMF_FieldBundleWrite
  use ESMF          , only : ESMF_FieldBundleRemove, ESMF_FieldBundleDestroy
  use ESMF          , only : ESMF_FieldBundleRedistStore, ESMF_FieldBundleRedist
  use ESMF          , only : ESMF_RouteHandleIsCreated
  use ESMF          , only : ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet, ESMF_FieldWriteVTK
  use ESMF          , only : ESMF_ArraySpec, ESMF_ArraySpecSet
  use ESMF          , only : ESMF_LogWrite, ESMF_FieldStatus_Flag, ESMF_GeomType_Flag
  use ESMF          , only : ESMF_Mesh, ESMF_MeshGet
  use ESMF          , only : ESMF_KIND_R8, ESMF_TYPEKIND_R8
  use ESMF          , only : ESMF_STAGGERLOC_CENTER, ESMF_MESHLOC_ELEMENT
  use ESMF          , only : ESMF_INDEX_DELOCAL, ESMF_INDEX_GLOBAL
  use ESMF          , only : ESMF_MAXSTR, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF          , only : ESMF_LOGMSG_ERROR, ESMF_LOGMSG_INFO
  use ESMF          , only : ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_COMPLETE

  use lnd_comp_types, only : noahmp_type, field_type
  use lnd_comp_kind , only : cl => shr_kind_cl
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

  integer, parameter           :: dbug = 1
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
       allocate(flds(4))
       flds(1)%short_name = 'sheleg'; flds(1)%ptr1r8 => noahmp%init%snow_water_equivalent
       flds(2)%short_name = 'snwdph'; flds(2)%ptr1r8 => noahmp%init%snow_depth
       flds(3)%short_name = 'canopy'; flds(3)%ptr1r8 => noahmp%init%canopy_water
       flds(4)%short_name = 'tsea'  ; flds(4)%ptr1r8 => noahmp%init%skin_temperature
       !flds(5)%short_name = 'stc'   ; flds(5)%ptr1r8 => noahmp%init%soil_temperature
       !flds(6)%short_name = 'smc'   ; flds(6)%ptr1r8 => noahmp%init%soil_moisture
       !flds(7)%short_name = 'slc'   ; flds(7)%ptr1r8 => noahmp%init%soil_liquid
    else 
       ! input file name
       write(filename, fmt="(A,I0,A)") trim(noahmp%nmlist%input_dir)//'C', maxval(noahmp%domain%nit), '.initial.tile#.nc'

       ! create field list
       allocate(flds(4))
       flds(1)%short_name = 'snow_water_equivalent'; flds(1)%ptr1r8 => noahmp%init%snow_water_equivalent
       flds(2)%short_name = 'snow_depth'           ; flds(2)%ptr1r8 => noahmp%init%snow_depth
       flds(3)%short_name = 'canopy_water'         ; flds(3)%ptr1r8 => noahmp%init%canopy_water
       flds(4)%short_name = 'skin_temperature'     ; flds(4)%ptr1r8 => noahmp%init%skin_temperature
       !flds(5)%short_name = 'soil_temperature'     ; flds(5)%ptr1r8 => noahmp%init%soil_temperature
       !flds(6)%short_name = 'soil_moisture'        ; flds(6)%ptr1r8 => noahmp%init%soil_moisture
       !flds(7)%short_name = 'soil_liquid'          ; flds(7)%ptr1r8 => noahmp%init%soil_liquid
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
    type(noahmp_type), intent(inout) :: noahmp
    integer          , intent(inout) :: rc

    ! local variables
    integer                          :: nt
    character(len=cl)                :: filename
    real(ESMF_KIND_R8), pointer      :: ptr(:,:,:)
    type(ESMF_Field)                 :: field
    character(len=*), parameter      :: subname=trim(modName)//':(read_restart) '
    !---------------------------------------------------------------------------

  end subroutine read_restart

  !===============================================================================
  !=============================================================================

  subroutine read_static(noahmp, rc)

    ! input/output variables
    type(noahmp_type), intent(inout) :: noahmp
    integer          , intent(inout) :: rc

  end subroutine read_static

  !===============================================================================
  !=============================================================================

  subroutine read_tiled_file(noahmp, filename, flds, rc)

    ! input/output variables
    type(noahmp_type), intent(inout) :: noahmp
    character(len=*),  intent(in)    :: filename
    type(field_type),  intent(in)    :: flds(:)
    integer,           intent(inout) :: rc

    ! local variables
    integer                     :: i
    type(ESMF_FieldBundle)      :: FBgrid, FBmesh
    type(ESMF_ArraySpec)        :: arraySpecR8
    type(ESMF_Field)            :: fgrid, fmesh
    character(len=*), parameter :: subname = trim(modName)//': (read_tiled_file) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called for '//trim(filename), ESMF_LOGMSG_INFO)

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
       ! set field type
       call ESMF_ArraySpecSet(arraySpecR8, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! create field on grid
       fgrid = ESMF_FieldCreate(noahmp%domain%grid, arraySpecR8, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name=trim(flds(i)%short_name), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! create field on mesh
       fmesh = ESMF_FieldCreate(noahmp%domain%mesh, flds(i)%ptr1r8, meshloc=ESMF_MESHLOC_ELEMENT, &
         name=trim(flds(i)%short_name), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
    ! Create routehandle to transfer data from grid to mesh
    !----------------------

    if (.not. ESMF_RouteHandleIsCreated(noahmp%domain%rh_grid2mesh, rc=rc)) then
       call ESMF_FieldBundleRedistStore(FBgrid, FBmesh, routehandle=noahmp%domain%rh_grid2mesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------
    ! Move data from ESMF grid to mesh
    !----------------------

    call ESMF_FieldBundleRedist(FBgrid, FBmesh, noahmp%domain%rh_grid2mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_diagnose(FBmesh, trim(subname), rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Clean memory
    !----------------------

    if (dbug > 0) then
       do i = 1, size(flds)
          ! get field from FB
          call ESMF_FieldBundleGet(FBmesh, fieldName=trim(flds(i)%short_name), field=fmesh, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! write to VTK file before removing field
          ! TODO: this needs to be extended for 3d fields
          call ESMF_FieldWriteVTK(fmesh, trim(flds(i)%short_name), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
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

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine read_tiled_file

  !===============================================================================
  !=============================================================================

  subroutine write_mosaic_output(filename, noahmp, now_time, rc)

    ! input/output variables
    character(len=*)  , intent(in)    :: filename
    type(noahmp_type) , intent(inout) :: noahmp
    real(ESMF_KIND_R8), intent(in)    :: now_time
    integer           , intent(inout) :: rc

  end subroutine write_mosaic_output

  !===============================================================================

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
    !-----------------------------------------------------------------------------

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
