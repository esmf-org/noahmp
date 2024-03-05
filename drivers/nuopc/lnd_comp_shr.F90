module lnd_comp_shr

  ! Coupling related shared routines
  use ESMF,  only : ESMF_GridComp, ESMF_FieldGet
  use ESMF,  only : ESMF_Mesh, ESMF_Field, ESMF_FieldCreate
  use ESMF,  only : ESMF_MeshGet, ESMF_LogWrite, ESMF_SUCCESS
  use ESMF,  only : ESMF_LOGERR_PASSTHRU, ESMF_TYPEKIND_I4
  use ESMF,  only : ESMF_INDEX_DELOCAL, ESMF_MESHLOC_ELEMENT
  use ESMF,  only : ESMF_LogFoundError, ESMF_LogFoundNetCDFError, ESMF_LOGMSG_INFO
  use ESMF,  only : ESMF_FieldRead, ESMF_FAILURE
  use NUOPC, only : NUOPC_CompAttributeGet

  use lnd_comp_types, only: model_type
  use lnd_comp_types, only: iglobal, iregional, imesh
  use lnd_comp_kind , only: cl => shr_kind_cl

  implicit none
  private

  interface FieldRead
     module procedure FieldRead1DI4
  end interface

  interface ReadConfig
     module procedure ReadConfigCS
     module procedure ReadConfigIS
     module procedure ReadConfigIV
  end interface

  public :: ChkErr
  public :: ChkErrNc
  public :: ReadNamelist
  public :: FieldRead

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer :: master_task = 0
  integer :: dbug = 1
  character(len=1), save  :: listDel  = ":"
  character(*), parameter :: modName =  "(lnd_comp_shr)"

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  logical function ChkErr(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc

    ChkErr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       ChkErr = .true.
    endif
  end function ChkErr

  !===============================================================================

  logical function ChkErrNc(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc

    ChkErrNc = .false.
    lrc = rc
    if (ESMF_LogFoundNetCDFError(lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       ChkErrNc = .true.
    endif
  end function ChkErrNc

  !=============================================================================

  subroutine ReadNamelist(gcomp, model, rc)

    ! ----------------------------------------------
    ! Reads coupling specific namelist options
    ! ----------------------------------------------

    use netcdf

    ! input/output variables
    type(ESMF_GridComp), intent(in)   :: gcomp
    type(model_type)  , intent(inout) :: model
    integer, intent(out) :: rc

    ! local variables
    integer :: n
    integer :: ncid, dimid, ncerr
    character(len=cl) :: cname
    character(len=cl) :: cvalue
    character(len=cl) :: msg 
    character(len=cl), allocatable :: valueList(:)
    logical :: isPresent, isSet
    character(len=*),parameter :: subname=trim(modName)//':(ReadNamelist) '
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! Generic options 
    !----------------------

    call ReadConfig(gcomp, 'DebugLevel', model%nmlist%debug_level)

    !----------------------
    ! Domain options (NUOPC cap)
    !----------------------

    call ReadConfig(gcomp, 'InputDir'  , model%nmlist%input_dir, dval='INPUT/')
    call ReadConfig(gcomp, 'MosaicFile', model%nmlist%mosaic_file)

    ! if mosaic file is provided, check number of tiles
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

       ! check if it is global or regional
       if (model%domain%ntiles > 1) then
          model%domain%gtype = iglobal
       else
          model%domain%gtype = iregional
       end if
    end if

    call ReadConfig(gcomp, 'MeshFile', model%nmlist%mesh_file)

    ! extra check for configuration
    if (trim(model%nmlist%mosaic_file) == '' .and. trim(model%nmlist%mesh_file) == '') then
       call ESMF_LogWrite(trim(subname)//': Both mosaic_file and mesh_file are empty. Please define one of them to create land domain. Exiting ...', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    if (trim(model%nmlist%mosaic_file) /= '' .and. trim(model%nmlist%mesh_file) /= '') then
       call ESMF_LogWrite(trim(subname)//': Both mosaic_file and mesh_file are defined. Please define only one of them to create land domain. Exiting ... ', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    call ReadConfig(gcomp, 'DomainFile', model%nmlist%domain_file)

    ! extra check for domain file
    if (trim(model%nmlist%mesh_file) /= '' .and. trim(model%nmlist%domain_file) == '') then
        call ESMF_LogWrite(trim(subname)//': domain_file needs to be given if mesh_file is provided. Exiting ...', ESMF_LOGMSG_INFO)
        rc = ESMF_FAILURE
        return
    end if

    call ReadConfig(gcomp, 'Layout', model%domain%layout)
    call ReadConfig(gcomp, 'Dims'  , model%domain%dims)

    !----------------------
    ! Domain options (NoahMP)
    !----------------------

    call ReadConfig(gcomp, 'NumSoilLayer', model%noahmp%config%domain%NumSoilLayer)
    
    !----------------------
    ! Model configuration
    !----------------------

    call ReadConfig(gcomp, 'OptDynamicVeg'             , model%noahmp%config%nmlist%OptDynamicVeg)
    call ReadConfig(gcomp, 'OptRainSnowPartition'      , model%noahmp%config%nmlist%OptRainSnowPartition)
    call ReadConfig(gcomp, 'OptSoilWaterTranspiration' , model%noahmp%config%nmlist%OptSoilWaterTranspiration)
    call ReadConfig(gcomp, 'OptGroundResistanceEvap'   , model%noahmp%config%nmlist%OptGroundResistanceEvap)
    call ReadConfig(gcomp, 'OptSurfaceDrag'            , model%noahmp%config%nmlist%OptSurfaceDrag)
    call ReadConfig(gcomp, 'OptStomataResistance'      , model%noahmp%config%nmlist%OptStomataResistance)
    call ReadConfig(gcomp, 'OptSnowAlbedo'             , model%noahmp%config%nmlist%OptSnowAlbedo)
    call ReadConfig(gcomp, 'OptCanopyRadiationTransfer', model%noahmp%config%nmlist%OptCanopyRadiationTransfer) 
    call ReadConfig(gcomp, 'OptSnowSoilTempTime'       , model%noahmp%config%nmlist%OptSnowSoilTempTime)
    call ReadConfig(gcomp, 'OptSnowThermConduct'       , model%noahmp%config%nmlist%OptSnowThermConduct)
    call ReadConfig(gcomp, 'OptSoilTemperatureBottom'  , model%noahmp%config%nmlist%OptSoilTemperatureBottom)
    call ReadConfig(gcomp, 'OptSoilSupercoolWater'     , model%noahmp%config%nmlist%OptSoilSupercoolWater)
    call ReadConfig(gcomp, 'OptRunoffSurface'          , model%noahmp%config%nmlist%OptRunoffSurface)
    call ReadConfig(gcomp, 'OptRunoffSubsurface'       , model%noahmp%config%nmlist%OptRunoffSubsurface)
    call ReadConfig(gcomp, 'OptSoilPermeabilityFrozen' , model%noahmp%config%nmlist%OptSoilPermeabilityFrozen)
    call ReadConfig(gcomp, 'OptDynVicInfiltration'     , model%noahmp%config%nmlist%OptDynVicInfiltration)
    call ReadConfig(gcomp, 'OptTileDrainage'           , model%noahmp%config%nmlist%OptTileDrainage)
    call ReadConfig(gcomp, 'OptIrrigation'             , model%noahmp%config%nmlist%OptIrrigation)
    call ReadConfig(gcomp, 'OptIrrigationMethod'       , model%noahmp%config%nmlist%OptIrrigationMethod)
    call ReadConfig(gcomp, 'OptCropModel'              , model%noahmp%config%nmlist%OptCropModel)
    call ReadConfig(gcomp, 'OptSoilProperty'           , model%noahmp%config%nmlist%OptSoilProperty)
    call ReadConfig(gcomp, 'OptPedotransfer'           , model%noahmp%config%nmlist%OptPedotransfer)
    call ReadConfig(gcomp, 'OptGlacierTreatment'       , model%noahmp%config%nmlist%OptGlacierTreatment)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ReadNamelist

  !=============================================================================

  subroutine ReadConfigIS(gcomp, cname, val, dval)

    ! ----------------------------------------------
    ! Reads namelist option (scalar, integer) 
    ! ----------------------------------------------

    use Machine

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    character(len=*)   , intent(in)    :: cname
    integer            , intent(inout) :: val
    integer, optional  , intent(in)    :: dval

    ! local variables
    integer           :: rc
    character(len=cl) :: cvalue
    logical           :: isPresent, isSet
    character(len=*),parameter :: subname=trim(modName)//':(ReadConfigIS) '
    ! ----------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name=trim(cname), value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) val
    else
       if (present(dval)) then
          val = dval
       else
          val = undefined_int
       end if
    end if
    call ESMF_LogWrite(trim(subname)//' '//trim(cname)//' = '//trim(cvalue), ESMF_LOGMSG_INFO)

  end subroutine ReadConfigIS

  !=============================================================================

  subroutine ReadConfigIV(gcomp, cname, val, dval)

    ! ----------------------------------------------
    ! Reads namelist option (vector, integer) 
    ! ----------------------------------------------

    use Machine

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    character(len=*)   , intent(in)    :: cname
    integer            , intent(inout) :: val(:)
    integer, optional  , intent(in)    :: dval

    ! local variables
    integer           :: rc, lsize, n
    character(len=cl) :: cvalue
    character(len=cl) :: ctmp
    character(len=cl) :: msg
    logical           :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(ReadConfigIV) '
    ! ----------------------------------------------

    lsize = size(val)

    call NUOPC_CompAttributeGet(gcomp, name=trim(cname), value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       do n = 1, lsize
          call shr_string_listGetName(cvalue, n, ctmp, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          read(ctmp,*) val(n)
       end do
    else
       if (present(dval)) then
          val(:) = dval
       else
          val(:) = undefined_int
       end if
    end if
    do n = 1, lsize
       write(msg, fmt='(A,I2,A,I4)') trim(subname)//': '//trim(cname)//'(',n,') = ', val(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    end do

  end subroutine ReadConfigIV

  !=============================================================================

  subroutine ReadConfigRV(gcomp, cname, val, dval)

    ! ----------------------------------------------
    ! Reads namelist option (vector, integer) 
    ! ----------------------------------------------

    use Machine

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    character(len=*)   , intent(in)    :: cname
    real               , intent(inout) :: val(:)
    real, optional     , intent(in)    :: dval

    ! local variables
    integer           :: rc, lsize, n
    character(len=cl) :: cvalue
    character(len=cl) :: ctmp
    character(len=cl) :: msg
    logical           :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(ReadConfigIV) '
    ! ----------------------------------------------

    lsize = size(val)

    !call NUOPC_CompAttributeGet(gcomp, name=trim(cname), value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    !if (chkerr(rc,__LINE__,u_FILE_u)) return
    !if (isPresent .and. isSet) then
    !   do n = 1, lsize
    !      call shr_string_listGetName(cvalue, n, ctmp, rc)
    !      if (chkerr(rc,__LINE__,u_FILE_u)) return
    !      read(ctmp,*) val(n)
    !   end do
    !else
    !   if (present(dval)) then
    !      val(:) = dval
    !   else
    !      val(:) = undefined_int
    !   end if
    !end if
    !do n = 1, lsize
    !   write(msg, fmt='(A,I2,A,I4)') trim(subname)//': '//trim(cname)//'(',n,') = ', val(n)
    !   call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    !end do

  end subroutine ReadConfigRV

  !=============================================================================

  subroutine ReadConfigCS(gcomp, cname, val, dval)

    ! ----------------------------------------------
    ! Reads namelist option (scalar, char) 
    ! ----------------------------------------------

    use Machine

    ! input/output variables
    type(ESMF_GridComp)       , intent(in)    :: gcomp
    character(len=*)          , intent(in)    :: cname
    character(len=*)          , intent(inout) :: val
    character(len=*), optional, intent(in)    :: dval

    ! local variables
    integer           :: rc
    character(len=cl) :: cvalue
    logical           :: isPresent, isSet
    character(len=*),parameter :: subname=trim(modName)//':(ReadConfigCS) '
    ! ----------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name=trim(cname), value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       val = trim(cvalue)
    else
       if (present(dval)) then
          val = trim(dval)
       else
          val = '' 
       end if
    end if
    call ESMF_LogWrite(trim(subname)//': '//trim(cname)//' = '//trim(val), ESMF_LOGMSG_INFO)

  end subroutine ReadConfigCS

  !===============================================================================
  subroutine shr_string_listGetName(list, k, name, rc)

    ! ----------------------------------------------
    ! Get name of k-th field in list
    ! It is adapted from CDEPS, shr_string_listGetName
    ! ----------------------------------------------

    implicit none

    ! input/output variables
    character(*)     , intent(in)  :: list    ! list/string
    integer          , intent(in)  :: k       ! index of field
    character(*)     , intent(out) :: name    ! k-th name in list
    integer          , intent(out) :: rc

    ! local variables
    integer :: i,n     ! generic indecies
    integer :: kFlds   ! number of fields in list
    integer :: i0,i1   ! name = list(i0:i1)
    character(*), parameter :: subName = '(shr_string_listGetName)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    !--- check that this is a valid index ---
    kFlds = shr_string_listGetNum(list)
    if (k < 1 .or. kFlds < k) then
      call ESMF_LogWrite(trim(subname)//": ERROR invalid index ", ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
    end if

    !--- start with whole list, then remove fields before and after desired
    !field ---
    i0 = 1
    i1 = len_trim(list)

    !--- remove field names before desired field ---
    do n=2,k
       i = index(list(i0:i1),listDel)
       i0 = i0 + i
    end do

    !--- remove field names after desired field ---
    if ( k < kFlds ) then
       i = index(list(i0:i1),listDel)
       i1 = i0 + i - 2
    end if

    !--- copy result into output variable ---
    name = list(i0:i1)//"   "

  end subroutine shr_string_listGetName

  !===============================================================================

  integer function shr_string_listGetNum(str)

    ! ----------------------------------------------
    ! Get number of fields in a string list
    ! It is adapted from CDEPS, shr_string_listGetNum
    ! ----------------------------------------------

    implicit none

    ! input/output variables
    character(*), intent(in) :: str   ! string to search

    ! local variables
    integer :: count ! counts occurances of char
    character(*), parameter :: subName = '(shr_string_listGetNum'
    ! ----------------------------------------------

    shr_string_listGetNum = 0

    if (len_trim(str) > 0) then
       count = shr_string_countChar(str,listDel)
       shr_string_listGetNum = count + 1
    endif

  end function shr_string_listGetNum

  !===============================================================================

  integer function shr_string_countChar(str,char,rc)

    ! ----------------------------------------------
    ! Count number of occurances of a character
    ! It is adapted from CDEPS, shr_string_countChar
    ! ----------------------------------------------

    implicit none

    ! input/output variables
    character(*), intent(in)       :: str   ! string to search
    character(1), intent(in)       :: char  ! char to search for
    integer, intent(out), optional :: rc    ! return code

    ! local variables
    integer :: count    ! counts occurances of char
    integer :: n        ! generic index
    character(*), parameter :: subName = '(shr_string_countChar)'
    ! ----------------------------------------------

    count = 0
    do n = 1, len_trim(str)
      if (str(n:n) == char) count = count + 1
    end do
    shr_string_countChar = count

  end function shr_string_countChar

  !===============================================================================

  subroutine FieldRead1DI4(mesh, filename, varname, var1d, rc)

    ! ----------------------------------------------
    ! Reads coupling specific namelist options
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_Mesh)     , intent(in)    :: mesh
    character(len=*)    , intent(in)    :: filename
    character(len=*)    , intent(in)    :: varname
    integer, allocatable, intent(inout) :: var1d(:)
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: lsize
    integer, pointer            :: ptr1d(:)
    type(ESMF_Field)            :: field
    character(len=*), parameter :: subname = trim(modName)//':(FieldRead) '
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query number of elements 
    call ESMF_MeshGet(mesh, numOwnedElements=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create temporary field
    field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_I4, indexflag=ESMF_INDEX_DELOCAL, &
       meshloc=ESMF_MESHLOC_ELEMENT, name=trim(varname), rc=rc) 
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! read data from file
    call ESMF_FieldRead(field, trim(filename), variableName=trim(varname), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get pointer and fill variable
    call ESMF_FieldGet(field, farrayptr=ptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (.not. allocated(var1d)) allocate(var1d(lsize))
    var1d(:) = ptr1d(:)
    if (associated(ptr1d)) nullify(ptr1d)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine FieldRead1DI4

end module lnd_comp_shr
