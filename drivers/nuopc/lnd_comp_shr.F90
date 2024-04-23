module lnd_comp_shr

  ! Coupling related shared routines
  use ESMF, only: ESMF_GridComp, ESMF_FieldGet
  use ESMF, only: ESMF_Mesh, ESMF_Field, ESMF_FieldCreate
  use ESMF, only: ESMF_MeshGet, ESMF_LogWrite, ESMF_SUCCESS
  use ESMF, only: ESMF_LOGERR_PASSTHRU, ESMF_TYPEKIND_I4
  use ESMF, only: ESMF_INDEX_DELOCAL, ESMF_MESHLOC_ELEMENT
  use ESMF, only: ESMF_LogFoundError, ESMF_LogFoundNetCDFError
  use ESMF, only: ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
  use ESMF, only: ESMF_FieldRead, ESMF_FAILURE
  use NUOPC,only: NUOPC_CompAttributeGet

  use lnd_comp_types, only: model_type
  use lnd_comp_types, only: iMosaic, iScrip
  use lnd_comp_kind , only: cl => shr_kind_cl
  use lnd_comp_kind , only: r8 => shr_kind_r8

  implicit none
  private

  interface ReadConfig
     module procedure ReadConfigCharScalar
     module procedure ReadConfigIntScalar
     module procedure ReadConfigIntVector
     module procedure ReadConfigLogScalar
     module procedure ReadConfigRealScalar
     module procedure ReadConfigRealVector
  end interface

  public :: ChkErr
  public :: ChkErrNc
  public :: ReadNamelist

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer :: master_task = 0
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
    character(len=*), parameter :: subname=trim(modName)//':(ReadNamelist) '
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! Generic run options 
    !----------------------

    call ReadConfig(gcomp, 'case_name' , model%nmlist%CaseName, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadConfig(gcomp, 'start_type', model%nmlist%RestartType, dval='startup', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(model%nmlist%RestartType) == 'startup') then
       model%nmlist%IsRestart = .false.
    else if (trim(model%nmlist%RestartType) == 'continue') then
       model%nmlist%IsRestart = .true.
    else
       call ESMF_LogWrite(trim(subname)//": ERROR in start_type. It can be only 'startup' and 'continue'.", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    call ReadConfig(gcomp, 'DebugLevel', model%nmlist%debug_level)

    !----------------------
    ! Soil options
    !----------------------

    associate( &
       NumSoilLayer => model%noahmp%config%domain%NumSoilLayer, &
       NumSnowLayerMax => model%noahmp%config%domain%NumSnowLayerMax)

       call ReadConfig(gcomp, 'NumSoilLayer', model%noahmp%config%domain%NumSoilLayer, dval=4)

       if (.not. allocated(model%noahmp%config%domain%DepthSoilLayer)) then
          allocate(model%noahmp%config%domain%DepthSoilLayer(1:NumSoilLayer))
       end if
       call ReadConfig(gcomp, 'DepthSoilLayer', model%noahmp%config%domain%DepthSoilLayer)

       if (.not. allocated(model%noahmp%config%domain%ThicknessSoilLayer)) then
          allocate(model%noahmp%config%domain%ThicknessSoilLayer(1:NumSoilLayer))
       end if
       call ReadConfig(gcomp, 'ThicknessSoilLayer', model%noahmp%config%domain%ThicknessSoilLayer)

       if (.not. allocated(model%noahmp%config%domain%SoilType)) then 
          allocate(model%noahmp%config%domain%SoilType(1:NumSoilLayer))
       end if

       call ReadConfig(gcomp, 'NumSnowLayerMax', model%noahmp%config%domain%NumSnowLayerMax, dval=3)

       if (.not. allocated(model%noahmp%config%domain%ThicknessSnowSoilLayer)) then
          allocate(model%noahmp%config%domain%ThicknessSnowSoilLayer(-NumSnowLayerMax+1:NumSoilLayer))
       end if
       call ReadConfig(gcomp, 'ThicknessSnowSoilLayer', model%noahmp%config%domain%ThicknessSnowSoilLayer)

       if (.not. allocated(model%noahmp%config%domain%DepthSnowSoilLayer)) then
          allocate(model%noahmp%config%domain%DepthSnowSoilLayer(-NumSnowLayerMax+1:NumSoilLayer))
       end if
       call ReadConfig(gcomp, 'DepthSnowSoilLayer', model%noahmp%config%domain%DepthSnowSoilLayer)

    end associate

    !----------------------
    ! Output options 
    !----------------------

    call ReadConfig(gcomp, 'OutputMode', model%nmlist%OutputMode, dval='all', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(model%nmlist%OutputMode) == 'all' .or. &
        trim(model%nmlist%OutputMode) == 'mid' .or. &
        trim(model%nmlist%OutputMode) == 'low') then
    else
       call ESMF_LogWrite(trim(subname)//": ERROR in OutputMode. It can be only 'all', 'mid' and 'low'.", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return 
    end if

    call ReadConfig(gcomp, 'OutputFreq', model%nmlist%OutputFreq, dval=3600)

    if (model%nmlist%OutputFreq == 0) then
       call ESMF_LogWrite(trim(subname)//": ERROR in OutputFreq. It can not be set to zero!", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return 
    end if

    !----------------------
    ! Domain options
    !----------------------

    call ReadConfig(gcomp, 'InputDir'  , model%nmlist%input_dir, dval='INPUT/', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadConfig(gcomp, 'MosaicFile', model%nmlist%mosaic_file, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
          model%domain%dtype = iMosaic
          model%nmlist%isGlobal = .true.
       else
          call ESMF_LogWrite(trim(subname)//': The mosaic_file is defined but number of tiles are less than 6! &
             Please define scrip_file for regional applications. Exiting ...', ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
    end if

    ! scrip file
    call ReadConfig(gcomp, 'ScripFile', model%nmlist%scrip_file, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(model%nmlist%scrip_file) /= '') then
       model%domain%dtype = iScrip
    end if

    ! regional vs. global domain
    call ReadConfig(gcomp, 'IsGlobal', model%nmlist%isGlobal, dval=.true.)

    ! extra check for configuration
    if (trim(model%nmlist%mosaic_file) == '' .and. trim(model%nmlist%scrip_file) == '') then
       call ESMF_LogWrite(trim(subname)//': Both mosaic_file and scrip_file options are empty. &
          Please define one of them to create land domain. Exiting ...', ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    if (trim(model%nmlist%mosaic_file) /= '' .and. trim(model%nmlist%scrip_file) /= '') then
       call ESMF_LogWrite(trim(subname)//': Both mosaic_file and scrip_file are defined. &
          Please define only one of them to create land domain. Exiting ... ', ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    call ReadConfig(gcomp, 'Layout', model%domain%layout, dval=(/1,1/))
    call ReadConfig(gcomp, 'Dims'  , model%domain%dims, dval=(/1,1/))

    !----------------------
    ! Input
    !----------------------

    call ReadConfig(gcomp, 'InputType', model%nmlist%InputType, dval='sfc', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(model%nmlist%InputType) /= 'sfc' .and. &
        trim(model%nmlist%InputType) /= 'custom') then
       call ESMF_LogWrite(trim(subname)//": ERROR in InputType. It can be only 'sfc' and 'custom'.", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    call ReadConfig(gcomp, 'InputIC', model%nmlist%InputIC, required=.true., rc=rc) 
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadConfig(gcomp, 'InputSoilType', model%nmlist%InputSoilType, required=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadConfig(gcomp, 'InputVegType', model%nmlist%InputVegType, required=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadConfig(gcomp, 'InputSlopeType', model%nmlist%InputSlopeType, required=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadConfig(gcomp, 'InputBottomTemp', model%nmlist%InputBottomTemp, required=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadConfig(gcomp, 'InputVegFracAnnMax', model%nmlist%InputVegFracAnnMax, required=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Model configuration
    !----------------------

    ! config
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

    if (model%noahmp%config%nmlist%OptSoilProperty == 2) then
       call ESMF_LogWrite(trim(subname)//": ERROR in OptSoilProperty. Only option 1 and 3 are supported. Exiting ...", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    call ReadConfig(gcomp, 'OptPedotransfer'           , model%noahmp%config%nmlist%OptPedotransfer)
    call ReadConfig(gcomp, 'OptGlacierTreatment'       , model%noahmp%config%nmlist%OptGlacierTreatment)

    if (model%noahmp%config%nmlist%OptSoilProperty /= 1) then
       call ESMF_LogWrite(trim(subname)//': ERROR in OptSoilProperty. It can be set only to 1 at this point. Exiting ... ', ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    ! domain
    call ReadConfig(gcomp, 'SoilTimeStep'   , model%noahmp%config%domain%SoilTimeStep, dval=-999.0_r8)
    call ReadConfig(gcomp, 'NumSoilLayer'   , model%noahmp%config%domain%NumSoilLayer)
    call ReadConfig(gcomp, 'FlagSoilProcess', model%noahmp%config%domain%FlagSoilProcess, dval=.false.)
    call ReadConfig(gcomp, 'LandUseDataName', model%noahmp%config%domain%LandUseDataName, dval='MODIFIED_IGBP_MODIS_NOAH', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !call ReadConfig(gcomp, 'LandUseDataName', model%noahmp%config%domain%LandUseDataName, dval='MODIFIED_IGBP_MODIS_NOAH')
       !dval='MODIFIED_IGBP_MODIS_NOAH', &
       !opts=(/'USGS                    ', &
       !       'MODIFIED_IGBP_MODIS_NOAH'/))

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ReadNamelist

  !=============================================================================

  subroutine ReadConfigIntScalar(gcomp, cname, val, dval)

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
    character(len=*),parameter :: subname=trim(modName)//':(ReadConfigIntScalar) '
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
    write(cvalue,*) val
    call ESMF_LogWrite(trim(subname)//' '//trim(cname)//' = '//trim(cvalue), ESMF_LOGMSG_INFO)

  end subroutine ReadConfigIntScalar

  !=============================================================================

  subroutine ReadConfigIntVector(gcomp, cname, val, dval)

    ! ----------------------------------------------
    ! Reads namelist option (vector, integer) 
    ! ----------------------------------------------

    use Machine

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    character(len=*)   , intent(in)    :: cname
    integer            , intent(inout) :: val(:)
    integer, optional  , intent(in)    :: dval(:)

    ! local variables
    integer           :: rc, lsize, n
    character(len=cl) :: cvalue
    character(len=cl) :: ctmp
    character(len=cl) :: msg
    logical           :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(ReadConfigIntVector) '
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
          val(:) = dval(:)
       else
          val(:) = undefined_int
       end if
    end if
    do n = 1, lsize
       write(msg, fmt='(A,I2,A,I4)') trim(subname)//': '//trim(cname)//'(',n,') = ', val(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    end do

  end subroutine ReadConfigIntVector

  !=============================================================================

  subroutine ReadConfigRealScalar(gcomp, cname, val, dval)

    ! ----------------------------------------------
    ! Reads namelist option (scalar, real) 
    ! ----------------------------------------------

    use Machine

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    character(len=*)   , intent(in)    :: cname
    real(r8)           , intent(inout) :: val
    real(r8), optional , intent(in)    :: dval

    ! local variables
    integer           :: rc
    character(len=cl) :: cvalue
    logical           :: isPresent, isSet
    character(len=*),parameter :: subname=trim(modName)//':(ReadConfigRealScalar) '
    ! ----------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name=trim(cname), value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) val
    else
       if (present(dval)) then
          val = dval
       else
          val = undefined_real
       end if
    end if
    write(cvalue,*) val
    call ESMF_LogWrite(trim(subname)//' '//trim(cname)//' = '//trim(cvalue), ESMF_LOGMSG_INFO)

  end subroutine ReadConfigRealScalar

  !=============================================================================

  subroutine ReadConfigRealVector(gcomp, cname, val, dval)

    ! ----------------------------------------------
    ! Reads namelist option (vector, real) 
    ! ----------------------------------------------

    use Machine

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    character(len=*)   , intent(in)    :: cname
    real(r8)           , intent(inout) :: val(:)
    real(r8), optional , intent(in)    :: dval(:)

    ! local variables
    integer           :: rc, lsize, n
    character(len=cl) :: cvalue
    character(len=cl) :: ctmp
    character(len=cl) :: msg
    logical           :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(ReadConfigRealVector) '
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
          val(:) = dval(:)
       else
          val(:) = undefined_real
       end if
    end if
    do n = 1, lsize
       write(msg, fmt='(A,I2,A,F10.3)') trim(subname)//': '//trim(cname)//'(',n,') = ', val(n)
       call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    end do

  end subroutine ReadConfigRealVector

  !=============================================================================

  !subroutine ReadConfigCharScalar(gcomp, cname, val, dval, required, opts, rc)
  subroutine ReadConfigCharScalar(gcomp, cname, val, dval, required, rc)

    ! ----------------------------------------------
    ! Reads namelist option (scalar, char) 
    ! ----------------------------------------------

    use Machine

    ! input/output variables
    type(ESMF_GridComp)       , intent(in)    :: gcomp
    character(len=*)          , intent(in)    :: cname
    character(len=*)          , intent(inout) :: val
    character(len=*), optional, intent(in)    :: dval
    logical         , optional, intent(in)    :: required
    !character(len=*), allocatable, optional, intent(in)    :: opts(:) 
    integer         ,           intent(out)   :: rc

    ! local variables
    integer           :: i
    character(len=cl) :: cvalue
    logical           :: isPresent, isSet
    logical           :: noDefault, validOpt
    character(len=*),parameter :: subname=trim(modName)//':(ReadConfigCharScalar) '
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! check the option required or not
    noDefault = .false.
    if (present(required)) noDefault = .true.

    ! get attribute
    call NUOPC_CompAttributeGet(gcomp, name=trim(cname), value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! check it is given or not
    if (isPresent .and. isSet) then
       val = trim(cvalue)
    else
       ! option is not provided, set to default if it is given
       if (present(dval)) then
          val = trim(dval)
       else
          ! throw error if it is required
          if (noDefault) then
             call ESMF_LogWrite(trim(subname)//": "//trim(cname)//" needs to be provided!", ESMF_LOGMSG_ERROR)
             rc = ESMF_FAILURE
             return
          else
             ! assign to empty string
             val = ''
          end if 
       end if
    end if

    !if (present(opts)) then
    !   ! check it is valid or not
    !   validOpt = .false.
    !   do i = 1, size(opts, dim=1)
    !      if (trim(val) == trim(opts(i))) then
    !         validOpt = .true.
    !         exit
    !      end if
    !   end do 

    !   ! throw error if it is not valid option
    !   if (.not. validOpt) then
    !      call ESMF_LogWrite(trim(subname)//": Invalid option given for "//trim(cname), ESMF_LOGMSG_ERROR)
    !      cvalue = ""
    !      do i = 1, size(opts, dim=1)
    !         write(cvalue, *) trim(cvalue)//" "//trim(opts(i))
    !      end do
    !      call ESMF_LogWrite(trim(subname)//": Valid options are "//trim(cvalue), ESMF_LOGMSG_ERROR)
    !      rc = ESMF_FAILURE
    !      return
    !   end if
    !end if

    call ESMF_LogWrite(trim(subname)//': '//trim(cname)//' = '//trim(val), ESMF_LOGMSG_INFO)

  end subroutine ReadConfigCharScalar

  !=============================================================================

  subroutine ReadConfigLogScalar(gcomp, cname, val, dval)

    ! ----------------------------------------------
    ! Reads namelist option (scalar, logical) 
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    character(len=*)   , intent(in)    :: cname
    logical            , intent(inout) :: val
    logical, optional  , intent(in)    :: dval

    ! local variables
    integer           :: rc
    character(len=cl) :: cvalue
    logical           :: isPresent, isSet
    character(len=cl) :: msg
    character(len=*),parameter :: subname=trim(modName)//':(ReadConfigLogScalar) '
    ! ----------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name=trim(cname), value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    val = .false. 
    if (isPresent .and. isSet) then
       if (trim(cvalue) .eq. '.true.' .or. trim(cvalue) .eq. 'true' .or. &
           trim(cvalue) .eq. '.TRUE.' .or. trim(cvalue) .eq. 'TRUE' .or. &
           trim(cvalue) .eq. '.True.' .or. trim(cvalue) .eq. 'True') then
          val = .true.
       end if
    else
       if (present(dval)) then
          val = dval
       end if
    end if

    write(msg, fmt='(A,I2,A,L)') trim(subname)//': '//trim(cname)//' = ', val
    call ESMF_LogWrite(trim(msg)//' = ', ESMF_LOGMSG_INFO)

  end subroutine ReadConfigLogScalar

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

end module lnd_comp_shr
