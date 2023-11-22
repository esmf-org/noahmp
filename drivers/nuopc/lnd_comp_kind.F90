module lnd_comp_kind 

  !----------------------------------------------------------------------------
  ! precision/kind constants add data public
  !----------------------------------------------------------------------------
  public

  integer, parameter :: shr_kind_r8 = selected_real_kind(12) ! 8 byte real
  integer, parameter :: shr_kind_r4 = selected_real_kind( 6) ! 4 byte real
  integer, parameter :: shr_kind_rn = kind(1.0)              ! native real
  integer, parameter :: shr_kind_i8 = selected_int_kind (13) ! 8 byte integer
  integer, parameter :: shr_kind_i4 = selected_int_kind ( 6) ! 4 byte integer
  integer, parameter :: shr_kind_iN = kind(1)                ! native integer
  integer, parameter :: shr_kind_cs = 80                     ! short char
  integer, parameter :: shr_kind_cm = 160                    ! mid-sized char
  integer, parameter :: shr_kind_cl = 256                    ! long char
  integer, parameter :: shr_kind_cx = 512                    ! extra-long char
  integer, parameter :: shr_kind_cxx= 4096                   ! extra-extra-long char

end module lnd_comp_kind
