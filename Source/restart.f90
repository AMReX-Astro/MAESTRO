module restart_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none
  
  private

  public :: fill_restart_data

contains

  subroutine fill_restart_data(restart_int,mba,chkdata,chk_p,chk_dsdt,chk_src_old, &
                               chk_src_new,chk_rho_omegadot2,chk_rho_Hnuc2, &
                               chk_rho_Hext,chk_thermal2,dt)

    use parallel
    use bl_prof_module
    use checkpoint_module
    use probin_module, only: check_base_name, max_levs, dm_in

    integer          , intent(in   ) :: restart_int
    real(dp_t)       , intent(  out) :: dt
    type(ml_boxarray), intent(  out) :: mba

    type(multifab)   , pointer        :: chkdata(:)
    type(multifab)   , pointer        :: chk_p(:)
    type(multifab)   , pointer        :: chk_dSdt(:)
    type(multifab)   , pointer        :: chk_src_old(:)
    type(multifab)   , pointer        :: chk_src_new(:)
    type(multifab)   , pointer        :: chk_rho_omegadot2(:)
    type(multifab)   , pointer        :: chk_rho_Hnuc2(:)
    type(multifab)   , pointer        :: chk_rho_Hext(:)
    type(multifab)   , pointer        :: chk_thermal2(:)
    character(len=5)                  :: check_index
    character(len=6)                  :: check_index6
    character(len=256)                :: check_file_name
    integer                           :: n,nlevs_local

    type(bl_prof_timer), save :: bpt

    call build(bpt, "fill_restart_data")

    if (restart_int <= 99999) then
       write(unit=check_index,fmt='(i5.5)') restart_int
       check_file_name = trim(check_base_name) // check_index
    else
       write(unit=check_index6,fmt='(i6.6)') restart_int
       check_file_name = trim(check_base_name) // check_index6
    endif

    if ( parallel_IOProcessor()) &
      print *,'Reading ', trim(check_file_name), ' to get state data for restart'
    call checkpoint_read(chkdata, chk_p, chk_dsdt, chk_src_old, chk_src_new, &
                         chk_rho_omegadot2, chk_rho_Hnuc2, chk_rho_Hext, &
                         chk_thermal2, check_file_name, &
                         dt, nlevs_local)

    call build(mba,nlevs_local,dm_in)
    mba%pd(1) =  bbox(get_boxarray(chkdata(1)))
    do n = 2,nlevs_local
      mba%pd(n) = refine(mba%pd(n-1),2)
      mba%rr(n-1,:) = 2
    end do
    do n = 1,nlevs_local
      call boxarray_build_copy(mba%bas(n), get_boxarray(chkdata(n))) 
    end do

    call destroy(bpt)

  end subroutine fill_restart_data

end module restart_module
