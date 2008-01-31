module restart_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none
  
  private

  public :: fill_restart_data

contains

  subroutine fill_restart_data(restart_int,mba,chkdata,chk_p,chk_dsdt,chk_src_old, &
                               chk_rho_omegadot2,time,dt)

    use parallel
    use bl_prof_module
    use checkpoint_module

    integer          , intent(in   ) :: restart_int
    real(dp_t)       , intent(  out) :: time,dt
    type(ml_boxarray), intent(  out) :: mba

    type(multifab)   , pointer        :: chkdata(:)
    type(multifab)   , pointer        :: chk_p(:)
    type(multifab)   , pointer        :: chk_dSdt(:)
    type(multifab)   , pointer        :: chk_src_old(:)
    type(multifab)   , pointer        :: chk_rho_omegadot2(:)
    character(len=8)                  :: sd_name
    integer                           :: n,nlevs,dm

    type(bl_prof_timer), save :: bpt

    call build(bpt, "fill_restart_data")

    write(unit=sd_name,fmt='("chk",i5.5)') restart_int
    if ( parallel_IOProcessor()) &
      print *,'Reading ',sd_name,' to get state data for restart'
    call checkpoint_read(chkdata, chk_p, chk_dsdt, chk_src_old, &
         chk_rho_omegadot2, sd_name, time, dt, nlevs)

    dm = chkdata(1)%dim

    call build(mba,nlevs,dm)
    mba%pd(1) =  bbox(get_boxarray(chkdata(1)))
    do n = 2,nlevs
      mba%pd(n) = refine(mba%pd(n-1),2)
      mba%rr(n-1,:) = 2
    end do
    do n = 1,nlevs
      call boxarray_build_copy(mba%bas(n), get_boxarray(chkdata(n))) 
    end do

    call destroy(bpt)

  end subroutine fill_restart_data

end module restart_module
