! Compute the timestep.  Here we use several different estimations and take


module estdt_module
  
  use bl_types
  use network,  only: nspec

  implicit none

  private

  public :: estdt

contains

  subroutine estdt(s,dx,dt,mla,the_bc_tower)

    use multifab_module
    use ml_layout_module
    use define_bc_module
    use geometry, only: dm, nlevs
    use probin_module, only: dt_mult_factor

    type(multifab),  intent(in   ) :: s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(  out) :: dt
    type(ml_layout), intent(inout) :: mla
    type(bc_tower),  intent(in   ) :: the_bc_tower

    ! local
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer :: lo(dm), hi(dm)
    integer :: n, i, ng_s

    real(kind=dp_t) :: dt_grid, dt_proc, dt_lev

    ng_s = s(1)%ng

    ! calculate the timestep
    dt     = 1.d20
    dt_lev = 1.d20

    do n = 1, nlevs

       dt_proc = 1.d20

       do i = 1, s(1)%nboxes

          if (multifab_remote(s(1),i)) cycle
          sp => dataptr(s(1),i)
          lo = lwb(get_box(s(1),i))
          hi = upb(get_box(s(1),i))
          
          dt_grid = 1.d20

          call estdt_2d(lo,hi,sp(:,:,1,:),ng_s,dx(1,1),dt_grid)

          dt_proc = min(dt_proc, dt_grid)

       enddo

       call parallel_reduce(dt_lev, dt_proc, MPI_MIN)

       dt = min(dt,dt_lev)
       dt = dt * dt_mult_factor

    enddo

  end subroutine estdt


  subroutine estdt_2d(lo, hi, s, ng_s, dx, dt)
  
    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use bl_constants_module

    integer        , intent(in   ) :: lo(:), hi(:), ng_s
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(inout) :: dt

    integer :: i, j
    real(kind=dp_t) :: diffusion_coeff

    do_diag = .false.

    do i = lo(1)-1,hi(1)+1
       do j = lo(2)-1,hi(2)+1

          den_eos(1)  = s(i,j,rho_comp)
          temp_eos(1) = s(i,j,temp_comp)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          call conducteos(eos_input_rt, den_eos, temp_eos, &
                          npts, nspec, &
                          xn_eos, &
                          p_eos, h_eos, e_eos, &
                          cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                          dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                          dpdX_eos, dhdX_eos, &
                          gam1_eos, cs_eos, s_eos, &
                          dsdt_eos, dsdr_eos, &
                          do_diag, conduct_eos)

          diffusion_coeff = abs(conduct_eos(1) / (den_eos(1) * cp_eos(1)))

          dt = min(dt, HALF*dx*dx / diffusion_coeff)

       enddo
    enddo

  end subroutine estdt_2d

  
end module estdt_module
