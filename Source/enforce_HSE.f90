module enforce_HSE_module

  use bl_types

  implicit none

  private

  public :: enforce_HSE, compute_th_from_rp, compute_rh_from_tp

contains

  subroutine enforce_HSE(nlevs,rho0,p0,grav_cell)

    use geometry, only: dr, r_start_coord, r_end_coord, numdisjointchunks, spherical, &
         base_cutoff_density_coord
    use restrict_base_module, only: fill_ghost_base

    integer,         intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(inout) ::   p0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell(:,0:)

    integer         :: n,i,r
    real(kind=dp_t) :: grav

    if (spherical .eq. 0) then

       ! gravity is constant
       grav = grav_cell(1,0)

       ! do level 1 first
       ! we start at r=1 since the pressure at r=0 is assumed correct
       do r=1,min(r_end_coord(1,1),base_cutoff_density_coord(1))
          p0(1,r) = p0(1,r-1) + (dr(1)/2.d0)*(rho0(1,r)+rho0(1,r-1))*grav
       end do
       do r=base_cutoff_density_coord(1)+1,r_end_coord(1,1)
          p0(1,r) = p0(1,r-1)
       end do

       do n=2,nlevs

          do i=1,numdisjointchunks(n)

             ! use a special stencil for the first point
             if (r_start_coord(n,i) .le. base_cutoff_density_coord(n)) then
                p0(n,r_start_coord(n,i)) = p0(n-1,r_start_coord(n,i)/2) &
                     + (2.d0/3.d0)*(rho0(n-1,r_start_coord(n,i)/2))*grav &
                     + (1.d0/3.d0)*(rho0(n,r_start_coord(n,i)))*grav
             else
                p0(n,r_start_coord(n,i)) = p0(n-1,r_start_coord(n,i)/2)
             end if

             ! iterate normally over the rest          
             do r=r_start_coord(n,i)+1,min(r_end_coord(n,i),base_cutoff_density_coord(n))
                p0(n,r) = p0(n,r-1) + (dr(n)/2.d0)*(rho0(n,r)+rho0(n,r-1))*grav
             end do
             do r=base_cutoff_density_coord(n)+1,r_end_coord(n,i)
                p0(n,r) = p0(n,r-1)
             end do

          end do

       end do

    else

    end if

    call fill_ghost_base(nlevs,p0,.true.)

  end subroutine enforce_HSE

  subroutine compute_th_from_rp(nlevs,s,p0,bc,mla)

    use multifab_module
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module
    use multifab_fill_ghost_module
    use variables, only: rhoh_comp, temp_comp
    use multifab_physbc_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    integer                  :: i,n,ng_s,dm
    integer                  :: lo(s(1)%dim),hi(s(1)%dim)

    ng_s = s(1)%ng
    dm = s(1)%dim

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call compute_th_from_rp_2d(sop(:,:,1,:), ng_s, lo, hi, p0(n,:))
          case (3)
             call compute_th_from_rp_3d(sop(:,:,:,:), ng_s, lo, hi, p0(n,:))             
          end select
       end do
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s(nlevs),rhoh_comp,1)
       call multifab_fill_boundary_c(s(nlevs),temp_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rhoh_comp,dm+rhoh_comp,1,bc(nlevs))
       call multifab_physbc(s(nlevs),temp_comp,dm+temp_comp,1,bc(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(s(n-1),rhoh_comp,s(n),rhoh_comp,mla%mba%rr(n-1,:))
          call ml_cc_restriction_c(s(n-1),temp_comp,s(n),temp_comp,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1),ng_s,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,dm+rhoh_comp,1)
          call multifab_fill_ghost_cells(s(n),s(n-1),ng_s,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,dm+temp_comp,1)

       enddo

    end if

  end subroutine compute_th_from_rp

  subroutine compute_th_from_rp_2d(s,ng_s,lo,hi,p0)

    use eos_module
    use network
    use variables

    integer           , intent(in   ) :: lo(:),hi(:),ng_s
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,:)  
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

    ! local
    integer    :: i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          den_eos(1) = s(i,j,rho_comp)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/s(i,j,rho_comp)
          temp_eos(1) = s(i,j,temp_comp)
          p_eos(1) = p0(j)

          call eos(eos_input_rp, den_eos, temp_eos, &
                   npts, nspec, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   do_diag)

          s(i,j,rhoh_comp) = den_eos(1)*h_eos(1)
          s(i,j,temp_comp) = temp_eos(1)

       end do
    end do

  end subroutine compute_th_from_rp_2d
  
  subroutine compute_th_from_rp_3d(s,ng_s,lo,hi,p0)

    use eos_module
    use network
    use variables

    integer           , intent(in   ) :: lo(:),hi(:),ng_s
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)  
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

    ! local
    integer    :: i,j,k

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             
             den_eos(1) = s(i,j,k,rho_comp)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/s(i,j,K,rho_comp)
             temp_eos(1) = s(i,j,k,temp_comp)
             p_eos(1) = p0(k)
             
             call eos(eos_input_rp, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)
             
             s(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             s(i,j,k,temp_comp) = temp_eos(1)
             
          end do
       end do
    end do

  end subroutine compute_th_from_rp_3d

  subroutine compute_rh_from_tp(nlevs,s,p0,bc,mla)

    use multifab_module
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module
    use multifab_fill_ghost_module
    use variables, only: rhoh_comp, temp_comp
    use multifab_physbc_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    integer                  :: i,n,ng_s,dm
    integer                  :: lo(s(1)%dim),hi(s(1)%dim)

    ng_s = s(1)%ng
    dm = s(1)%dim

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call compute_rh_from_tp_2d(sop(:,:,1,:), ng_s, lo, hi, p0(n,:))
          case (3)
             call compute_rh_from_tp_3d(sop(:,:,:,:), ng_s, lo, hi, p0(n,:))             
          end select
       end do
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s(nlevs),rhoh_comp,1)
       call multifab_fill_boundary_c(s(nlevs), rho_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rhoh_comp,dm+rhoh_comp,1,bc(nlevs))
       call multifab_physbc(s(nlevs),temp_comp,dm+rho_comp,1,bc(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(s(n-1),rhoh_comp,s(n),rhoh_comp,mla%mba%rr(n-1,:))
          call ml_cc_restriction_c(s(n-1),temp_comp,s(n), rho_comp,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1),ng_s,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,dm+rhoh_comp,1)
          call multifab_fill_ghost_cells(s(n),s(n-1),ng_s,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,dm+rho_comp,1)

       enddo

    end if
    
  end subroutine compute_rh_from_tp

  subroutine compute_rh_from_tp_2d(s,ng_s,lo,hi,p0)

    use eos_module
    use network
    use variables

    integer           , intent(in   ) :: lo(:),hi(:),ng_s
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,:)  
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

    ! local
    integer    :: i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          den_eos(1) = s(i,j,rho_comp)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/s(i,j,rho_comp)
          temp_eos(1) = s(i,j,temp_comp)
          p_eos(1) = p0(j)

          call eos(eos_input_tp, den_eos, temp_eos, &
                   npts, nspec, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   do_diag)

          s(i,j,rhoh_comp) = den_eos(1)*h_eos(1)
          s(i,j,rho_comp) =  den_eos(1)

       end do
    end do

  end subroutine compute_rh_from_tp_2d
  
  subroutine compute_rh_from_tp_3d(s,ng_s,lo,hi,p0)

    use eos_module
    use network
    use variables

    integer           , intent(in   ) :: lo(:),hi(:),ng_s
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)  
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

    ! local
    integer    :: i,j,k

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             
             den_eos(1) = s(i,j,k,rho_comp)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/s(i,j,K,rho_comp)
             temp_eos(1) = s(i,j,k,temp_comp)
             p_eos(1) = p0(k)
             
             call eos(eos_input_tp, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)
             
             s(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             s(i,j,k,rho_comp)  = den_eos(1)
             
          end do
       end do
    end do

  end subroutine compute_rh_from_tp_3d

end module enforce_HSE_module
