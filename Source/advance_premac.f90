! advance_premac is the driver routine that orchestrates the creation
! of the edge-centered, 1/2-time advective velocities, without enforcing
! the divergence constraint.

module pre_advance_module

  use bl_types, only: dp_t
  use multifab_module, only: multifab, multifab_build, multifab_build_edge, &
                             multifab_plus_plus_c, &
                             get_layout, nghost, destroy
  use ml_layout_module, only: ml_layout
  use define_bc_module, only: bc_level

  implicit none

  private
  public :: advance_premac

contains

  subroutine advance_premac(uold,sold,umac,gpi,normal,w0,w0mac, &
                            w0_force,w0_force_cart_vec,rho0,grav_cell,dx,dt, &
                            the_bc_level,mla)

    use bl_prof_module, only: bl_prof_timer, build, destroy
    use velpred_module, only: velpred
    use mkutrans_module, only: mkutrans
    use mk_vel_force_module, only: mk_vel_force
    use addw0_module, only: addw0
    use bl_constants_module, only: ONE
    use variables, only: rho_comp
    use fill_3d_module, only: put_1d_array_on_cart
    use probin_module, only: ppm_trace_forces

    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: gpi(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: w0_force(:,0:)
    type(multifab) , intent(in   ) :: w0_force_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    ! local
    type(multifab) ::  force(mla%nlevel)
    type(multifab) :: utrans(mla%nlevel,mla%dim)
    type(multifab) ::  ufull(mla%nlevel)
    integer        :: n,comp,dm,nlevs
    logical        :: is_final_update

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_premac")

    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       if (ppm_trace_forces == 0) then
          call multifab_build(force(n),get_layout(uold(n)),dm,1)
       else
          ! tracing needs more ghost cells
          call multifab_build(force(n),get_layout(uold(n)),dm,uold(n)%ng)
       endif

       call multifab_build(ufull(n),get_layout(uold(n)),dm,nghost(uold(n)))
    end do

    ! create fullu = uold + w0
    call put_1d_array_on_cart(w0,ufull,1,.true.,.true.,dx,the_bc_level,mla)
    do n=1,nlevs
       call multifab_plus_plus_c(ufull(n),1,uold(n),1,dm,nghost(uold(n)))
    end do    

    !*************************************************************
    !     Create utrans.
    !*************************************************************

    do n=1,nlevs
       do comp=1,dm
          call multifab_build_edge(utrans(n,comp), mla%la(n),1,1,comp)
       end do
    end do

    call mkutrans(uold,ufull,utrans,w0,w0mac,dx,dt,the_bc_level,mla)


    !*************************************************************
    !     Create force, initializing with pressure gradient and buoyancy terms, and
    !     the utilde . gradw0 force
    !*************************************************************
    is_final_update = .false.
    call mk_vel_force(force,is_final_update, &
                      uold,utrans,w0,w0mac,gpi,sold,rho_comp,normal, &
                      rho0(:,:),grav_cell,dx,w0_force,w0_force_cart_vec, &
                      the_bc_level,mla,.true.)


    !*************************************************************
    !     Add w0 to trans velocities.
    !*************************************************************
    
    if (dm > 1) then
       call addw0(utrans,the_bc_level,mla,w0,w0mac,mult=ONE)
    end if

    !*************************************************************
    !     Create the edge states to be used for the MAC velocity 
    !*************************************************************

    call velpred(uold,ufull,umac,utrans,force,w0,w0mac,dx,dt,the_bc_level,mla)

    do n = 1,nlevs
       call destroy(force(n))
       call destroy(ufull(n))
       do comp=1,dm
          call destroy(utrans(n,comp))
       end do
    end do

    call destroy(bpt)

  end subroutine advance_premac

end module pre_advance_module
