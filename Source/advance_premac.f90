module pre_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private
  public :: advance_premac

contains

  subroutine advance_premac(uold,sold,umac,gpi,normal,w0,w0mac, &
                            w0_force,w0_force_cart_vec,rho0,grav_cell,dx,dt,the_bc_level,mla)

    use bl_prof_module
    use velpred_module
    use mkutrans_module
    use mk_vel_force_module
    use geometry, only: dm, nlevs
    use addw0_module
    use bl_constants_module
    use variables, only: rho_comp
    use fill_3d_module, only: put_1d_array_on_cart

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
    type(multifab) ::  force(nlevs)
    type(multifab) :: utrans(nlevs,dm)
    type(multifab) :: ufull(nlevs)
    integer        :: n,comp
    logical        :: is_final_update

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_premac")

    do n=1,nlevs
       call multifab_build(force(n),get_layout(uold(n)),dm,1)
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

    if (dm > 1) then
       call mkutrans(uold,ufull,utrans,w0,w0mac,dx,dt,the_bc_level)
    end if

    !*************************************************************
    !     Create force, initializing with pressure gradient and buoyancy terms.
    !*************************************************************
    is_final_update = .false.
    call mk_vel_force(force,is_final_update, &
                      uold,utrans,w0,w0mac,gpi,sold,rho_comp, &
                      rho0(:,:),grav_cell,dx,w0_force,w0_force_cart_vec, &
                      the_bc_level,mla)

    call add_utilde_force(force,normal,utrans,w0,dx,the_bc_level,mla)

    !*************************************************************
    !     Add w0 to trans velocities.
    !*************************************************************
    
    if (dm > 1) then
       call addw0(utrans,w0,w0mac,mult=ONE)
    end if

    !*************************************************************
    !     Create the edge states to be used for the MAC velocity 
    !*************************************************************

    call velpred(uold,ufull,umac,utrans,force,normal,w0,w0mac,dx,dt,the_bc_level,mla)

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
