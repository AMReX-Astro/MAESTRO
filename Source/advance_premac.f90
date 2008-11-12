module pre_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private
  public :: advance_premac

contains

  subroutine advance_premac(uold,sold,umac,gpres,normal,w0,w0mac, &
                            w0_force,w0_force_cart_vec,rho0,grav_cell,dx,dt,the_bc_level,mla)

    use bl_prof_module
    use velpred_module
    use mkutrans_module
    use mk_vel_force_module
    use probin_module, only: edge_nodal_flag
    use geometry, only: dm, nlevs
    use addw0_module
    use bl_constants_module

    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: gpres(:)
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
    integer        :: n,comp

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_premac")

    do n=1,nlevs
       call multifab_build(force(n),uold(n)%la,dm,1)
    end do

    !*************************************************************
    !     Create force, initializing with pressure gradient and buoyancy terms.
    !*************************************************************

    call mk_vel_force(force,uold,gpres,sold,normal,rho0(:,:),grav_cell,dx,the_bc_level,mla)

    call add_w0_force(force,w0_force,w0_force_cart_vec,the_bc_level,mla)

    !*************************************************************
    !     Create utrans.
    !*************************************************************

    do n=1,nlevs
       do comp=1,dm
          call multifab_build(utrans(n,comp), mla%la(n),1,1,nodal=edge_nodal_flag(comp,:))
       end do
    end do

    call mkutrans(uold,utrans,w0,w0mac,dx,dt,the_bc_level)

    !*************************************************************
    !     Add w0 to trans velocities.
    !*************************************************************
    
    call addw0(utrans,w0,w0mac,mult=ONE)

    !*************************************************************
    !     Create the edge states to be used for the MAC velocity 
    !*************************************************************

    call velpred(uold,umac,utrans,force,normal,w0,w0mac,dx,dt,the_bc_level,mla)

    do n = 1,nlevs
       call destroy(force(n))
       do comp=1,dm
          call destroy(utrans(n,comp))
       end do
    end do

    call destroy(bpt)

  end subroutine advance_premac

end module pre_advance_module
