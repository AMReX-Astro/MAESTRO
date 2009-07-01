module velocity_advance_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: velocity_advance

contains

  subroutine velocity_advance(mla,uold,unew,sold,rhohalf,umac,gpres, &
                              normal,w0,w0mac,w0_force,w0_force_cart_vec, &
                              rho0_old,rho0_nph,grav_cell_old,grav_cell_nph,dx,dt, &
                              the_bc_level,sponge)
    use addw0_module
    use bl_prof_module
    use update_vel_module
    use bl_constants_module
    use mk_vel_force_module
    use make_edge_scal_module
    use probin_module, only: verbose, edge_nodal_flag
    use variables, only: rho_comp
    use geometry, only: dm, nlevs

    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(inout) :: unew(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(in   ) :: rhohalf(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: gpres(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: w0_force(:,0:)
    type(multifab) , intent(in   ) :: w0_force_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_nph(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell_old(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell_nph(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(in   ) :: sponge(:)

    type(multifab)  :: force(nlevs)
    type(multifab)  :: uedge(nlevs,dm)
    logical         :: is_vel
    logical         :: is_final_update
    integer         :: velpred,n,comp
    real(kind=dp_t) :: smin,smax

    type(bl_prof_timer), save :: bpt

    call build(bpt, "velocity_advance")

    is_vel  = .true.
    velpred = 0

    do n = 1, nlevs
       call multifab_build(force(n),uold(n)%la,dm,1)
    end do

    !********************************************************
    !     Create the velocity forcing term at time n using rho 
    !********************************************************

    is_final_update = .false.
    call mk_vel_force(force,is_final_update, &
                      uold,umac,gpres,sold,normal, &
                      rho0_old,grav_cell_old,dx,the_bc_level,mla)

    call add_w0_force(force,w0_force,w0_force_cart_vec,the_bc_level,mla)

    !********************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !********************************************************
    
    call addw0(umac,w0,w0mac,mult=ONE)
    
    !********************************************************
    !     Create the edge states of velocity using the MAC velocity plus w0 on edges. 
    !********************************************************
    
    do n=1,nlevs
       do comp=1,dm
          call multifab_build(uedge(n,comp),mla%la(n),dm,0,nodal=edge_nodal_flag(comp,:))
       end do
    end do

    call make_edge_scal(uold,uedge,umac,force,normal,w0,w0mac, &
                        dx,dt,is_vel,the_bc_level,1,1,dm,.false.,mla)

    !********************************************************
    !     Subtract w0 from MAC velocities.
    !********************************************************

    call addw0(umac,w0,w0mac,mult=-ONE)

    !********************************************************
    !     Now create the force at half-time using rhohalf 
    !********************************************************

    is_final_update = .true.
    call mk_vel_force(force,is_final_update, &
                      uold,umac,gpres,rhohalf,normal, &
                      rho0_nph,grav_cell_nph,dx,the_bc_level,mla)

    call add_w0_force(force,w0_force,w0_force_cart_vec,the_bc_level,mla)

    !********************************************************
    !     Update the velocity with convective differencing
    !********************************************************
    
    call update_velocity(uold,unew,umac,uedge,force,normal,w0,w0mac, &
                         dx,dt,sponge,mla,the_bc_level)

    do n = 1, nlevs
       call destroy(force(n))
       do comp=1,dm
          call destroy(uedge(n,comp))
       end do
    end do

    if ( verbose .ge. 1 ) then
       do n = 1, nlevs
          smin = multifab_min_c(unew(n),1)
          smax = multifab_max_c(unew(n),1)
          if (parallel_IOProcessor()) write(6,1001) smin,smax
          smin = multifab_min_c(unew(n),2) 
          smax = multifab_max_c(unew(n),2)
          if (parallel_IOProcessor()) write(6,1002) smin,smax
          if (dm .eq. 3) then
             smin = multifab_min_c(unew(n),3) 
             smax = multifab_max_c(unew(n),3)
             if (parallel_IOProcessor()) write(6,1003) smin,smax
          end if
          if (parallel_IOProcessor()) write(6,1004)
          
1001      format('... new min/max : x-velocity       ',e17.10,2x,e17.10)
1002      format('... new min/max : y-velocity       ',e17.10,2x,e17.10)
1003      format('... new min/max : z-velocity       ',e17.10,2x,e17.10)
1004      format(' ')
          
       end do
    end if

    call destroy(bpt)

  end subroutine velocity_advance

end module velocity_advance_module
