module velocity_advance_module

  use bl_types
  use multifab_module
  use mkflux_module
  use mk_vel_force_module
  use addw0_module
  use make_grav_module
  use update_module
  use define_bc_module

  implicit none

contains

   subroutine velocity_advance(uold,unew,sold,rhohalf,&
                               umac,uedge,utrans,gp,p,&
                               normal,w0,w0_cart_vec,s0,s0_nph,&
                               dx,time,dt, &
                               the_bc_level, &
                               verbose)
 
      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: unew
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: umac(:)
      type(multifab) , intent(inout) :: rhohalf
      type(multifab) , intent(inout) :: utrans(:)
      type(multifab) , intent(inout) :: uedge(:)
      type(multifab) , intent(inout) :: gp
      type(multifab) , intent(inout) :: p
      type(multifab) , intent(in   ) :: normal
      real(kind=dp_t), intent(in   ) :: w0(:)
      type(multifab) , intent(in   ) :: w0_cart_vec
      real(kind=dp_t), intent(in   ) :: s0(:,:),s0_nph(:,:)

      real(kind=dp_t), intent(in   ) :: dx(:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level
 
      integer        , intent(in   ) :: verbose
 
      type(multifab) :: force
 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: unp(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: vmp(:,:,:,:)
      real(kind=dp_t), pointer:: wmp(:,:,:,:)
      real(kind=dp_t), pointer:: w0p(:,:,:,:)
      real(kind=dp_t), pointer:: utp(:,:,:,:)
      real(kind=dp_t), pointer:: vtp(:,:,:,:)
      real(kind=dp_t), pointer:: wtp(:,:,:,:)
      real(kind=dp_t), pointer::  np(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
      real(kind=dp_t), pointer:: uepx(:,:,:,:)
      real(kind=dp_t), pointer:: uepy(:,:,:,:)
      real(kind=dp_t), pointer:: uepz(:,:,:,:)
!
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:) 
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)

      real(kind=dp_t), allocatable:: grav_cell(:)
!
      integer :: nr,velpred,edge_based
      integer :: lo(uold%dim),hi(uold%dim)
      integer :: i,n,comp,dm,ng_cell,ng_rho
      logical :: is_vel,is_conservative(uold%dim)
      logical :: dummy_false
      real(kind=dp_t) :: mult

      ng_cell = uold%ng
      ng_rho  = rhohalf%ng
      dm      = uold%dim

      is_conservative = .false.
      dummy_false     = .false.

      call multifab_build(force,uold%la,dm,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the velocity forcing term at time n using rho 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nr = size(s0,dim=1)
      allocate(grav_cell(nr))
      call make_grav_cell(grav_cell,s0(:,rho_comp))
      call mk_vel_force(force,gp,sold,normal,s0,grav_cell,dx,rho_comp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Add w0 to MAC velocities (trans velocities already have w0).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mult = ONE
      call addw0(umac,w0,w0_cart_vec,dx,mult)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of velocity using the MAC velocity plus w0 on edges. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      velpred = 0
      is_vel = .true.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop  => dataptr(uold, i)
         sop  => dataptr(sold, i)
         uepx => dataptr(uedge(1), i)
         uepy => dataptr(uedge(2), i)
         ump  => dataptr(umac(1), i)
         vmp  => dataptr(umac(2), i)
         utp  => dataptr(utrans(1), i)
         vtp  => dataptr(utrans(2), i)
          fp  => dataptr(force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              do n = 1,dm
                call mkflux_2d(uop(:,:,1,:), uop(:,:,1,:), &
                               uepx(:,:,1,:), uepy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), w0, &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%ell_bc_level_array(i,:,:,n:), &
                               velpred, ng_cell, s0(:,rho_comp), dummy_false, n)
              end do
            case (3)
              do n = 1,dm
                uepz => dataptr(uedge(3), i)
                wmp  => dataptr(umac(3), i)
                wtp  => dataptr(utrans(3), i)
                w0p  => dataptr(w0_cart_vec, i)
                call mkflux_3d(uop(:,:,:,:), uop(:,:,:,:), &
                               uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), w0, w0p(:,:,:,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%ell_bc_level_array(i,:,:,n:), &
                               velpred, ng_cell, s0(:,rho_comp), dummy_false, n)
              end do
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Subtract w0 from MAC velocities.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mult = -ONE
      call addw0(umac,w0,w0_cart_vec,dx,mult)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Now create the force at half-time using rhohalf 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call make_grav_cell(grav_cell,s0_nph(:,rho_comp))
      call mk_vel_force(force,gp,rhohalf,normal,s0_nph,grav_cell,dx,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Update the velocity with convective differencing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop => dataptr(uold, i)
         unp => dataptr(unew, i)
         ump => dataptr(umac(1), i)
         vmp => dataptr(umac(2), i)
         uepx => dataptr(uedge(1), i)
         uepy => dataptr(uedge(2), i)
          fp => dataptr(force, i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call update_velocity_2d( &
                             uop(:,:,1,:), unp(:,:,1,:), &
                             ump(:,:,1,1), vmp(:,:,1,1), &
                             uepx(:,:,1,:), uepy(:,:,1,:), &
                             fp(:,:,1,:), s0, w0, &
                             lo, hi, ng_cell, dx, time, dt, verbose)
            case (3)
              wmp => dataptr(umac(3), i)
              uepz => dataptr(uedge(3), i)
              w0p  => dataptr(w0_cart_vec, i)
              call update_velocity_3d( &
                             uop(:,:,:,:), unp(:,:,:,:), &
                              ump(:,:,:,1),  vmp(:,:,:,1),  wmp(:,:,:,1), &
                             uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                             fp(:,:,:,:), w0, w0p(:,:,:,:), &
                             lo, hi, ng_cell, dx, time, dt, verbose)
         end select
      end do

      call multifab_destroy(force)

      deallocate(grav_cell)

   end subroutine velocity_advance

end module velocity_advance_module
