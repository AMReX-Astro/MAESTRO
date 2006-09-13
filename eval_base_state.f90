module base_state_module

  ! adjust the base state quantities in response to the heating.
  ! This is step 3 of ABRZ2.

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_module
  use heating_module
  use mkflux_module
  use make_div_coeff_module
  use eos_module

  implicit none

contains

   subroutine eval_base_state (vel,p0_old,p0_new,rho0_old,rho0_nph,rho0_new,temp0,rhoh0_old,rhoh0_new, &
                               gam1,div_coeff_n,div_coeff_nph,div_coeff_half, &
                               rhohalf,grav,dx,dt,time,div_coef_type,anelastic_cutoff)

      real(kind=dp_t), intent(  out) :: vel(:)
      real(kind=dp_t), intent(in   ) :: p0_old(:), rho0_old(:), rhoh0_old(:)
      real(kind=dp_t), intent(  out) :: p0_new(:), rho0_new(:), rho0_nph(:), rhoh0_new(:)
      real(kind=dp_t), intent(inout) :: temp0(:),gam1(:),div_coeff_n(:),div_coeff_nph(:),div_coeff_half(:)
      type(multifab) , intent(in   ) :: rhohalf
      real(kind=dp_t), intent(in   ) :: dx(:),dt,grav,time,anelastic_cutoff
      integer        , intent(in   ) :: div_coef_type

      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer:: sep(:,:,:,:)
      real(kind=dp_t), pointer:: ufp(:,:,:,:)
      real(kind=dp_t), pointer:: rhp(:,:,:,:)
      real(kind=dp_t), pointer:: uap(:,:,:,:)
      integer :: lo(rhohalf%dim),hi(rhohalf%dim),dm
      integer :: i

      dm = rhohalf%dim

      do i = 1, rhohalf%nboxes
         if ( multifab_remote(rhohalf, i) ) cycle
         rhp => dataptr(rhohalf , i)
         lo =  lwb(get_box(rhohalf, i))
         hi =  upb(get_box(rhohalf, i))
         select case (dm)
            case (2)
              call eval_base_state_2d(vel,p0_old,p0_new,rho0_old,rho0_nph,rho0_new,temp0, &
                                      rhoh0_old,rhoh0_new, &
                                      gam1,div_coeff_n,div_coeff_nph, &
                                      div_coeff_half,rhp(:,:,1,1),grav,lo,hi, &
                                      dx,dt,time,div_coef_type,anelastic_cutoff)
            case (3)
              call eval_base_state_3d(vel,p0_old,p0_new,rho0_old,rho0_nph,rho0_new,temp0, &
                                      rhoh0_old,rhoh0_new, &
                                      gam1,div_coeff_n,div_coeff_nph, &
                                      div_coeff_half,rhp(:,:,:,1),grav,lo,hi, &
                                      dx,dt,time,div_coef_type,anelastic_cutoff)
         end select
      end do

   end subroutine eval_base_state

   subroutine eval_base_state_2d (vel,p0_old,p0_new,rho0_old,rho0_nph,rho0_new,temp0, &
                                  rhoh0_old,rhoh0_new, &
                                  gam1,div_coeff_n,div_coeff_nph, &
                                  div_coeff_half,rhohalf,grav,lo,hi, & 
                                  dx,dt,time,div_coef_type,anelastic_cutoff)

      implicit none
      integer, intent(in) :: lo(:), hi(:)
      real(kind=dp_t), intent(  out) :: vel(lo(2):)
      real(kind=dp_t), intent(in   ) :: p0_old(lo(2):), rho0_old(lo(2):), rhoh0_old(lo(2):)
      real(kind=dp_t), intent(  out) :: p0_new(lo(2):), rho0_new(lo(2):), rhoh0_new(lo(2):)
      real(kind=dp_t), intent(  out) :: rho0_nph(lo(2):)
      real(kind=dp_t), intent(inout) :: temp0(lo(2):)
      real(kind=dp_t), intent(inout) :: gam1(lo(2):)
      real(kind=dp_t), intent(inout) :: div_coeff_n   (lo(2):)
      real(kind=dp_t), intent(inout) :: div_coeff_nph (lo(2):)
      real(kind=dp_t), intent(inout) :: div_coeff_half(lo(2):)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) :: dx(:),dt,grav,time,anelastic_cutoff
      integer        , intent(in   ) :: div_coef_type

!     Local variables
      integer :: i, j
      real(kind=dp_t) :: coeff
      real(kind=dp_t) :: max_vel, denom
      real(kind=dp_t) :: sigma_H
      real(kind=dp_t) :: half_time

      real (kind = dp_t), allocatable :: H(:,:)
      real (kind = dp_t), allocatable :: force(:)
      real (kind = dp_t), allocatable :: edge(:)
      real (kind = dp_t), allocatable :: temp_array(:)
      real (kind = dp_t), allocatable :: temp_array_half(:)

      do_diag = .false.

      denom = ONE / dble(hi(1)-lo(1)+1)

      allocate(     force(lo(2):hi(2)))
      allocate(      edge(lo(2):hi(2)+1))
      allocate(temp_array(lo(2):hi(2)))
      allocate(temp_array_half(lo(2):hi(2)+1))

      allocate(H(lo(1):hi(1),lo(2):hi(2)))
      half_time = time + HALF * dt
      call get_H_2d(H,lo,hi,dx,half_time)

      max_vel = ZERO

      ! Initialize velocity to zero.
      vel = ZERO

      do j = lo(2),hi(2)
        sigma_H = ZERO
        do i = lo(1), hi(1)
        ! Compute the coefficient of heating in the divu expression
           den_row(1) = rhohalf(i,j)
          temp_row(1) = temp0(j)
             p_row(1) = p0_old(j)
           ! (rho,P) --> h, etc
           input_flag = 4

           call eos(input_flag, den_row, temp_row, npts, nspecies, &
                xmass, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, &
                pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
                s_row, do_diag)

           ! coeff = p_T / (rho * c_p * p_rho)
           coeff = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))

           sigma_H = sigma_H + coeff * H(i,j)
        end do
        sigma_H = sigma_H * denom
        vel(j+1) = vel(j) + sigma_H * dx(2)
        max_vel = max(max_vel, abs(vel(j+1)))

      end do

      ! HACK HACK HACK 
!     vel     = ZERO
!     max_vel = ZERO

      print *,'MAX CFL FRAC OF DISPL ',max_vel * dt / dx(2)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE P0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      force = ZERO
      call mkflux_1d(p0_old,edge,vel,force,lo(2),dx(2),dt)
      do j = lo(2), hi(2)
        p0_new(j) = p0_old(j) - dt / dx(2) * HALF * (vel(j) + vel(j+1)) * (edge(j+1) - edge(j))
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHO0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = lo(2),hi(2)
        force(j) = rho0_old(j) * (vel(j+1) - vel(j)) / dx(2)
      end do
      call mkflux_1d(rho0_old,edge,vel,force,lo(2),dx(2),dt)
      do j = lo(2), hi(2)

        rho0_new(j) = rho0_old(j) - dt / dx(2) * (edge(j+1) * vel(j+1) - edge(j) * vel(j))
        rho0_new(j) = max(rho0_new(j), rho0_old(hi(2)))

        rho0_nph(j) = HALF * (rho0_old(j) + rho0_new(j))
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     MAKE TEMP0, RHOH0 AND GAM1 FROM P0 AND RHO0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = lo(2), hi(2)

            den_row(1) = rho0_new(j)
           temp_row(1) = temp0(j)
              p_row(1) = p0_new(j)

!        ! (rho,P) --> T, h
         input_flag = 4

         call eos(input_flag, den_row, temp_row, npts, nspecies, &
                  xmass, aion, zion, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, &
                  pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
                  s_row, do_diag)

           temp0(j) = temp_row(1)
            gam1(j) = gam1_row(1)

           rhoh0_new(j) = rho0_new(j) * h_row(1)

      end do

      if (div_coef_type .eq. 2) then
        div_coeff_n   = rho0_new
        div_coeff_nph = rho0_nph
      else 
        call make_div_coeff(temp_array,temp_array_half,rho0_new,p0_new,gam1,grav,dx(2),anelastic_cutoff)
        div_coeff_nph = HALF * (temp_array + div_coeff_n)
        div_coeff_n   =         temp_array
        div_coeff_half=         temp_array_half
      end if

      deallocate(H)
      deallocate(force)
      deallocate(edge)
      deallocate(temp_array)
      deallocate(temp_array_half)

   end subroutine eval_base_state_2d

   subroutine eval_base_state_3d (vel,p0_old,p0_new,rho0_old,rho0_nph,rho0_new,temp0, &
                                  rhoh0_old,rhoh0_new, &
                                  gam1,div_coeff_n,div_coeff_nph, &
                                  div_coeff_half,rhohalf,grav,lo,hi, & 
                                  dx,dt,time,div_coef_type,anelastic_cutoff)

      implicit none
      integer, intent(in) :: lo(:), hi(:)
      real(kind=dp_t), intent(  out) :: vel(lo(3):)
      real(kind=dp_t), intent(in   ) :: p0_old(lo(3):), rho0_old(lo(3):), rhoh0_old(lo(3):)
      real(kind=dp_t), intent(  out) :: p0_new(lo(3):), rho0_new(lo(3):), rhoh0_new(lo(3):)
      real(kind=dp_t), intent(  out) :: rho0_nph(lo(3):)
      real(kind=dp_t), intent(inout) :: temp0(lo(3):)
      real(kind=dp_t), intent(inout) :: gam1(lo(3):)
      real(kind=dp_t), intent(inout) :: div_coeff_n   (lo(3):)
      real(kind=dp_t), intent(inout) :: div_coeff_nph (lo(3):)
      real(kind=dp_t), intent(inout) :: div_coeff_half(lo(3):)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(in   ) :: dx(:),dt,grav,time, anelastic_cutoff
      integer        , intent(in   ) :: div_coef_type

!     Local variables
      integer :: i, j, k
      real(kind=dp_t) :: coeff
      real(kind=dp_t) :: max_vel, denom
      real(kind=dp_t) :: sigma_H
      real(kind=dp_t) :: half_time

      real (kind = dp_t), allocatable :: H(:,:,:)
      real (kind = dp_t), allocatable :: force(:)
      real (kind = dp_t), allocatable :: edge(:)
      real (kind = dp_t), allocatable :: temp_array(:)
      real (kind = dp_t), allocatable :: temp_array_half(:)

      do_diag = .false.

      denom = ONE / dble(hi(1)-lo(1)+1)

      allocate(     force(lo(3):hi(3)))
      allocate(      edge(lo(3):hi(3)+1))
      allocate(temp_array(lo(3):hi(3)))
      allocate(temp_array_half(lo(3):hi(3)+1))

      allocate(H(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      half_time = time + HALF * dt
      call get_H_3d(H,lo,hi,dx,half_time)

      max_vel = ZERO

      ! Initialize velocity to zero.
      vel = ZERO

      do k = lo(3),hi(3)
        sigma_H = ZERO
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
        ! Compute the coefficient of heating in the divu expression
           den_row(1) = rhohalf(i,j,k)
          temp_row(1) = temp0(k)
             p_row(1) = p0_old(k)
           ! (rho,P) --> h, etc
           input_flag = 4
           call eos(input_flag, den_row, temp_row, npts, nspecies, &
                xmass, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, &
                pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
                s_row, do_diag)

           ! coeff = p_T / (rho * c_p * p_rho)
           coeff = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))

           sigma_H = sigma_H + coeff * H(i,j,k)
        end do
        end do
        sigma_H = sigma_H * denom
        vel(k+1) = vel(k) + sigma_H * dx(3)
        max_vel = max(max_vel, abs(vel(k+1)))

      end do

      ! HACK HACK HACK 
!     vel     = ZERO
!     max_vel = ZERO

      print *,'MAX CFL FRAC OF DISPL ',max_vel * dt / dx(3)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE P0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      force = ZERO
      call mkflux_1d(p0_old,edge,vel,force,lo(3),dx(3),dt)
      do k = lo(3), hi(3)
        p0_new(k) = p0_old(k) - dt / dx(3) * HALF * (vel(k) + vel(k+1)) * (edge(k+1) - edge(k))
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHO0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k = lo(3),hi(3)
        force(k) = rho0_old(k) * (vel(k+1) - vel(k)) / dx(3)
      end do
      call mkflux_1d(rho0_old,edge,vel,force,lo(3),dx(3),dt)
      do k = lo(3), hi(3)

        rho0_new(k) = rho0_old(k) - dt / dx(3) * (edge(k+1) * vel(k+1) - edge(k) * vel(k))
        rho0_new(k) = max(rho0_new(k), rho0_old(hi(3)))

        rho0_nph(k) = HALF * (rho0_old(k) + rho0_new(k))
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     MAKE TEMP0, RHOH0 AND GAM1 FROM P0 AND RHO0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k = lo(3), hi(3)

            den_row(1) = rho0_new(k)
           temp_row(1) = temp0(k)
              p_row(1) = p0_new(k)

!        ! (rho,P) --> T, h
         input_flag = 4

         call eos(input_flag, den_row, temp_row, npts, nspecies, &
                  xmass, aion, zion, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, &
                  pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
                  s_row, do_diag)

           temp0(k) = temp_row(1)
            gam1(k) = gam1_row(1)

           rhoh0_new(k) = rho0_new(k) * h_row(1)

      end do

      if (div_coef_type .eq. 2) then
        div_coeff_n   = rho0_new
        div_coeff_nph = rho0_nph
      else 
        call make_div_coeff(temp_array,temp_array_half,rho0_new,p0_new,gam1,grav,dx(2),anelastic_cutoff)
        div_coeff_nph = HALF * (temp_array + div_coeff_n)
        div_coeff_n   =         temp_array
        div_coeff_half=         temp_array_half
      end if

      deallocate(H)
      deallocate(force)
      deallocate(edge)
      deallocate(temp_array)
      deallocate(temp_array_half)

   end subroutine eval_base_state_3d

end module base_state_module
