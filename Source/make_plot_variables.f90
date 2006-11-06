
module plot_variables_module

  use bl_types
  use bc_module
  use multifab_module
  use eos_module
  use fill_3d_module
  use network
  use variables
  use geometry

  implicit none

contains

  subroutine make_deltagamma (deltagamma,comp,s,rho0,p0,t0,rhoX0)

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: deltagamma
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: rho0(:), p0(:), t0(:), rhoX0(:,:)
    
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: mp(:,:,:,:)
    real(kind=dp_t), pointer:: dp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i

    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       dp => dataptr(deltagamma, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call makedeltagamma_2d(dp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng, rho0, p0, t0, rhoX0)
       case (3)
          call makedeltagamma_3d(dp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng, rho0, p0, t0, rhoX0)
       end select
    end do

  end subroutine make_deltagamma

  subroutine makedeltagamma_2d (deltagamma,s,lo,hi,ng,rho0,p0,t0,rhoX0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    rho0(lo(2):)
    real (kind = dp_t), intent(in   ) ::      p0(lo(2):)
    real (kind = dp_t), intent(in   ) ::      t0(lo(2):)
    real (kind = dp_t), intent(in   ) ::   rhoX0(lo(2):,:)

!     Local variables
    integer :: i, j

    integer :: input_flag
    real (kind = dp_t) :: gam10
      
    do_diag = .false.

    do j = lo(2), hi(2)

       ! first compute the base state gamma1
       den_row(1) = rho0(j)
       p_row(1) = p0(j)
       temp_row(1) = t0(j)   ! used as an initial guess

       xn_zone(:) = rhoX0(j,:)/rho0(j)

       input_flag = 4

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                do_diag)

       gam10 = gam1_row(1)

       do i = lo(1), hi(1)

          ! compute the thermodynamics of the full state
          den_row(1) = s(i,j,rho_comp)
          h_row(1) = s(i,j,rhoh_comp) / s(i,j,rho_comp)
          temp_row(1) = t0(j)
            
          xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          input_flag = 2

          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)

          deltagamma(i,j) = gam1_row(1) - gam10

       enddo
    enddo

  end subroutine makedeltagamma_2d
    
  subroutine makedeltagamma_3d (deltagamma,s,lo,hi,ng,rho0,p0,t0,rhoX0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):,lo(3):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    rho0(lo(3):)
    real (kind = dp_t), intent(in   ) ::      p0(lo(3):)
    real (kind = dp_t), intent(in   ) ::      t0(lo(3):)
    real (kind = dp_t), intent(in   ) ::   rhoX0(lo(3):,:)


!     Local variables
    integer :: i, j, k
    integer :: input_flag
    real (kind = dp_t) :: gam10
      
    do_diag = .false.

    do k = lo(3), hi(3)

       ! first compute the base state gamma1
       den_row(1) = rho0(k)
       p_row(1) = p0(k)
       temp_row(1) = t0(k)   ! used as an initial guess

       xn_zone(:) = rhoX0(k,:)/rho0(k)

       input_flag = 4

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                do_diag)

       gam10 = gam1_row(1)

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! compute the thermodynamics of the full state
             den_row(1) = s(i,j,k,rho_comp)
             h_row(1) = s(i,j,k,rhoh_comp) / s(i,j,k,rho_comp)
             temp_row(1) = t0(k)
            
             xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)
             
             input_flag = 2

             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, &
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             deltagamma(i,j,k) = gam1_row(1) - gam10

          enddo
       enddo
    enddo

  end subroutine makedeltagamma_3d

  subroutine make_deltap (deltap,comp,s,p0,temp0)

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: deltap
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: p0(:),temp0(:)
    
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: mp(:,:,:,:)
    real(kind=dp_t), pointer:: dp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i

    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       dp => dataptr(deltap, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call makedeltap_2d(dp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng, p0, temp0)
       case (3)
          call makedeltap_3d(dp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng, p0, temp0)
       end select
    end do

  end subroutine make_deltap

  subroutine makedeltap_2d (deltap,s,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: deltap(lo(1):,lo(2):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::      p0(lo(2):)
    real (kind = dp_t), intent(in   ) ::   temp0(lo(2):)
    
!     Local variables
    integer :: i, j
    
    integer :: input_flag
      
    do_diag = .false.
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          den_row(1) = s(i,j,rho_comp)
          p_row(1) = p0(j)
          h_row(1) = s(i,j,rhoh_comp) / s(i,j,rho_comp)
          xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          ! (rho,h) --> T,p, etc
          input_flag = 2

          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)

          deltap(i,j) = (p_row(1)-p0(j))/ p0(j)

       enddo
    enddo

  end subroutine makedeltap_2d

  subroutine makedeltap_3d (deltap,s,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: deltap(lo(1):,lo(2):,lo(3):)  
    real (kind = dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::      p0(lo(3):)
    real (kind = dp_t), intent(in   ) ::   temp0(lo(3):)

!     Local variables
    integer :: i, j, k

    integer :: input_flag
    
    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_row(1) = s(i,j,k,rho_comp)
             p_row(1) = p0(k)
             h_row(1) = s(i,j,k,rhoh_comp) / s(i,j,k,rho_comp)
             xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

             ! (rho,h) --> T,p, etc
             input_flag = 2

             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, &
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             deltap(i,j,k) = (p_row(1)-p0(k))/ p0(k)

          enddo
       enddo
    enddo
    
  end subroutine makedeltap_3d

  subroutine make_enthalpy (enthalpy,comp,s)

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: enthalpy
    type(multifab) , intent(in   ) :: s

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: rp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i

    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       rp => dataptr(enthalpy, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call make_enthalpy_2d(rp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng)
       case (3)
          call make_enthalpy_3d(rp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng)
       end select
    end do

  end subroutine make_enthalpy

  subroutine make_enthalpy_2d (enthalpy,s,lo,hi,ng)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: enthalpy(lo(1):,lo(2):)  
    real (kind = dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:,:)

!     Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          enthalpy(i,j) = s(i,j,rhoh_comp) / s(i,j,rho_comp)
       enddo
    enddo

  end subroutine make_enthalpy_2d

  subroutine make_enthalpy_3d (enthalpy,s,lo,hi,ng)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: enthalpy(lo(1):,lo(2):, lo(3):)  
    real (kind = dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:, lo(3)-ng:,:)
    
!     Local variables
    integer :: i, j, k
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             enthalpy(i,j,k) = s(i,j,k,rhoh_comp) / s(i,j,k,rho_comp)
          enddo
       enddo
    end do
    
  end subroutine make_enthalpy_3d

   subroutine make_machno (machno,comp,u,s,p0,temp0)

      integer        , intent(in   ) :: comp
      type(multifab) , intent(inout) :: machno
      type(multifab) , intent(in   ) :: u,s
      real(kind=dp_t), intent(in   ) :: p0(:),temp0(:)

      real(kind=dp_t), pointer:: up(:,:,:,:)
      real(kind=dp_t), pointer:: sp(:,:,:,:)
      real(kind=dp_t), pointer:: mp(:,:,:,:)
      real(kind=dp_t), pointer:: dp(:,:,:,:)
      integer :: lo(s%dim),hi(s%dim),ng,dm
      integer :: i

      ng = s%ng
      dm = s%dim

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
         up => dataptr(u, i)
         sp => dataptr(s, i)
         mp => dataptr(machno, i)
         lo =  lwb(get_box(s, i))
         hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call makemachno_2d(mp(:,:,1,comp),up(:,:,1,:),sp(:,:,1,:), lo, hi, ng, p0, temp0)
            case (3)
              call makemachno_3d(mp(:,:,:,comp),up(:,:,:,:),sp(:,:,:,:), lo, hi, ng, p0, temp0)
         end select
      end do

   end subroutine make_machno

   subroutine makemachno_2d (machno,u,s,lo,hi,ng,p0,temp0)

      implicit none
      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: machno(lo(1):,lo(2):)  
      real (kind = dp_t), intent(in   ) ::       u(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(in   ) ::      p0(lo(2):)
      real (kind = dp_t), intent(in   ) ::   temp0(lo(2):)

!     Local variables
      integer :: i, j
      real (kind = dp_t) :: vel

      integer :: input_flag
      integer :: nx
      
      do_diag = .false.

      do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_row(1) = s(i,j,rho_comp)
             p_row(1) = p0(j)

             xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

             ! (rho,P) --> T,h, etc
             input_flag = 4

             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, &
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             vel = sqrt(u(i,j,1)*u(i,j,1) + u(i,j,2)*u(i,j,2))
             machno(i,j) = vel / cs_row(1)

          enddo
      enddo

   end subroutine makemachno_2d

   subroutine makemachno_3d (machno,u,s,lo,hi,ng,p0,temp0)

      implicit none
      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: machno(lo(1):,lo(2):,lo(3):)  
      real (kind = dp_t), intent(in   ) ::       u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::      p0(lo(3):)
      real (kind = dp_t), intent(in   ) ::   temp0(lo(3):)

!     Local variables
      integer :: i, j, k
      real (kind = dp_t) :: vel

      integer :: input_flag
      integer :: nx
      
      do_diag = .false.

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               den_row(1) = s(i,j,k,rho_comp)
               p_row(1) = p0(k)

               xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

               ! (rho,P) --> T,h, etc
               input_flag = 4

               call eos(input_flag, den_row, temp_row, &
                        npts, nspec, &
                        xn_zone, aion, zion, &
                        p_row, h_row, e_row, &
                        cv_row, cp_row, xne_row, eta_row, pele_row, &
                        dpdt_row, dpdr_row, dedt_row, dedr_row, &
                        dpdX_row, dhdX_row, &
                        gam1_row, cs_row, s_row, &
                        do_diag)

               vel = sqrt(u(i,j,k,1)*u(i,j,k,1) + u(i,j,k,2)*u(i,j,k,2) + u(i,j,k,3)*u(i,j,k,3))
               machno(i,j,k) = vel / cs_row(1)

            enddo
         enddo
      enddo

   end subroutine makemachno_3d

  subroutine make_rhopert (rhopert,comp,s,s0,dx)

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: rhopert
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: s0(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: pp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i
    
    ng = s%ng
    dm = s%dim
    
    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       pp => dataptr(rhopert, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call makerhopert_2d(pp(:,:,1,comp), sp(:,:,1,rho_comp), &
                              lo, hi, ng, s0(:,rho_comp))
       case (3)
          if (spherical .eq. 0) then
            call makerhopert_3d_cart(pp(:,:,:,comp), sp(:,:,:,rho_comp), &
                                     lo, hi, ng, s0(:,rho_comp))
          else
            call makerhopert_3d_sphr(pp(:,:,:,comp), sp(:,:,:,rho_comp), &
                                     lo, hi, ng, s0, dx)
          end if
       end select
    end do

  end subroutine make_rhopert

  subroutine makerhopert_2d (rhopert,s,lo,hi,ng,rho0)

    implicit none
    
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: rhopert(lo(1):,lo(2):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:)
    real (kind = dp_t), intent(in   ) :: rho0(lo(2):)

!     Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rhopert(i,j) = s(i,j) - rho0(j)
       enddo
    enddo
    
  end subroutine makerhopert_2d
  
  subroutine makerhopert_3d_cart (rhopert,s,lo,hi,ng,rho0)

    implicit none
    
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: rhopert(lo(1):,lo(2):, lo(3):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:, lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: rho0(lo(3):)
    
    !     Local variables
    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhopert(i,j,k) = s(i,j,k) - rho0(k)
          enddo
       enddo
    end do
    
  end subroutine makerhopert_3d_cart
  
  subroutine makerhopert_3d_sphr (rhopert,s,lo,hi,ng,s0,dx)

    implicit none
    
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: rhopert(lo(1):,lo(2):, lo(3):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:, lo(3)-ng:)
    real (kind = dp_t), intent(in   ) ::      s0(lo(3):,:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    
    !     Local variables
    integer :: i, j, k
    real (kind = dp_t), allocatable :: rho0_cart(:,:,:)

    allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

    call fill_3d_data(rho0_cart,s0(:,rho_comp),dx,0)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhopert(i,j,k) = s(i,j,k) - rho0_cart(i,j,k)
          enddo
       enddo
    end do

    deallocate(rho0_cart)
    
  end subroutine makerhopert_3d_sphr

  subroutine make_tfromH (T,comp,state,p0,temp0)

    integer        , intent(inout) :: comp
    type(multifab) , intent(inout) :: T
    type(multifab) , intent(in   ) :: state
    real(kind=dp_t), intent(in   ) ::    p0(:)
    real(kind=dp_t), intent(in   ) :: temp0(:)

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(state%dim),hi(state%dim),ng,dm
    integer :: i

    ng = state%ng
    dm = state%dim

    do i = 1, state%nboxes
       if ( multifab_remote(state, i) ) cycle
       sp => dataptr(state, i)
       tp => dataptr(T    , i)
       lo =  lwb(get_box(state, i))
       hi =  upb(get_box(state, i))
       select case (dm)
       case (2)
          call maketfromH_2d(tp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng, p0, temp0)
       case (3)
          call maketfromH_3d(tp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng, p0, temp0)
       end select
    end do

  end subroutine make_tfromH

  subroutine maketfromH_2d (T,state,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::     T(lo(1)   :,lo(2):)  
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    p0(lo(2):)
    real (kind = dp_t), intent(in   ) :: temp0(lo(2):)

    !     Local variables
    integer :: i, j

    integer :: nx

    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T
            
          den_row(1)  = state(i,j,rho_comp)
          h_row(1)    = state(i,j,rhoh_comp) / state(i,j,rho_comp)
          p_row(1)    = p0(j)
          temp_row(1) = temp0(j)
          xn_zone(:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          input_flag = 2
  
          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)
          
          T(i,j) = log(temp_row(1))/log(10.)
          
       enddo
    enddo

  end subroutine maketfromH_2d

  subroutine maketfromH_3d (T,state,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::     T(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    p0(lo(3):)
    real (kind = dp_t), intent(in   ) :: temp0(lo(3):)

    !     Local variables
    integer :: i, j, k

    integer :: nx

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, H) --> T
            
             den_row(1)  = state(i,j,k,rho_comp)
             h_row(1)    = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             p_row(1)    = p0(k)
             temp_row(1) = temp0(k)
             xn_zone(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

             input_flag = 2
  
             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, &
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             T(i,j,k) = log(temp_row(1))/log(10.)

          enddo
       enddo
    enddo

  end subroutine maketfromH_3d

  subroutine make_tfromrho (t,comp,s,t0,p0,time,dx)

    integer        , intent(inout) :: comp
    type(multifab) , intent(inout) :: t
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: t0(:),p0(:)
    real(kind=dp_t), intent(in   ) :: time,dx(:)

    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i

    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       tp => dataptr(t, i)
       sp => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call maketfromrho_2d(tp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng, t0, p0, time, dx)
       case (3)
          call maketfromrho_3d(tp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng, t0, p0, time, dx)
       end select
    end do

  end subroutine make_tfromrho

  subroutine maketfromrho_2d (t,s,lo,hi,ng,t0,p0,time,dx)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out) :: t(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) :: t0(lo(2):)
    real (kind=dp_t), intent(in   ) :: p0(lo(2):)
    real (kind=dp_t), intent(in   ) :: time,dx(:)

    !     Local variables
    integer :: i, j

    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          den_row(1) = s(i,j,rho_comp)
          temp_row(1) = t0(j)
          p_row(1) = p0(j)
          xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          ! (rho,P) --> T,h
          input_flag = 4

          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, & 
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)
          
          t(i,j) = log(temp_row(1))/log(10.)
       enddo
    enddo

  end subroutine maketfromrho_2d

  subroutine maketfromrho_3d (t,s,lo,hi,ng,t0,p0,time,dx)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out) :: t(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: t0(lo(3):)
    real (kind=dp_t), intent(in   ) :: p0(lo(3):)
    real (kind=dp_t), intent(in   ) :: time,dx(:)

    !     Local variables
    integer :: i, j, k

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_row(1) = s(i,j,k,rho_comp)
             temp_row(1) = t0(k)
             p_row(1) = p0(k)
             xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

             ! (rho,P) --> T,h
             input_flag = 4

             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, & 
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             t(i,j,k) = log(temp_row(1))/log(10.)

          enddo
       enddo
    enddo

  end subroutine maketfromrho_3d

  subroutine make_tpert (T,comp,state,p0,temp0)

    integer        , intent(inout) :: comp
    type(multifab) , intent(inout) :: T
    type(multifab) , intent(in   ) :: state
    real(kind=dp_t), intent(in   ) ::    p0(:)
    real(kind=dp_t), intent(in   ) :: temp0(:)
    
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(state%dim),hi(state%dim),ng,dm
    integer :: i
    
    ng = state%ng
    dm = state%dim

    do i = 1, state%nboxes
       if ( multifab_remote(state, i) ) cycle
       sp => dataptr(state, i)
       tp => dataptr(T    , i)
       lo =  lwb(get_box(state, i))
       hi =  upb(get_box(state, i))
       select case (dm)
       case (2)
          call maketpert_2d(tp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng, p0, temp0)
       case (3)
          call maketpert_3d(tp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng, p0, temp0)
       end select
    end do

  end subroutine make_tpert

  subroutine maketpert_2d (T,state,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::     T(lo(1)   :,lo(2):)  
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    p0(lo(2):)
    real (kind = dp_t), intent(in   ) :: temp0(lo(2):)

    integer :: i, j

    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
            
          den_row(1) = state(i,j,rho_comp)
          temp_row(1) = temp0(j)
          p_row(1) = p0(j)
          xn_zone(:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          ! (rho, P) --> T
          input_flag = 4

          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)
          
          T(i,j) = temp_row(1) - temp0(j)
       enddo
    enddo

  end subroutine maketpert_2d

  subroutine maketpert_3d (T,state,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::     T(lo(1)   :,lo(2):   ,lo(3)   :  )  
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    p0(lo(3):)
    real (kind = dp_t), intent(in   ) :: temp0(lo(3):)

    integer :: i, j, k

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            
             den_row(1) = state(i,j,k,rho_comp)
             temp_row(1) = temp0(k)
             p_row(1) = p0(k)
             xn_zone(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)
                          
             ! (rho, P) --> T
             input_flag = 4

             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, &
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             T(i,j,k) = temp_row(1) - temp0(k)
          enddo
       enddo
    enddo

  end subroutine maketpert_3d

   subroutine make_XfromrhoX (plotdata,comp,s)

      integer        , intent(in   ) :: comp
      type(multifab) , intent(in   ) :: s
      type(multifab) , intent(inout) :: plotdata

      real(kind=dp_t), pointer:: sp(:,:,:,:)
      real(kind=dp_t), pointer:: pp(:,:,:,:)
      integer :: lo(s%dim),hi(s%dim),ng,dm
      integer :: i

      ng = s%ng
      dm = s%dim

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
         sp => dataptr(s, i)
         pp => dataptr(plotdata, i)
         lo =  lwb(get_box(s, i))
         hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call makeXfromrhoX_2d(pp(:,:,1,comp:),sp(:,:,1,1),sp(:,:,1,spec_comp:), lo, hi, ng)
            case (3)
              call makeXfromrhoX_3d(pp(:,:,:,comp:),sp(:,:,:,1),sp(:,:,:,spec_comp:), lo, hi, ng)
         end select
      end do

   end subroutine make_XfromrhoX

   subroutine makeXfromrhoX_2d (X,rho,rhoX,lo,hi,ng)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) ::    X(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng:,lo(2)-ng:  )

!     Local variables
      integer :: i, j, n

      do n = 1, nspec
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         X(i,j,n) = rhoX(i,j,n) / rho(i,j)
      enddo
      enddo
      enddo

   end subroutine makeXfromrhoX_2d

   subroutine makeXfromrhoX_3d (X,rho,rhoX,lo,hi,ng)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) ::    X(lo(1)   :,lo(2)   :,lo(3)   :,:)  
      real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:  )

!     Local variables
      integer :: i, j, k, n

      do n = 1, nspec
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         X(i,j,k,n) = rhoX(i,j,k,n) / rho(i,j,k)
      enddo
      enddo
      enddo
      enddo

   end subroutine makeXfromrhoX_3d

end module plot_variables_module
