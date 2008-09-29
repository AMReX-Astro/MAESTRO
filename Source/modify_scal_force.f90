module modify_scal_force_module

  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: modify_scal_force

contains

  subroutine modify_scal_force(force,s,umac,base,base_edge,w0,dx,base_cart, &
                               comp,mla,the_bc_level)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical, dm, nlevs
    use variables, only: foextrap_comp
    use multifab_physbc_module
    use ml_restriction_module
    use multifab_fill_ghost_module

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.
    
    type(multifab) , intent(inout) :: force(:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: base(:,0:), base_edge(:,0:), w0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(multifab) , intent(in   ) :: base_cart(:)
    integer        , intent(in   ) :: comp
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    
    ! local
    integer :: i,ng_f,ng_s,ng_um,ng_b,n
    integer :: lo(dm),hi(dm)
    integer :: domlo(dm),domhi(dm)    

    type(box) :: domain

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: bcp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "modify_scal_force")
    
    ng_s = s(1)%ng
    ng_f = force(1)%ng
    ng_um = umac(1,1)%ng
    ng_b = base_cart(1)%ng

    do n=1,nlevs

       domain = layout_get_pd(s(n)%la)
       domlo = lwb(domain)
       domhi = upb(domain)
       
       do i=1,force(n)%nboxes
          if ( multifab_remote(force(n),i) ) cycle
          fp => dataptr(force(n),i)
          sp => dataptr(s(n),i)
          ump => dataptr(umac(n,1),i)
          vmp => dataptr(umac(n,2),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call modify_scal_force_2d(fp(:,:,1,comp),ng_f,sp(:,:,1,comp),ng_s,lo,hi, &
                                       ump(:,:,1,1),vmp(:,:,1,1),ng_um,base(n,:), &
                                       base_edge(n,:),w0(n,:),dx(n,:))
          case(3)
             wmp  => dataptr(umac(n,3), i)
             if (spherical .eq. 1) then
                bcp => dataptr(base_cart(n), i)
                call modify_scal_force_3d_sphr(n,fp(:,:,:,comp),ng_f,sp(:,:,:,comp),ng_s, &
                                               lo,hi,domlo,domhi, &
                                               ump(:,:,:,1),vmp(:,:,:,1), &
                                               wmp(:,:,:,1),ng_um,bcp(:,:,:,1),ng_b, &
                                               w0(n,:),dx(n,:))
             else
                call modify_scal_force_3d_cart(fp(:,:,:,comp),ng_f,sp(:,:,:,comp),ng_s, &
                                               lo,hi,ump(:,:,:,1), &
                                               vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
                                               base(n,:),base_edge(n,:), &
                                               w0(n,:),dx(n,:))
             end if
          end select
       end do
       
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(force(nlevs),comp,1)
       
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(force(nlevs),comp,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(force(n-1),comp,force(n),comp,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(force(n),force(n-1), &
                                         force(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         comp,foextrap_comp,1)

       enddo

    end if

    call destroy(bpt)
    
  end subroutine modify_scal_force
  
  subroutine modify_scal_force_2d(force,ng_f,s,ng_s,lo,hi,umac,vmac,ng_um,base,base_edge,w0,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_s,ng_um
    real(kind=dp_t), intent(  out) :: force(lo(1)-ng_f :,lo(2)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  base(0:), base_edge(0:)
    real(kind=dp_t), intent(in   ) ::    w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer :: i,j
    real(kind=dp_t) :: divu,divbaseu
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          divu =  (umac(i+1,j) - umac(i,j)) / dx(1) &
               +( (vmac(i,j+1) + w0(j+1)) &
                 -(vmac(i,j)   + w0(j)  ) ) / dx(2)

          divbaseu = base(j)*(umac(i+1,j) - umac(i,j))/dx(1) &
                            +(vmac(i,j+1) * base_edge(j+1) - &
                              vmac(i,j  ) * base_edge(j  ) )/ dx(2)
          
          force(i,j) = force(i,j) - (s(i,j)-base(j))*divu - divbaseu 
       end do
    end do
          
  end subroutine modify_scal_force_2d
  
  subroutine modify_scal_force_3d_cart(force,ng_f,s,ng_s,lo,hi,umac,vmac,wmac,ng_um,base,base_edge,w0,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_s,ng_um
    real(kind=dp_t), intent(  out) :: force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) :: base(0:), base_edge(0:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer :: i,j,k
    real(kind=dp_t) :: divu,divbaseu
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                   +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                   +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)

             divu = divu + (w0(k+1)-w0(k))/dx(3)

             divbaseu = base(k)*( (umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                                 +(vmac(i,j+1,k) - vmac(i,j,k))/dx(2) ) &
                                 +(wmac(i,j,k+1) * base_edge(k+1) &
                                 - wmac(i,j,k  ) * base_edge(k  ))/ dx(3)

             force(i,j,k) = force(i,j,k) - (s(i,j,k)-base(k))*divu - divbaseu 
          end do
       end do
    end do
    
  end subroutine modify_scal_force_3d_cart
  
  subroutine modify_scal_force_3d_sphr(n,force,ng_f,s,ng_s,lo,hi,domlo,domhi, &
                                       umac,vmac,wmac,ng_um,base_cart,ng_b,w0,dx)

    use geometry, only: nr_fine, r_edge_loc, dr, r_cc_loc
    use fill_3d_module
    use bl_constants_module
    
    integer        , intent(in   ) :: n,lo(:),hi(:),domlo(:),domhi(:),ng_f,ng_s,ng_um,ng_b
    real(kind=dp_t), intent(  out) :: force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    
    real(kind=dp_t), intent(in   ) :: base_cart(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    ! Local variables
    integer :: i,j,k,r
    real(kind=dp_t) :: divumac,divbaseu
    real(kind=dp_t) :: base_xlo,base_xhi
    real(kind=dp_t) :: base_ylo,base_yhi
    real(kind=dp_t) :: base_zlo,base_zhi
    
    real(kind=dp_t) :: divu(0:nr_fine-1)
    real(kind=dp_t) :: divu_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)
    
    do r=0,nr_fine-1
       divu(r) = (r_edge_loc(n,r+1)**2 * w0(r+1) - &
                  r_edge_loc(n,r  )**2 * w0(r  ) ) / &
                 (dr(n)*r_cc_loc(n,r)**2)
    end do

    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,divu,divu_cart,lo,hi,dx,0,0)
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             divumac = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                      +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                      +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
             
             if (i.lt.domhi(1)) then
                base_xhi = HALF * (base_cart(i,j,k) + base_cart(i+1,j,k))
             else
                base_xhi = base_cart(i,j,k)
             end if
             if (i.gt.domlo(1)) then
                base_xlo = HALF * (base_cart(i,j,k) + base_cart(i-1,j,k))
             else
                base_xlo = base_cart(i,j,k)
             end if

             if (j.lt.domhi(2)) then
                base_yhi = HALF * (base_cart(i,j,k) + base_cart(i,j+1,k))
             else
                base_yhi = base_cart(i,j,k)
             end if
             if (j.gt.domlo(2)) then
                base_ylo = HALF * (base_cart(i,j,k) + base_cart(i,j-1,k))
             else
                base_ylo = base_cart(i,j,k)
             end if

             if (k.lt.domhi(3)) then
                base_zhi = HALF * (base_cart(i,j,k) + base_cart(i,j,k+1))
             else
                base_zhi = base_cart(i,j,k)
             end if
             if (k.gt.domlo(3)) then
                base_zlo = HALF * (base_cart(i,j,k) + base_cart(i,j,k-1))
             else
                base_zlo = base_cart(i,j,k)
             end if
             
             divbaseu = (umac(i+1,j,k)*base_xhi - umac(i,j,k)*base_xlo)/dx(1) + &
                        (vmac(i,j+1,k)*base_yhi - vmac(i,j,k)*base_ylo)/dx(2) + &
                        (wmac(i,j,k+1)*base_zhi - wmac(i,j,k)*base_zlo)/dx(3)
             
             force(i,j,k) = force(i,j,k) - divbaseu &
                  -(s(i,j,k)-base_cart(i,j,k))*(divumac+divu_cart(i,j,k,1)) 
          end do
       end do
    end do
    
  end subroutine modify_scal_force_3d_sphr
  
end module modify_scal_force_module
 
