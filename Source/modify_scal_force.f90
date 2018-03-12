module modify_scal_force_module

  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: modify_scal_force

contains

  subroutine modify_scal_force(force,s,umac,s0,s0_edge,w0,dx,s0_cart, &
                               comp,mla,the_bc_level,fullform)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical
    use variables, only: foextrap_comp
    use ml_restrict_fill_module

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.
    
    type(multifab) , intent(inout) :: force(:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:), s0_edge(:,0:), w0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(multifab) , intent(in   ) :: s0_cart(:)
    integer        , intent(in   ) :: comp
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    logical, optional, intent(in ) :: fullform

    ! local
    integer :: i,ng_f,ng_s,ng_um,ng_b,n
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: domlo(mla%dim),domhi(mla%dim),dm,nlevs

    type(box) :: domain

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: bcp(:,:,:,:)

    logical :: do_fullform

    type(bl_prof_timer), save :: bpt

    call build(bpt, "modify_scal_force")

    do_fullform = .false. ; if (present(fullform)    ) do_fullform = fullform

    
    dm = mla%dim
    nlevs = mla%nlevel

    ng_s  = nghost(s(1))
    ng_f  = nghost(force(1))
    ng_um = nghost(umac(1,1))
    ng_b  = nghost(s0_cart(1))

    do n=1,nlevs

       domain = get_pd(get_layout(s(n)))
       domlo  = lwb(domain)
       domhi  = upb(domain)
       
       do i=1, nfabs(force(n))
          fp => dataptr(force(n),i)
          sp => dataptr(s(n),i)
          ump => dataptr(umac(n,1),i)
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))
          select case (dm)
          case (1)
             call modify_scal_force_1d(fp(:,1,1,comp),ng_f,sp(:,1,1,comp),ng_s,lo,hi, &
                                       ump(:,1,1,1),ng_um,s0(n,:), &
                                       s0_edge(n,:),w0(n,:),dx(n,:),do_fullform)
          case (2)
             vmp => dataptr(umac(n,2),i)
             call modify_scal_force_2d(fp(:,:,1,comp),ng_f,sp(:,:,1,comp),ng_s,lo,hi, &
                                       ump(:,:,1,1),vmp(:,:,1,1),ng_um,s0(n,:), &
                                       s0_edge(n,:),w0(n,:),dx(n,:),do_fullform)
          case(3)
             vmp => dataptr(umac(n,2),i)
             wmp  => dataptr(umac(n,3), i)
             if (spherical .eq. 1) then
                bcp => dataptr(s0_cart(n), i)
                call modify_scal_force_3d_sphr(fp(:,:,:,comp),ng_f,sp(:,:,:,comp),ng_s, &
                                               lo,hi,domlo,domhi, &
                                               ump(:,:,:,1),vmp(:,:,:,1), &
                                               wmp(:,:,:,1),ng_um,bcp(:,:,:,1),ng_b, &
                                               w0(1,:),dx(n,:),do_fullform)
             else
                call modify_scal_force_3d_cart(fp(:,:,:,comp),ng_f,sp(:,:,:,comp),ng_s, &
                                               lo,hi,ump(:,:,:,1), &
                                               vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
                                               s0(n,:),s0_edge(n,:), &
                                               w0(n,:),dx(n,:),do_fullform)
             end if
          end select
       end do
       
    enddo

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,force,mla%mba%rr,the_bc_level, &
                              icomp=comp, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=force(1)%ng)

    call destroy(bpt)
    
  end subroutine modify_scal_force
  
  subroutine modify_scal_force_1d(force,ng_f,s,ng_s,lo,hi,umac,ng_um, &
                                  s0,s0_edge,w0,dx,do_fullform)

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_s,ng_um
    real(kind=dp_t), intent(  out) :: force(lo(1)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s :)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  s0(0:), s0_edge(0:)
    real(kind=dp_t), intent(in   ) ::    w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical,         intent(in   ) :: do_fullform
    
    integer :: i
    real(kind=dp_t) :: divu,divs0u
    
    do i = lo(1),hi(1)

       ! umac does not contain w0
       divu =  (umac(i+1) - umac(i)) / dx(1)
          
       ! add w0 contribution
       divu = divu + (w0(i+1)-w0(i))/dx(1)

       if (do_fullform) then
          force(i) = force(i) - s(i)*divu

       else
          
          divs0u = (umac(i+1) * s0_edge(i+1) - &
                      umac(i  ) * s0_edge(i  ) )/ dx(1)
       
          force(i) = force(i) - (s(i)-s0(i))*divu - divs0u 
       endif

    end do
          
  end subroutine modify_scal_force_1d
  
  subroutine modify_scal_force_2d(force,ng_f,s,ng_s,lo,hi,umac,vmac,ng_um, &
                                  s0,s0_edge,w0,dx,do_fullform)

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_s,ng_um
    real(kind=dp_t), intent(  out) :: force(lo(1)-ng_f :,lo(2)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  s0(0:), s0_edge(0:)
    real(kind=dp_t), intent(in   ) ::    w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical        , intent(in   ) :: do_fullform
    
    integer :: i,j
    real(kind=dp_t) :: divu,divs0u
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          ! umac does not contain w0
          divu =  (umac(i+1,j) - umac(i,j)) / dx(1) &
                 +(vmac(i,j+1) - vmac(i,j)) / dx(2)

          ! add w0 contribution
          divu = divu + (w0(j+1)-w0(j))/dx(2)


          if (do_fullform) then
             force(i,j) = force(i,j) - s(i,j)*divu 
          
          else
             divs0u = s0(j)*(umac(i+1,j) - umac(i,j))/dx(1) &
                               +(vmac(i,j+1) * s0_edge(j+1) - &
                                 vmac(i,j  ) * s0_edge(j  ) )/ dx(2)
          
             force(i,j) = force(i,j) - (s(i,j)-s0(j))*divu - divs0u 
          endif

       end do
    end do
          
  end subroutine modify_scal_force_2d
  
  subroutine modify_scal_force_3d_cart(force,ng_f,s,ng_s,lo,hi,umac,vmac,wmac,ng_um, &
                                       s0,s0_edge,w0,dx,do_fullform)

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_s,ng_um
    real(kind=dp_t), intent(  out) :: force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) :: s0(0:), s0_edge(0:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical        , intent(in   ) :: do_fullform
    
    integer :: i,j,k
    real(kind=dp_t) :: divu,divs0u
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,divu,divs0u)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             ! umac does not contain w0
             divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                   +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                   +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)

             ! add w0 contribution
             divu = divu + (w0(k+1)-w0(k))/dx(3)

             if (do_fullform) then
                force(i,j,k) = force(i,j,k) - s(i,j,k)*divu

             else

                divs0u = s0(k)*( (umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                                    +(vmac(i,j+1,k) - vmac(i,j,k))/dx(2) ) &
                                    +(wmac(i,j,k+1) * s0_edge(k+1) &
                                    - wmac(i,j,k  ) * s0_edge(k  ))/ dx(3)

                force(i,j,k) = force(i,j,k) - (s(i,j,k)-s0(k))*divu - divs0u 
             endif

          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine modify_scal_force_3d_cart
  
  subroutine modify_scal_force_3d_sphr(force,ng_f,s,ng_s,lo,hi,domlo,domhi, &
                                       umac,vmac,wmac,ng_um,s0_cart,ng_b,w0,dx,do_fullform)

    use geometry, only: nr_fine, r_edge_loc, dr, r_cc_loc
    use fill_3d_module
    use bl_constants_module
    
    integer        , intent(in   ) :: lo(:),hi(:),domlo(:),domhi(:),ng_f,ng_s,ng_um,ng_b
    real(kind=dp_t), intent(  out) :: force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    
    real(kind=dp_t), intent(in   ) :: s0_cart(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical        , intent(in   ) :: do_fullform
    
    ! Local variables
    integer :: i,j,k,r
    real(kind=dp_t) :: divumac,divs0u
    real(kind=dp_t) :: s0_xlo,s0_xhi
    real(kind=dp_t) :: s0_ylo,s0_yhi
    real(kind=dp_t) :: s0_zlo,s0_zhi
    
    real(kind=dp_t) :: divu(0:nr_fine-1)
    real(kind=dp_t), allocatable :: divu_cart(:,:,:,:)

    allocate(divu_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    !$OMP PARALLEL DO PRIVATE(r)    
    do r=0,nr_fine-1
       divu(r) = (r_edge_loc(1,r+1)**2 * w0(r+1) - &
                  r_edge_loc(1,r  )**2 * w0(r  ) ) / &
                 (dr(1)*r_cc_loc(1,r)**2)
    end do
    !$OMP END PARALLEL DO

    ! compute w0 contribution to divu
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,divu,divu_cart,lo,hi,dx,0)
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,divumac,s0_xhi,s0_xlo,s0_yhi,s0_ylo) &
    !$OMP PRIVATE(s0_zhi,s0_zlo,divs0u)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             ! umac does not contain w0
             divumac = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                      +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                      +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)


             if (do_fullform) then

                force(i,j,k) = force(i,j,k) - s(i,j,k)*(divumac+divu_cart(i,j,k,1))                 
             else
             
                if (i.lt.domhi(1)) then
                   s0_xhi = HALF * (s0_cart(i,j,k) + s0_cart(i+1,j,k))
                else
                   s0_xhi = s0_cart(i,j,k)
                end if
                if (i.gt.domlo(1)) then
                   s0_xlo = HALF * (s0_cart(i,j,k) + s0_cart(i-1,j,k))
                else
                   s0_xlo = s0_cart(i,j,k)
                end if
                
                if (j.lt.domhi(2)) then
                   s0_yhi = HALF * (s0_cart(i,j,k) + s0_cart(i,j+1,k))
                else
                   s0_yhi = s0_cart(i,j,k)
                end if
                if (j.gt.domlo(2)) then
                   s0_ylo = HALF * (s0_cart(i,j,k) + s0_cart(i,j-1,k))
                else
                   s0_ylo = s0_cart(i,j,k)
                end if
                
                if (k.lt.domhi(3)) then
                   s0_zhi = HALF * (s0_cart(i,j,k) + s0_cart(i,j,k+1))
                else
                   s0_zhi = s0_cart(i,j,k)
                end if
                if (k.gt.domlo(3)) then
                   s0_zlo = HALF * (s0_cart(i,j,k) + s0_cart(i,j,k-1))
                else
                   s0_zlo = s0_cart(i,j,k)
                end if
                
                divs0u = (umac(i+1,j,k)*s0_xhi - umac(i,j,k)*s0_xlo)/dx(1) + &
                           (vmac(i,j+1,k)*s0_yhi - vmac(i,j,k)*s0_ylo)/dx(2) + &
                           (wmac(i,j,k+1)*s0_zhi - wmac(i,j,k)*s0_zlo)/dx(3)
             
                force(i,j,k) = force(i,j,k) - divs0u &
                     -(s(i,j,k)-s0_cart(i,j,k))*(divumac+divu_cart(i,j,k,1)) 

             endif

          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    deallocate(divu_cart)

  end subroutine modify_scal_force_3d_sphr
  
end module modify_scal_force_module
 
