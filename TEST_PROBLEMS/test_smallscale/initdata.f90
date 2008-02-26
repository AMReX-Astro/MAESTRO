module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use fill_3d_module
  use eos_module
  use variables
  use network
  use geometry
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initscalardata, initveldata, scalar_diags

contains

  subroutine initscalardata(nlevs,s,s0,p0,dx,bc,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(inout) :: s0(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla


    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(s(1)%dim),hi(s(1)%dim),ng,dm,i,n

    ng = s(1)%ng
    dm = s(1)%dim

    do n=1,nlevs

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0(n,:,:))
          case (3)
             call initscalardata_3d(n,sop(:,:,:,:), lo, hi, ng, dx(n,:), s0(n,:,:))
          end select
       end do

       call multifab_fill_boundary(s(n))
       call multifab_physbc(s(n),rho_comp,dm+rho_comp,nscal,bc(n))

       ! set base state rho, rhoh, and species to zero
       ! we do not zero base state temperature because it is used as an initial
       ! guess to the eos later in the code
       do i=0,nr(n)-1
          s0(n,i,rho_comp) = 0.d0
          s0(n,i,rhoh_comp) = 0.d0
          s0(n,i,spec_comp:spec_comp+nspec-1) = 0.d0
       end do

    end do ! end loop over levels

    do n=nlevs,2,-1
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                      bc(n-1),bc(n),1,dm+rho_comp,nscal)
    enddo

  end subroutine initscalardata

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real(kind=dp_t), intent(inout) ::    s0(0:,:)

    !     Local variables
    integer :: i,j

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          s(i,j,rho_comp) = s0(j,rho_comp)
          s(i,j,rhoh_comp) = s0(j,rhoh_comp)
          s(i,j,temp_comp) = s0(j,temp_comp)
          s(i,j,spec_comp:spec_comp+nspec-1) = s0(j,spec_comp:spec_comp+nspec-1)
          s(i,j,trac_comp:trac_comp+ntrac-1) = s0(j,trac_comp:trac_comp+ntrac-1)
       enddo
    enddo

  end subroutine initscalardata_2d

  subroutine initscalardata_3d(n,s,lo,hi,ng,dx,s0)

    integer, intent(in) :: n, lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real(kind=dp_t), intent(inout) ::    s0(0:,:)

    !     Local variables
    integer :: i, j, k, comp

    ! initial the domain with the base state
    s = ZERO

    if (spherical .eq. 1) then
       ! initialize the scalars
       call fill_3d_data(n,s(:,:,:,rho_comp), s0(:,rho_comp), lo,hi,dx,ng)
       call fill_3d_data(n,s(:,:,:,rhoh_comp),s0(:,rhoh_comp),lo,hi,dx,ng)
       call fill_3d_data(n,s(:,:,:,temp_comp),s0(:,temp_comp),lo,hi,dx,ng)
       do comp = spec_comp, spec_comp+nspec-1
          call fill_3d_data(n,s(:,:,:,comp),s0(:,comp),lo,hi,dx,ng)
       end do
       do comp = trac_comp, trac_comp+ntrac-1
          call fill_3d_data(n,s(:,:,:,comp),s0(:,comp),lo,hi,dx,ng)
       end do
    else 
       ! initialize the scalars
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                s(i,j,k,rho_comp) = s0(k,rho_comp)
                s(i,j,k,rhoh_comp) = s0(k,rhoh_comp)
                s(i,j,k,temp_comp) = s0(k,temp_comp)
                s(i,j,k,spec_comp:spec_comp+nspec-1) = s0(k,spec_comp:spec_comp+nspec-1)
                s(i,j,k,trac_comp:trac_comp+ntrac-1) = s0(k,trac_comp:trac_comp+ntrac-1)
             enddo
          enddo
       enddo
    end if

  end subroutine initscalardata_3d

  subroutine initveldata(nlevs,u,s0,p0,dx,bc,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla


    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(u(1)%dim),hi(u(1)%dim),ng,dm
    integer :: i,n

    ng = u(1)%ng
    dm = u(1)%dim

    do n=1,nlevs

       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), s0(n,:,:))
          case (3)
             call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), s0(n,:,:))
          end select
       end do

       call multifab_fill_boundary(u(n))
       call multifab_physbc(u(n),1,1,dm,bc(n))

    enddo

    do n=nlevs,2,-1
       call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                      bc(n-1),bc(n),1,1,dm)
    enddo

  end subroutine initveldata

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0)

    use probin_module, only: prob_lo_y

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)

    ! local
    integer ndum, j, dm
    parameter (ndum = 31)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum)
    real(kind=dp_t) :: loloc,hiloc,flameloc

    dm = size(dx)

    lamsolfile = 'flame_4.e7_screen_left.out'

    flameloc = ONE

    do j=lo(2),hi(2)

       loloc = prob_lo_y +  dble(j)     *dx(dm) - flameloc
       hiloc = prob_lo_y + (dble(j)+ONE)*dx(dm) - flameloc

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       u(lo(1):hi(1),j,1) = 0.0d0
       u(lo(1):hi(1),j,2) = state1d(2)

    enddo

  end subroutine initveldata_2d

  subroutine initveldata_3d(u,lo,hi,ng,dx,s0)

    use probin_module, only: prob_lo_z

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)

    ! local
    integer ndum, k, dm
    parameter (ndum = 31)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum)
    real(kind=dp_t) :: loloc,hiloc

    dm = size(dx)

    lamsolfile = 'flame_4.e7_screen_left.out'

    do k=lo(3),hi(3)

       loloc = prob_lo_z +  dble(k)     *dx(dm)
       hiloc = prob_lo_z + (dble(k)+ONE)*dx(dm)

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       u(lo(1):hi(1),lo(2):hi(2),k,1:2) = 0.0d0
       u(lo(1):hi(1),lo(2):hi(2),k,3) = state1d(2)

    enddo

  end subroutine initveldata_3d

  subroutine scalar_diags(istep,s,s0,dx)

    integer        , intent(in   ) :: istep
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in)    :: s0(:,:)
    real(kind=dp_t), intent(in)    :: dx(:)

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i,n

    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sop => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))

       select case (dm)
       case (2)
          call scalar_diags_2d(istep, sop(:,:,1,:), lo, hi, ng, dx, s0)
       case (3)
          ! 3d case not written yet
       end select
    end do

  end subroutine scalar_diags

  subroutine scalar_diags_2d(istep, s,lo,hi,ng,dx,s0)

    integer, intent(in) :: istep, lo(:), hi(:), ng
    real (kind = dp_t), intent(in) ::  s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in) :: dx(:)
    real(kind=dp_t)   , intent(in) :: s0(0:,:)

    ! Local variables
    integer :: i, j, n
    real(kind=dp_t) :: fac, stot, smax
    character(len=11) :: file_name

    write(unit=file_name,fmt='("rhodiag",i4.4)') istep
    open(90,file=file_name)

    fac = ONE / dble(hi(1)-lo(1)+1)
    do j = lo(2), hi(2)
       stot = ZERO
       smax = ZERO
       do i = lo(1), hi(1)
          stot = stot + (s(i,j,rho_comp) - s0(j,rho_comp))
          smax = max(smax,abs(s(i,j,rho_comp) - s0(j,rho_comp)))
       enddo
       write(90,*) j,stot*fac/ s0(j,rho_comp), smax / s0(j,rho_comp)
    enddo

  end subroutine scalar_diags_2d

end module init_module
