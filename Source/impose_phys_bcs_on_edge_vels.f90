module impose_phys_bcs_on_edges_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module

  implicit none

  private

  public :: impose_phys_bcs_on_edges

contains

  subroutine impose_phys_bcs_on_edges(u,uedge,the_bc_level)

    use bl_prof_module
    use geometry, only: dm, nlevs

    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: uedge(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    real(kind=dp_t), pointer :: utp(:,:,:,:)
    real(kind=dp_t), pointer :: vtp(:,:,:,:)
    real(kind=dp_t), pointer :: wtp(:,:,:,:)
    integer                  :: lo(dm),hi(dm)
    integer                  :: i,n,ng_ut

    type(bl_prof_timer), save :: bpt

    call build(bpt, "impose_phys_bcs_on_edges")

    ng_ut = uedge(1,1)%ng

    do n=1,nlevs

       do i=1,u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          utp => dataptr(uedge(n,1),i)
          vtp => dataptr(uedge(n,2),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call impose_phys_bcs_2d( &
                              utp(:,:,1,1), vtp(:,:,1,1), ng_ut, lo,hi, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          case (3)
             wtp => dataptr(uedge(n,3), i)
             call impose_phys_bcs_3d( &
                              utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                              lo, hi, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          end select
       end do

    end do

    call destroy(bpt)

  end subroutine impose_phys_bcs_on_edges

  subroutine impose_phys_bcs_2d(uedge,vedge,ng_ut,lo,hi,phys_bc)

    use bc_module

    integer,         intent(in   ) :: lo(:),hi(:),ng_ut
    real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vedge(lo(1)-ng_ut:,lo(2)-ng_ut:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    ! Local variables
    integer :: is,ie,js,je,bc

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)
    
    ! impose lo j side bc's
    select case(phys_bc(2,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(is:ie+1,js-1) = ZERO
    case (OUTLET)
       uedge(is:ie+1,js-1) = uedge(is:ie+1,js)
       vedge(is:ie  ,js-1) = vedge(is:ie  ,js)
    case (SYMMETRY)
       uedge(is:ie+1,js-1) = uedge(is:ie+1,js  )
       vedge(is:ie  ,js-1) = vedge(is:ie  ,js+1)
    case (INTERIOR, PERIODIC)
    case  default 
       call bl_error("impose_phys_bcs_2d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi j side bc's
    select case(phys_bc(2,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(is:ie+1,je+1) = ZERO
    case (OUTLET)
       uedge(is:ie+1,je+1) = uedge(is:ie+1,je)
       vedge(is:ie  ,je+1) = vedge(is:ie  ,je)
    case (SYMMETRY)
       uedge(is:ie+1,je+1) = uedge(is:ie+1,je  )
       vedge(is:ie  ,je+1) = vedge(is:ie  ,je-1)
    case (INTERIOR, PERIODIC)
    case  default 
       call bl_error("impose_phys_bcs_2d: invalid boundary type phys_bc(2,2)")
    end select

    ! impose lo i side bc's
    select case(phys_bc(1,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       vedge(is-1,js:je+1) = ZERO
    case (OUTLET)
       uedge(is-1,js:je  ) = uedge(is,js:je  )
       vedge(is-1,js:je+1) = vedge(is,js:je+1)
    case (SYMMETRY)
       uedge(is-1,js:je  ) = uedge(is+1,js:je  )
       vedge(is-1,js:je+1) = vedge(is  ,js:je+1)
    case (INTERIOR, PERIODIC)
    case  default 
       call bl_error("impose_phys_bcs_2d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi i side bc's
    select case(phys_bc(1,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       vedge(ie+1,js:je+1) = ZERO
    case (OUTLET)
       uedge(ie+1,js:je  ) = uedge(ie,js:je  )
       vedge(ie+1,js:je+1) = vedge(ie,js:je+1)
    case (SYMMETRY)
       uedge(ie+1,js:je  ) = uedge(ie-1,js:je  )
       vedge(ie+1,js:je+1) = vedge(ie  ,js:je+1)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("impose_phys_bcs_2d: invalid boundary type phys_bc(1,2)")
    end select

  end subroutine impose_phys_bcs_2d

  subroutine impose_phys_bcs_3d(uedge,vedge,wedge,ng_ut,lo,hi,phys_bc)

    use bc_module

    integer,         intent(in   ) :: lo(:),hi(:),ng_ut
    real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vedge(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) :: wedge(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    ! Local variables
    integer :: is,ie,js,je,ks,ke

    is = lo(1)
    js = lo(2)
    ks = lo(3)
    ie = hi(1)
    je = hi(2)
    ke = hi(3)

    ! impose lo k side bc's
    select case(phys_bc(3,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,:,ks-1) = ZERO
       vedge(:,:,ks-1) = ZERO
    case (OUTLET)
       uedge(:,:,ks-1) = uedge(:,:,ks)
       vedge(:,:,ks-1) = vedge(:,:,ks)
       wedge(:,:,ks-1) = wedge(:,:,ks)
    case (SYMMETRY)
       uedge(:,:,ks-1) = uedge(:,:,ks)
       vedge(:,:,ks-1) = vedge(:,:,ks)
       wedge(:,:,ks-1) = wedge(:,:,ks+1)
    case (INTERIOR, PERIODIC)
    case  default 
       call bl_error("impose_phys_bcs_3d: invalid boundary type phys_bc(3,1)")
    end select

    ! impose hi k side bc's
    select case(phys_bc(3,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,:,ke+1) = ZERO
       vedge(:,:,ke+1) = ZERO
    case (OUTLET)
       uedge(:,:,ke+1) = uedge(:,:,ke)
       vedge(:,:,ke+1) = vedge(:,:,ke)
       wedge(:,:,ke+1) = wedge(:,:,ke)
    case (SYMMETRY)
       uedge(:,:,ke+1) = uedge(:,:,ke)
       vedge(:,:,ke+1) = vedge(:,:,ke)
       wedge(:,:,ke+1) = wedge(:,:,ke-1)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("impose_phys_bcs_3d: invalid boundary type phys_bc(3,2)")
    end select

    ! impose lo j side bc's
    select case(phys_bc(2,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,js-1,:) = ZERO
       wedge(:,js-1,:) = ZERO
    case (OUTLET)
       uedge(:,js-1,:) = uedge(:,js,:)
       vedge(:,js-1,:) = vedge(:,js,:)
       wedge(:,js-1,:) = wedge(:,js,:)
    case (SYMMETRY)
       uedge(:,js-1,:) = uedge(:,js  ,:)
       vedge(:,js-1,:) = vedge(:,js+1,:)
       wedge(:,js-1,:) = wedge(:,js  ,:)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("impose_phys_bcs_3d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi j side bc's
    select case(phys_bc(2,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       uedge(:,je+1,:) = ZERO
       wedge(:,je+1,:) = ZERO
    case (OUTLET)
       uedge(:,je+1,:) = uedge(:,je,:)
       vedge(:,je+1,:) = vedge(:,je,:)
       wedge(:,je+1,:) = wedge(:,je,:)
    case (SYMMETRY)
       uedge(:,je+1,:) = uedge(:,je  ,:)
       vedge(:,je+1,:) = vedge(:,je-1,:)
       wedge(:,je+1,:) = wedge(:,je  ,:)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("impose_phys_bcs_3d: invalid boundary type phys_bc(2,2)")
    end select

    ! impose lo i side bc's
    select case(phys_bc(1,1))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       vedge(is-1,:,:) = ZERO
       wedge(is-1,:,:) = ZERO
    case (OUTLET)
       uedge(is-1,:,:) = uedge(is,:,:)
       vedge(is-1,:,:) = vedge(is,:,:)
       wedge(is-1,:,:) = wedge(is,:,:)
    case (SYMMETRY)
       uedge(is-1,:,:) = uedge(is+1,:,:)
       vedge(is-1,:,:) = vedge(is  ,:,:)
       wedge(is-1,:,:) = wedge(is  ,:,:)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("impose_phys_bcs_3d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi i side bc's
    select case(phys_bc(1,2))
    case (INLET, SLIP_WALL, NO_SLIP_WALL)
       vedge(ie+1,:,:) = ZERO
       wedge(ie+1,:,:) = ZERO
    case (OUTLET)
       print *,'SETTING UEDGE ',ie+1,lo(3)-ng_ut
       uedge(ie+1,:,:) = uedge(ie,:,:)
       vedge(ie+1,:,:) = vedge(ie,:,:)
       wedge(ie+1,:,:) = wedge(ie,:,:)
    case (SYMMETRY)
       uedge(ie+1,:,:) = uedge(ie-1,:,:)
       vedge(ie+1,:,:) = vedge(ie  ,:,:)
       wedge(ie+1,:,:) = wedge(ie  ,:,:)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("impose_phys_bcs_3d: invalid boundary type phys_bc(1,2)")
    end select

  end subroutine impose_phys_bcs_3d
  
end module impose_phys_bcs_on_edges_module
