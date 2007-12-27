module phihalf_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none
  
  private

  public :: make_S_at_halftime, make_at_halftime

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_S_at_halftime (nlevs,shalf,sold,snew)

     use bl_prof_module

     integer        , intent(in   ) :: nlevs
     type(multifab) , intent(inout) :: shalf(:)
     type(multifab) , intent(in   ) :: sold(:)
     type(multifab) , intent(in   ) :: snew(:)

     real(kind=dp_t), pointer:: shp(:,:,:,:)
     real(kind=dp_t), pointer:: sop(:,:,:,:)
     real(kind=dp_t), pointer:: snp(:,:,:,:)
     integer :: lo(shalf(1)%dim),hi(shalf(1)%dim),ng_h,ng_o,dm
     integer :: i,in_comp,out_comp,n

     type(bl_prof_timer), save :: bpt

     call build(bpt, "make_S_at_halftime")
     
     dm = shalf(1)%dim
     ng_h = shalf(1)%ng
     ng_o = sold(1)%ng
     
     in_comp = 1
     out_comp = 1
     
     do n = 1, nlevs
        
        do i = 1, shalf(n)%nboxes
           if ( multifab_remote(shalf(n), i) ) cycle
           shp => dataptr(shalf(n), i)
           sop => dataptr(sold(n), i)
           snp => dataptr(snew(n), i)
           lo =  lwb(get_box(shalf(n), i))
           hi =  upb(get_box(shalf(n), i))
           select case (dm)
           case (2)
              call make_at_halftime_2d(shp(:,:,1,out_comp),sop(:,:,1,in_comp), &
                                       snp(:,:,1,in_comp),lo,hi,ng_h,ng_o)
           case (3)
              call make_at_halftime_3d(shp(:,:,:,out_comp),sop(:,:,:,in_comp), &
                                       snp(:,:,:,in_comp),lo,hi,ng_h,ng_o)
           end select
        end do
        
     enddo

     call destroy(bpt)

   end subroutine make_S_at_halftime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_at_halftime(nlevs,phihalf,sold,snew,in_comp,out_comp,dx,the_bc_level,mla)

     use multifab_physbc_module
     use ml_restriction_module, only : ml_cc_restriction
     use multifab_fill_ghost_module
     
     integer        , intent(in   ) :: nlevs
     type(multifab) , intent(inout) :: phihalf(:)
     type(multifab) , intent(in   ) :: sold(:)
     type(multifab) , intent(in   ) :: snew(:)
     integer        , intent(in   ) :: in_comp,out_comp
     real(kind=dp_t), intent(in   ) :: dx(:,:)
     type(bc_level) , intent(in   ) :: the_bc_level(:)
     type(ml_layout), intent(inout) :: mla
     
     real(kind=dp_t), pointer:: rhp(:,:,:,:)
     real(kind=dp_t), pointer:: rop(:,:,:,:)
     real(kind=dp_t), pointer:: rnp(:,:,:,:)
     integer   :: lo(phihalf(1)%dim),hi(phihalf(1)%dim)
     integer   :: ng_h,ng_o,dm,i,n

     dm = phihalf(1)%dim
     ng_h = phihalf(1)%ng
     ng_o = sold(1)%ng

     do n = 1, nlevs

        do i = 1, phihalf(n)%nboxes
           if ( multifab_remote(phihalf(n), i) ) cycle
           rhp => dataptr(phihalf(n), i)
           rop => dataptr(sold(n), i)
           rnp => dataptr(snew(n), i)
           lo =  lwb(get_box(phihalf(n), i))
           hi =  upb(get_box(phihalf(n), i))
           select case (dm)
           case (2)
              call make_at_halftime_2d(rhp(:,:,1,out_comp),rop(:,:,1,in_comp), &
                                       rnp(:,:,1,in_comp),lo,hi,ng_h,ng_o)
           case (3)
              call make_at_halftime_3d(rhp(:,:,:,out_comp),rop(:,:,:,in_comp), &
                                       rnp(:,:,:,in_comp),lo,hi,ng_h,ng_o)
           end select
        end do
        
        if(ng_h .gt. 0) then
           call multifab_fill_boundary(phihalf(n))
           call multifab_physbc(phihalf(n),out_comp,dm+in_comp,1,dx(n,:),the_bc_level(n))
        endif
        
     enddo ! end loop over nlevs

     if(ng_h .gt. 0) then
        do n = nlevs, 2, -1
           call ml_cc_restriction(phihalf(n-1),phihalf(n),mla%mba%rr(n-1,:))

           call multifab_fill_ghost_cells(phihalf(n),phihalf(n-1), &
                                          ng_h,mla%mba%rr(n-1,:), &
                                          the_bc_level(n-1), the_bc_level(n  ), &
                                          1,dm+in_comp,1)
        enddo
     endif

   end subroutine make_at_halftime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_at_halftime_2d(phihalf,phiold,phinew,lo,hi,ng_half,ng_old)

     use bl_constants_module

     implicit none
     
     integer         , intent(in   ) :: lo(:),hi(:),ng_half,ng_old
     real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:,lo(2)-ng_half:)
     real (kind=dp_t), intent(in   ) :: phiold(lo(1)-ng_old:,lo(2)-ng_old:)
     real (kind=dp_t), intent(in   ) :: phinew(lo(1)-ng_old:,lo(2)-ng_old:)
     
     !  Local variables
     integer :: i, j
     
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           phihalf(i,j) = HALF * (phiold(i,j) + phinew(i,j))
        end do
     end do
     
   end subroutine make_at_halftime_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   subroutine make_at_halftime_3d(phihalf,phiold,phinew,lo,hi,ng_half,ng_old)

     use bl_constants_module
     
     implicit none
     
     integer         , intent(in   ) :: lo(:),hi(:),ng_half,ng_old
     real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:,lo(2)-ng_half:,lo(3)-ng_half:)
     real (kind=dp_t), intent(in   ) :: phiold(lo(1)-ng_old:,lo(2)-ng_old:,lo(3)-ng_old:)
     real (kind=dp_t), intent(in   ) :: phinew(lo(1)-ng_old:,lo(2)-ng_old:,lo(3)-ng_old:)

     ! Local variables
     integer :: i, j, k
     
     do k = lo(3),hi(3)
        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
              phihalf(i,j,k) = HALF * (phiold(i,j,k) + phinew(i,j,k))
           end do
        end do
     end do
     
   end subroutine make_at_halftime_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
 end module phihalf_module
