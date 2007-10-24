module convert_h

  use bl_types
  use bc_module
  use define_bc_module
  use multifab_module
  use boxarray_module
  use stencil_module
  use macproject_module
  use fill_3d_module
  use probin_module
  use network

  implicit none

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine convert_bigH_to_h(mla,s)

  type(ml_layout), intent(inout) :: mla
  type(multifab) , intent(inout) :: s(:)

! Local
  real(kind=dp_t), pointer    :: sp(:,:,:,:)
  integer                     :: lo(s(1)%dim),hi(s(1)%dim)
  integer                     :: ng, n, dm, nlevs, i

  ng = s(1)%ng
  dm = mla%dim
  nlevs = mla%nlevel
  
  do n=1,nlevs
     do i=1,s(n)%nboxes
        if (multifab_remote(s(n),i)) cycle
        sp         => dataptr(s(n),i)
        lo = lwb(get_box(s(n), i))
        hi = upb(get_box(s(n), i))
        select case (dm)
        case (2)
           call convert_bigH_to_h_2d(lo,hi,sp(:,:,1,:),ng)
        case (3)
        end select
     end do
  enddo

end subroutine convert_bigH_to_h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine convert_h_to_bigH(mla,s)

  type(ml_layout), intent(inout) :: mla
  type(multifab) , intent(inout) :: s(:)

! Local
  real(kind=dp_t), pointer    :: sp(:,:,:,:)
  integer                     :: lo(s(1)%dim),hi(s(1)%dim)
  integer                     :: ng, n, dm, nlevs, i

  ng = s(1)%ng
  dm = mla%dim
  nlevs = mla%nlevel
  
  do n=1,nlevs
     do i=1,s(n)%nboxes
        if (multifab_remote(s(n),i)) cycle
        sp         => dataptr(s(n),i)
        lo = lwb(get_box(s(n), i))
        hi = upb(get_box(s(n), i))
        select case (dm)
        case (2)
           call convert_h_to_bigH_2d(lo,hi,sp(:,:,1,:),ng)
        case (3)
           call convert_h_to_bigH_3d(lo,hi,sp(:,:,:,:),ng)
        end select
     end do
  enddo

end subroutine convert_h_to_bigH


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine convert_bigH_to_h_2d(lo,hi,s2,ng)

 integer        , intent(in   ) :: lo(:),hi(:),ng
 real(kind=dp_t), intent(inout) :: s2(lo(1)-ng:,lo(2)-ng:,:)

! Local
  integer :: i,j,n
  real(kind=dp_t) :: qreact

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        s2(i,j,rhoh_comp) = s2(i,j,rhoh_comp)/s2(i,j,rho_comp)

        qreact = 0.0d0
        do n=1,nspec
           qreact = qreact + s2(i,j,spec_comp+n-1)*ebin(n)/s2(i,j,rho_comp)
        enddo
        s2(i,j,rhoh_comp) = s2(i,j,rhoh_comp) - qreact

        s2(i,j,rhoh_comp) = s2(i,j,rhoh_comp)*s2(i,j,rho_comp)

     enddo
  enddo

end subroutine convert_bigH_to_h_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine convert_bigH_to_h_3d(lo,hi,s2,ng)

 integer        , intent(in   ) :: lo(:),hi(:),ng
 real(kind=dp_t), intent(inout) :: s2(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

! Local
  integer :: i,j,k,n
  real(kind=dp_t) :: qreact

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           
           s2(i,j,k,rhoh_comp) = s2(i,j,k,rhoh_comp)/s2(i,j,k,rho_comp)
           
           qreact = 0.0d0
           do n=1,nspec
              qreact = qreact + s2(i,j,k,spec_comp+n-1)*ebin(n)/s2(i,j,k,rho_comp)
           enddo
           s2(i,j,k,rhoh_comp) = s2(i,j,k,rhoh_comp) - qreact
           
           s2(i,j,k,rhoh_comp) = s2(i,j,k,rhoh_comp)*s2(i,j,k,rho_comp)

        enddo
     enddo
  enddo

end subroutine convert_bigH_to_h_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine convert_h_to_bigH_2d(lo,hi,s2,ng)

 integer        , intent(in   ) :: lo(:),hi(:),ng
 real(kind=dp_t), intent(inout) :: s2(lo(1)-ng:,lo(2)-ng:,:)

! Local
  integer :: i,j,n
  real(kind=dp_t) :: qreact

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        s2(i,j,rhoh_comp) = s2(i,j,rhoh_comp)/s2(i,j,rho_comp)

        qreact = 0.0d0
        do n=1,nspec
           qreact = qreact + s2(i,j,spec_comp+n-1)*ebin(n)/s2(i,j,rho_comp)
        enddo
        s2(i,j,rhoh_comp) = s2(i,j,rhoh_comp) + qreact

        s2(i,j,rhoh_comp) = s2(i,j,rhoh_comp)*s2(i,j,rho_comp)

     enddo
  enddo

end subroutine convert_h_to_bigH_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine convert_h_to_bigH_3d(lo,hi,s2,ng)

 integer        , intent(in   ) :: lo(:),hi(:),ng
 real(kind=dp_t), intent(inout) :: s2(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

! Local
  integer :: i,j,k,n
  real(kind=dp_t) :: qreact

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           
           s2(i,j,k,rhoh_comp) = s2(i,j,k,rhoh_comp)/s2(i,j,k,rho_comp)
           
           qreact = 0.0d0
           do n=1,nspec
              qreact = qreact + s2(i,j,k,spec_comp+n-1)*ebin(n)/s2(i,j,k,rho_comp)
           enddo
           s2(i,j,k,rhoh_comp) = s2(i,j,k,rhoh_comp) + qreact
           
           s2(i,j,k,rhoh_comp) = s2(i,j,k,rhoh_comp)*s2(i,j,k,rho_comp)
           
        enddo
     enddo
  enddo

end subroutine convert_h_to_bigH_3d

end module convert_h
