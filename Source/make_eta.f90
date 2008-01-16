module make_eta_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: make_eta

contains

  subroutine make_eta(nlevs,eta,sold,umac,sflux,dx,mla)

    use bl_constants_module
    use geometry, only: spherical

    integer           , intent(in   ) :: nlevs
    real(kind=dp_t)   , intent(inout) :: eta(:,0:,:)
    type(multifab)    , intent(in   ) :: sold(:)
    type(multifab)    , intent(in   ) :: umac(:,:)
    type(multifab)    , intent(in   ) :: sflux(:,:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)
    type(ml_layout)   , intent(inout) :: mla

    ! local
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)

    type(box) :: domain

    integer :: domlo(sold(1)%dim),domhi(sold(1)%dim)
    integer :: lo(sold(1)%dim),hi(sold(1)%dim)
    integer :: i,n,dm

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_eta")

    dm = sold(1)%dim

    do n=1,nlevs

       domain = layout_get_pd(sold(n)%la)
       domlo = lwb(domain)
       domhi = upb(domain)

       do i = 1, sold(n)%nboxes

          select case (dm)
          case (2)

          case (3)

             if (spherical .eq. 0) then

             else

             end if

          end select
       end do

    end do

    call destroy(bpt)

  end subroutine make_eta

end module make_eta_module
