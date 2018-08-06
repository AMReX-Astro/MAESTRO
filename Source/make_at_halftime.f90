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

  subroutine make_S_at_halftime(mla,shalf,sold,snew,dt,derivative,the_bc_level)

    use bl_prof_module
    use ml_restrict_fill_module
    use variables, only: foextrap_comp

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: shalf(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(in   ) :: snew(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    logical        , intent(in   ) :: derivative
    real(kind=dp_t), intent(in   ) :: dt

    real(kind=dp_t), pointer:: shp(:,:,:,:)
    real(kind=dp_t), pointer:: sop(:,:,:,:)
    real(kind=dp_t), pointer:: snp(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng_h,ng_o,ng_n
    integer :: i,in_comp,out_comp,n,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_S_at_halftime")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_h = nghost(shalf(1))
    ng_o = nghost(sold(1))
    ng_n = nghost(snew(1))

    in_comp = 1
    out_comp = 1

    do n = 1, nlevs

       do i = 1, nfabs(shalf(n))
          shp => dataptr(shalf(n), i)
          sop => dataptr(sold(n), i)
          snp => dataptr(snew(n), i)
          lo =  lwb(get_box(shalf(n), i))
          hi =  upb(get_box(shalf(n), i))
          select case (dm)
          case (1)
             call make_at_halftime_1d(shp(:,1,1,out_comp),sop(:,1,1,in_comp), &
                                      snp(:,1,1,in_comp),lo,hi,ng_h,ng_o,ng_n, &
                                      dt,derivative)
          case (2)
             call make_at_halftime_2d(shp(:,:,1,out_comp),sop(:,:,1,in_comp), &
                                      snp(:,:,1,in_comp),lo,hi,ng_h,ng_o,ng_n, &
                                      dt,derivative)
          case (3)
             call make_at_halftime_3d(shp(:,:,:,out_comp),sop(:,:,:,in_comp), &
                                      snp(:,:,:,in_comp),lo,hi,ng_h,ng_o,ng_n, &
                                      dt,derivative)
          end select
       end do

    end do

    call ml_restrict_and_fill(nlevs, shalf, mla%mba%rr, the_bc_level, &
         icomp=1, bcomp=foextrap_comp, nc=1, ng=ng_h)

    call destroy(bpt)

  end subroutine make_S_at_halftime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime(phihalf,sold,snew,in_comp,out_comp,dt,derivative,the_bc_level,mla)

    use ml_restrict_fill_module

    type(multifab) , intent(inout) :: phihalf(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(in   ) :: snew(:)
    integer        , intent(in   ) :: in_comp,out_comp
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    logical        , intent(in   ) :: derivative
    real(kind=dp_t), intent(in   ) :: dt
    
    real(kind=dp_t), pointer:: rhp(:,:,:,:)
    real(kind=dp_t), pointer:: rop(:,:,:,:)
    real(kind=dp_t), pointer:: rnp(:,:,:,:)
    integer   :: lo(mla%dim),hi(mla%dim)
    integer   :: ng_h,ng_o,ng_n,i,n,dm,nlevs

    dm = mla%dim
    nlevs = mla%nlevel

    ng_h = nghost(phihalf(1))
    ng_o = nghost(sold(1))
    ng_n = nghost(snew(1))

    do n = 1, nlevs
       do i = 1, nfabs(phihalf(n))
          rhp => dataptr(phihalf(n), i)
          rop => dataptr(sold(n), i)
          rnp => dataptr(snew(n), i)
          lo =  lwb(get_box(phihalf(n), i))
          hi =  upb(get_box(phihalf(n), i))
          select case (dm)
          case (1)
             call make_at_halftime_1d(rhp(:,1,1,out_comp),rop(:,1,1,in_comp), &
                                      rnp(:,1,1,in_comp),lo,hi,ng_h,ng_o,ng_n,&
                                      dt,derivative)
          case (2)
             call make_at_halftime_2d(rhp(:,:,1,out_comp),rop(:,:,1,in_comp), &
                                      rnp(:,:,1,in_comp),lo,hi,ng_h,ng_o,ng_n,&
                                      dt,derivative)
          case (3)
             call make_at_halftime_3d(rhp(:,:,:,out_comp),rop(:,:,:,in_comp), &
                                      rnp(:,:,:,in_comp),lo,hi,ng_h,ng_o,ng_n,&
                                      dt,derivative)
          end select
       end do
    end do

    call ml_restrict_and_fill(nlevs, phihalf, mla%mba%rr, the_bc_level, &
         icomp=1, bcomp=dm+in_comp, nc=1, ng=ng_h)

  end subroutine make_at_halftime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime_1d(phihalf,phiold,phinew,lo,hi,ng_half,ng_old,ng_new,&
                                 dt,derivative)

    use bl_constants_module
    use probin_module, only: derivative_mode
    integer         , intent(in   ) :: lo(:),hi(:),ng_half,ng_old,ng_new
    real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:)
    real (kind=dp_t), intent(in   ) ::  phiold(lo(1)-ng_old :)
    real (kind=dp_t), intent(in   ) ::  phinew(lo(1)-ng_new :)
    logical         , intent(in   ) :: derivative
    real (kind=dp_t), intent(in   ) :: dt
    
    !  Local variables
    integer :: i
    if (derivative) then
            do i = lo(1),hi(1)
               phihalf(i) = phiold(i) + HALF * dt * phinew(i)
            end do
    else
            do i = lo(1),hi(1)
               phihalf(i) = HALF * (phiold(i) + phinew(i))
            end do
    endif
  end subroutine make_at_halftime_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime_2d(phihalf,phiold,phinew,lo,hi,ng_half,ng_old,ng_new,&
                                 dt,derivative)

    use bl_constants_module
    
    integer         , intent(in   ) :: lo(:),hi(:),ng_half,ng_old,ng_new
    real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:,lo(2)-ng_half:)
    real (kind=dp_t), intent(in   ) ::  phiold(lo(1)-ng_old :,lo(2)-ng_old :)
    real (kind=dp_t), intent(in   ) ::  phinew(lo(1)-ng_new :,lo(2)-ng_new :)
    logical         , intent(in   ) :: derivative
    real (kind=dp_t), intent(in   ) :: dt

    !  Local variables
    integer :: i, j
    if (derivative) then
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  phihalf(i,j) = phiold(i,j) + HALF * dt * phinew(i,j)
               end do
            end do
    else
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  phihalf(i,j) = HALF * (phiold(i,j) + phinew(i,j))
               end do
            end do
    endif
  end subroutine make_at_halftime_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime_3d(phihalf,phiold,phinew,lo,hi,ng_half,ng_old,ng_new,&
                                 dt,derivative)

    use bl_constants_module
 
    integer         , intent(in   ) :: lo(:),hi(:),ng_half,ng_old,ng_new
    real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:,lo(2)-ng_half:,lo(3)-ng_half:)
    real (kind=dp_t), intent(in   ) ::  phiold(lo(1)-ng_old :,lo(2)-ng_old :,lo(3)-ng_old :)
    real (kind=dp_t), intent(in   ) ::  phinew(lo(1)-ng_new :,lo(2)-ng_new :,lo(3)-ng_new :)
    logical         , intent(in   ) :: derivative
    real (kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: i, j, k
    if (derivative) then
        !$OMP PARALLEL DO PRIVATE(i,j,k)
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     phihalf(i,j,k) = phiold(i,j,k) + HALF * dt * phinew(i,j,k)
                  end do
               end do
            end do
            !$OMP END PARALLEL DO

    else
            !$OMP PARALLEL DO PRIVATE(i,j,k)
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     phihalf(i,j,k) = HALF * (phiold(i,j,k) + phinew(i,j,k))
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
    endif
  end subroutine make_at_halftime_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module phihalf_module
