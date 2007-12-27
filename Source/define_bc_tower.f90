module define_bc_module

  use bl_types
  use ml_layout_module
  use bc_module

  type bc_level

     integer   :: dim    = 0
     integer   :: ngrids = 0
     type(box) :: domain 
     integer, pointer :: phys_bc_level_array(:,:,:) => Null()
     integer, pointer ::  adv_bc_level_array(:,:,:,:) => Null()
     integer, pointer ::  ell_bc_level_array(:,:,:,:) => Null()

  end type bc_level

  type bc_tower

     integer :: dim     = 0
     integer :: nlevels = 0
     type(bc_level), pointer :: bc_tower_array(:) => Null()

  end type bc_tower

  private

  public :: bc_level, bc_tower, bc_tower_build, bc_tower_destroy

contains

  subroutine bc_tower_build(bct,mla,domain_bc,domain_box,nspec)

    use bl_prof_module
    use variables, only: nscal

    type(bc_tower ), intent(  out) :: bct
    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: domain_bc(:,:)
    type(box)      , intent(in   ) :: domain_box(:)
    integer        , intent(in   ) :: nspec

    integer :: i,ngrids
    integer :: default_value

    type(bl_prof_timer), save :: bpt

    call build(bpt, "bc_tower_build")

    bct%nlevels = mla%nlevel
    bct%dim     = mla%dim

    allocate(bct%bc_tower_array(bct%nlevels))
    do i = 1,bct%nlevels
       ngrids = layout_nboxes(mla%la(i))
       bct%bc_tower_array(i)%dim    = bct%dim
       bct%bc_tower_array(i)%ngrids = ngrids
       bct%bc_tower_array(i)%domain = domain_box(i)

       allocate(bct%bc_tower_array(i)%phys_bc_level_array(0:ngrids,bct%dim,2))
       default_value = INTERIOR
       call phys_bc_level_build(bct%bc_tower_array(i)%phys_bc_level_array,mla%la(i), &
                                domain_bc,default_value)

       allocate(bct%bc_tower_array(i)%adv_bc_level_array(0:ngrids,bct%dim,2,bct%dim+nscal+2))
       default_value = INTERIOR
       call adv_bc_level_build(bct%bc_tower_array(i)%adv_bc_level_array, &
                               bct%bc_tower_array(i)%phys_bc_level_array,default_value,nspec)

       allocate(bct%bc_tower_array(i)%ell_bc_level_array(0:ngrids,bct%dim,2,bct%dim+nscal+2))
       default_value = BC_INT
       call ell_bc_level_build(bct%bc_tower_array(i)%ell_bc_level_array, &
                               bct%bc_tower_array(i)%phys_bc_level_array,default_value,nspec)
    end do

    call destroy(bpt)

  end subroutine bc_tower_build

  subroutine bc_tower_destroy(bct)

    type(bc_tower), intent(inout) :: bct

    integer :: i

    do i = 1,bct%nlevels
       deallocate(bct%bc_tower_array(i)%phys_bc_level_array)
       deallocate(bct%bc_tower_array(i)%adv_bc_level_array)
       deallocate(bct%bc_tower_array(i)%ell_bc_level_array)
    end do
    deallocate(bct%bc_tower_array)

  end subroutine bc_tower_destroy

  subroutine phys_bc_level_build(phys_bc_level,la_level,domain_bc,default_value)

    integer     , intent(inout) :: phys_bc_level(0:,:,:)
    integer     , intent(in   ) :: domain_bc(:,:)
    type(layout), intent(in   ) :: la_level
    integer     , intent(in   ) :: default_value
    type(box) :: bx,pd
    integer :: i,dm

    pd = layout_get_pd(la_level) 

    phys_bc_level = default_value

    i = 0
    do dm = 1,layout_dim(la_level)
       phys_bc_level(i,dm,1) = domain_bc(dm,1)
       phys_bc_level(i,dm,2) = domain_bc(dm,2)
    end do

    do i = 1,layout_nboxes(la_level)
       bx = layout_get_box(la_level,i)
       do dm = 1,layout_dim(la_level)
          if (bx%lo(dm) == pd%lo(dm)) phys_bc_level(i,dm,1) = domain_bc(dm,1)
          if (bx%hi(dm) == pd%hi(dm)) phys_bc_level(i,dm,2) = domain_bc(dm,2)
       end do
    end do

  end subroutine phys_bc_level_build

  subroutine adv_bc_level_build(adv_bc_level,phys_bc_level,default_value,nspec)

    use variables, only: rho_comp, rhoh_comp, spec_comp, temp_comp, trac_comp, press_comp, foextrap_comp, ntrac

    ! define boundary conditions for the advection problem

    integer  , intent(inout) ::  adv_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value
    integer  , intent(in   ) :: nspec

    integer :: dm
    integer :: n,d,i

    adv_bc_level = default_value

!    *** 2-D ***
!   COMP = 1     : x-velocity
!   COMP = 2     : y-velocity
!   COMP = 3...  : density, (rho h), (rho X)_i, temp, tracers, pressure

!    *** 3-D ***
!   COMP = 1     : x-velocity
!   COMP = 2     : y-velocity
!   COMP = 3     : z-velocity
!   COMP = 4...  : density, (rho h), (rho X)_i, temp, tracers, pressure

    dm = size(adv_bc_level,dim=2)

    do n  = 0,size(adv_bc_level,dim=1)-1
    do d  = 1,dm
    do i  = 1,2
       if (phys_bc_level(n,d,i) == SLIP_WALL) then
          adv_bc_level(n,d,i,     1:dm)                        = HOEXTRAP      ! tangential vel
          adv_bc_level(n,d,i,        d)                        = EXT_DIR       ! normal vel
          adv_bc_level(n,d,i, rho_comp+dm)                     = HOEXTRAP      ! density
          adv_bc_level(n,d,i,rhoh_comp+dm)                     = HOEXTRAP      ! (rho h)
          adv_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= HOEXTRAP      ! (rho X)_i
          adv_bc_level(n,d,i,temp_comp+dm)                     = HOEXTRAP      ! temperature
          adv_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= HOEXTRAP      ! tracers
          adv_bc_level(n,d,i,press_comp)                       = FOEXTRAP      ! pressure
          adv_bc_level(n,d,i,foextrap_comp)                    = FOEXTRAP      ! first order extrap

       else if (phys_bc_level(n,d,i) == NO_SLIP_WALL) then
          adv_bc_level(n,d,i,     1:dm)                        = EXT_DIR       ! velocity
          adv_bc_level(n,d,i, rho_comp+dm)                     = HOEXTRAP      ! density
          adv_bc_level(n,d,i,rhoh_comp+dm)                     = HOEXTRAP      ! (rho h)
          adv_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= HOEXTRAP      ! (rho X)_i
          adv_bc_level(n,d,i,temp_comp+dm)                     = HOEXTRAP      ! temperature
          adv_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= HOEXTRAP      ! tracers
          adv_bc_level(n,d,i,press_comp)                       = FOEXTRAP      ! pressure
          adv_bc_level(n,d,i,foextrap_comp)                    = FOEXTRAP      ! first order extrap
   
       else if (phys_bc_level(n,d,i) == INLET) then
          adv_bc_level(n,d,i,     1:dm)                        = EXT_DIR       ! velocity
          adv_bc_level(n,d,i, rho_comp+dm)                     = EXT_DIR       ! density
          adv_bc_level(n,d,i,rhoh_comp+dm)                     = EXT_DIR       ! (rho h)
          adv_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= EXT_DIR       ! (rho X)_i
          adv_bc_level(n,d,i,temp_comp+dm)                     = EXT_DIR       ! temperature
          adv_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= EXT_DIR       ! tracers
          adv_bc_level(n,d,i,press_comp)                       = FOEXTRAP      ! pressure
          adv_bc_level(n,d,i,foextrap_comp)                    = FOEXTRAP      ! first order extrap

       else if (phys_bc_level(n,d,i) == OUTLET) then
          adv_bc_level(n,d,i,     1:dm)                        = FOEXTRAP      ! velocity
          adv_bc_level(n,d,i, rho_comp+dm)                     = FOEXTRAP      ! density
          adv_bc_level(n,d,i,rhoh_comp+dm)                     = FOEXTRAP      ! (rho h)
          adv_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= FOEXTRAP      ! (rho X)_i
          adv_bc_level(n,d,i,temp_comp+dm)                     = FOEXTRAP      ! temperature
          adv_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= FOEXTRAP      ! tracers
          adv_bc_level(n,d,i,press_comp)                       = EXT_DIR       ! pressure
          adv_bc_level(n,d,i,foextrap_comp)                    = FOEXTRAP      ! first order extrap

       else if (phys_bc_level(n,d,i) == SYMMETRY) then
          adv_bc_level(n,d,i,     1:dm)                        = REFLECT_EVEN  ! tangential vel
          adv_bc_level(n,d,i,        d)                        = REFLECT_ODD   ! normal vel
          adv_bc_level(n,d,i, rho_comp+dm)                     = REFLECT_EVEN  ! density
          adv_bc_level(n,d,i,rhoh_comp+dm)                     = REFLECT_EVEN  ! (rho h)
          adv_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= REFLECT_EVEN  ! (rho X)_i
          adv_bc_level(n,d,i,temp_comp+dm)                     = REFLECT_EVEN  ! temperature
          adv_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= REFLECT_EVEN  ! tracers
          adv_bc_level(n,d,i,press_comp)                       = REFLECT_EVEN  ! pressure
          adv_bc_level(n,d,i,foextrap_comp)                    = FOEXTRAP      ! first order extrap

       end if
    end do
    end do
    end do

  end subroutine adv_bc_level_build

  subroutine ell_bc_level_build(ell_bc_level,phys_bc_level,default_value,nspec)

    use variables, only: rho_comp, rhoh_comp, spec_comp, temp_comp, trac_comp, press_comp, foextrap_comp, ntrac

    ! define boundary conditions for the elliptic problem

    implicit none

    integer  , intent(inout) ::  ell_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value
    integer  , intent(in   ) :: nspec

    integer :: dm
    integer :: n,d,i

    ell_bc_level = default_value

!    *** 2-D ***
!   COMP = 1     : x-velocity
!   COMP = 2     : y-velocity
!   COMP = 3...  : density, (rho h), (rho X)_i, temp, tracers, pressure

!    *** 3-D ***
!   COMP = 1     : x-velocity
!   COMP = 2     : y-velocity
!   COMP = 3     : z-velocity
!   COMP = 4...  : density, (rho h), (rho X)_i, temp, tracers, pressure

    dm = size(ell_bc_level,dim=2)

    do n = 0,size(ell_bc_level,dim=1)-1
    do d = 1,dm
    do i = 1,2
       if (phys_bc_level(n,d,i) == SLIP_WALL) then
          ell_bc_level(n,d,i,                      1:dm)       = BC_NEU   ! tangential vel.
          ell_bc_level(n,d,i,                         d)       = BC_DIR   ! normal vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_NEU   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_NEU   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_NEU   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_NEU   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_NEU   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_NEU   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_NEU   ! first order extrap
       else if (phys_bc_level(n,d,i) == NO_SLIP_WALL) then
          ell_bc_level(n,d,i,                      1:dm)       = BC_DIR   ! vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_NEU   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_NEU   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_NEU   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_NEU   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_NEU   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_NEU   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_NEU   ! first order extrap
       else if (phys_bc_level(n,d,i) == INLET) then
          ell_bc_level(n,d,i,                      1:dm)       = BC_DIR   ! vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_DIR   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_DIR   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_DIR   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_DIR   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_DIR   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_NEU   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_NEU   ! first order extrap
       else if (phys_bc_level(n,d,i) == OUTLET) then
          ell_bc_level(n,d,i,                      1:dm)       = BC_NEU   ! tangential vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_NEU   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_NEU   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_NEU   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_NEU   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_NEU   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_DIR   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_NEU   ! first order extrap
       else if (phys_bc_level(n,d,i) == SYMMETRY) then
          ell_bc_level(n,d,i,                      1:dm)       = BC_NEU   ! tangential vel.
          ell_bc_level(n,d,i,                         d)       = BC_DIR   ! normal vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_NEU   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_NEU   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_NEU   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_NEU   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_NEU   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_NEU   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_NEU   ! first order extrap
       else if (phys_bc_level(n,d,i) == PERIODIC) then
          ell_bc_level(n,d,i,                      1:dm      ) = BC_PER   ! vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_PER   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_PER   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_PER   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_PER   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_PER   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_PER   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_PER   ! first order extrap
       end if
    end do
    end do
    end do

  end subroutine ell_bc_level_build

end module define_bc_module
