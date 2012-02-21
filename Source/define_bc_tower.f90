module define_bc_module

  use box_module, only: box
  use layout_module, only: layout, layout_get_pd, layout_get_box, layout_dim, layout_nboxes
  use bc_module, only: EXT_DIR, FOEXTRAP, HOEXTRAP, REFLECT_EVEN, REFLECT_ODD, BC_DIR, &
                       BC_NEU, BC_PER, BC_INT, INTERIOR, &
                       INLET, NO_SLIP_WALL, OUTLET, PERIODIC, SLIP_WALL, SYMMETRY

  implicit none

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
     integer :: max_level_built = 0
     type(bc_level), pointer :: bc_tower_array(:) => Null()
     integer       , pointer :: domain_bc(:,:) => Null()

  end type bc_tower

  private

  public :: bc_level, bc_tower, bc_tower_init, bc_tower_level_build, bc_tower_destroy

contains

  subroutine bc_tower_init(bct,num_levs,dm,phys_bc_in)

    type(bc_tower ), intent(  out) :: bct
    integer        , intent(in   ) :: num_levs
    integer        , intent(in   ) :: dm
    integer        , intent(in   ) :: phys_bc_in(:,:)

    integer :: n

    bct%nlevels = num_levs
    bct%dim     = dm
    allocate(bct%bc_tower_array(bct%nlevels))
    allocate(bct%domain_bc(dm,2))

    do n = 1, num_levs
      bct%bc_tower_array(n)%ngrids = -1
    end do

    bct%domain_bc(:,:) = phys_bc_in(:,:)

  end subroutine bc_tower_init

  subroutine bc_tower_level_build(bct,n,la)

    use variables, only : nscal

    type(bc_tower ), intent(inout) :: bct
    integer        , intent(in   ) :: n
    type(layout)   , intent(in   ) :: la

    integer :: ngrids
    integer :: default_value

    if (bct%bc_tower_array(n)%ngrids > 0) then
      deallocate(bct%bc_tower_array(n)%phys_bc_level_array)
      deallocate(bct%bc_tower_array(n)%adv_bc_level_array)
      deallocate(bct%bc_tower_array(n)%ell_bc_level_array)
    end if

    ngrids = layout_nboxes(la)
    bct%bc_tower_array(n)%dim    = bct%dim
    bct%bc_tower_array(n)%ngrids = ngrids
    bct%bc_tower_array(n)%domain = layout_get_pd(la)

    allocate(bct%bc_tower_array(n)%phys_bc_level_array(0:ngrids,bct%dim,2))
    default_value = INTERIOR
    call phys_bc_level_build(bct%bc_tower_array(n)%phys_bc_level_array,la, &
                             bct%domain_bc,default_value)

    allocate(bct%bc_tower_array(n)%adv_bc_level_array(0:ngrids,bct%dim,2,bct%dim+nscal+3))
    default_value = INTERIOR
    call adv_bc_level_build(bct%bc_tower_array(n)%adv_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array,default_value)

    allocate(bct%bc_tower_array(n)%ell_bc_level_array(0:ngrids,bct%dim,2,bct%dim+nscal+3))
    default_value = BC_INT
    call ell_bc_level_build(bct%bc_tower_array(n)%ell_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array,default_value)

     bct%max_level_built = n

  end subroutine bc_tower_level_build

  subroutine bc_tower_destroy(bct)

    type(bc_tower), intent(inout) :: bct

    integer :: i

    do i = 1,bct%nlevels
       deallocate(bct%bc_tower_array(i)%phys_bc_level_array)
       deallocate(bct%bc_tower_array(i)%adv_bc_level_array)
       deallocate(bct%bc_tower_array(i)%ell_bc_level_array)
    end do
    deallocate(bct%bc_tower_array)
    deallocate(bct%domain_bc)

  end subroutine bc_tower_destroy

  subroutine phys_bc_level_build(phys_bc_level,la_level,domain_bc,default_value)

    integer     , intent(inout) :: phys_bc_level(0:,:,:)
    integer     , intent(in   ) :: domain_bc(:,:)
    type(layout), intent(in   ) :: la_level
    integer     , intent(in   ) :: default_value
    type(box) :: bx,pd
    integer :: i,d

    pd = layout_get_pd(la_level) 

    phys_bc_level = default_value

    i = 0
    do d = 1,layout_dim(la_level)
       phys_bc_level(i,d,1) = domain_bc(d,1)
       phys_bc_level(i,d,2) = domain_bc(d,2)
    end do

    do i = 1,layout_nboxes(la_level)
       bx = layout_get_box(la_level,i)
       do d = 1,layout_dim(la_level)
          if (bx%lo(d) == pd%lo(d)) phys_bc_level(i,d,1) = domain_bc(d,1)
          if (bx%hi(d) == pd%hi(d)) phys_bc_level(i,d,2) = domain_bc(d,2)
       end do
    end do

  end subroutine phys_bc_level_build

  subroutine adv_bc_level_build(adv_bc_level,phys_bc_level,default_value)

    use variables, only: rho_comp, rhoh_comp, spec_comp, temp_comp, trac_comp, press_comp, &
         foextrap_comp, hoextrap_comp, ntrac
    use network, only: nspec

    ! define boundary conditions for the advection problem

    integer  , intent(inout) ::  adv_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: n,d,i,dm

    adv_bc_level = default_value

    dm = size(adv_bc_level,dim=2)

!    *** 2-D ***
!   COMP = 1     : x-velocity
!   COMP = 2     : y-velocity
!   COMP = 3...  : density, (rho h), (rho X)_i, temp, tracers, pressure

!    *** 3-D ***
!   COMP = 1     : x-velocity
!   COMP = 2     : y-velocity
!   COMP = 3     : z-velocity
!   COMP = 4...  : density, (rho h), (rho X)_i, temp, tracers, pressure

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
          adv_bc_level(n,d,i,hoextrap_comp)                    = HOEXTRAP      ! higher order extrap
       else if (phys_bc_level(n,d,i) == NO_SLIP_WALL) then
          adv_bc_level(n,d,i,     1:dm)                        = EXT_DIR       ! velocity
          adv_bc_level(n,d,i, rho_comp+dm)                     = HOEXTRAP      ! density
          adv_bc_level(n,d,i,rhoh_comp+dm)                     = HOEXTRAP      ! (rho h)
          adv_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= HOEXTRAP      ! (rho X)_i
          adv_bc_level(n,d,i,temp_comp+dm)                     = HOEXTRAP      ! temperature
          adv_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= HOEXTRAP      ! tracers
          adv_bc_level(n,d,i,press_comp)                       = FOEXTRAP      ! pressure
          adv_bc_level(n,d,i,foextrap_comp)                    = FOEXTRAP      ! first order extrap
          adv_bc_level(n,d,i,hoextrap_comp)                    = HOEXTRAP      ! higher order extrap
       else if (phys_bc_level(n,d,i) == INLET) then
          adv_bc_level(n,d,i,     1:dm)                        = EXT_DIR       ! velocity
          adv_bc_level(n,d,i, rho_comp+dm)                     = EXT_DIR       ! density
          adv_bc_level(n,d,i,rhoh_comp+dm)                     = EXT_DIR       ! (rho h)
          adv_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= EXT_DIR       ! (rho X)_i
          adv_bc_level(n,d,i,temp_comp+dm)                     = EXT_DIR       ! temperature
          adv_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= EXT_DIR       ! tracers
          adv_bc_level(n,d,i,press_comp)                       = FOEXTRAP      ! pressure
          adv_bc_level(n,d,i,foextrap_comp)                    = FOEXTRAP      ! first order extrap
          adv_bc_level(n,d,i,hoextrap_comp)                    = HOEXTRAP      ! higher order extrap
       else if (phys_bc_level(n,d,i) == OUTLET) then
          adv_bc_level(n,d,i,     1:dm)                        = FOEXTRAP      ! velocity
          adv_bc_level(n,d,i, rho_comp+dm)                     = FOEXTRAP      ! density
          adv_bc_level(n,d,i,rhoh_comp+dm)                     = FOEXTRAP      ! (rho h)
          adv_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= FOEXTRAP      ! (rho X)_i
          adv_bc_level(n,d,i,temp_comp+dm)                     = FOEXTRAP      ! temperature
          adv_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= FOEXTRAP      ! tracers
          adv_bc_level(n,d,i,press_comp)                       = EXT_DIR       ! pressure
          adv_bc_level(n,d,i,foextrap_comp)                    = FOEXTRAP      ! first order extrap
          adv_bc_level(n,d,i,hoextrap_comp)                    = HOEXTRAP      ! higher order extrap
       else if (phys_bc_level(n,d,i) == SYMMETRY) then
          adv_bc_level(n,d,i,     1:dm)                        = REFLECT_EVEN  ! tangential vel
          adv_bc_level(n,d,i,        d)                        = REFLECT_ODD   ! normal vel
          adv_bc_level(n,d,i, rho_comp+dm)                     = REFLECT_EVEN  ! density
          adv_bc_level(n,d,i,rhoh_comp+dm)                     = REFLECT_EVEN  ! (rho h)
          adv_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= REFLECT_EVEN  ! (rho X)_i
          adv_bc_level(n,d,i,temp_comp+dm)                     = REFLECT_EVEN  ! temperature
          adv_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= REFLECT_EVEN  ! tracers
          adv_bc_level(n,d,i,press_comp)                       = REFLECT_EVEN  ! pressure
          adv_bc_level(n,d,i,foextrap_comp)                    = REFLECT_EVEN  ! first order extrap -- overridden by symmetry
          adv_bc_level(n,d,i,hoextrap_comp)                    = REFLECT_EVEN  ! higher order extrap -- overridden by symmetry
       end if
    end do
    end do
    end do

  end subroutine adv_bc_level_build

  subroutine ell_bc_level_build(ell_bc_level,phys_bc_level,default_value)

    use variables, only: rho_comp, rhoh_comp, spec_comp, temp_comp, trac_comp, press_comp, &
         foextrap_comp, hoextrap_comp, ntrac
    use network, only: nspec

    ! define boundary conditions for the elliptic problem

    integer  , intent(inout) ::  ell_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: n,d,i,dm

    dm = size(ell_bc_level,dim=2)

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
          ell_bc_level(n,d,i,hoextrap_comp)                    = BC_NEU   ! higher order extrap
       else if (phys_bc_level(n,d,i) == NO_SLIP_WALL) then
          ell_bc_level(n,d,i,                      1:dm)       = BC_DIR   ! vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_NEU   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_NEU   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_NEU   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_NEU   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_NEU   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_NEU   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_NEU   ! first order extrap
          ell_bc_level(n,d,i,hoextrap_comp)                    = BC_NEU   ! higher order extrap
       else if (phys_bc_level(n,d,i) == INLET) then
          ell_bc_level(n,d,i,                      1:dm)       = BC_DIR   ! vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_DIR   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_DIR   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_DIR   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_DIR   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_DIR   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_NEU   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_NEU   ! first order extrap
          ell_bc_level(n,d,i,hoextrap_comp)                    = BC_NEU   ! higher order extrap
       else if (phys_bc_level(n,d,i) == OUTLET) then
          ell_bc_level(n,d,i,                      1:dm)       = BC_NEU   ! tangential vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_NEU   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_NEU   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_NEU   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_NEU   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_NEU   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_DIR   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_NEU   ! first order extrap
          ell_bc_level(n,d,i,hoextrap_comp)                    = BC_NEU   ! higher order extrap
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
          ell_bc_level(n,d,i,hoextrap_comp)                    = BC_NEU   ! higher order extrap
       else if (phys_bc_level(n,d,i) == PERIODIC) then
          ell_bc_level(n,d,i,                      1:dm      ) = BC_PER   ! vel.
          ell_bc_level(n,d,i,rho_comp+dm)                      = BC_PER   ! density
          ell_bc_level(n,d,i,rhoh_comp+dm)                     = BC_PER   ! (rho h)
          ell_bc_level(n,d,i,spec_comp+dm:spec_comp+dm+nspec-1)= BC_PER   ! (rho X)_i
          ell_bc_level(n,d,i,temp_comp+dm)                     = BC_PER   ! temperature
          ell_bc_level(n,d,i,trac_comp+dm:trac_comp+dm+ntrac-1)= BC_PER   ! tracers
          ell_bc_level(n,d,i,press_comp)                       = BC_PER   ! pressure
          ell_bc_level(n,d,i,foextrap_comp)                    = BC_PER   ! first order extrap
          ell_bc_level(n,d,i,hoextrap_comp)                    = BC_PER   ! higher order extrap
       end if
    end do
    end do
    end do

  end subroutine ell_bc_level_build

end module define_bc_module
