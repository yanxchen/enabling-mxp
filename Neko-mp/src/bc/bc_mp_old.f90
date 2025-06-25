! Copyright (c) 2020-2024, The Neko Authors
! All rights reserved.
!
!> Defines mixed-precision boundary conditions
module bc_mp
  use neko_config
  use num_types
  use device
  use dofmap, only : dofmap_t
  use dofmap_mp, only : dofmap_sp_t, dofmap_dp_t
  use coef_mp, only : coef_dp_t, coef_sp_t
!   use space, only : space_t
  use space_mp, only : space_dp_t
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use facet_zone, only : facet_zone_t
  use stack, only : stack_i4t2_t
  use tuple, only : tuple_i4_t
  use utils, only : neko_error, linear_index, split_string
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private

  !> Base type for a mixed-precision boundary condition
  type, public, abstract :: bc_mp_t
     !> The linear index of each node in each boundary facet
     integer, allocatable :: msk(:)
     !> A list of facet ids (1 to 6), one for each element in msk
     integer, allocatable :: facet(:)
     !> Map of degrees of freedom
     type(dofmap_t), pointer :: dof
     type(dofmap_sp_t), pointer :: dof_sp => null()
     type(dofmap_dp_t), pointer :: dof_dp => null()
     !> Mixed-precision SEM coefficients
     type(coef_dp_t), pointer :: coef
     !> The mesh
     type(mesh_t), pointer :: msh
     !> The function space
     type(space_dp_t), pointer :: Xh
     !> Index tuples (facet, element) marked as part of the boundary condition
     type(stack_i4t2_t) :: marked_facet
     !> Device pointer for msk
     type(c_ptr) :: msk_d = C_NULL_PTR
     !> Device pointer for facet
     type(c_ptr) :: facet_d = C_NULL_PTR
   contains
     !> Constructor
     procedure, pass(this) :: init_base => bc_mp_init_base
     !> Destructor
     procedure, pass(this) :: free_base => bc_mp_free_base
     !> Mark a facet on an element as part of the boundary condition
     procedure, pass(this) :: mark_facet => bc_mp_mark_facet
     !> Finalize the construction of the bc by populating the msk and facet arrays
     procedure, pass(this) :: finalize => bc_mp_finalize
     !> Apply the boundary condition to a scalar field (SP)
     procedure(bc_mp_apply_scalar_sp), pass(this), deferred :: apply_scalar_sp
     !> Apply the boundary condition to a vector field (SP)
     procedure(bc_mp_apply_vector_sp), pass(this), deferred :: apply_vector_sp
     !> Apply the boundary condition to a scalar field (DP)
     procedure(bc_mp_apply_scalar_dp), pass(this), deferred :: apply_scalar_dp
     !> Apply the boundary condition to a vector field (DP)
     procedure(bc_mp_apply_vector_dp), pass(this), deferred :: apply_vector_dp
     !> Destructor
     procedure(bc_mp_destructor), pass(this), deferred :: free
  end type bc_mp_t

  !> Pointer to boundary condition [mp]
  type, private :: bcp_mp_t
     class(bc_mp_t), pointer :: bcp
  end type bcp_mp_t

  !> List of boundary conditions [mp]
  type, public :: bc_list_mp_t
     private
     !> Number of boundary conditions in the list
     integer :: n
     !> Maximum number of boundary conditions in the list
     integer :: size
     !> List of boundary conditions
     type(bcp_mp_t), allocatable :: bc(:)
  end type bc_list_mp_t

  abstract interface
     !> Apply the boundary condition to a scalar field (SP)
     subroutine bc_mp_apply_scalar_sp(this, x, n, t, tstep)
       import :: bc_mp_t, sp
       class(bc_mp_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=sp), intent(inout), dimension(n) :: x
       real(kind=sp), intent(in), dimension(n) :: t
       integer, intent(in), optional :: tstep
     end subroutine bc_mp_apply_scalar_sp

     !> Apply the boundary condition to a vector field (SP)
     subroutine bc_mp_apply_vector_sp(this, x, y, z, n, t, tstep)
       import :: bc_mp_t, sp
       class(bc_mp_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=sp), intent(inout), dimension(n) :: x
       real(kind=sp), intent(inout), dimension(n) :: y
       real(kind=sp), intent(inout), dimension(n) :: z
       real(kind=sp), intent(in), dimension(n) :: t
       integer, intent(in), optional :: tstep
     end subroutine bc_mp_apply_vector_sp

     !> Apply the boundary condition to a scalar field (DP)
     subroutine bc_mp_apply_scalar_dp(this, x, n, t, tstep)
       import :: bc_mp_t, dp
       class(bc_mp_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=dp), intent(inout), dimension(n) :: x
       real(kind=dp), intent(in), dimension(n) :: t
       integer, intent(in), optional :: tstep
     end subroutine bc_mp_apply_scalar_dp

     !> Apply the boundary condition to a vector field (DP)
     subroutine bc_mp_apply_vector_dp(this, x, y, z, n, t, tstep)
       import :: bc_mp_t, dp
       class(bc_mp_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=dp), intent(inout), dimension(n) :: x
       real(kind=dp), intent(inout), dimension(n) :: y
       real(kind=dp), intent(inout), dimension(n) :: z
       real(kind=dp), intent(in), dimension(n) :: t
       integer, intent(in), optional :: tstep
     end subroutine bc_mp_apply_vector_dp

     !> Destructor
     subroutine bc_mp_destructor(this)
       import :: bc_mp_t
       class(bc_mp_t), intent(inout) :: this
     end subroutine bc_mp_destructor
  end interface

  interface bc_list_apply_mp
     module procedure bc_list_apply_scalar_mp, bc_list_apply_vector_mp
  end interface bc_list_apply_mp

  public :: bc_list_init_mp, bc_list_free_mp, bc_list_add_mp, &
       bc_list_apply_scalar_mp, bc_list_apply_vector_mp, bc_list_apply_mp

contains

  !> Constructor
  subroutine bc_mp_init_base(this, coef)
    class(bc_mp_t), intent(inout) :: this
    type(coef_dp_t), target, intent(in) :: coef

    call this%free_base

    this%coef => coef
    this%dof_dp => coef%dof
    this%msh => coef%msh
    this%Xh => coef%Xh
  end subroutine bc_mp_init_base

  !> Free base class data
  subroutine bc_mp_free_base(this)
    class(bc_mp_t), intent(inout) :: this

    call this%marked_facet%free()

    nullify(this%Xh)
    nullify(this%msh)
    nullify(this%dof)
    nullify(this%coef)

    if (allocated(this%msk)) deallocate(this%msk)
    if (allocated(this%facet)) deallocate(this%facet)
    ! device...
  end subroutine bc_mp_free_base

  !> Mark a facet on an element as part of the boundary condition
  subroutine bc_mp_mark_facet(this, facet, el)
    class(bc_mp_t), intent(inout) :: this
    integer, intent(in) :: facet
    integer, intent(in) :: el
    type(tuple_i4_t) :: t

    t%x = (/facet, el/)
    call this%marked_facet%push(t)
  end subroutine bc_mp_mark_facet

  !> Finalize the construction of the bc by populating the msk and facet arrays
  subroutine bc_mp_finalize(this)
    class(bc_mp_t), target, intent(inout) :: this
    type(tuple_i4_t), pointer :: bfp(:)
    type(tuple_i4_t) :: bc_facet
    integer :: facet_size, facet, el
    integer :: i, j, k, l, msk_c
    integer :: lx, ly, lz, n

    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz

    !>@todo add 2D case

    ! Note we assume that lx = ly = lz
    facet_size = lx**2
    allocate(this%msk(0:facet_size * this%marked_facet%size()))
    allocate(this%facet(0:facet_size * this%marked_facet%size()))

    msk_c = 0
    bfp => this%marked_facet%array()

    ! Loop through each (facet, element) id tuple
    ! Then loop over all the nodes of the face and compute their linear index
    ! This index goes into this%msk, whereas the corresponding face id goes into
    ! this%facet
    do i = 1, this%marked_facet%size()
       bc_facet = bfp(i)
       facet = bc_facet%x(1)
       el = bc_facet%x(2)
       select case (facet)
       case (1)
          do l = 1, lz
             do k = 1, ly
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(1,k,l,el,lx,ly,lz)
                this%facet(msk_c) = 1
             end do
          end do
       case (2)
          do l = 1, lz
             do k = 1, ly
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(lx,k,l,el,lx,ly,lz)
                this%facet(msk_c) = 2
             end do
          end do
       case(3)
          do l = 1, lz
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j,1,l,el,lx,ly,lz)
                this%facet(msk_c) = 3
             end do
          end do
       case(4)
          do l = 1, lz
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j,ly,l,el,lx,ly,lz)
                this%facet(msk_c) = 4
             end do
          end do
       case(5)
          do k = 1, ly
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j,k,1,el,lx,ly,lz)
                this%facet(msk_c) = 5
             end do
          end do
       case(6)
          do k = 1, ly
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j,k,lz,el,lx,ly,lz)
                this%facet(msk_c) = 6
             end do
          end do
       end select
    end do

    this%msk(0) = msk_c
    this%facet(0) = msk_c

    if (NEKO_BCKND_DEVICE .eq. 1) then
       !!!
    end if
  end subroutine bc_mp_finalize

  !> Constructor for a list of boundary conditions
  !! @param size The size of the list to allocate.
  subroutine bc_list_init_mp(bclst, size)
    type(bc_list_mp_t), intent(inout), target :: bclst
    integer, optional :: size
    integer :: n, i

    call bc_list_free_mp(bclst)

    if (present(size)) then
       n = size
    else
       n = 1
    end if

    allocate(bclst%bc(n))

    do i = 1, n
       bclst%bc(i)%bcp => null()
    end do

    bclst%n = 0
    bclst%size = n
  end subroutine bc_list_init_mp

  !> Destructor for a list of boundary conditions
  !! @note This will only nullify all pointers, not deallocate any
  !! conditions pointed to by the list
  subroutine bc_list_free_mp(bclst)
    type(bc_list_mp_t), intent(inout) :: bclst

    if (allocated(bclst%bc)) then
       deallocate(bclst%bc)
    end if

    bclst%n = 0
    bclst%size = 0

  end subroutine bc_list_free_mp

!> Add a condition to a list of boundary conditions
  !! @param bc The boundary condition to add.
  subroutine bc_list_add_mp(bclst, bc)
    type(bc_list_mp_t), intent(inout) :: bclst
    class(bc_mp_t), intent(inout), target :: bc
    type(bcp_mp_t), allocatable :: tmp(:)

    !> Do not add if bc is empty
    if(bc%marked_facet%size() .eq. 0) return

    if (bclst%n .ge. bclst%size) then
       bclst%size = bclst%size * 2
       allocate(tmp(bclst%size))
       tmp(1:bclst%n) = bclst%bc
       call move_alloc(tmp, bclst%bc)
    end if

    bclst%n = bclst%n + 1
    bclst%bc(bclst%n)%bcp => bc

  end subroutine bc_list_add_mp


  !> Apply a list of boundary conditions to a scalar field
  subroutine bc_list_apply_scalar_mp(bclst, x, n, t, tstep)
    type(bc_list_mp_t), intent(inout) :: bclst
    integer, intent(in) :: n
    class(*), intent(inout), dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i
    real(kind=sp), allocatable :: t_sp(:)
    real(kind=dp), allocatable :: t_dp(:)

    ! Allocate and fill time arrays if needed
    if (present(t)) then
      allocate(t_sp(n))
      allocate(t_dp(n))
      t_sp = real(t, sp)
      t_dp = real(t, dp)
    end if

    ! Apply boundary conditions
    select type (x)
    type is (real(kind=sp))
       do i = 1, bclst%n
          select type (bcp => bclst%bc(i)%bcp)
          class is (bc_mp_t)
             if (present(t) .and. present(tstep)) then
                call bcp%apply_scalar_sp(x, n, t_sp, tstep)
             else if (present(t)) then
                call bcp%apply_scalar_sp(x, n, t_sp)
             else if (present(tstep)) then
                call bcp%apply_scalar_sp(x, n, [(real(0.0, sp), i=1,n)], tstep)
             else
                call bcp%apply_scalar_sp(x, n, [(real(0.0, sp), i=1,n)])
             end if
          end select
       end do
    type is (real(kind=dp))
       do i = 1, bclst%n
          select type (bcp => bclst%bc(i)%bcp)
          class is (bc_mp_t)
             if (present(t) .and. present(tstep)) then
                call bcp%apply_scalar_dp(x, n, t_dp, tstep)
             else if (present(t)) then
                call bcp%apply_scalar_dp(x, n, t_dp)
             else if (present(tstep)) then
                call bcp%apply_scalar_dp(x, n, [(real(0.0_dp, dp), i=1,n)], tstep)
             else
                call bcp%apply_scalar_dp(x, n, [(real(0.0_dp, dp), i=1,n)])
             end if
          end select
       end do
    end select

    ! Cleanup
    if (allocated(t_sp)) deallocate(t_sp)
    if (allocated(t_dp)) deallocate(t_dp)
  end subroutine bc_list_apply_scalar_mp

  subroutine bc_list_apply_vector_mp(bclst, x, y, z, n, t, tstep)
    type(bc_list_mp_t), intent(inout) :: bclst
    integer, intent(in) :: n
    class(*), intent(inout), dimension(n) :: x, y, z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i
    real(kind=sp), allocatable :: t_sp(:)
    real(kind=dp), allocatable :: t_dp(:)

    ! Allocate and fill time arrays if needed
    if (present(t)) then
      allocate(t_sp(n))
      allocate(t_dp(n))
      t_sp = real(t, sp)
      t_dp = real(t, dp)
    end if

    ! Apply boundary conditions
    select type (x)
    type is (real(kind=sp))
       select type (y)
       type is (real(kind=sp))
          select type (z)
          type is (real(kind=sp))
             do i = 1, bclst%n
                select type (bcp => bclst%bc(i)%bcp)
                class is (bc_mp_t)
                   if (present(t) .and. present(tstep)) then
                      call bcp%apply_vector_sp(x, y, z, n, t_sp, tstep)
                   else if (present(t)) then
                      call bcp%apply_vector_sp(x, y, z, n, t_sp)
                   else if (present(tstep)) then
                      call bcp%apply_vector_sp(x, y, z, n, [(real(0.0, sp), i=1,n)], tstep)
                   else
                      call bcp%apply_vector_sp(x, y, z, n, [(real(0.0, sp), i=1,n)])
                   end if
                end select
             end do
          end select
       end select
    type is (real(kind=dp))
       select type (y)
       type is (real(kind=dp))
          select type (z)
          type is (real(kind=dp))
             do i = 1, bclst%n
                select type (bcp => bclst%bc(i)%bcp)
                class is (bc_mp_t)
                   if (present(t) .and. present(tstep)) then
                      call bcp%apply_vector_dp(x, y, z, n, t_dp, tstep)
                   else if (present(t)) then
                      call bcp%apply_vector_dp(x, y, z, n, t_dp)
                   else if (present(tstep)) then
                      call bcp%apply_vector_dp(x, y, z, n, [(real(0.0_dp, dp), i=1,n)], tstep)
                   else
                      call bcp%apply_vector_dp(x, y, z, n, [(real(0.0_dp, dp), i=1,n)])
                   end if
                end select
             end do
          end select
       end select
    end select

    ! Cleanup
    if (allocated(t_sp)) deallocate(t_sp)
    if (allocated(t_dp)) deallocate(t_dp)
  end subroutine bc_list_apply_vector_mp

end module bc_mp
