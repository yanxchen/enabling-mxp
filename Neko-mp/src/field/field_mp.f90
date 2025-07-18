!> Mixed-precision field module
module field_mp
  use neko_config, only : NEKO_BCKND_DEVICE
  use device_math
  use num_types, only : rp, sp, dp
  use math, only : add2, copy, cadd, add2_sp, copy_sp, cadd_sp
  use mesh, only : mesh_t
  use space_mp, only : space_dp_t, space_sp_t, space_ne_sp, space_ne_dp
  use dofmap_mp, only : dofmap_sp_t, dofmap_dp_t
  use field
  implicit none

  private

  type, public :: field_sp_t
     real(kind=sp), allocatable :: x(:,:,:,:) !< Field data

     type(space_sp_t), pointer :: Xh   !< Function space \f$ X_h \f$
     type(mesh_t), pointer :: msh   !< Mesh
     type(dofmap_sp_t), pointer :: dof !< Dofmap

     logical :: internal_dofmap = .false. !< Does the field have an own dofmap
     character(len=80) :: name            !< Name of the field
   !   type(c_ptr) :: x_d = C_NULL_PTR
   contains
     procedure, private, pass(this) :: init_common => field_init_common_sp
     procedure, private, pass(this) :: init_external_dof => &
          field_init_external_dof_sp
     procedure, private, pass(this) :: init_internal_dof => &
          field_init_internal_dof_sp
     procedure, private, pass(this) :: assign_field => field_assign_field_sp
     procedure, private, pass(this) :: assign_scalar => field_assign_scalar_sp
     procedure, private, pass(this) :: add_field => field_add_field_sp
     procedure, private, pass(this) :: add_scalar => field_add_scalar_sp
     procedure, pass(this) :: free => field_free_sp
     !> Return the size of the field.
     procedure, pass(this) :: size => field_size_sp
     !> Initialise a field
     generic :: init => init_external_dof, init_internal_dof
     !> Assignemnt to current field
     generic :: assignment(=) => assign_field, assign_scalar
     !> Add to current field
     !! @note We don't overload operator(+), to avoid
     !! the extra assignemnt operator
     generic :: add => add_field, add_scalar
  end type field_sp_t

  !> field_ptr_sp_t, To easily obtain a pointer to a field
  type, public ::  field_ptr_sp_t
     type(field_sp_t), pointer :: ptr => null()
  end type field_ptr_sp_t

   type, public :: field_dp_t
      real(kind=dp), allocatable :: x(:,:,:,:) !< Field data

      type(space_dp_t), pointer :: Xh   !< Function space \f$ X_h \f$
      type(mesh_t), pointer :: msh   !< Mesh
      type(dofmap_dp_t), pointer :: dof !< Dofmap

      logical :: internal_dofmap = .false. !< Does the field have an own dofmap
      character(len=80) :: name            !< Name of the field
   !  type(c_ptr) :: x_d = C_NULL_PTR
   contains
      procedure, private, pass(this) :: init_common => field_init_common_dp
      procedure, private, pass(this) :: init_external_dof => &
         field_init_external_dof_dp
      procedure, private, pass(this) :: init_internal_dof => &
         field_init_internal_dof_dp
      procedure, private, pass(this) :: assign_field => field_assign_field_dp
      procedure, private, pass(this) :: assign_scalar => field_assign_scalar_dp
      procedure, private, pass(this) :: add_field => field_add_field_dp
      procedure, private, pass(this) :: add_scalar => field_add_scalar_dp
      procedure, pass(this) :: free => field_free_dp
      !> Return the size of the field.
      procedure, pass(this) :: size => field_size_dp
      !> Initialise a field
      generic :: init => init_external_dof, init_internal_dof
      !> Assignemnt to current field
      generic :: assignment(=) => assign_field, assign_scalar
      !> Add to current field
      !! @note We don't overload operator(+), to avoid
      !! the extra assignemnt operator
      generic :: add => add_field, add_scalar
   end type field_dp_t

   !> field_ptr_dp_t, To easily obtain a pointer to a field
   type, public ::  field_ptr_dp_t
      type(field_dp_t), pointer :: ptr => null()
   end type field_ptr_dp_t
    
contains

  !> Initialize a field @a this on the mesh @a msh using an internal dofmap
  subroutine field_init_internal_dof_sp(this, msh, space, fld_name)
    class(field_sp_t), intent(inout) :: this      !< Field to be initialized
    type(mesh_t), target, intent(in) :: msh    !< underlying mesh of the field
    type(space_sp_t), target, intent(in) :: space !< Function space for the field
    character(len=*), optional :: fld_name     !< Name of the field

    call this%free()

    this%Xh => space
    this%msh => msh

    allocate(this%dof)
    call this%dof%init(this%msh, this%Xh)
    this%internal_dofmap = .true.

    if (present(fld_name)) then
       call this%init_common(fld_name)
    else
       call this%init_common()
    end if

  end subroutine field_init_internal_dof_sp

  !> Initialize a field @a this on the mesh @a msh using an internal dofmap
  subroutine field_init_external_dof_sp(this, dof, fld_name)
    class(field_sp_t), intent(inout) :: this      !< Field to be initialized
    type(dofmap_sp_t), target, intent(in) :: dof  !< External dofmap for the field
    character(len=*), optional :: fld_name     !< Name of the field

    call this%free()

    this%dof => dof
    this%Xh => dof%Xh
    this%msh => dof%msh

    if (present(fld_name)) then
       call this%init_common(fld_name)
    else
       call this%init_common()
    end if

  end subroutine field_init_external_dof_sp

  !> Initialize a field @a this
  subroutine field_init_common_sp(this, fld_name)
    class(field_sp_t), intent(inout) :: this  !< Field to be initialized
    character(len=*), optional :: fld_name !< Name of the field
    integer :: ierr
    integer :: n

    associate(lx => this%Xh%lx, ly => this%Xh%ly, &
         lz => this%Xh%lz, nelv => this%msh%nelv)

      if (.not. allocated(this%x)) then
         allocate(this%x(lx, ly, lz, nelv), stat = ierr)
         this%x = 0d0
      end if

      if (present(fld_name)) then
         this%name = fld_name
      else
         this%name = "Field"
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         ! n = lx * ly * lz * nelv
         ! call device_map(this%x, this%x_d, n)
      end if
    end associate

  end subroutine field_init_common_sp

  !> Deallocate a field @a f
  subroutine field_free_sp(this)
    class(field_sp_t), intent(inout) :: this

    if (allocated(this%x)) then
       deallocate(this%x)
    end if

    if (this%internal_dofmap) then
       deallocate(this%dof)
       this%internal_dofmap = .false.
    end if

    nullify(this%msh)
    nullify(this%Xh)
    nullify(this%dof)

   !  if (c_associated(this%x_d)) then
   !     call device_free(this%x_d)
   !  end if

  end subroutine field_free_sp

  !> Assignment \f$ this = G \f$
  !! @note @a this will be initialized if it has a different size than
  !! @a G or it's not allocated
  subroutine field_assign_field_sp(this, g)
    class(field_sp_t), intent(inout) :: this
    type(field_sp_t), intent(in) :: g

    if (allocated(this%x)) then
       if (space_ne_sp(this%Xh, g%Xh)) then
          call this%free()
       end if
    end if

    this%Xh => g%Xh
    this%msh => g%msh
    this%dof => g%dof


    this%Xh%lx = g%Xh%lx
    this%Xh%ly = g%Xh%ly
    this%Xh%lz = g%Xh%lz

    if (.not. allocated(this%x)) then

       allocate(this%x(this%Xh%lx, this%Xh%ly, this%Xh%lz, this%msh%nelv))

       if (NEKO_BCKND_DEVICE .eq. 1) then
         !  call device_map(this%x, this%x_d, this%size())
       end if

    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
      !  call device_copy(this%x_d, g%x_d, this%size())
    else
       call copy_sp(this%x, g%x, this%dof%size())
    end if

  end subroutine field_assign_field_sp

  !> Assignment \f$ this = a \f$
  subroutine field_assign_scalar_sp(this, a)
    class(field_sp_t), intent(inout) :: this
    real(kind=sp), intent(in) :: a
    integer :: i, j, k, l

    if (NEKO_BCKND_DEVICE .eq. 1) then
      !  call device_cfill(this%x_d, a, this%size())
    else
       do i = 1, this%msh%nelv
          do l = 1, this%Xh%lz
             do k = 1, this%Xh%ly
                do j = 1, this%Xh%lx
                   this%x(j, k, l, i) = a
                end do
             end do
          end do
       end do
    end if

  end subroutine field_assign_scalar_sp

  !> Add \f$ this(u_1, u_2, ... , u_n) =
  !! this(u_1, u_2, ... , u_n) + G(u_1, u_2, ... , u_n) \f$
  !! @note Component wise
  subroutine field_add_field_sp(this, g)
    class(field_sp_t), intent(inout) :: this
    type(field_sp_t), intent(in) :: g

    if (NEKO_BCKND_DEVICE .eq. 1) then
      !  call device_add2(this%x_d, g%x_d, this%size())
    else
       call add2_sp(this%x, g%x, this%size())
    end if

  end subroutine field_add_field_sp


  !> Add \f$ this(u_1, u_2, ... , u_n) =
  !! this(u_1, u_2, ... , u_n) + a \f$
  subroutine field_add_scalar_sp(this, a)
    class(field_sp_t), intent(inout) :: this
    real(kind=sp), intent(in) :: a

    if (NEKO_BCKND_DEVICE .eq. 1) then
      !  call device_cadd(this%x_d, a, this%size())
    else
       call cadd_sp(this%x, a, this%size())
    end if

  end subroutine field_add_scalar_sp

  !> Return the size of the field based on the underlying dofmap.
  pure function field_size_sp(this) result(size)
    class(field_sp_t), intent(in) :: this
    integer :: size

    size = this%dof%size()
  end function field_size_sp

    subroutine field_init_internal_dof_dp(this, msh, space, fld_name)
      class(field_dp_t), intent(inout) :: this      !< Field to be initialized
      type(mesh_t), target, intent(in) :: msh    !< underlying mesh of the field
      type(space_dp_t), target, intent(in) :: space !< Function space for the field
      character(len=*), optional :: fld_name     !< Name of the field
  
      call this%free()
  
      this%Xh => space
      this%msh => msh
  
      allocate(this%dof)
      call this%dof%init(this%msh, this%Xh)
      this%internal_dofmap = .true.
  
      if (present(fld_name)) then
         call this%init_common(fld_name)
      else
         call this%init_common()
      end if
  
    end subroutine field_init_internal_dof_dp
  
    !> Initialize a field @a this on the mesh @a msh using an internal dofmap
    subroutine field_init_external_dof_dp(this, dof, fld_name)
      class(field_dp_t), intent(inout) :: this      !< Field to be initialized
      type(dofmap_dp_t), target, intent(in) :: dof  !< External dofmap for the field
      character(len=*), optional :: fld_name     !< Name of the field
  
      call this%free()
  
      this%dof => dof
      this%Xh => dof%Xh
      this%msh => dof%msh
  
      if (present(fld_name)) then
         call this%init_common(fld_name)
      else
         call this%init_common()
      end if
  
    end subroutine field_init_external_dof_dp
  
    !> Initialize a field @a this
    subroutine field_init_common_dp(this, fld_name)
      class(field_dp_t), intent(inout) :: this  !< Field to be initialized
      character(len=*), optional :: fld_name !< Name of the field
      integer :: ierr
      integer :: n
  
      associate(lx => this%Xh%lx, ly => this%Xh%ly, &
           lz => this%Xh%lz, nelv => this%msh%nelv)
  
        if (.not. allocated(this%x)) then
           allocate(this%x(lx, ly, lz, nelv), stat = ierr)
           this%x = 0d0
        end if
  
        if (present(fld_name)) then
           this%name = fld_name
        else
           this%name = "Field"
        end if
  
        if (NEKO_BCKND_DEVICE .eq. 1) then
         !   n = lx * ly * lz * nelv
         !   call device_map(this%x, this%x_d, n)
        end if
      end associate
  
    end subroutine field_init_common_dp
  
    !> Deallocate a field @a f
    subroutine field_free_dp(this)
      class(field_dp_t), intent(inout) :: this
  
      if (allocated(this%x)) then
         deallocate(this%x)
      end if
  
      if (this%internal_dofmap) then
         deallocate(this%dof)
         this%internal_dofmap = .false.
      end if
  
      nullify(this%msh)
      nullify(this%Xh)
      nullify(this%dof)
  
      ! if (c_associated(this%x_d)) then
      !    call device_free(this%x_d)
      ! end if
  
    end subroutine field_free_dp
  
    !> Assignment \f$ this = G \f$
    !! @note @a this will be initialized if it has a different size than
    !! @a G or it's not allocated
    subroutine field_assign_field_dp(this, g)
      class(field_dp_t), intent(inout) :: this
      type(field_dp_t), intent(in) :: g
  
      if (allocated(this%x)) then
         if (space_ne_dp(this%Xh, g%Xh)) then
            call this%free()
         end if
      end if
  
      this%Xh => g%Xh
      this%msh => g%msh
      this%dof => g%dof
  
  
      this%Xh%lx = g%Xh%lx
      this%Xh%ly = g%Xh%ly
      this%Xh%lz = g%Xh%lz
  
      if (.not. allocated(this%x)) then
  
         allocate(this%x(this%Xh%lx, this%Xh%ly, this%Xh%lz, this%msh%nelv))
  
         if (NEKO_BCKND_DEVICE .eq. 1) then
            ! call device_map(this%x, this%x_d, this%size())
         end if
  
      end if
  
      if (NEKO_BCKND_DEVICE .eq. 1) then
         ! call device_copy(this%x_d, g%x_d, this%size())
      else
         call copy(this%x, g%x, this%dof%size())
      end if
  
    end subroutine field_assign_field_dp
  
    !> Assignment \f$ this = a \f$
    subroutine field_assign_scalar_dp(this, a)
      class(field_dp_t), intent(inout) :: this
      real(kind=dp), intent(in) :: a
      integer :: i, j, k, l
  
      if (NEKO_BCKND_DEVICE .eq. 1) then
         ! call device_cfill(this%x_d, a, this%size())
      else
         do i = 1, this%msh%nelv
            do l = 1, this%Xh%lz
               do k = 1, this%Xh%ly
                  do j = 1, this%Xh%lx
                     this%x(j, k, l, i) = a
                  end do
               end do
            end do
         end do
      end if
  
    end subroutine field_assign_scalar_dp
  
    !> Add \f$ this(u_1, u_2, ... , u_n) =
    !! this(u_1, u_2, ... , u_n) + G(u_1, u_2, ... , u_n) \f$
    !! @note Component wise
    subroutine field_add_field_dp(this, g)
      class(field_dp_t), intent(inout) :: this
      type(field_dp_t), intent(in) :: g
  
      if (NEKO_BCKND_DEVICE .eq. 1) then
         ! call device_add2(this%x_d, g%x_d, this%size())
      else
         call add2(this%x, g%x, this%size())
      end if
  
    end subroutine field_add_field_dp
  
  
    !> Add \f$ this(u_1, u_2, ... , u_n) =
    !! this(u_1, u_2, ... , u_n) + a \f$
    subroutine field_add_scalar_dp(this, a)
      class(field_dp_t), intent(inout) :: this
      real(kind=dp), intent(in) :: a
  
      if (NEKO_BCKND_DEVICE .eq. 1) then
         ! call device_cadd(this%x_d, a, this%size())
      else
         call cadd(this%x, a, this%size())
      end if
  
    end subroutine field_add_scalar_dp
  
    !> Return the size of the field based on the underlying dofmap.
    pure function field_size_dp(this) result(size)
      class(field_dp_t), intent(in) :: this
      integer :: size
  
      size = this%dof%size()
    end function field_size_dp

end module field_mp
