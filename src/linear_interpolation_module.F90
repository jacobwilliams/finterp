!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!  Multidimensional linear interpolation/extrapolation.
!
!  Uses repeated linear interpolation to evaluate
!  functions \(f(x), f(x,y), f(x,y,z), f(x,y,z,q), f(x,y,z,q,r), f(x,y,z,q,r,s) \)
!  which have been tabulated at the nodes of an n-dimensional rectangular grid.
!  If any coordinate \( (x_i, y_i, ...) \) lies outside the range of the corresponding
!  variable, then extrapolation is performed using the two nearest points.
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:
#ifdef REAL32
!      `real(kind=real32)` [4 bytes]
#elif REAL64
!      `real(kind=real64)` [8 bytes]
#elif REAL128
!      `real(kind=real128)` [16 bytes]
#else
!      `real(kind=real64)` [8 bytes]
#endif

    module linear_interpolation_module

    use iso_fortran_env

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: finterp_rk = real32   !! real kind used by this module [4 bytes]
#elif REAL64
    integer,parameter,public :: finterp_rk = real64   !! real kind used by this module [8 bytes]
#elif REAL128
    integer,parameter,public :: finterp_rk = real128  !! real kind used by this module [16 bytes]
#else
    integer,parameter,public :: finterp_rk = real64   !! real kind used by this module [8 bytes]
#endif

    integer,parameter :: wp = finterp_rk  !! local copy of `finterp_rk` with a shorter name

    real(wp),parameter,private :: zero = 0.0_wp  !! numeric constant
    real(wp),parameter,private :: one  = 1.0_wp  !! numeric constant

    type,public,abstract :: linear_interp_class
        !! Base class for the linear interpolation types
        private
        logical :: initialized = .false. !! if the class was properly initialized
    contains
        private
        procedure(destroy_func),deferred,public :: destroy  !! destructor
        procedure :: check_inputs
    end type linear_interp_class

    abstract interface
        pure elemental subroutine destroy_func(me)  !! interface for bspline destructor routines
        import :: linear_interp_class
        implicit none
        class(linear_interp_class),intent(inout) :: me
        end subroutine destroy_func
    end interface

    type,extends(linear_interp_class),public :: linear_interp_1d
        !! Class for 1d linear interpolation.
        private
        real(wp),dimension(:),allocatable :: f
        real(wp),dimension(:),allocatable :: x
        integer :: ilox = 1
        contains
        private
        procedure,public :: initialize => initialize_1d
        procedure,public :: evaluate   => interp_1d
        procedure,public :: destroy    => destroy_1d
        final :: finalize_1d
    end type linear_interp_1d

    type,extends(linear_interp_class),public :: linear_interp_2d
        !! Class for 2d linear interpolation.
        private
        real(wp),dimension(:,:),allocatable :: f
        real(wp),dimension(:),allocatable :: x
        real(wp),dimension(:),allocatable :: y
        integer :: ilox = 1
        integer :: iloy = 1
        contains
        private
        procedure,public :: initialize => initialize_2d
        procedure,public :: evaluate   => interp_2d
        procedure,public :: destroy    => destroy_2d
        final :: finalize_2d
    end type linear_interp_2d

    type,extends(linear_interp_class),public :: linear_interp_3d
        !! Class for 3d linear interpolation.
        private
        real(wp),dimension(:,:,:),allocatable :: f
        real(wp),dimension(:),allocatable :: x
        real(wp),dimension(:),allocatable :: y
        real(wp),dimension(:),allocatable :: z
        integer :: ilox  = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        contains
        private
        procedure,public :: initialize => initialize_3d
        procedure,public :: evaluate   => interp_3d
        procedure,public :: destroy    => destroy_3d
        final :: finalize_3d
    end type linear_interp_3d

    type,extends(linear_interp_class),public :: linear_interp_4d
        !! Class for 4d linear interpolation.
        private
        real(wp),dimension(:,:,:,:),allocatable :: f
        real(wp),dimension(:),allocatable :: x
        real(wp),dimension(:),allocatable :: y
        real(wp),dimension(:),allocatable :: z
        real(wp),dimension(:),allocatable :: q
        integer :: ilox  = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        integer :: iloq  = 1
        contains
        private
        procedure,public :: initialize => initialize_4d
        procedure,public :: evaluate   => interp_4d
        procedure,public :: destroy    => destroy_4d
        final :: finalize_4d
    end type linear_interp_4d

    type,extends(linear_interp_class),public :: linear_interp_5d
        !! Class for 5d linear interpolation.
        private
        real(wp),dimension(:,:,:,:,:),allocatable :: f
        real(wp),dimension(:),allocatable :: x
        real(wp),dimension(:),allocatable :: y
        real(wp),dimension(:),allocatable :: z
        real(wp),dimension(:),allocatable :: q
        real(wp),dimension(:),allocatable :: r
        integer :: ilox  = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        integer :: iloq  = 1
        integer :: ilor  = 1
        contains
        private
        procedure,public :: initialize => initialize_5d
        procedure,public :: evaluate   => interp_5d
        procedure,public :: destroy    => destroy_5d
        final :: finalize_5d
    end type linear_interp_5d

    type,extends(linear_interp_class),public :: linear_interp_6d
        !! Class for 6d linear interpolation.
        private
        real(wp),dimension(:,:,:,:,:,:),allocatable :: f
        real(wp),dimension(:),allocatable :: x
        real(wp),dimension(:),allocatable :: y
        real(wp),dimension(:),allocatable :: z
        real(wp),dimension(:),allocatable :: q
        real(wp),dimension(:),allocatable :: r
        real(wp),dimension(:),allocatable :: s
        integer :: ilox  = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        integer :: iloq  = 1
        integer :: ilor  = 1
        integer :: ilos  = 1
        contains
        private
        procedure,public :: initialize => initialize_6d
        procedure,public :: evaluate   => interp_6d
        procedure,public :: destroy    => destroy_6d
        final :: finalize_6d
    end type linear_interp_6d

    type,extends(linear_interp_1d),public :: nearest_interp_1d
        !! Class for 1d nearest neighbor interpolation.
        contains
        procedure,public :: evaluate => nearest_1d
    end type nearest_interp_1d

    type,extends(linear_interp_2d),public :: nearest_interp_2d
        !! Class for 2d nearest neighbor interpolation.
        contains
        procedure,public :: evaluate => nearest_2d
    end type nearest_interp_2d

    type,extends(linear_interp_3d),public :: nearest_interp_3d
        !! Class for 3d nearest neighbor interpolation.
        contains
        procedure,public :: evaluate => nearest_3d
    end type nearest_interp_3d

    type,extends(linear_interp_4d),public :: nearest_interp_4d
        !! Class for 4d nearest neighbor interpolation.
        contains
        procedure,public :: evaluate => nearest_4d
    end type nearest_interp_4d

    type,extends(linear_interp_5d),public :: nearest_interp_5d
        !! Class for 5d nearest neighbor interpolation.
        contains
        procedure,public :: evaluate => nearest_5d
    end type nearest_interp_5d

    type,extends(linear_interp_6d),public :: nearest_interp_6d
        !! Class for 6d nearest neighbor interpolation.
        contains
        procedure,public :: evaluate => nearest_6d
    end type nearest_interp_6d

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_1d]] type.

    pure elemental subroutine finalize_1d(me)

    implicit none

    type(linear_interp_1d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_2d]] type.

    pure elemental subroutine finalize_2d(me)

    implicit none

    type(linear_interp_2d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_3d]] type.

    pure elemental subroutine finalize_3d(me)

    implicit none

    type(linear_interp_3d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_4d]] type.

    pure elemental subroutine finalize_4d(me)

    implicit none

    type(linear_interp_4d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_5d]] type.

    pure elemental subroutine finalize_5d(me)

    implicit none

    type(linear_interp_5d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_6d]] type.

    pure elemental subroutine finalize_6d(me)

    implicit none

    type(linear_interp_6d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_1d]] class.

    pure elemental subroutine destroy_1d(me)

    implicit none

    class(linear_interp_1d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    me%ilox = 1
    me%initialized = .false.

    end subroutine destroy_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_2d]] class.

    pure elemental subroutine destroy_2d(me)

    implicit none

    class(linear_interp_2d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    me%ilox = 1
    me%iloy = 1
    me%initialized = .false.

    end subroutine destroy_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_3d]] class.

    pure elemental subroutine destroy_3d(me)

    implicit none

    class(linear_interp_3d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    if (allocated(me%z)) deallocate(me%z)
    me%ilox = 1
    me%iloy = 1
    me%iloz = 1
    me%initialized = .false.

    end subroutine destroy_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_4d]] class.

    pure elemental subroutine destroy_4d(me)

    implicit none

    class(linear_interp_4d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    if (allocated(me%z)) deallocate(me%z)
    if (allocated(me%q)) deallocate(me%q)
    me%ilox = 1
    me%iloy = 1
    me%iloz = 1
    me%iloq = 1
    me%initialized = .false.

    end subroutine destroy_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_5d]] class.

    pure elemental subroutine destroy_5d(me)

    implicit none

    class(linear_interp_5d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    if (allocated(me%z)) deallocate(me%z)
    if (allocated(me%q)) deallocate(me%q)
    if (allocated(me%r)) deallocate(me%r)
    me%ilox = 1
    me%iloy = 1
    me%iloz = 1
    me%iloq = 1
    me%ilor = 1
    me%initialized = .false.

    end subroutine destroy_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_6d]] class.

    pure elemental subroutine destroy_6d(me)

    implicit none

    class(linear_interp_6d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    if (allocated(me%z)) deallocate(me%z)
    if (allocated(me%q)) deallocate(me%q)
    if (allocated(me%r)) deallocate(me%r)
    if (allocated(me%s)) deallocate(me%s)
    me%ilox = 1
    me%iloy = 1
    me%iloz = 1
    me%iloq = 1
    me%ilor = 1
    me%ilos = 1
    me%initialized = .false.

    end subroutine destroy_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_1d]] class.

    pure subroutine initialize_1d(me,x,f,istat)

    implicit none

    class(linear_interp_1d),intent(inout) :: me
    real(wp),dimension(:),intent(in)      :: x
    real(wp),dimension(:),intent(in)      :: f
    integer,intent(out)                   :: istat  !! `0`   : no problems,
                                                    !! `1`   : `x` is not strictly increasing,
                                                    !! `10`  : `x` is not equal to size(f,1),
                                                    !! `100` : cannot use linear interpolation for only one point.

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10

    if (istat==0) then
        call me%check_inputs(x=x,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x))); me%f = f
            allocate(me%x(size(x))); me%x = x
            me%initialized = .true.
        end if
    end if

    end subroutine initialize_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_2d]] class.

    pure subroutine initialize_2d(me,x,y,f,istat)

    implicit none

    class(linear_interp_2d),intent(inout) :: me
    real(wp),dimension(:),intent(in)      :: x
    real(wp),dimension(:),intent(in)      :: y
    real(wp),dimension(:,:),intent(in)    :: f
    integer,intent(out)                   :: istat  !! `0`   : no problems,
                                                    !! `1`   : `x` is not strictly increasing,
                                                    !! `2`   : `y` is not strictly increasing,
                                                    !! `10`  : `x` is not equal to size(f,1),
                                                    !! `20`  : `y` is not equal to size(f,2),
                                                    !! `100` : cannot use linear interpolation for only one point.

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20

    if (istat==0) then
        call me%check_inputs(x=x,y=y,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            me%initialized = .true.
        end if
    end if

    end subroutine initialize_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_3d]] class.

    pure subroutine initialize_3d(me,x,y,z,f,istat)

    implicit none

    class(linear_interp_3d),intent(inout) :: me
    real(wp),dimension(:),intent(in)      :: x
    real(wp),dimension(:),intent(in)      :: y
    real(wp),dimension(:),intent(in)      :: z
    real(wp),dimension(:,:,:),intent(in)  :: f
    integer,intent(out)                   :: istat  !! `0`   : no problems,
                                                    !! `1`   : `x` is not strictly increasing,
                                                    !! `2`   : `y` is not strictly increasing,
                                                    !! `3`   : `z` is not strictly increasing,
                                                    !! `10`  : `x` is not equal to size(f,1),
                                                    !! `20`  : `y` is not equal to size(f,2),
                                                    !! `30`  : `z` is not equal to size(f,3),
                                                    !! `100` : cannot use linear interpolation for only one point.

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20
    if (istat==0 .and. size(z)/=size(f,3)) istat = 30

    if (istat==0) then
        call me%check_inputs(x=x,y=y,z=z,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y),size(z))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            allocate(me%z(size(z))); me%z = z
            me%initialized = .true.
        end if
    end if

    end subroutine initialize_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_4d]] class.

    pure subroutine initialize_4d(me,x,y,z,q,f,istat)

    implicit none

    class(linear_interp_4d),intent(inout)  :: me
    real(wp),dimension(:),intent(in)       :: x
    real(wp),dimension(:),intent(in)       :: y
    real(wp),dimension(:),intent(in)       :: z
    real(wp),dimension(:),intent(in)       :: q
    real(wp),dimension(:,:,:,:),intent(in) :: f
    integer,intent(out)                    :: istat !! `0`   : no problems,
                                                    !! `1`   : `x` is not strictly increasing,
                                                    !! `2`   : `y` is not strictly increasing,
                                                    !! `3`   : `z` is not strictly increasing,
                                                    !! `4`   : `q` is not strictly increasing,
                                                    !! `10`  : `x` is not equal to size(f,1),
                                                    !! `20`  : `y` is not equal to size(f,2),
                                                    !! `30`  : `z` is not equal to size(f,3),
                                                    !! `40`  : `q` is not equal to size(f,4),
                                                    !! `100` : cannot use linear interpolation for only one point.

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20
    if (istat==0 .and. size(z)/=size(f,3)) istat = 30
    if (istat==0 .and. size(q)/=size(f,4)) istat = 40

    if (istat==0) then
        call me%check_inputs(x=x,y=y,z=z,q=q,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y),size(z),size(q))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            allocate(me%z(size(z))); me%z = z
            allocate(me%q(size(q))); me%q = q
            me%initialized = .true.
        end if
    end if

    end subroutine initialize_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_5d]] class.

    pure subroutine initialize_5d(me,x,y,z,q,r,f,istat)

    implicit none

    class(linear_interp_5d),intent(inout)    :: me
    real(wp),dimension(:),intent(in)         :: x
    real(wp),dimension(:),intent(in)         :: y
    real(wp),dimension(:),intent(in)         :: z
    real(wp),dimension(:),intent(in)         :: q
    real(wp),dimension(:),intent(in)         :: r
    real(wp),dimension(:,:,:,:,:),intent(in) :: f
    integer,intent(out)                      :: istat   !! `0`   : no problems,
                                                        !! `1`   : `x` is not strictly increasing,
                                                        !! `2`   : `y` is not strictly increasing,
                                                        !! `3`   : `z` is not strictly increasing,
                                                        !! `4`   : `q` is not strictly increasing,
                                                        !! `5`   : `r` is not strictly increasing,
                                                        !! `10`  : `x` is not equal to size(f,1),
                                                        !! `20`  : `y` is not equal to size(f,2),
                                                        !! `30`  : `z` is not equal to size(f,3),
                                                        !! `40`  : `q` is not equal to size(f,4),
                                                        !! `50`  : `r` is not equal to size(f,5),
                                                        !! `100` : cannot use linear interpolation for only one point.

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20
    if (istat==0 .and. size(z)/=size(f,3)) istat = 30
    if (istat==0 .and. size(q)/=size(f,4)) istat = 40
    if (istat==0 .and. size(r)/=size(f,5)) istat = 50

    if (istat==0) then
        call me%check_inputs(x=x,y=y,z=z,q=q,r=r,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y),size(z),size(q),size(r))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            allocate(me%z(size(z))); me%z = z
            allocate(me%q(size(q))); me%q = q
            allocate(me%r(size(r))); me%r = r
            me%initialized = .true.
        end if
    end if

    end subroutine initialize_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_6d]] class.

    pure subroutine initialize_6d(me,x,y,z,q,r,s,f,istat)

    implicit none

    class(linear_interp_6d),intent(inout)      :: me
    real(wp),dimension(:),intent(in)           :: x
    real(wp),dimension(:),intent(in)           :: y
    real(wp),dimension(:),intent(in)           :: z
    real(wp),dimension(:),intent(in)           :: q
    real(wp),dimension(:),intent(in)           :: r
    real(wp),dimension(:),intent(in)           :: s
    real(wp),dimension(:,:,:,:,:,:),intent(in) :: f
    integer,intent(out)                        :: istat !! `0`   : no problems,
                                                        !! `1`   : `x` is not strictly increasing,
                                                        !! `2`   : `y` is not strictly increasing,
                                                        !! `3`   : `z` is not strictly increasing,
                                                        !! `4`   : `q` is not strictly increasing,
                                                        !! `5`   : `r` is not strictly increasing,
                                                        !! `6`   : `s` is not strictly increasing,
                                                        !! `10`  : `x` is not equal to size(f,1),
                                                        !! `20`  : `y` is not equal to size(f,2),
                                                        !! `30`  : `z` is not equal to size(f,3),
                                                        !! `40`  : `q` is not equal to size(f,4),
                                                        !! `50`  : `r` is not equal to size(f,5),
                                                        !! `60`  : `s` is not equal to size(f,6),
                                                        !! `100` : cannot use linear interpolation for only one point.

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20
    if (istat==0 .and. size(z)/=size(f,3)) istat = 30
    if (istat==0 .and. size(q)/=size(f,4)) istat = 40
    if (istat==0 .and. size(r)/=size(f,5)) istat = 50
    if (istat==0 .and. size(s)/=size(f,6)) istat = 60

    if (istat==0) then
        call me%check_inputs(x=x,y=y,z=z,q=q,r=r,s=s,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y),size(z),size(q),size(r),size(s))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            allocate(me%z(size(z))); me%z = z
            allocate(me%q(size(q))); me%q = q
            allocate(me%r(size(r))); me%r = r
            allocate(me%s(size(s))); me%s = s
            me%initialized = .true.
        end if
    end if

    end subroutine initialize_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  1D linear interpolation routine.

    pure subroutine interp_1d(me,x,f,istat)

    implicit none

    class(linear_interp_1d),intent(inout) :: me
    real(wp),intent(in)                   :: x
    real(wp),intent(out)                  :: f     !! Interpolated \( f(x) \)
    integer,intent(out),optional          :: istat !! `0`  : no problems,
                                                   !! `-1` : class has not been initialized

    integer,dimension(2) :: ix
    real(wp) :: p1
    real(wp) :: q1
    integer :: mflag

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)

        q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
        p1 = one-q1

        f = p1*me%f(ix(1)) + q1*me%f(ix(2))
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine interp_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  2D linear interpolation routine.

    pure subroutine interp_2d(me,x,y,f,istat)

    implicit none

    class(linear_interp_2d),intent(inout) :: me
    real(wp),intent(in)                   :: x
    real(wp),intent(in)                   :: y
    real(wp),intent(out)                  :: f     !! Interpolated \( f(x,y) \)
    integer,intent(out),optional          :: istat !! `0`  : no problems,
                                                   !! `-1` : class has not been initialized

    integer,dimension(2) :: ix, iy
    real(wp) :: p1, p2
    real(wp) :: q1, q2
    integer :: mflag
    real(wp) :: fx1, fx2

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)

        q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
        q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
        p1 = one-q1
        p2 = one-q2

        fx1 = p1*me%f(ix(1),iy(1)) + q1*me%f(ix(2),iy(1))
        fx2 = p1*me%f(ix(1),iy(2)) + q1*me%f(ix(2),iy(2))

        f = p2*fx1 + q2*fx2
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine interp_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  3D linear interpolation routine.

    pure subroutine interp_3d(me,x,y,z,f,istat)

    implicit none

    class(linear_interp_3d),intent(inout) :: me
    real(wp),intent(in)                   :: x
    real(wp),intent(in)                   :: y
    real(wp),intent(in)                   :: z
    real(wp),intent(out)                  :: f     !! Interpolated \( f(x,y,z) \)
    integer,intent(out),optional          :: istat !! `0`  : no problems,
                                                   !! `-1` : class has not been initialized

    integer,dimension(2) :: ix, iy, iz
    real(wp) :: p1, p2, p3
    real(wp) :: q1, q2, q3
    integer :: mflag
    real(wp) :: fx11, fx21, fx12, fx22, fxy1, fxy2

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)
        call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag)

        q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
        q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
        q3 = (z-me%z(iz(1)))/(me%z(iz(2))-me%z(iz(1)))
        p1 = one-q1
        p2 = one-q2
        p3 = one-q3

        fx11 = p1*me%f(ix(1),iy(1),iz(1)) + q1*me%f(ix(2),iy(1),iz(1))
        fx21 = p1*me%f(ix(1),iy(2),iz(1)) + q1*me%f(ix(2),iy(2),iz(1))
        fx12 = p1*me%f(ix(1),iy(1),iz(2)) + q1*me%f(ix(2),iy(1),iz(2))
        fx22 = p1*me%f(ix(1),iy(2),iz(2)) + q1*me%f(ix(2),iy(2),iz(2))
        fxy1 = p2*fx11 + q2*fx21
        fxy2 = p2*fx12 + q2*fx22

        f = p3*fxy1 + q3*fxy2
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine interp_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  4D linear interpolation routine.

    pure subroutine interp_4d(me,x,y,z,q,f,istat)

    implicit none

    class(linear_interp_4d),intent(inout) :: me
    real(wp),intent(in)                   :: x
    real(wp),intent(in)                   :: y
    real(wp),intent(in)                   :: z
    real(wp),intent(in)                   :: q
    real(wp),intent(out)                  :: f     !! Interpolated \( f(x,y,z,q) \)
    integer,intent(out),optional          :: istat  !! `0`  : no problems,
                                                    !! `-1` : class has not been initialized

    integer,dimension(2) :: ix, iy, iz, iq
    real(wp) :: p1, p2, p3, p4
    real(wp) :: q1, q2, q3, q4
    integer :: mflag
    real(wp) :: fx111,fx211,fx121,fx221,fxy11,fxy21,fxyz1,&
                fx112,fx212,fx122,fx222,fxy12,fxy22,fxyz2

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)
        call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag)
        call dintrv(me%q,q,me%iloq,iq(1),iq(2),mflag)

        q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
        q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
        q3 = (z-me%z(iz(1)))/(me%z(iz(2))-me%z(iz(1)))
        q4 = (q-me%q(iq(1)))/(me%q(iq(2))-me%q(iq(1)))
        p1 = one-q1
        p2 = one-q2
        p3 = one-q3
        p4 = one-q4

        fx111 = p1*me%f(ix(1),iy(1),iz(1),iq(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(1))
        fx211 = p1*me%f(ix(1),iy(2),iz(1),iq(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(1))
        fx121 = p1*me%f(ix(1),iy(1),iz(2),iq(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(1))
        fx221 = p1*me%f(ix(1),iy(2),iz(2),iq(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(1))
        fx112 = p1*me%f(ix(1),iy(1),iz(1),iq(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(2))
        fx212 = p1*me%f(ix(1),iy(2),iz(1),iq(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(2))
        fx122 = p1*me%f(ix(1),iy(1),iz(2),iq(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(2))
        fx222 = p1*me%f(ix(1),iy(2),iz(2),iq(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(2))

        fxy11 = p2*fx111 + q2*fx211
        fxy21 = p2*fx121 + q2*fx221
        fxy12 = p2*fx112 + q2*fx212
        fxy22 = p2*fx122 + q2*fx222

        fxyz1 = p3*fxy11 + q3*fxy21
        fxyz2 = p3*fxy12 + q3*fxy22

        f = p4*fxyz1 + q4*fxyz2
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine interp_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  5D linear interpolation routine.

    pure subroutine interp_5d(me,x,y,z,q,r,f,istat)

    implicit none

    class(linear_interp_5d),intent(inout) :: me
    real(wp),intent(in)                   :: x
    real(wp),intent(in)                   :: y
    real(wp),intent(in)                   :: z
    real(wp),intent(in)                   :: q
    real(wp),intent(in)                   :: r
    real(wp),intent(out)                  :: f       !! Interpolated \( f(x,y,z,q,r) \)
    integer,intent(out),optional          :: istat   !! `0`  : no problems,
                                                     !! `-1` : class has not been initialized

    integer,dimension(2) :: ix, iy, iz, iq, ir
    real(wp) :: p1, p2, p3, p4, p5
    real(wp) :: q1, q2, q3, q4, q5
    integer :: mflag
    real(wp) :: fx1111, fx2111, fx1211, fx2211, fx1121, fx2121, fx1221, fx2221, &
                fxy111, fxy211, fxy121, fxy221, fxyz11, fxyz21, fxyzq1, fx1112, &
                fx2112, fx1212, fx2212, fx1122, fx2122, fx1222, fx2222, fxy112, &
                fxy212, fxy122, fxy222, fxyz12, fxyz22, fxyzq2

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)
        call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag)
        call dintrv(me%q,q,me%iloq,iq(1),iq(2),mflag)
        call dintrv(me%r,r,me%ilor,ir(1),ir(2),mflag)

        q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
        q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
        q3 = (z-me%z(iz(1)))/(me%z(iz(2))-me%z(iz(1)))
        q4 = (q-me%q(iq(1)))/(me%q(iq(2))-me%q(iq(1)))
        q5 = (r-me%r(ir(1)))/(me%r(ir(2))-me%r(ir(1)))
        p1 = one-q1
        p2 = one-q2
        p3 = one-q3
        p4 = one-q4
        p5 = one-q5

        fx1111 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(1))
        fx2111 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(1))
        fx1211 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(1))
        fx2211 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(1))
        fx1121 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(1))
        fx2121 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(1))
        fx1221 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(1))
        fx2221 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(1))
        fx1112 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(2))
        fx2112 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(2))
        fx1212 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(2))
        fx2212 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(2))
        fx1122 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(2))
        fx2122 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(2))
        fx1222 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(2))
        fx2222 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(2))

        fxy111 = p2*fx1111 + q2*fx2111
        fxy211 = p2*fx1211 + q2*fx2211
        fxy121 = p2*fx1121 + q2*fx2121
        fxy221 = p2*fx1221 + q2*fx2221
        fxy112 = p2*fx1112 + q2*fx2112
        fxy212 = p2*fx1212 + q2*fx2212
        fxy122 = p2*fx1122 + q2*fx2122
        fxy222 = p2*fx1222 + q2*fx2222

        fxyz11 = p3*fxy111 + q3*fxy211
        fxyz21 = p3*fxy121 + q3*fxy221
        fxyz12 = p3*fxy112 + q3*fxy212
        fxyz22 = p3*fxy122 + q3*fxy222

        fxyzq1 = p4*fxyz11 + q4*fxyz21
        fxyzq2 = p4*fxyz12 + q4*fxyz22

        f = p5*fxyzq1 + q5*fxyzq2
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine interp_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  6D linear interpolation routine.

    pure subroutine interp_6d(me,x,y,z,q,r,s,f,istat)

    implicit none

    class(linear_interp_6d),intent(inout) :: me
    real(wp),intent(in)                   :: x
    real(wp),intent(in)                   :: y
    real(wp),intent(in)                   :: z
    real(wp),intent(in)                   :: q
    real(wp),intent(in)                   :: r
    real(wp),intent(in)                   :: s
    real(wp),intent(out)                  :: f        !! Interpolated \( f(x,y,z,q,r,s) \)
    integer,intent(out),optional          :: istat    !! `0`  : no problems,
                                                      !! `-1` : class has not been initialized

    integer,dimension(2) :: ix, iy, iz, iq, ir, is
    real(wp) :: p1, p2, p3, p4, p5, p6
    real(wp) :: q1, q2, q3, q4, q5, q6
    integer :: mflag
    real(wp) :: fx11111, fx21111, fx12111, fx22111, fx11211, fx21211, fx12211, &
                fx22211, fxy1111, fxy2111, fxy1211, fxy2211, fxyz111, fxyz211, &
                fxyzq11, fx11121, fx21121, fx12121, fx22121, fx11221, fx21221, &
                fx12221, fx22221, fxy1121, fxy2121, fxy1221, fxy2221, fxyz121, &
                fxyz221, fxyzq21, fx11112, fx21112, fx12112, fx22112, fx11212, &
                fx21212, fx12212, fx22212, fxy1112, fxy2112, fxy1212, fxy2212, &
                fxyz112, fxyz212, fxyzq12, fx11122, fx21122, fx12122, fx22122, &
                fx11222, fx21222, fx12222, fx22222, fxy1122, fxy2122, fxy1222, &
                fxy2222, fxyz122, fxyz222, fxyzq22, fxyzqr1, fxyzqr2

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)
        call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag)
        call dintrv(me%q,q,me%iloq,iq(1),iq(2),mflag)
        call dintrv(me%r,r,me%ilor,ir(1),ir(2),mflag)
        call dintrv(me%s,s,me%ilos,is(1),is(2),mflag)

        q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
        q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
        q3 = (z-me%z(iz(1)))/(me%z(iz(2))-me%z(iz(1)))
        q4 = (q-me%q(iq(1)))/(me%q(iq(2))-me%q(iq(1)))
        q5 = (r-me%r(ir(1)))/(me%r(ir(2))-me%r(ir(1)))
        q6 = (s-me%s(is(1)))/(me%s(is(2))-me%s(is(1)))
        p1 = one-q1
        p2 = one-q2
        p3 = one-q3
        p4 = one-q4
        p5 = one-q5
        p6 = one-q6

        fx11111 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(1),is(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(1),is(1))
        fx21111 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(1),is(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(1),is(1))
        fx12111 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(1),is(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(1),is(1))
        fx22111 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(1),is(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(1),is(1))
        fx11211 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(1),is(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(1),is(1))
        fx21211 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(1),is(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(1),is(1))
        fx12211 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(1),is(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(1),is(1))
        fx22211 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(1),is(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(1),is(1))
        fx11121 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(2),is(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(2),is(1))
        fx21121 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(2),is(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(2),is(1))
        fx12121 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(2),is(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(2),is(1))
        fx22121 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(2),is(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(2),is(1))
        fx11221 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(2),is(1)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(2),is(1))
        fx21221 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(2),is(1)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(2),is(1))
        fx12221 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(2),is(1)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(2),is(1))
        fx22221 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(2),is(1)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(2),is(1))
        fx11112 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(1),is(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(1),is(2))
        fx21112 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(1),is(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(1),is(2))
        fx12112 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(1),is(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(1),is(2))
        fx22112 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(1),is(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(1),is(2))
        fx11212 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(1),is(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(1),is(2))
        fx21212 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(1),is(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(1),is(2))
        fx12212 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(1),is(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(1),is(2))
        fx22212 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(1),is(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(1),is(2))
        fx11122 = p1*me%f(ix(1),iy(1),iz(1),iq(1),ir(2),is(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(1),ir(2),is(2))
        fx21122 = p1*me%f(ix(1),iy(2),iz(1),iq(1),ir(2),is(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(1),ir(2),is(2))
        fx12122 = p1*me%f(ix(1),iy(1),iz(2),iq(1),ir(2),is(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(1),ir(2),is(2))
        fx22122 = p1*me%f(ix(1),iy(2),iz(2),iq(1),ir(2),is(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(1),ir(2),is(2))
        fx11222 = p1*me%f(ix(1),iy(1),iz(1),iq(2),ir(2),is(2)) + q1*me%f(ix(2),iy(1),iz(1),iq(2),ir(2),is(2))
        fx21222 = p1*me%f(ix(1),iy(2),iz(1),iq(2),ir(2),is(2)) + q1*me%f(ix(2),iy(2),iz(1),iq(2),ir(2),is(2))
        fx12222 = p1*me%f(ix(1),iy(1),iz(2),iq(2),ir(2),is(2)) + q1*me%f(ix(2),iy(1),iz(2),iq(2),ir(2),is(2))
        fx22222 = p1*me%f(ix(1),iy(2),iz(2),iq(2),ir(2),is(2)) + q1*me%f(ix(2),iy(2),iz(2),iq(2),ir(2),is(2))

        fxy1111 = p2*fx11111 + q2*fx21111
        fxy2111 = p2*fx12111 + q2*fx22111
        fxy1211 = p2*fx11211 + q2*fx21211
        fxy2211 = p2*fx12211 + q2*fx22211
        fxy1121 = p2*fx11121 + q2*fx21121
        fxy2121 = p2*fx12121 + q2*fx22121
        fxy1221 = p2*fx11221 + q2*fx21221
        fxy2221 = p2*fx12221 + q2*fx22221
        fxy1112 = p2*fx11112 + q2*fx21112
        fxy2112 = p2*fx12112 + q2*fx22112
        fxy1212 = p2*fx11212 + q2*fx21212
        fxy2212 = p2*fx12212 + q2*fx22212
        fxy1122 = p2*fx11122 + q2*fx21122
        fxy2122 = p2*fx12122 + q2*fx22122
        fxy1222 = p2*fx11222 + q2*fx21222
        fxy2222 = p2*fx12222 + q2*fx22222

        fxyz111 = p3*fxy1111 + q3*fxy2111
        fxyz211 = p3*fxy1211 + q3*fxy2211
        fxyz121 = p3*fxy1121 + q3*fxy2121
        fxyz221 = p3*fxy1221 + q3*fxy2221
        fxyz112 = p3*fxy1112 + q3*fxy2112
        fxyz212 = p3*fxy1212 + q3*fxy2212
        fxyz122 = p3*fxy1122 + q3*fxy2122
        fxyz222 = p3*fxy1222 + q3*fxy2222

        fxyzq11 = p4*fxyz111 + q4*fxyz211
        fxyzq21 = p4*fxyz121 + q4*fxyz221
        fxyzq12 = p4*fxyz112 + q4*fxyz212
        fxyzq22 = p4*fxyz122 + q4*fxyz222

        fxyzqr1 = p5*fxyzq11 + q5*fxyzq21
        fxyzqr2 = p5*fxyzq12 + q5*fxyzq22

        f = p6*fxyzqr1 + q6*fxyzqr2
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine interp_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the indices in `xt` that bound `x`, to use for interpolation.
!  If outside the range, then the indices are returned that can
!  be used for extrapolation.
!  Precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   iright=2,    mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   iright=i+1,  mflag=0
!         if   xt(n) <= x           then ileft=n-1, iright=n,    mflag=1
!```
!
!### History
!
!  * interv written by carl de boor [5]
!  * dintrv author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * Jacob Williams, 2/24/2015 : updated to free-form Fortran.
!  * Jacob Williams, 2/17/2016 : additional refactoring (eliminated GOTOs).
!  * Jacob Williams, 2/22/2016 : modified bspline-fortran `dintrv` routine for
!    linear interpolation/extrapolation use.
!  * Jacob Williams, 10/9/2019 : added optional `inearest` output.

    pure subroutine dintrv(xt,x,ilo,ileft,iright,mflag,inearest)

    implicit none

    real(wp),dimension(:),intent(in) :: xt       !! a knot or break point vector
    real(wp),intent(in)              :: x        !! argument
    integer,intent(inout)            :: ilo      !! an initialization parameter which must be set
                                                 !! to 1 the first time the array `xt` is
                                                 !! processed by dintrv. `ilo` contains information for
                                                 !! efficient processing after the initial call and `ilo`
                                                 !! must not be changed by the user.  each dimension
                                                 !! requires a distinct `ilo` parameter.
    integer,intent(out)              :: ileft    !! left index
    integer,intent(out)              :: iright   !! right index
    integer,intent(out)              :: mflag    !! signals when `x` lies out of bounds
    integer,intent(out),optional     :: inearest !! nearest index

    integer :: ihi, istep, imid, n

    n = size(xt)

    if (n==1) then
        ! this is only allowed for nearest interpolation
        if (present(inearest)) then
            inearest = 1
            return
        end if
    end if

    ihi = ilo + 1
    if ( ihi>=n ) then
        if ( x>=xt(n) ) then
            mflag = 1
            ileft = n-1
            iright= n
            if (present(inearest)) inearest = n
            return
        end if
        if ( n<=1 ) then
            mflag = -1
            ileft = 1
            iright= 2
            if (present(inearest)) inearest = 1
            return
        end if
        ilo = n - 1
        ihi = n
    endif

    if ( x>=xt(ihi) ) then

        ! now x >= xt(ilo). find upper bound
        istep = 1
        do
            ilo = ihi
            ihi = ilo + istep
            if ( ihi>=n ) then
                if ( x>=xt(n) ) then
                    mflag = 1
                    ileft = n-1
                    iright= n
                    if (present(inearest)) inearest = n
                    return
                end if
                ihi = n
            elseif ( x>=xt(ihi) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    else

        if ( x>=xt(ilo) ) then
            mflag = 0
            ileft = ilo
            iright= ilo+1
            if (present(inearest)) then
                if ( abs(x-xt(ileft)) <= abs(x-xt(iright)) ) then
                    inearest = ileft
                else
                    inearest = iright
                end if
            end if
            return
        end if
        ! now x <= xt(ihi). find lower bound
        istep = 1
        do
            ihi = ilo
            ilo = ihi - istep
            if ( ilo<=1 ) then
                ilo = 1
                if ( x<xt(1) ) then
                    mflag = -1
                    ileft = 1
                    iright= 2
                    if (present(inearest)) inearest = 1
                    return
                end if
            elseif ( x<xt(ilo) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    endif

    ! now xt(ilo) <= x < xt(ihi). narrow the interval
    do
        imid = (ilo+ihi)/2
        if ( imid==ilo ) then
            mflag = 0
            ileft = ilo
            iright= ilo+1
            if (present(inearest)) then
                if ( abs(x-xt(ileft)) <= abs(x-xt(iright)) ) then
                    inearest = ileft
                else
                    inearest = iright
                end if
            end if
            return
        end if
        ! note. it is assumed that imid = ilo in case ihi = ilo+1
        if ( x<xt(imid) ) then
            ihi = imid
        else
            ilo = imid
        endif
    end do

    end subroutine dintrv
!*****************************************************************************************

!*****************************************************************************************
!>
!  1D nearest neighbor interpolation routine.

    pure subroutine nearest_1d(me,x,f,istat)

    implicit none

    class(nearest_interp_1d),intent(inout) :: me
    real(wp),intent(in)                    :: x
    real(wp),intent(out)                   :: f     !! Nearest \( f(x) \)
    integer,intent(out),optional           :: istat !! `0`  : no problems,
                                                    !! `-1` : class has not been initialized

    integer :: mflag
    integer,dimension(2) :: ix
    integer :: i

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag,i)

        f = me%f(i)
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine nearest_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  2D nearest neighbor interpolation routine.

    pure subroutine nearest_2d(me,x,y,f,istat)

    implicit none

    class(nearest_interp_2d),intent(inout) :: me
    real(wp),intent(in)                    :: x
    real(wp),intent(in)                    :: y
    real(wp),intent(out)                   :: f     !! Nearest \( f(x,y) \)
    integer,intent(out),optional           :: istat !! `0`  : no problems,
                                                    !! `-1` : class has not been initialized

    integer :: mflag
    integer,dimension(2) :: ix, iy
    integer :: i, j

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag,i)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag,j)

        f = me%f(i,j)
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine nearest_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  3D nearest neighbor interpolation routine.

    pure subroutine nearest_3d(me,x,y,z,f,istat)

    implicit none

    class(nearest_interp_3d),intent(inout) :: me
    real(wp),intent(in)                    :: x
    real(wp),intent(in)                    :: y
    real(wp),intent(in)                    :: z
    real(wp),intent(out)                   :: f     !! Nearest \( f(x,y,z) \)
    integer,intent(out),optional           :: istat !! `0`  : no problems,
                                                    !! `-1` : class has not been initialized

    integer :: mflag
    integer,dimension(2) :: ix, iy, iz
    integer :: i, j, k

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag,i)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag,j)
        call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag,k)

        f = me%f(i,j,k)
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine nearest_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  4D nearest neighbor interpolation routine.

    pure subroutine nearest_4d(me,x,y,z,q,f,istat)

    implicit none

    class(nearest_interp_4d),intent(inout) :: me
    real(wp),intent(in)                    :: x
    real(wp),intent(in)                    :: y
    real(wp),intent(in)                    :: z
    real(wp),intent(in)                    :: q
    real(wp),intent(out)                   :: f     !! Nearest \( f(x,y,z,q) \)
    integer,intent(out),optional           :: istat !! `0`  : no problems,
                                                    !! `-1` : class has not been initialized

    integer :: mflag
    integer,dimension(2) :: ix, iy, iz, iq
    integer :: i, j, k, l

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag,i)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag,j)
        call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag,k)
        call dintrv(me%q,q,me%iloq,iq(1),iq(2),mflag,l)

        f = me%f(i,j,k,l)
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine nearest_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  5D nearest neighbor interpolation routine.

    pure subroutine nearest_5d(me,x,y,z,q,r,f,istat)

    implicit none

    class(nearest_interp_5d),intent(inout) :: me
    real(wp),intent(in)                    :: x
    real(wp),intent(in)                    :: y
    real(wp),intent(in)                    :: z
    real(wp),intent(in)                    :: q
    real(wp),intent(in)                    :: r
    real(wp),intent(out)                   :: f     !! Nearest \( f(x,y,z,q,r) \)
    integer,intent(out),optional           :: istat !! `0`  : no problems,
                                                    !! `-1` : class has not been initialized

    integer :: mflag
    integer,dimension(2) :: ix, iy, iz, iq, ir
    integer :: i, j, k, l, m

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag,i)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag,j)
        call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag,k)
        call dintrv(me%q,q,me%iloq,iq(1),iq(2),mflag,l)
        call dintrv(me%r,r,me%ilor,ir(1),ir(2),mflag,m)

        f = me%f(i,j,k,l,m)
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine nearest_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  6D nearest neighbor interpolation routine.

    pure subroutine nearest_6d(me,x,y,z,q,r,s,f,istat)

    implicit none

    class(nearest_interp_6d),intent(inout) :: me
    real(wp),intent(in)                    :: x
    real(wp),intent(in)                    :: y
    real(wp),intent(in)                    :: z
    real(wp),intent(in)                    :: q
    real(wp),intent(in)                    :: r
    real(wp),intent(in)                    :: s
    real(wp),intent(out)                   :: f     !! Nearest \( f(x,y,z,q,r,s) \)
    integer,intent(out),optional           :: istat !! `0`  : no problems,
                                                    !! `-1` : class has not been initialized

    integer :: mflag
    integer,dimension(2) :: ix, iy, iz, iq, ir, is
    integer :: i, j, k, l, m, n

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag,i)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag,j)
        call dintrv(me%z,z,me%iloz,iz(1),iz(2),mflag,k)
        call dintrv(me%q,q,me%iloq,iq(1),iq(2),mflag,l)
        call dintrv(me%r,r,me%ilor,ir(1),ir(2),mflag,m)
        call dintrv(me%s,s,me%ilos,is(1),is(2),mflag,n)

        f = me%f(i,j,k,l,m,n)
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine nearest_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the validity of the inputs to the initialize routines.
!  Prints warning message if there is an error,
!  and also sets `ierr` (/=0 if there were any errors).
!
!  Supports up to 6D: x,y,z,q,r,s
!
!# History
!  * Jacob Williams, 2/24/2015 : Created this routine.
!  * Jacob Williams, 2/23/2016 : modified for linear interp module.

    pure subroutine check_inputs(me,x,y,z,q,r,s,ierr)

    implicit none

    class(linear_interp_class),intent(in)      :: me
    real(wp),dimension(:),intent(in),optional  :: x     !! `x` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: y     !! `y` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: z     !! `z` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: q     !! `q` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: r     !! `r` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: s     !! `s` abscissa vector
    integer,intent(out)                        :: ierr  !! `0`   : no problems,
                                                        !! `1`   : `x` is not strictly increasing,
                                                        !! `2`   : `y` is not strictly increasing,
                                                        !! `3`   : `z` is not strictly increasing,
                                                        !! `4`   : `q` is not strictly increasing,
                                                        !! `5`   : `r` is not strictly increasing,
                                                        !! `6`   : `s` is not strictly increasing,
                                                        !! `100` : cannot use linear interpolation for only one point.

    ierr = 0  ! initialize

    if (present(x)) call check(x,1,ierr); if (ierr/=0) return
    if (present(y)) call check(y,2,ierr); if (ierr/=0) return
    if (present(z)) call check(z,3,ierr); if (ierr/=0) return
    if (present(q)) call check(q,4,ierr); if (ierr/=0) return
    if (present(r)) call check(r,5,ierr); if (ierr/=0) return
    if (present(s)) call check(s,6,ierr); if (ierr/=0) return

    if (ierr == 0) then
        select type (me)
        class is (nearest_interp_1d)
        class is (nearest_interp_2d)
        class is (nearest_interp_3d)
        class is (nearest_interp_4d)
        class is (nearest_interp_5d)
        class is (nearest_interp_6d)
        class default
            ! need at least two points for linear interpolation:
            if (size(x)==1) ierr = 100
        end select
    end if

    contains
!*****************************************************************************************

        pure subroutine check(v,error_code,ierr)

        implicit none

        real(wp),dimension(:),intent(in) :: v          !! abcissae vector
        integer,intent(in)               :: error_code !! error code for check
        integer,intent(inout)            :: ierr       !! will be set to `error_code` if there is a problem

        integer :: i  !! counter
        integer :: n  !! size of the input `v` array

        n = size(v)
        do i=2,n
            if (v(i) <= v(i-1)) then
                ierr = error_code
                exit
            end if
        end do

        end subroutine check

    end subroutine check_inputs
!*****************************************************************************************

!*****************************************************************************************
    end module linear_interpolation_module
!*****************************************************************************************
