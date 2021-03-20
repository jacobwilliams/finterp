!*****************************************************************************************
!
!> Units test for 1d nearest neighbor interpolation.

    program nearest_neighbor_test

    use linear_interpolation_module
    use,intrinsic :: iso_fortran_env, only: wp => real64

    implicit none

    integer,parameter :: nx = 6  !! number of points in x
    real(wp),dimension(nx),parameter :: x         = [1.0_wp,2.0_wp,3.0_wp,4.0_wp,5.0_wp,6.0_wp]
    real(wp),dimension(nx),parameter :: fcn_1d    = [1.0_wp,2.0_wp,3.0_wp,4.0_wp,5.0_wp,6.0_wp]
    real(wp),dimension(nx),parameter :: x_to_test = [0.1_wp,2.1_wp,3.4_wp,4.0_wp,4.9_wp,9.0_wp]
    real(wp),parameter :: tol = 1.0e-14_wp  !! error tolerance

    type(nearest_interp_1d) :: s1
    real(wp) :: val,tru,err,errmax
    integer :: i,idx,iflag

    ! initialize
    call s1%initialize(x,fcn_1d,iflag)

    if (iflag/=0) then
        write(*,*) 'Error initializing 1D nearest neighbor interpolator. iflag=',iflag
        stop 1
    end if

    ! compute max error at interpolation points
     errmax = 0.0_wp
     do i=1,nx
        call s1%evaluate(x_to_test(i),val)
        tru     = x(i)
        write(*,*) val, tru
        err     = abs(tru-val)
        errmax  = max(err,errmax)
     end do

    ! check max error against tolerance
    write(*,*) '1D: max error:', errmax
    if (errmax >= tol) then
        write(*,*)  ' ** test failed ** '
    else
        write(*,*)  ' ** test passed ** '
    end if
    write(*,*) ''

    end program nearest_neighbor_test