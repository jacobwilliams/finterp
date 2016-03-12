!*****************************************************************************************
!
!> Units test for 2d-6d linear interpolation.

    program linear_interpolation_test

    use linear_interpolation_module
    use,intrinsic :: iso_fortran_env, only: wp => real64

    implicit none

    integer,parameter :: nx = 6     !! number of points in x
    integer,parameter :: ny = 6     !! number of points in y
    integer,parameter :: nz = 6     !! number of points in z
    integer,parameter :: nq = 6     !! number of points in q
    integer,parameter :: nr = 6     !! number of points in r
    integer,parameter :: ns = 6     !! number of points in s

    real(wp) :: x(nx),y(ny),z(nz),q(nq),r(nr),s(ns)
    real(wp) :: fcn_1d(nx)
    real(wp) :: fcn_2d(nx,ny)
    real(wp) :: fcn_3d(nx,ny,nz)
    real(wp) :: fcn_4d(nx,ny,nz,nq)
    real(wp) :: fcn_5d(nx,ny,nz,nq,nr)
    real(wp) :: fcn_6d(nx,ny,nz,nq,nr,ns)

    type(linear_interp_1d) :: s1
    type(linear_interp_2d) :: s2
    type(linear_interp_3d) :: s3
    type(linear_interp_4d) :: s4
    type(linear_interp_5d) :: s5
    type(linear_interp_6d) :: s6

    real(wp) :: tol,rnd
    real(wp),dimension(6) :: val,tru,err,errmax
    logical :: fail
    integer :: i,j,k,l,m,n,idx,idy,idz,idq,idr,ids
    integer,dimension(6) :: iflag

    fail = .false.
    tol = 1.0e-14_wp
    do i=1,nx
        x(i) = dble(i-1)/dble(nx-1)
    end do
    do j=1,ny
        y(j) = dble(j-1)/dble(ny-1)
    end do
    do k=1,nz
        z(k) = dble(k-1)/dble(nz-1)
    end do
    do l=1,nq
        q(l) = dble(l-1)/dble(nq-1)
    end do
    do m=1,nr
        r(m) = dble(m-1)/dble(nr-1)
    end do
    do n=1,ns
        s(n) = dble(n-1)/dble(ns-1)
    end do
    do i=1,nx
                        fcn_1d(i) = f1(x(i))
        do j=1,ny
                        fcn_2d(i,j) = f2(x(i),y(j))
           do k=1,nz
                        fcn_3d(i,j,k) = f3(x(i),y(j),z(k))
              do l=1,nq
                        fcn_4d(i,j,k,l) = f4(x(i),y(j),z(k),q(l))
                 do m=1,nr
                        fcn_5d(i,j,k,l,m) = f5(x(i),y(j),z(k),q(l),r(m))
                     do n=1,ns
                        fcn_6d(i,j,k,l,m,n) = f6(x(i),y(j),z(k),q(l),r(m),s(n))
                     end do
                 end do
              end do
           end do
        end do
    end do

    ! initialize
    call s1%initialize(x,fcn_1d,            iflag(1))
    call s2%initialize(x,y,fcn_2d,          iflag(2))
    call s3%initialize(x,y,z,fcn_3d,        iflag(3))
    call s4%initialize(x,y,z,q,fcn_4d,      iflag(4))
    call s5%initialize(x,y,z,q,r,fcn_5d,    iflag(5))
    call s6%initialize(x,y,z,q,r,s,fcn_6d,  iflag(6))

    if (any(iflag/=0)) then
        do i=1,6
            if (iflag(i)/=0) then
                write(*,*) 'Error initializing ',i,'D interpolator. iflag=',iflag(i)
            end if
        end do
        stop 1
    end if

    ! compute max error at interpolation points

     errmax = 0.0_wp
     do i=1,nx
                        call RANDOM_NUMBER(rnd); rnd = (rnd - 0.5_wp) / 100.0_wp
                        call s1%evaluate(x(i)+rnd,val(1))
                        tru(1)    = f1(x(i)+rnd)
                        err(1)    = abs(tru(1)-val(1))
                        errmax(1) = max(err(1),errmax(1))
        do j=1,ny
                        call RANDOM_NUMBER(rnd); rnd = rnd / 100.0_wp
                        call s2%evaluate(x(i)+rnd,y(j)+rnd,val(2))
                        tru(2)    = f2(x(i)+rnd,y(j)+rnd)
                        err(2)    = abs(tru(2)-val(2))
                        errmax(2) = max(err(2),errmax(2))
           do k=1,nz
                        call RANDOM_NUMBER(rnd); rnd = rnd / 100.0_wp
                        call s3%evaluate(x(i)+rnd,y(j)+rnd,z(k)+rnd,val(3))
                        tru(3)    = f3(x(i)+rnd,y(j)+rnd,z(k)+rnd)
                        err(3)    = abs(tru(3)-val(3))
                        errmax(3) = max(err(3),errmax(3))
              do l=1,nq
                        call RANDOM_NUMBER(rnd); rnd = rnd / 100.0_wp
                        call s4%evaluate(x(i)+rnd,y(j)+rnd,z(k)+rnd,q(l)+rnd,val(4))
                        tru(4)    = f4(x(i)+rnd,y(j)+rnd,z(k)+rnd,q(l)+rnd)
                        err(4)    = abs(tru(4)-val(4))
                        errmax(4) = max(err(4),errmax(4))
                do m=1,nr
                        call RANDOM_NUMBER(rnd); rnd = rnd / 100.0_wp
                        call s5%evaluate(x(i)+rnd,y(j)+rnd,z(k)+rnd,q(l)+rnd,r(m)+rnd,val(5))
                        tru(5)    = f5(x(i)+rnd,y(j)+rnd,z(k)+rnd,q(l)+rnd,r(m)+rnd)
                        err(5)    = abs(tru(5)-val(5))
                        errmax(5) = max(err(5),errmax(5))
                    do n=1,ns
                        call RANDOM_NUMBER(rnd); rnd = rnd / 100.0_wp
                        call s6%evaluate(x(i)+rnd,y(j)+rnd,z(k)+rnd,q(l)+rnd,r(m)+rnd,s(n)+rnd,val(6))
                        tru(6)    = f6(x(i)+rnd,y(j)+rnd,z(k)+rnd,q(l)+rnd,r(m)+rnd,s(n)+rnd)
                        err(6)    = abs(tru(6)-val(6))
                        errmax(6) = max(err(6),errmax(6))
                    end do
                end do
              end do
           end do
        end do
     end do

    ! check max error against tolerance
    do i=1,6
        write(*,*) i,'D: max error:', errmax(i)
        if (errmax(i) >= tol) then
            write(*,*)  ' ** test failed ** '
        else
            write(*,*)  ' ** test passed ** '
        end if
        write(*,*) ''
    end do

    contains

        real(wp) function f1(x) !! 1d test function
        implicit none
        real(wp),intent(in) :: x
        f1 = x + x + x + x + x + x
        end function f1

        real(wp) function f2(x,y) !! 2d test function
        implicit none
        real(wp),intent(in) :: x,y
        f2 = x + y + x + y + x + y
        end function f2

        real(wp) function f3 (x,y,z) !! 3d test function
        implicit none
        real(wp),intent(in) :: x,y,z
        f3 = x + y + z + x + y + z
        end function f3

        real(wp) function f4 (x,y,z,q) !! 4d test function
        implicit none
        real(wp),intent(in) ::  x,y,z,q
        f4 = x + y + z + q + x + y
        end function f4

        real(wp) function f5 (x,y,z,q,r) !! 5d test function
        implicit none
        real(wp),intent(in) ::  x,y,z,q,r
        f5 = x + y + z + q + r + x
        end function f5

        real(wp) function f6 (x,y,z,q,r,s) !! 6d test function
        implicit none
        real(wp),intent(in) ::  x,y,z,q,r,s
        f6 = x + y + z + q + r + s
        end function f6

    end program linear_interpolation_test
