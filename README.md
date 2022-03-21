![finterp](media/logo.png)
============

FINTERP: Multidimensional Linear Interpolation with Modern Fortran

## Status

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/finterp.svg?style=plastic)](https://github.com/jacobwilliams/finterp/releases/latest)
![Build Status](https://github.com/jacobwilliams/finterp/actions/workflows/CI.yml/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/jacobwilliams/finterp/branch/master/graph/badge.svg?token=BHtd51oUTE)](https://codecov.io/gh/jacobwilliams/finterp)

## Description

Can be used to perform multidimensional (1D-6D) linear interpolation of data on a regular grid. The code is written in modern Fortran (2003/2008) and is object-oriented and thread safe.

## Usage

There are six classes (`linear_interp_1d`, `linear_interp_2d`, `linear_interp_3d`, `linear_interp_4d`, `linear_interp_5d`, and `linear_interp_6d`). Each has three methods: `initialize`, `evaluate`, and `destroy`.

```fortran
real(wp) :: x(nx),y(ny),z(nz),q(nq),r(nr),s(ns)
real(wp) :: fcn_1d(nx)
real(wp) :: fcn_2d(nx,ny)
real(wp) :: fcn_3d(nx,ny,nz)
real(wp) :: fcn_4d(nx,ny,nz,nq)
real(wp) :: fcn_5d(nx,ny,nz,nq,nr)
real(wp) :: fcn_6d(nx,ny,nz,nq,nr,ns)
real(wp) :: xval,yval,zval,qval,rval,sval,fval
integer :: iflag

type(linear_interp_1d) :: s1
type(linear_interp_2d) :: s2
type(linear_interp_3d) :: s3
type(linear_interp_4d) :: s4
type(linear_interp_5d) :: s5
type(linear_interp_6d) :: s6

!populate the arrays
!...

!initialize the class:
call s1%initialize(x,fcn_1d,iflag)
call s2%initialize(x,y,fcn_2d,iflag)
call s3%initialize(x,y,z,fcn_3d,iflag)
call s4%initialize(x,y,z,q,fcn_4d,iflag)
call s5%initialize(x,y,z,q,r,fcn_5d,iflag)
call s6%initialize(x,y,z,q,r,s,fcn_6d,iflag)

!interpolate:
call s1%evaluate(xval,fval)
call s2%evaluate(xval,yval,fval)
call s3%evaluate(xval,yval,zval,fval)
call s4%evaluate(xval,yval,zval,qval,fval)
call s5%evaluate(xval,yval,zval,qval,rval,fval)
call s6%evaluate(xval,yval,zval,qval,rval,sval,fval)

!free memory:
call s1%destroy()
call s2%destroy()
call s3%destroy()
call s4%destroy()
call s5%destroy()
call s6%destroy()
```

## Nearest Neighbor Interpolation

The library also includes classes for nearest neighbor interpolation (`nearest_interp_1d`, `nearest_interp_2d`, ...). The interfaces are the same as for the linear classes.

## Compiling

A `fmp.toml` file is provided for compiling finterp with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

```
fpm build --profile release
```

By default, the library is built with double precision (`real64`) real values. Explicitly specifying the real kind can be done using the following preprocessor flags:

Preprocessor flag | Kind  | Number of bytes
----------------- | ----- | ---------------
`REAL32`  | `real(kind=real32)`  | 4
`REAL64`  | `real(kind=real64)`  | 8
`REAL128` | `real(kind=real128)` | 16

For example, to build a single precision version of the library, use:

```
fpm build --profile release --flag "-DREAL32"
```

To run the unit tests:

```
fpm test
```

To use `finterp` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
finterp = { git="https://github.com/jacobwilliams/finterp.git" }
```

or, to use a specific version:
```toml
[dependencies]
finterp = { git="https://github.com/jacobwilliams/finterp.git", tag = "1.3.0"  }
```

To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```ford finterp.md```

## Documentation

The latest API documentation can be found [here](https://jacobwilliams.github.io/finterp/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

## License

The finterp source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/finterp/blob/master/LICENSE) (BSD-style).

## See also

 * [bspline-fortran](https://github.com/jacobwilliams/bspline-fortran), B-spline interpolation.
 * [PCHIP](https://github.com/jacobwilliams/PCHIP), Piecewise Cubic Hermite Interpolation.
 * [Regridpack](https://github.com/jacobwilliams/regridpack), Linear or cubic interpolation for 1D-4D grids.