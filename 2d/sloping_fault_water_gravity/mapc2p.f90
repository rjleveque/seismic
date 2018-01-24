  !=====================================================
  subroutine mapc2p(xc,zc,xp,zp)
  !=====================================================
    ! Maps for sloping fault
    ! on input,  (xc,zc) is a computational grid point
    ! on output, (xp,zp) is corresponding point in physical space

    use fault_module, only: center, theta, xcb

    implicit none
    real (kind=8), intent(in) :: xc, zc
    real (kind=8), intent(out) :: xp, zp

    ! Variables from setprob:
    real(kind=8) :: zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf
    common /topography/ zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf

    real(kind=8) :: fault_zshift
    common /mapping/ fault_zshift

    ! Local variables
    real (kind=8) :: ls, tol, x_rot, z_rot, x_floor, z_floor, zc_tmp
    real (kind=8) :: slope, floor_scale, floor_shift

    ! compute location in grid scaled or shifted to line up with ocean floor
    if (xc > xlower_shelf) then
      floor_scale = zlower_shelf/zlower_ocean
      floor_shift = zlower_shelf-zlower_ocean
    else if (xc > xlower_slope) then
      slope = (zlower_shelf - zlower_ocean)/(xlower_shelf - xlower_slope)
      floor_scale = (zlower_ocean + (xc-xlower_slope)*slope)/zlower_ocean
      floor_shift = (xc - xlower_slope)*slope
    else
      floor_scale = 1.d0
      floor_shift = 0.d0
    end if
    x_floor = xc
    if (zc > zlower_ocean) then
      z_floor = zc*floor_scale
    elseif (zc < center(2)) then
       z_floor = zc
    else
      z_floor = zc + (zc - center(2))/(zlower_ocean - center(2))*floor_shift
    end if

    ! compute location in grid rotated to line up with fault
    zc_tmp = zc + fault_zshift
    x_rot = center(1) + dcos(theta)*(xc-center(1)) + dsin(theta)*(zc_tmp-center(2))
    z_rot = center(2) - dsin(theta)*(xc-center(1)) + dcos(theta)*(zc_tmp-center(2))

    ! use distance function to interpolate between rotated grid and current grid
    if (xc < xcb(1)) then
      ls = dsqrt((xc-xcb(1))**2 + (zc_tmp-center(2))**2)
    elseif (xc > xcb(2)) then
      ls = dsqrt((xc-xcb(2))**2 + (zc_tmp-center(2))**2)
    else
      ls = dabs(zc_tmp - center(2))
    end if
    tol = -center(2) + zlower_ocean
    ! ensure that there is a neighborhood around the fault that is aligned to
    ! the fault
    ls = ls + zlower_ocean
    tol = tol + zlower_ocean
    if (ls < 0.0) then
      xp = x_rot
      zp = z_rot
    elseif (ls < tol) then
      xp = (tol-ls)/tol*x_rot + ls/tol*x_floor
      zp = (tol-ls)/tol*z_rot + ls/tol*z_floor
    else
      xp = x_floor
      zp = z_floor
    end if

    return
    end
