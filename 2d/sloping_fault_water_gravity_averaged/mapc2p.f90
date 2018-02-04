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
    real (kind=8) :: ls, tol, x_rot, z_rot, zc_tmp

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
      xp = (tol-ls)/tol*x_rot + ls/tol*xc
      zp = (tol-ls)/tol*z_rot + ls/tol*zc
    else
      xp = xc
      zp = zc
    end if

    return
    end
