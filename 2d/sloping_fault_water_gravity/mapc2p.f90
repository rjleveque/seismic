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
    real(kind=8) :: zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf, scale
    common /topography/ zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf, scale

    ! Local variables
    real (kind=8) :: ls, tol, xrot, zrot, factor, slope

    factor = 1.d0
    ! ! compute stretch factor in the z direction for computational grid to match
    ! ! ocean floor
    ! if (xc < xlower_slope) then
    !   factor = 1.d0
    ! else if (xc < xlower_shelf) then
    !   slope = (zlower_shelf - zlower_ocean)/(xlower_shelf - xlower_slope)
    !   factor = (zlower_ocean + (xc - xlower_slope)*slope)/zlower_ocean
    ! end if
    ! compute stretch factor in the z direction for computational grid to match
    ! fault depth and interpolate between that value and existing stretch factor
    if (zc < zlower_ocean) then
      if (zc > center(2)) then
        factor = (zc-center(2))/(zlower_ocean-center(2))*factor + &
                 (zlower_ocean-zc)/(zlower_ocean-center(2))*scale
      else
        factor = scale
      end if
    end if
    xp = xc
    zp = zc*factor

    ! compute seperate physical grid that is the current grid, but rotated to
    ! match the fault
    xrot = center(1) + dcos(theta)*(xp-center(1)) + dsin(theta)*(zp-center(2))
    zrot = center(2) - dsin(theta)*(xp-center(1)) + dcos(theta)*(zp-center(2))

    ! use distance function to interpolate between rotated grid and current grid
    if (xp < xcb(1)) then
      ls = dsqrt((xp-xcb(1))**2 + (zp-center(2))**2)
    elseif (xc > xcb(2)) then
      ls = dsqrt((xp-xcb(2))**2 + (zp-center(2))**2)
    else
      ls = dabs(zp - center(2))
    end if
    tol = -center(2) + zlower_ocean
    ls = ls + zlower_ocean ! this adds a buffer around the fault
    if (ls < 0.0) then
      xp = xrot
      zp = zrot
    elseif (ls < tol) then
      xp = (tol-ls)/tol*xrot + ls/tol*xp
      zp = (tol-ls)/tol*zrot + ls/tol*zp
    end if

    return
    end
