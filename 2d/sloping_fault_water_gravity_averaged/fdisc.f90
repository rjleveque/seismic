function fdisc(x,z)
  implicit none
  real(kind=8), intent(in) :: x, z

!     # for computing cell averages for initial data that has a
!     # discontinuity along some curve.  fdisc should be negative to the
!     # left of the curve and positive to the right
!     # fdisc is positive in the lower layer, negative in upper

  real(kind=8) :: zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf
  common /topography/ zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf

  real(kind=8) :: fdisc, slope_of_slope, zfloor

  ! flat ocean with slope up to flat shelf, piecewise linear:
  slope_of_slope = (zlower_ocean - zlower_shelf) / (xlower_slope - xlower_shelf)
  zfloor = zlower_ocean
  if (x > xlower_slope) zfloor = zlower_ocean + slope_of_slope*(x-xlower_slope)
  if (x > xlower_shelf) zfloor = zlower_shelf

  fdisc = zfloor - z

  return
end function fdisc
