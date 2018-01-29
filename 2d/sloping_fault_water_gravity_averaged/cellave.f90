subroutine cellave(xlow,zlow,dx,dz,w)
  ! compute the fraction of the cell that is water
  implicit none

  real(kind=8), intent(in) :: xlow, zlow, dx, dz
  real(kind=8), intent(out) :: w

  real(kind=8) :: zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf
  common /topography/ zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf

  integer :: i
  logical :: all_water, all_rock
  real(kind=8) :: xx(4), zz(4), f(4), b, h1, h2
  real(kind=8) :: fdisc

  external fdisc

  xx(1) = xlow
  xx(2) = xlow
  xx(3) = xlow+dx
  xx(4) = xlow+dx
  zz(1) = zlow
  zz(2) = zlow+dz
  zz(3) = zlow+dz
  zz(4) = zlow

  all_water = .true.
  all_rock = .true.
  do i=1,4
    f(i) = fdisc(xx(i),zz(i))
    all_water = all_water .and. (f(i) <= 0.d0)
    all_rock = all_rock .and. (f(i) >= 0.d0)
  end do

  if (all_water) then
    w = 1.d0
    return
  end if
  if (all_rock) then
    w = 0.d0
    return
  end if

  if (f(1) >= 0.d0 .and. f(2) < 0.d0 .and. f(3) < 0.d0 .and. f(4) < 0.d0) then
    ! compute area based on bottom left triangle
    w = 0.5*f(1)*dx/(-f(1)+f(4))*(-f(1))
    w = dx*dz-w
  elseif (f(1) >= 0.d0 .and. f(2) < 0.d0 .and. f(3) >= 0.d0 .and. f(4) >= 0.d0) then
    ! compute area based on top left triangle
    w = 0.5*abs(f(2))*dx/(f(3)-f(2))*(dz-f(1))
  elseif (f(1) >= 0.d0 .and. f(2) < 0.d0 .and. f(3) < 0.d0 .and. f(4) >= 0.d0) then
    ! compute area based on trapezoid (horiztonal split)
    w = 0.5d0*(abs(f(2)) + abs(f(3)))*dx
  elseif (f(1) >= 0.d0 .and. f(2) >= 0.d0 .and. f(3) < 0.d0 .and. f(4) >= 0.d0) then
    ! compute area based on top right triangle
    w = 0.5d0*abs(f(3))*dx*(1.d0 - 1.d0/(-f(2)+f(3))*(dz-f(1)))
  elseif (f(1) < 0.d0 .and. f(2) < 0.d0 .and. f(3) >= 0.d0 .and. f(4) >= 0.d0) then
    ! compute area based on trapezoid (vertical split)
    w = 0.5*(dx/(f(3)-f(1))*(-f(1)) + dx/(f(3)-f(1))*(dz-f(1)))*dz
  elseif (f(1) < 0.d0 .and. f(2) < 0.d0 .and. f(3) < 0.d0 .and. f(4) >= 0.d0) then
    ! compute area based on bottom right triangle
    w = 0.5d0*f(4)*dx*(1.d0 - 1.d0/(f(4)-f(1))*(-f(1)))
    w = dx*dz - w
  else
    print *, 'f', f
    stop
  end if

  w = w/(dx*dz)
  if (w < 0.d0 .or. w > 1.d0) then
    print *, 'w', w
    stop
  end if

end subroutine cellave
