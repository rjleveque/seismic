!   ==================
    subroutine setprob
!   ==================

    use fault_module, only: load_fault, center

    implicit none

    character*12 fname
    integer iunit

    real(kind=8) :: ABLdepth, ABLxpos(2), ABLypos
    common /ablparam/ ABLdepth, ABLxpos, ABLypos

    real(kind=8) :: lambda_plate, mu_plate, rho_plate, lambda_water, mu_water, rho_water, g
    common /material/ lambda_plate, mu_plate, rho_plate, lambda_water, mu_water, rho_water, g

    real(kind=8) :: zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf
    common /topography/ zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf

    real(kind=8) :: fault_zshift
    common /mapping/ fault_zshift

    integer :: nx, ny
    real(kind=8) :: xlower, xupper, ylower, yupper, dy
!

    fname = 'fault.data'

    call load_fault(fname)
!
    iunit = 7
    fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(iunit,*) zlower_ocean
    read(iunit,*) xlower_slope
    read(iunit,*) xlower_shelf
    read(iunit,*) zlower_shelf
    read(iunit,*) ABLdepth ! xlower_beach
    read(iunit,*) ABLdepth ! xlower_shore
    read(iunit,*) ABLdepth ! zlower_shore
    read(iunit,*) ABLdepth
    close(iunit)


!     # Read in grid parameters
    fname = 'claw.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(iunit,*) xlower ! num_dim
    read(iunit,*) xlower, ylower ! lower[0], lower[1]
    read(iunit,*) xupper, yupper ! upper[0], upper[1]
    read(iunit,*) nx, ny ! nx, ny
    close(iunit)

!     # Compute ABL position
    ABLxpos(1) = xlower + ABLdepth
    ABLxpos(2) = xupper - ABLdepth
    ABLypos = ylower + ABLdepth

!    # Compute vertical shift for computational grid to line up with fault depth
    dy = (yupper-ylower)/ny
    fault_zshift = center(2) + ceiling(-center(2)/dy)*dy

!   # Set material parameters
    lambda_plate = 60.d9  ! Pa
    mu_plate = 30.d9      ! Pa
    rho_plate = 2500.d0   ! kg/m**3

    lambda_water = 2.202256d9  ! Pa
    mu_water = 0.d0      ! Pa
    rho_water = 1000.d0   ! kg/m**3

    g = 9.8 ! m/s**2

    return
    end
