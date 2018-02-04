!   ==================
    subroutine setprob
!   ==================

    use fault_module, only: load_fault

    implicit none

    character*12 fname
    integer iunit

    real(kind=8) :: ABLdepth, ABLxpos(2), ABLypos
    common /ablparam/ ABLdepth, ABLxpos, ABLypos

    real(kind=8) :: lambda_plate, mu_plate, rho_plate, lambda_water, mu_water, rho_water
    common /material/ lambda_plate, mu_plate, rho_plate, lambda_water, mu_water, rho_water

    real (kind=8) :: scaling
    common /water/  scaling

    real (kind=8) :: xlower, xupper, ylower
!
!
    fname = 'fault.data'

    call load_fault(fname)

    iunit = 7
    fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(iunit,*) scaling
    read(iunit,*) ABLdepth
    close(iunit)

!     # Read in grid parameters
    fname = 'claw.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(iunit,*) xlower ! num_dim
    read(iunit,*) xlower, ylower ! lower[0], lower[1]
    read(iunit,*) xupper ! upper[0]
    close(iunit)

!     # Compute ABL position
    ABLxpos(1) = xlower + ABLdepth
    ABLxpos(2) = xupper - ABLdepth
    ABLypos = ylower + ABLdepth

!   # Set material parameters
    lambda_plate = 60.d9  ! Pa
    mu_plate = 30.d9      ! Pa
    rho_plate = 2500.d0   ! kg/m**3

    lambda_water = 2.202256d9  ! Pa
    mu_water = 0.d0      ! Pa
    rho_water = 1000.d0   ! kg/m**3

    return
    end
