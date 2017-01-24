!   ==================
    subroutine setprob
!   ==================

    use fault_module, only: load_fault

    implicit none

    real(kind=8) :: ABLdepth, ABLxpos(2), ABLypos
    common /ablparam/ ABLdepth, ABLxpos, ABLypos

    character*12 :: fname
    integer :: iunit
    real(kind=8) :: xlower, xupper, ylower

    fname = 'fault.data'

    call load_fault(fname)


!     # Set the material parameters for the acoustic equations

    iunit = 7
    fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
    call opendatafile(iunit, fname)

!     # Read in parameters:

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


    return
    end
