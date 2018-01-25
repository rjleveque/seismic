!   ==================
    subroutine setprob
!   ==================

    use fault_module, only: load_fault, center

    implicit none

    character*12 fname
    integer iunit

    real(kind=8) :: rho1,amu1,alam1,rho2,amu2,alam2

      common /comaux/ rho1,amu1,alam1,rho2,amu2,alam2

!
!     # fault parameters:
!
      fname = 'fault.data'
      call load_fault(fname)

!     # Set the material parameters for the elasticity equations
!
!
      iunit = 7
      fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
      call opendatafile(iunit, fname)

!
!     # Piecewise constant medium 
!     # Material parameters

      read(7,*) rho1
      read(7,*) alam1
      read(7,*) amu1

      read(7,*) rho2
      read(7,*) alam2
      read(7,*) amu2

      return

    end subroutine setprob
