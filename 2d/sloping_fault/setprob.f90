!   ==================
    subroutine setprob
!   ==================

    implicit none

    character*12 fname
    integer iunit

    REAL (kind=8) :: fcenter(2), theta, xcb(2), mindepth
    common /fault/  fcenter, theta, xcb, mindepth

    real (kind=8) :: bacenter, bawidth, baelev, slcenter, slwidth, &
      shecenter, shewidth, sheelev, becenter, bewidth, shoelev
    common /surface/  bacenter, bawidth, baelev, slcenter, slwidth, &
      shecenter, shewidth, sheelev, becenter, bewidth, shoelev

    REAL (kind=8) :: bcenter(2), bwidth, belev
    common /basin/  bcenter, bwidth, belev

    real (kind=8) :: fwidth, xc, yc, xp, yp, ls
    integer :: i
!
!
      iunit = 7
      fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
      call opendatafile(iunit, fname)


!
      read(7,*) fcenter(1)
      read(7,*) fwidth
      read(7,*) theta
      read(7,*) fcenter(2)
      read(7,*) bacenter
      read(7,*) bawidth
      read(7,*) baelev
      read(7,*) slcenter
      read(7,*) slwidth
      read(7,*) shecenter
      read(7,*) shewidth
      read(7,*) sheelev
      read(7,*) becenter
      read(7,*) bewidth
      read(7,*) shoelev

      fcenter(2) = -fcenter(2)
      xcb(1) = fcenter(1) - 0.5*fwidth
      xcb(2) = fcenter(1) + 0.5*fwidth

      ! will need to do better than this once topography is entered via
      ! topo files
      mindepth = 1.0d9

      do i = 1,10000
        xc = xcb(1)-0.5d0*fwidth + (i-1.d0)*2.d0*fwidth/(9999.d0)
        yc = 0.d0
        call mapc2p(xc,yc,xp,yp)
        if (xc < xcb(1)) then
          ls = dsqrt((xc-xcb(1))**2 + (yp-fcenter(2))**2)
        elseif (xc > xcb(2)) then
          ls = dsqrt((xc-xcb(2))**2 + (yp-fcenter(2))**2)
        else
          ls = dabs(yp-fcenter(2))
        end if
        mindepth = dmin1(mindepth,ls)
      end do

      write(6,*) mindepth

    return
    end
