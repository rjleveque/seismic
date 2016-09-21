  !=====================================================
  subroutine mapc2p(xc,yc,xp,yp)
  !=====================================================
    ! Maps for sloping fault
    ! on input,  (xc,yc) is a computational grid point
    ! on output, (xp,yp) is corresponding point in physical space

    implicit none
    real (kind=8), intent(in) :: xc,yc
    real (kind=8), intent(out) :: xp,yp

    ! Variables from setprob:
    real (kind=8) :: fcenter(2), theta, xcb(2), mindepth
    common /fault/  fcenter, theta, xcb, mindepth

    real (kind=8) :: bacenter, bawidth, baelev, slcenter, slwidth, &
      shecenter, shewidth, sheelev, becenter, bewidth, shoelev
    common /surface/  bacenter, bawidth, baelev, slcenter, slwidth, &
      shecenter, shewidth, sheelev, becenter, bewidth, shoelev


    ! Local variables
    real (kind=8) :: ls, alpha, xrot, yrot, xsf, ysf, tol

    tol = mindepth

    ! Grid aligned with surface (xsf,ysf)
    xsf = xc
    if (xc < bacenter + 0.5d0*bawidth) then
      ysf = yc + baelev
    elseif (xc < slcenter + 0.5d0*slwidth) then
      ysf = yc + baelev + ((xc - slcenter)/slwidth + 0.5d0)*(sheelev - baelev)
    elseif (xc < shecenter + 0.5d0*shewidth) then
      ysf = yc + sheelev
    elseif (xc < becenter + 0.5d0*bewidth) then
      ysf = yc + sheelev + ((xc - becenter)/bewidth + 0.5d0)*(shoelev - sheelev)
    else
      ysf = yc + shoelev
    end if

    ! Grid aligned with fault (xrot,yrot)
    if (xc < xcb(1)) then
      ls = dsqrt((xc-xcb(1))**2 + (yc-fcenter(2))**2)
    elseif (xc > xcb(2)) then
      ls = dsqrt((xc-xcb(2))**2 + (yc-fcenter(2))**2)
    else
      ls = dabs(yc - fcenter(2))
    end if

    alpha = ls/tol
    xrot = fcenter(1) + dcos(theta)*(xc-fcenter(1)) + dsin(theta)*(yc-fcenter(2))
    yrot = fcenter(2) - dsin(theta)*(xc-fcenter(1)) + dcos(theta)*(yc-fcenter(2))

    ! Linear interpolation between two grids
    if (alpha < 1.d0) then
      xp = (1.d0-alpha)*xrot + alpha*xsf
      !xp = xc
      yp = (1.d0-alpha)*yrot + alpha*ysf
    else
      xp = xsf
      yp = ysf
    end if

    return
    end
