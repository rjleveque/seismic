  !=====================================================
  subroutine mapc2p(xc,yc,xp,yp)
  !=====================================================
    ! Maps for sloping sea floor
    ! on input,  (xc,yc) is a computational grid point
    ! on output, (xp,yp) is corresponding point in physical space

    implicit none
    real (kind=8), intent(in) :: xc,yc
    real (kind=8), intent(out) :: xp,yp

    ! Local variables
    real (kind=8) :: xp1,yp0,ypi1,ypi2,spi,yci,ypi,s1,s2

    xp1 = 40d3
    yp0 = -40d3
    ypi1 = -10d3
    ypi2 = -4d3
    yci = 0.8d0
    spi = (ypi2 - ypi1)/(2.d0*xp1)
    xp = xp1 * xc
    ypi = ypi1 + spi*(xp + xp1)
    s1 = (ypi - yp0)/yci
    s2 = -ypi/(1.d0-yci)
    if (yc < yci) then
        yp = yp0 + s1*yc
      else
        yp = -s2*(1.d0 - yc)
      endif

    return
    end

