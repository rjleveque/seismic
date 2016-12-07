subroutine setaux1d(ixyz,mbc,mx1,x1lower,x2val,x3val,dx1,dx2,dx3,t,maux,aux,iout)

    ! set auxiliary variables in current row with direction specified by ixyz
    !
    ! ixyz = 1 => row in x direction, e2 is y direction, e3 is z direction
    ! ixyz = 2 => row in y direction, e2 is z direction, e3 is x direction
    ! ixyz = 3 => row in z direction, e2 is x direction, e3 is y direction
    !
    ! Auxiliary variables:
    !       1 rho
    !       2 lambda
    !       3 mu
    !       4 cp
    !       5 cs
    !       6 slip
    !       7 nx at lower wall in e1 direction
    !       8 ny at lower wall in e1 direction
    !       9 nz at lower wall in e1 direction
    !       10 area ratio of lower wall in e1 direction
    !       11 nx at lower wall in e2 direction
    !       12 ny at lower wall in e2 direction
    !       13 nz at lower wall in e2 direction
    !       14 area ratio of lower wall in e2 direction
    !       15 nx at lower wall in e3 direction
    !       16 ny at lower wall in e3 direction
    !       17 nz at lower wall in e3 direction
    !       18 area ratio of lower wall in e3 direction

    implicit none
    integer, intent(in) :: ixyz,mbc,mx1,maux,iout
    real(kind=8), intent(in) :: x1lower,x2val,x3val,dx1,dx2,dx3,t
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx1+mbc,3)

    real(kind=8) :: x1cell, cp, cs, rho_cell, lambda_cell, mu_cell
    real(kind=8) :: x1ccorn(4), x2ccorn(4), x3ccorn(4)
    real(kind=8) :: x1pcorn(4), x2pcorn(4), x3pcorn(4), mag
    real(kind=8) :: nx1, nx2, nx3
    integer :: i

    real (kind=8) :: center(3), theta, xcb(2), ycb(2), mindepth
    common /fault/  center, theta, xcb, ycb, mindepth

    lambda_cell = 60.d9  ! Pa
    mu_cell = 30.d9      ! Pa
    rho_cell = 2500.d0   ! kg/m**3
    cp = dsqrt((lambda_cell + 2.d0*mu_cell)/rho_cell)
    cs = dsqrt(mu_cell/rho_cell)

    do i=1-mbc,mx1 + mbc
      x1cell = x1lower + (i-0.5d0)*dx1
      aux(1,i,iout) = rho_cell
      aux(2,i,iout) = lambda_cell
      aux(3,i,iout) = mu_cell
      aux(4,i,iout) = cp
      aux(5,i,iout) = cs

      ! set slip
      if ((ixyz == 3) .and. (t < 1.d0) .and. &
          (abs(x1cell-0.5d0*dx1-center(3)) < 0.5d0*dx1) .and. &
          (xcb(1) <= x2val) .and. (x2val <= xcb(2)) .and. &
          (ycb(1) <= x3val) .and. (x3val <= ycb(2))) then

        aux(6,i,iout) = exp(-( &
                           ((x2val-center(1))/(center(1)-xcb(1)))**2 &
                           + ((x3val-center(2))/(center(2)-ycb(1)))**2 &
                           ))
      else
        aux(6,i,iout) = 0.d0
      end if

      ! compute mapping info for lower face in e1 direction
      x1ccorn(1) = x1cell - 0.5d0*dx1 ! p1 dot e1
      x2ccorn(1) = x2val - 0.5d0*dx2  ! p1 dot e2
      x3ccorn(1) = x3val - 0.5d0*dx3  ! p1 dot e3

      x1ccorn(2) = x1cell - 0.5d0*dx1 ! p2 dot e1
      x2ccorn(2) = x2val - 0.5d0*dx2  ! p2 dot e2
      x3ccorn(2) = x3val + 0.5d0*dx3  ! p2 dot e3

      x1ccorn(3) = x1cell - 0.5d0*dx1 ! p3 dot e1
      x2ccorn(3) = x2val + 0.5d0*dx2  ! p3 dot e2
      x3ccorn(3) = x3val + 0.5d0*dx3  ! p3 dot e3

      x1ccorn(4) = x1cell - 0.5d0*dx1 ! p4 dot e1
      x2ccorn(4) = x2val + 0.5d0*dx2  ! p4 dot e2
      x3ccorn(4) = x3val - 0.5d0*dx3  ! p4 dot e3
      if (ixyz == 1) then
        call mapc2p(x1ccorn(1), x2ccorn(1), x3ccorn(1), x1pcorn(1), x2pcorn(1), x3pcorn(1))
        call mapc2p(x1ccorn(2), x2ccorn(2), x3ccorn(2), x1pcorn(2), x2pcorn(2), x3pcorn(2))
        call mapc2p(x1ccorn(3), x2ccorn(3), x3ccorn(3), x1pcorn(3), x2pcorn(3), x3pcorn(3))
        call mapc2p(x1ccorn(4), x2ccorn(4), x3ccorn(4), x1pcorn(4), x2pcorn(4), x3pcorn(4))
      else if (ixyz == 2) then
        call mapc2p(x3ccorn(1), x1ccorn(1), x2ccorn(1), x3pcorn(1), x1pcorn(1), x2pcorn(1))
        call mapc2p(x3ccorn(2), x1ccorn(2), x2ccorn(2), x3pcorn(2), x1pcorn(2), x2pcorn(2))
        call mapc2p(x3ccorn(3), x1ccorn(3), x2ccorn(3), x3pcorn(3), x1pcorn(3), x2pcorn(3))
        call mapc2p(x3ccorn(4), x1ccorn(4), x2ccorn(4), x3pcorn(4), x1pcorn(4), x2pcorn(4))
      else if (ixyz == 3) then
        call mapc2p(x2ccorn(1), x3ccorn(1), x1ccorn(1), x2pcorn(1), x3pcorn(1), x1pcorn(1))
        call mapc2p(x2ccorn(2), x3ccorn(2), x1ccorn(2), x2pcorn(2), x3pcorn(2), x1pcorn(2))
        call mapc2p(x2ccorn(3), x3ccorn(3), x1ccorn(3), x2pcorn(3), x3pcorn(3), x1pcorn(3))
        call mapc2p(x2ccorn(4), x3ccorn(4), x1ccorn(4), x2pcorn(4), x3pcorn(4), x1pcorn(4))
      end if

      nx1 = (x2pcorn(3) - x2pcorn(1))*(x3pcorn(2) - x3pcorn(4)) - (x2pcorn(2) - x2pcorn(4))*(x3pcorn(3) - x3pcorn(1))
      nx2 = (x1pcorn(2) - x1pcorn(4))*(x3pcorn(3) - x3pcorn(1)) - (x1pcorn(3) - x1pcorn(1))*(x3pcorn(2) - x3pcorn(4))
      nx3 = (x1pcorn(3) - x1pcorn(1))*(x2pcorn(2) - x2pcorn(4)) - (x1pcorn(2) - x1pcorn(4))*(x2pcorn(3) - x2pcorn(1))
      mag = dsqrt(nx1*nx1 + nx2*nx2 + nx3*nx3)
      nx1 = nx1/mag
      nx2 = nx2/mag
      nx3 = nx3/mag

      if (ixyz == 1) then
        aux(7,i,iout) = nx1
        aux(8,i,iout) = nx2
        aux(9,i,iout) = nx3
      else if (ixyz == 2) then
        aux(7,i,iout) = nx3
        aux(8,i,iout) = nx1
        aux(9,i,iout) = nx2
      else if (ixyz == 3) then
        aux(7,i,iout) = nx2
        aux(8,i,iout) = nx3
        aux(9,i,iout) = nx1
      end if

      aux(10,i,iout) = 0.5d0*mag/(dx2*dx3)

      ! compute mapping info for lower face in e2 direction
      x1ccorn(1) = x1cell - 0.5d0*dx1 ! p1 dot e1
      x2ccorn(1) = x2val - 0.5d0*dx2  ! p1 dot e2
      x3ccorn(1) = x3val - 0.5d0*dx3  ! p1 dot e3

      x1ccorn(2) = x1cell + 0.5d0*dx1 ! p2 dot e1
      x2ccorn(2) = x2val - 0.5d0*dx2  ! p2 dot e2
      x3ccorn(2) = x3val - 0.5d0*dx3  ! p2 dot e3

      x1ccorn(3) = x1cell + 0.5d0*dx1 ! p3 dot e1
      x2ccorn(3) = x2val - 0.5d0*dx2  ! p3 dot e2
      x3ccorn(3) = x3val + 0.5d0*dx3  ! p3 dot e3

      x1ccorn(4) = x1cell - 0.5d0*dx1 ! p4 dot e1
      x2ccorn(4) = x2val - 0.5d0*dx2  ! p4 dot e2
      x3ccorn(4) = x3val + 0.5d0*dx3  ! p4 dot e3
      if (ixyz == 1) then
        call mapc2p(x1ccorn(1), x2ccorn(1), x3ccorn(1), x1pcorn(1), x2pcorn(1), x3pcorn(1))
        call mapc2p(x1ccorn(2), x2ccorn(2), x3ccorn(2), x1pcorn(2), x2pcorn(2), x3pcorn(2))
        call mapc2p(x1ccorn(3), x2ccorn(3), x3ccorn(3), x1pcorn(3), x2pcorn(3), x3pcorn(3))
        call mapc2p(x1ccorn(4), x2ccorn(4), x3ccorn(4), x1pcorn(4), x2pcorn(4), x3pcorn(4))
      else if (ixyz == 2) then
        call mapc2p(x3ccorn(1), x1ccorn(1), x2ccorn(1), x3pcorn(1), x1pcorn(1), x2pcorn(1))
        call mapc2p(x3ccorn(2), x1ccorn(2), x2ccorn(2), x3pcorn(2), x1pcorn(2), x2pcorn(2))
        call mapc2p(x3ccorn(3), x1ccorn(3), x2ccorn(3), x3pcorn(3), x1pcorn(3), x2pcorn(3))
        call mapc2p(x3ccorn(4), x1ccorn(4), x2ccorn(4), x3pcorn(4), x1pcorn(4), x2pcorn(4))
      else if (ixyz == 3) then
        call mapc2p(x2ccorn(1), x3ccorn(1), x1ccorn(1), x2pcorn(1), x3pcorn(1), x1pcorn(1))
        call mapc2p(x2ccorn(2), x3ccorn(2), x1ccorn(2), x2pcorn(2), x3pcorn(2), x1pcorn(2))
        call mapc2p(x2ccorn(3), x3ccorn(3), x1ccorn(3), x2pcorn(3), x3pcorn(3), x1pcorn(3))
        call mapc2p(x2ccorn(4), x3ccorn(4), x1ccorn(4), x2pcorn(4), x3pcorn(4), x1pcorn(4))
      end if

      nx1 = (x2pcorn(3) - x2pcorn(1))*(x3pcorn(2) - x3pcorn(4)) - (x2pcorn(2) - x2pcorn(4))*(x3pcorn(3) - x3pcorn(1))
      nx2 = (x1pcorn(2) - x1pcorn(4))*(x3pcorn(3) - x3pcorn(1)) - (x1pcorn(3) - x1pcorn(1))*(x3pcorn(2) - x3pcorn(4))
      nx3 = (x1pcorn(3) - x1pcorn(1))*(x2pcorn(2) - x2pcorn(4)) - (x1pcorn(2) - x1pcorn(4))*(x2pcorn(3) - x2pcorn(1))
      mag = dsqrt(nx1*nx1 + nx2*nx2 + nx3*nx3)
      nx1 = nx1/mag
      nx2 = nx2/mag
      nx3 = nx3/mag

      if (ixyz == 1) then
        aux(11,i,iout) = nx1
        aux(12,i,iout) = nx2
        aux(13,i,iout) = nx3
      else if (ixyz == 2) then
        aux(11,i,iout) = nx3
        aux(12,i,iout) = nx1
        aux(13,i,iout) = nx2
      else if (ixyz == 3) then
        aux(11,i,iout) = nx2
        aux(12,i,iout) = nx3
        aux(13,i,iout) = nx1
      end if

      aux(14,i,iout) = 0.5d0*mag/(dx1*dx3)

      ! compute mapping info for lower face in e3 direction
      x1ccorn(1) = x1cell - 0.5d0*dx1 ! p1 dot e1
      x2ccorn(1) = x2val - 0.5d0*dx2  ! p1 dot e2
      x3ccorn(1) = x3val - 0.5d0*dx3  ! p1 dot e3

      x1ccorn(2) = x1cell - 0.5d0*dx1 ! p2 dot e1
      x2ccorn(2) = x2val + 0.5d0*dx2  ! p2 dot e2
      x3ccorn(2) = x3val - 0.5d0*dx3  ! p2 dot e3

      x1ccorn(3) = x1cell + 0.5d0*dx1 ! p3 dot e1
      x2ccorn(3) = x2val + 0.5d0*dx2  ! p3 dot e2
      x3ccorn(3) = x3val - 0.5d0*dx3  ! p3 dot e3

      x1ccorn(4) = x1cell + 0.5d0*dx1 ! p4 dot e1
      x2ccorn(4) = x2val - 0.5d0*dx2  ! p4 dot e2
      x3ccorn(4) = x3val - 0.5d0*dx3  ! p4 dot e3
      if (ixyz == 1) then
        call mapc2p(x1ccorn(1), x2ccorn(1), x3ccorn(1), x1pcorn(1), x2pcorn(1), x3pcorn(1))
        call mapc2p(x1ccorn(2), x2ccorn(2), x3ccorn(2), x1pcorn(2), x2pcorn(2), x3pcorn(2))
        call mapc2p(x1ccorn(3), x2ccorn(3), x3ccorn(3), x1pcorn(3), x2pcorn(3), x3pcorn(3))
        call mapc2p(x1ccorn(4), x2ccorn(4), x3ccorn(4), x1pcorn(4), x2pcorn(4), x3pcorn(4))
      else if (ixyz == 2) then
        call mapc2p(x3ccorn(1), x1ccorn(1), x2ccorn(1), x3pcorn(1), x1pcorn(1), x2pcorn(1))
        call mapc2p(x3ccorn(2), x1ccorn(2), x2ccorn(2), x3pcorn(2), x1pcorn(2), x2pcorn(2))
        call mapc2p(x3ccorn(3), x1ccorn(3), x2ccorn(3), x3pcorn(3), x1pcorn(3), x2pcorn(3))
        call mapc2p(x3ccorn(4), x1ccorn(4), x2ccorn(4), x3pcorn(4), x1pcorn(4), x2pcorn(4))
      else if (ixyz == 3) then
        call mapc2p(x2ccorn(1), x3ccorn(1), x1ccorn(1), x2pcorn(1), x3pcorn(1), x1pcorn(1))
        call mapc2p(x2ccorn(2), x3ccorn(2), x1ccorn(2), x2pcorn(2), x3pcorn(2), x1pcorn(2))
        call mapc2p(x2ccorn(3), x3ccorn(3), x1ccorn(3), x2pcorn(3), x3pcorn(3), x1pcorn(3))
        call mapc2p(x2ccorn(4), x3ccorn(4), x1ccorn(4), x2pcorn(4), x3pcorn(4), x1pcorn(4))
      end if

      nx1 = (x2pcorn(3) - x2pcorn(1))*(x3pcorn(2) - x3pcorn(4)) - (x2pcorn(2) - x2pcorn(4))*(x3pcorn(3) - x3pcorn(1))
      nx2 = (x1pcorn(2) - x1pcorn(4))*(x3pcorn(3) - x3pcorn(1)) - (x1pcorn(3) - x1pcorn(1))*(x3pcorn(2) - x3pcorn(4))
      nx3 = (x1pcorn(3) - x1pcorn(1))*(x2pcorn(2) - x2pcorn(4)) - (x1pcorn(2) - x1pcorn(4))*(x2pcorn(3) - x2pcorn(1))
      mag = dsqrt(nx1*nx1 + nx2*nx2 + nx3*nx3)
      nx1 = nx1/mag
      nx2 = nx2/mag
      nx3 = nx3/mag

      if (ixyz == 1) then
        aux(15,i,iout) = nx1
        aux(16,i,iout) = nx2
        aux(17,i,iout) = nx3
      else if (ixyz == 2) then
        aux(15,i,iout) = nx3
        aux(16,i,iout) = nx1
        aux(17,i,iout) = nx2
      else if (ixyz == 3) then
        aux(15,i,iout) = nx2
        aux(16,i,iout) = nx3
        aux(17,i,iout) = nx1
      end if

      aux(18,i,iout) = 0.5d0*mag/(dx1*dx2)

    end do

end subroutine setaux1d
