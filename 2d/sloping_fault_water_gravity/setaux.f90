subroutine setaux(mbc,mx,mz,xlower,zlower,dx,dz,maux,aux)

    implicit none
    integer, intent(in) :: mbc,mx,mz,maux
    real(kind=8), intent(in) :: xlower,zlower,dx,dz
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc,1-mbc:mz+mbc)
    real(kind=8) :: xcell, zcell, cp, cs
    real(kind=8) :: rho_cell, lambda_cell, mu_cell
    integer :: i,j

    ! Arrays to temporarily store computational and physical corners of grid cells
    real(kind=8) :: xccorn(4),zccorn(4),xpcorn(4),zpcorn(4)
    real(kind=8) :: norm, xn, yn, areap, pi2

    real(kind=8) :: ABLdepth, ABLxpos(2), ABLypos
    common /ablparam/ ABLdepth, ABLxpos, ABLypos

    real(kind=8) :: lambda_plate, mu_plate, rho_plate, lambda_water, mu_water, rho_water, g
    common /material/ lambda_plate, mu_plate, rho_plate, lambda_water, mu_water, rho_water, g

    real(kind=8) :: zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf, scale
    common /topography/ zlower_ocean, xlower_slope, xlower_shelf, zlower_shelf, scale



! c     #   (lambda = nu*E/((1+nu)(1-2nu))), E=young modulus, nu=poisson ration
! c     #   aux(1,i,j) is the density of the elastic material
! c     #   aux(2,i,j) is the lambda Lame parameter for elasticity
! c     #   aux(3,i,j) is the mu Lame parameter for elasticity (shear)
!
! c     #   aux(4,i,j) is cp is the P-wave speed (pressure wave)
! c     #   aux(5,i,j) is cs is the S-wave speed (shear wave)
!
! c     #   aux(6,i,j) is x component of normal at "left" boundary of grid point (i,j)
! c     #   aux(7,i,j) is y component of normal at "left" boundary of grid point (i,j)
! c     #   aux(8,i,j) =  is ratio of left edge to dz
!
! c     #   aux(9,i,j) is x component of normal at "bottom" boundary of grid point (i,j)
! c     #   aux(10,i,j) is y component of normal at "bottom" boundary of grid point (i,j)
! c     #   aux(11,i,j) =  is ratio of bottom edge to dx
!
! c     #   aux(12,i,j) = kappa  is ratio of cell area to dx*dz
!
! c     #   aux(13,i,j) = slip:
    !
    pi2 = 2.d0*datan(1.d0)

    ! Loop over all cells
      do j=1-mbc,mz + mbc
        zcell = zlower + (j-0.5d0)*dz
        do i=1-mbc,mx + mbc

          if (zcell .gt. zlower_ocean) then
            lambda_cell = lambda_water
            mu_cell = mu_water
            rho_cell = rho_water
          else
            lambda_cell = lambda_plate
            mu_cell = mu_plate
            rho_cell = rho_plate
          endif

          aux(1,i,j) = rho_cell
          aux(2,i,j) = lambda_cell
          aux(3,i,j) = mu_cell

          ! Calculate pressure and shear wave speeds
          cp = dsqrt((lambda_cell + 2*mu_cell)/rho_cell)
          cs = dsqrt(mu_cell/rho_cell)
          aux(4,i,j) = cp
          aux(5,i,j) = cs

          ! ============================================
          !BEGINS CALCULATING AND SAVING NORMALS, LENGTH RATIOS AND AREA RATIO FOR MAPPED GRID

          ! Computes corners of grid cell
! c           # lower left corner:
          xccorn(1) = xlower + float(i-1)*dx
          zccorn(1) = zlower + float(j-1)*dz
          call mapc2p(xccorn(1),zccorn(1),xpcorn(1),zpcorn(1))

! c           # upper left corner:
          xccorn(2) = xccorn(1)
          zccorn(2) = zccorn(1) + dz
          call mapc2p(xccorn(2),zccorn(2),xpcorn(2),zpcorn(2))

! c           # upper right corner:
          xccorn(3) = xccorn(1) + dx
          zccorn(3) = zccorn(1) + dz
          call mapc2p(xccorn(3),zccorn(3),xpcorn(3),zpcorn(3))

! c           # lower right corner:
          xccorn(4) = xccorn(1) + dx
          zccorn(4) = zccorn(1)
          call mapc2p(xccorn(4),zccorn(4),xpcorn(4),zpcorn(4))

          ! Compute inner normals
          !left
          xn = zpcorn(2) - zpcorn(1)
          yn = -(xpcorn(2) - xpcorn(1))
          if (dabs(yn) < 1e-10) then
            yn = 0.0
          end if
          norm = dsqrt(xn*xn + yn*yn)
          aux(6,i,j) = xn/norm
          aux(7,i,j) = yn/norm
          aux(8,i,j) = norm/dz

          !bottom
          xn = -(zpcorn(4) - zpcorn(1))
          yn = xpcorn(4) - xpcorn(1)
          if (dabs(xn) < 1e-10) then
            xn = 0.0
          end if
          norm = dsqrt(xn*xn + yn*yn)
          aux(9,i,j) = xn/norm
          aux(10,i,j) = yn/norm
          aux(11,i,j) = norm/dx

          ! Computes area of grid cell using cross product of diagonal
          areap = (xpcorn(3) - xpcorn(1))*(zpcorn(2) - zpcorn(4))
          areap = areap - (zpcorn(3)-zpcorn(1))*(xpcorn(2)-xpcorn(4))
          areap = 0.5*dabs(areap)
          aux(12,i,j) = areap/(dx*dz)

          ! Initialize slip to zero
          aux(13,i,j) = 0.d0

          ! set absorbing layer factor in x direction
          if (ABLdepth > 1.d-10) then
            if (xcell .le. ABLxpos(1)) then
              aux(14,i,j) = 1.d0/(1.d0 + dtan(pi2*(ABLxpos(1) - xcell)/ABLdepth)**2)
            elseif (xcell .ge. ABLxpos(2)) then
              aux(14,i,j) = 1.d0/(1.d0 + dtan(pi2*(xcell - ABLxpos(2))/ABLdepth)**2)
            else
              aux(14,i,j) = 1.d0
            end if
          end if

          ! set absorbing layer factor in y direction
          if (ABLdepth > 1.d-10) then
            if (zcell .le. ABLypos) then
              aux(15,i,j) = 1.d0/(1.d0 + dtan(pi2*(ABLypos - zcell)/ABLdepth)**2)
            else
              aux(15,i,j) = 1.d0
            end if
          end if

        end do
      end do

end subroutine setaux
