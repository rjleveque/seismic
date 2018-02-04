
c
c ------------------------------------------------------------------
c
      subroutine bc2amr(val,aux,nrow,ncol,meqn,naux,
     &                  hx, hy, level, time,
     &                  xlower_patch, xupper_patch,
     &                  ylower_patch,yupper_patch)

c  ---------------------------------------------------------------
c  Modified for elasticity with stress imposed at top boundary.
c  Assume extrapolation will be used at other boundaries.
c  ---------------------------------------------------------------
c
c
c :::::::::: bc2amr ::::::::::::::::::::::::::::::::::::::::::::::;
c
c     Take a grid patch with mesh widths hx,hy, of dimensions nrow by
c     ncol,  and set the values of any piece of
c     of the patch which extends outside the physical domain
c     using the boundary conditions.
c
c     ------------------------------------------------
c     # Standard boundary condition choices for amr2ez in clawpack
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (bottom), 4 (top):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary conditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     #            =  4  sphere bcs (left half maps to right half of same
c     #                  side, and vice versa), as if domain folded in half
c     ------------------------------------------------
c
c     The corners of the grid patch are at
c        (xlower_patch,ylower_patch)  --  lower left corner
c        (xupper_patch,yupper_patch) --  upper right corner
c
c     The physical domain itself is a rectangle bounded by
c        (xlower,ylower)  -- lower left corner
c        (xupper,yupper)  -- upper right corner
c
c     the picture is the following:
c
c               _____________________ (xupper,yupper)
c              |                     |
c          _________ (xupper_patch,yupper_patch)   |
c          |   |    |                |
c          |   |    |                |
c          |   |    |                |
c          |___|____|                |
c (xlower_patch,ylower_patch) |                     |
c              |                     |
c              |_____________________|
c   (xlower,ylower)
c
c
c     Any cells that lie outside the physical domain are ghost cells whose
c     values should be set in this routine.  This is tested for by comparing
c     xlower_patch with xlower to see if values need to be set at the left, as in
c     the figure above, and similarly at the other boundaries.
c
c     Patches are guaranteed to have at least 1 row of cells filled
c     with interior values so it is possible to  extrapolate.
c     Fix trimbd if you want more than 1 row pre-set.
c
c     Make sure the order the boundaries are specified is correct
c     so that diagonal corner cells are also properly taken care of.
c
c     Periodic boundaries are set before calling this routine, so if the
c     domain is periodic in one direction only you
c     can safely extrapolate in the other direction.
c
c     Don't overwrite ghost cells in periodic directions!
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      use amr_module, only: mthbc,xlower,ylower,xupper,yupper
      use amr_module, only: xperdom,yperdom,spheredom

      implicit double precision (a-h,o-z)

      real*8  val(meqn,nrow,ncol), aux(naux,nrow,ncol)
      integer nrow,ncol,meqn,naux,level
      real*8  hx,hy,time, hxmarg, hymarg
      real*8  xlo_patch,xhi_patch,ylo_patch,yhi_patch
      integer nxl,nxr,ibeg,nyb,nyt,jbeg,i,j,m

      real(kind=8) :: lambda_plate, mu_plate, rho_plate, lambda_water,
     &     mu_water, rho_water, g
      common /material/ lambda_plate, mu_plate, rho_plate, lambda_water,
     &     mu_water, rho_water, g

      hxmarg = hx*.01
      hymarg = hy*.01

      if (xperdom .and. (yperdom .or. spheredom)) go to 499
c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      if (xlower_patch .ge. xlower-hxmarg) then
c        # not a physical boundary -- no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 199
         endif
c
c     # number of grid cells from this patch lying outside physical domain:
      nxl = (xlower+hxmarg-xlower_patch)/hx
c
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*)
     &   '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2amr'
      stop
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do 115 j = 1,ncol
         do 115 i=1,nxl
            do 115 m=1,meqn
               val(m,i,j) = val(m,nxl+1,j)
  115       continue
      go to 199

  120 continue
c     # periodic:   handled elsewhere in amr
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 j = 1,ncol
         do 135 i=1,nxl
            do 135 m=1,meqn
               val(m,i,j) = val(m,2*nxl+1-i,j)
  135       continue
c     # negate the normal velocity:
      do 136 j = 1,ncol
         do 136 i=1,nxl
            val(2,i,j) = -val(2,i,j)
  136    continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      if (xupper_patch .le. xupper+hxmarg) then
c        # not a physical boundary --  no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 299
         endif
c
c     # number of grid cells lying outside physical domain:
      nxr = (xupper_patch - xupper + hxmarg)/hx
      ibeg = max0(nrow-nxr+1, 1)
c
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*)
     &   '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2amr'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do 215 j = 1,ncol
         do 215 i=ibeg,nrow
            do 215 m=1,meqn
               val(m,i,j) = val(m,ibeg-1,j)
  215       continue
      go to 299

  220 continue
c     # periodic:   handled elsewhere in amr
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 j = 1,ncol
         do 235 i=ibeg,nrow
            do 235 m=1,meqn
               val(m,i,j) = val(m,2*ibeg-1-i,j)
  235       continue
c     # negate the normal velocity:
      do 236 j = 1,ncol
         do 236 i=ibeg,nrow
            val(2,i,j) = -val(2,i,j)
  236    continue
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      if (ylower_patch .ge. ylower-hymarg) then
c        # not a physical boundary -- no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 399
         endif
c
c     # number of grid cells lying outside physical domain:
      nyb = (ylower+hymarg-ylower_patch)/hy
c
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # user specified - tangential motion
      do 305 j=1,nyb
         do 305 i=1,nrow
            do 305 m=1,meqn
               val(m,i,j) =  val(m,i,2*nyb+1-j)
  305       continue
c     # negate the normal velocity and impose tangential velocity:
      pi = acos(-1.d0)
      do 306 j=1,nyb
         do 306 i=1,nrow
            xcell = xlower_patch + (i-0.5d0)*hx
            if (xcell.gt.0.95d0 .and. xcell.lt.1.05d0) then
                if (time < 0.1d0) then
                    s = 1.d0 - cos(2*pi*time/0.1d0)
                  else
                    s = 0.d0
                  endif
                val(4,i,j) = s - val(4,i,j)
                val(5,i,j) = -val(5,i,j)
            else
                val(4,i,j) = -val(4,i,j)
                val(5,i,j) = -val(5,i,j)
            endif
  306    continue
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 j=1,nyb
         do 315 i=1,nrow
            do 315 m=1,meqn
                val(m,i,j) = val(m,i,nyb+1)
  315       continue
      go to 399

  320 continue
c     # periodic:   handled elsewhere in amr
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 j=1,nyb
         do 335 i=1,nrow
            do 335 m=1,meqn
               val(m,i,j) =  val(m,i,2*nyb+1-j)
  335       continue
c     # negate the normal velocity:
      do 336 j=1,nyb
         do 336 i=1,nrow
            val(2,i,j) = -val(2,i,j)
            val(3,i,j) = -val(3,i,j)
  336    continue
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      if (yupper_patch .le. yupper+hymarg) then
c        # not a physical boundary --  no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 499
         endif
c
c     # number of grid cells lying outside physical domain:
      nyt = (yupper_patch - yupper + hymarg)/hy
      jbeg = max0(ncol-nyt+1, 1)
c
      go to (400,410,420,430) mthbc(4)+1
c

  400 continue
c     # First-order extrapolate u, v, h
c     # Negate sigma_xy, although it should already be zero
c     # Set sigma_xx, sigma_yy to specified pressure
      do 405 j=jbeg,ncol
        do 405 i=1,nrow
            s = -rho_water*g*val(6,i,jbeg-1)
            val(1,i,j) = 2.d0*s - val(1,i,2*jbeg-1-j)
            val(2,i,j) = 2.d0*s - val(2,i,2*jbeg-1-j)
            val(3,i,j) = -val(3,i,2*jbeg-1-j)
            val(4,i,j) = val(4,i,2*jbeg-1-j)
            val(5,i,j) = val(5,i,2*jbeg-1-j)
            val(6,i,j) = val(6,i,2*jbeg-1-j)
  405 continue
      go to 499

  410 continue
c     # zero-order extrapolation:
      do 415 j=jbeg,ncol
         do 415 i=1,nrow
            do 415 m=1,meqn
               val(m,i,j) =  val(m,i,jbeg-1)
  415       continue
      go to 499

  420 continue
c     # periodic:   handled elsewhere in amr
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 j=jbeg,ncol
         do 435 i=1,nrow
            do 435 m=1,meqn
               val(m,i,j) =  val(m,i,2*jbeg-1-j)
  435       continue
c     # negate the normal velocity:
      do 436 j=jbeg,ncol
         do 436 i=1,nrow
            val(3,i,j) = -val(3,i,j)
  436    continue
      go to 499

  499 continue

      return
      end