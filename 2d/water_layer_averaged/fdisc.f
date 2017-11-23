c
c
c
c     =================================================
      function fdisc(x,y)
c     =================================================

      implicit double precision (a-h,o-z)

c
c     # for computing cell averages for initial data that has a
c     # discontinuity along some curve.  fdisc should be negative to the 
c     # left of the curve and positive to the right

c     # three layers with one linear and one quadratic interface
c     # fdisc is positive in the middle layer, negative in the other two

      fdisc = 0.4d0 + 0.2d0*x - y
c
      return
      end

