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

c     # fdisc is positive in the lower layer, negative in upper

      fdisc = -10000. + 6000.*(x+40.d3)/80.d3  - y
c
      return
      end

