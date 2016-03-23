from pylab import *
from clawpack.geoclaw import dtopotools
fault = dtopotools.Fault()
fault.subfaults = []

width = sqrt(50e3**2 + 8.6e3**2)

ndip = 1

for i in range(ndip):
    subfault = dtopotools.SubFault()
    subfault.mu = 3e10
    subfault.dip = 9.9
    subfault.width = width / ndip
    x = i*50e3/ndip
    subfault.depth = -15e3 - x*sin(subfault.dip * pi/180.)
    subfault.slip = 1.  #exp(-((x-25e3)/12e3)**2)
    subfault.rake = 90
    subfault.strike = 0
    subfault.length = 1000e3
    subfault.longitude = x
    subfault.latitude = 0.
    subfault.coordinate_specification = 'top center'
    subfault.calculate_geometry()

    fault.subfaults.append(subfault)

print "Mw = ", fault.Mw()

x = linspace(-1.5,1.5,1000)
y = array([0.])
dtopo = fault.create_dtopography(x,y,[1.])

figure()
plot(dtopo.x*111e3,dtopo.dZ[0,0,:])
show()
