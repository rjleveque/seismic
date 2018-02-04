import numpy
from pylab import *
from clawpack.seismic.mappings import Mapping2D

def test(mfault):

    from clawpack.clawutil.data import ClawData

    probdata = ClawData()
    probdata.read('setprob.data',force=True)

    fault = dtopotools.Fault()
    fault.read('fault.data')

    mapping = Mapping(fault)

    domain_depth = probdata.domain_depth
    domain_width = probdata.domain_width

    # num of cells here determined in an identical fashion to that in setrun.py
    # additional comments can be found there
    dx = mapping.fault_width/mfault
    num_cells_above = numpy.rint(mapping.fault_depth/dx)
    dy = mapping.fault_depth/num_cells_above
    mx = int(numpy.ceil(domain_width/dx)) # mx
    my = int(numpy.ceil(domain_depth/dy)) # my
    mr = mx - mfault

    x = linspace(mapping.xcenter-0.5*mapping.fault_width - numpy.floor(mr/2.0)*dx, mapping.xcenter+0.5*mapping.fault_width + numpy.ceil(mr/2.0)*dx, mx+1)
    y = linspace(-my*dy, 0.0, my+1)
    xc,yc = meshgrid(x,y)
    xp,yp = mapping.mapc2p(xc,yc)
    figure()
    plot(xp,yp,'k-')
    plot(xp.T,yp.T,'k-')
    plot((mapping.xp1,mapping.xp2),(mapping.yp1,mapping.yp2),'-g')
    axis('scaled')


class Mapping(Mapping2D):

    def __init__(self, fault, water_scaling=0):
        super(Mapping,self).__init__(fault)
        self.water_scaling = water_scaling

    def mapc2p(self,xc,yc):
        """
        map computational grid to physical grid that rotates near the fault
        so cell edges match the fault line.  Linear interpolation is used to
        adjust the rotation angle based on distance from fault in computational space.
        The variable tol ensures the physical grid also lines up with a horizontal sea floor
        """

        # constucted signed distance function in computational domain
        ls = numpy.abs(yc - self.ycenter)
        ls = numpy.where(xc < self.xcl, numpy.sqrt((xc-self.xcl)**2 + (yc-self.ycenter)**2), ls)
        ls = numpy.where(xc > self.xcr, numpy.sqrt((xc-self.xcr)**2 + (yc-self.ycenter)**2), ls)

        # define grid that is rotated to line up with fault
        xrot = self.xcenter + numpy.cos(self.theta)*(xc-self.xcenter) + numpy.sin(self.theta)*(yc-self.ycenter)
        yrot = self.ycenter - numpy.sin(self.theta)*(xc-self.xcenter) + numpy.cos(self.theta)*(yc-self.ycenter)

        # Interpolate between rotated grid and cartesian grid near the fault,
        # using cartesian grid far away from fault.
        tol = self.fault_depth
        xp = xc
        yp = yc
        yp = numpy.where(yc > 0, yc*self.water_scaling, yc)
        xp = numpy.where(ls < tol, (tol-ls)/tol*xrot + ls/tol*xc, xp)
        yp = numpy.where(ls < tol, (tol-ls)/tol*yrot + ls/tol*yc, yp)

        return xp,yp
