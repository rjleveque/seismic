import numpy as np
from pylab import *
from clawpack.seismic.mappings import Mapping2D
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
from make_topo_and_grid import get_oceanfloor_parameters

def test(mfault):

    from clawpack.clawutil.data import ClawData

    probdata = ClawData()
    probdata.read('setprob.data',force=True)
    fault = dtopotools.Fault()
    fault.read('fault.data')
    mapping = Mapping(fault, probdata)

    xlower_slope = mapping.xlower_slope
    xlower_shelf = mapping.xlower_shelf
    zlower_ocean = mapping.zlower_ocean
    fault_center = mapping.xcenter
    fault_depth = mapping.fault_depth
    fault_width = mapping.fault_width

    # num of cells here determined in an identical fashion to that in setrun.py
    # additional comments can be found there
    xupper_domain = mapping.xlower_shore
    # determine cell number and set computational boundaries
    target_dx = fault_width/mfault
    # x direction
    num_cells_across_slope = np.ceil((xlower_shelf-xlower_slope)/target_dx)
    dx = (xlower_shelf-xlower_slope)/num_cells_across_slope
    num_cells_above_slope = np.ceil((xupper_domain - xlower_shelf)/dx)
    target_num_cells = np.rint(probdata.domain_width/dx)
    num_cells_below_slope = target_num_cells - num_cells_above_slope - num_cells_across_slope
    mx = int(target_num_cells)
    xlower_domain = xlower_slope - num_cells_below_slope*dx
    xupper_domain = xlower_shelf + num_cells_above_slope*dx

    # z direction
    num_cells_across_ocean = np.ceil(-zlower_ocean/dx)
    dz = -zlower_ocean/num_cells_across_ocean
    mz = int(np.rint(probdata.domain_depth/dz))
    zlower_domain = -mz*dz
    zupper_domain = 0.0

    mapping.fault_zshift = -fault_depth + np.ceil(fault_depth/dz)*dz

    x = linspace(xlower_domain, xupper_domain, mx+1)
    z = linspace(zlower_domain, zupper_domain, mz+1)
    xc,zc = meshgrid(x,z)
    xp,zp = mapping.mapc2p(xc,zc)
    figure()
    plot(xp,zp,'k-')
    plot(xp.T,zp.T,'k-')
    plot((mapping.xp1,mapping.xp2),(mapping.zp1,mapping.zp2),'-g')
    axis('scaled')


class Mapping(Mapping2D):

    def __init__(self, fault, probdata, clawdata=None):
        super(Mapping,self).__init__(fault)
        # Obtain topography parameters
        self.xlower_slope = probdata.xlower_slope
        self.xlower_shelf = probdata.xlower_shelf
        self.xlower_beach = probdata.xlower_beach
        self.xlower_shore = probdata.xlower_shore
        self.zlower_ocean = probdata.zlower_ocean
        self.zlower_shelf = probdata.zlower_shelf
        self.zlower_shore = probdata.zlower_shore

        self.fault_zshift = 0.0
        if (clawdata is not None):
            dz = (clawdata.upper[1]-clawdata.lower[1])/clawdata.num_cells[1]
            self.fault_zshift = self.zcenter + np.ceil(-self.zcenter/dz)*dz

    def mapc2p(self,xc,zc):
        """
        map computational grid to physical grid that rotates near the fault
        so cell edges match the fault line.  Linear interpolation is used to
        adjust the rotation angle based on distance from fault in computational space.
        The variable tol ensures the physical grid also lines up with a horizontal sea floor
        """

        # define grid that is scaled or shifted to line up with ocean floor
        slope = (self.zlower_shelf - self.zlower_ocean)/(self.xlower_shelf - self.xlower_slope)
        scale = np.ones(np.shape(zc))
        scale = where(xc > self.xlower_slope,
             (self.zlower_ocean + (xc-self.xlower_slope)*slope)/self.zlower_ocean, scale)
        scale = where(xc > self.xlower_shelf, self.zlower_shelf/self.zlower_ocean, scale)
        shift = np.zeros(np.shape(zc))
        shift = where(xc > self.xlower_slope, (xc-self.xlower_slope)*slope, shift)
        shift = where(xc > self.xlower_shelf, self.zlower_shelf-self.zlower_ocean, shift)
        shift = where(zc < self.zlower_ocean,
            (zc + self.fault_depth)/(self.zlower_ocean+self.fault_depth)*shift, shift)
        shift = where(zc < -self.fault_depth, 0.0, shift)
        x_floor = xc
        z_floor = where(zc > self.zlower_ocean, zc*scale, zc+shift)

        # define grid that is rotated to line up with fault
        zc = zc + self.fault_zshift
        x_rot = self.xcenter + np.cos(self.theta)*(xc-self.xcenter) + np.sin(self.theta)*(zc-self.zcenter)
        z_rot = self.zcenter - np.sin(self.theta)*(xc-self.xcenter) + np.cos(self.theta)*(zc-self.zcenter)
        # construct level set function in computational domain
        tol = self.fault_depth + self.zlower_ocean
        ls = np.abs(zc - self.zcenter)
        ls = np.where(xc < self.xcl, np.sqrt((xc-self.xcl)**2 + (zc-self.zcenter)**2), ls)
        ls = np.where(xc > self.xcr, np.sqrt((xc-self.xcr)**2 + (zc-self.zcenter)**2), ls)
        ls += self.zlower_ocean
        tol += self.zlower_ocean

        # interpolate between grids
        xp = np.where(ls < tol, (tol-ls)/tol*x_rot + ls/tol*x_floor, x_floor)
        zp = np.where(ls < tol, (tol-ls)/tol*z_rot + ls/tol*z_floor, z_floor)
        xp = np.where(ls < 0.0, x_rot, xp)
        zp = np.where(ls < 0.0, z_rot, zp)

        return xp,zp
