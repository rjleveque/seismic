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

    zlower_ocean = mapping.zlower_ocean
    fault_center = mapping.xcenter
    fault_depth = mapping.fault_depth
    fault_width = mapping.fault_width

    # num of cells here determined in an identical fashion to that in setrun.py
    # additional comments can be found there
    xupper_domain = mapping.xlower_shore
    # determine cell number and set computational boundaries
    dx = fault_width/mfault
    # x direction
    num_cells_above_fault = np.ceil((xupper_domain - (fault_center+0.5*fault_width))/dx)
    target_num_cells = np.rint(probdata.domain_width/dx)
    num_cells_below_fault = target_num_cells - num_cells_above_fault - mfault
    mx = int(target_num_cells)
    xlower_domain = fault_center-0.5*fault_width - num_cells_below_fault*dx
    xupper_domain = fault_center+0.5*fault_width + num_cells_above_fault*dx

    # z direction
    dz = -mapping.zlower_ocean
    mz = int(np.rint(probdata.domain_depth/dz))
    zlower_domain = -mz*dz
    zupper_domain = 0.0

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

    def __init__(self, fault, probdata):
        super(Mapping,self).__init__(fault, probdata)

    def mapc2p(self,xc,zc):
        """
        map computational grid to physical grid that rotates near the fault
        so cell edges match the fault line.  Linear interpolation is used to
        adjust the rotation angle based on distance from fault in computational space.
        The variable tol ensures the physical grid also lines up with a horizontal sea floor
        """

        mult = np.ones(np.shape(zc))
        slope = (self.zlower_shelf - self.zlower_ocean)/(self.xlower_shelf - self.xlower_slope)
        mult = where(xc > self.xlower_slope,
            (self.zlower_ocean + (xc-self.xlower_slope)*slope)/self.zlower_ocean, mult)
        mult = where(xc > self.xlower_shelf, self.zlower_shelf/self.zlower_ocean, mult)
        mult = where(zc < self.zlower_ocean,
            (-self.fault_depth - zc)/(-self.fault_depth-self.zlower_ocean)*mult
            + (zc-self.zlower_ocean)/(-self.fault_depth-self.zlower_ocean)*self.factor, mult)
        mult = where(zc < -self.fault_depth, self.factor, mult)
        xp = xc
        zp = zc*mult

        # define grid that is rotated to line up with fault
        x_rot = self.xcenter + np.cos(self.theta)*(xp-self.xcenter) + np.sin(self.theta)*(zp-self.zcenter)
        z_rot = self.zcenter - np.sin(self.theta)*(xp-self.xcenter) + np.cos(self.theta)*(zp-self.zcenter)
        # construct level set function in computational domain
        tol = self.fault_depth + self.zlower_ocean
        ls = np.abs(zp - self.zcenter)
        ls = np.where(xp < self.xcl, np.sqrt((xp-self.xcl)**2 + (zp-self.zcenter)**2), ls)
        ls = np.where(xp > self.xcr, np.sqrt((xp-self.xcr)**2 + (zp-self.zcenter)**2), ls)
        ls += self.zlower_ocean

        # interpolate between grids
        xp = np.where(ls < tol, (tol-ls)/tol*x_rot + ls/tol*xp, xp)
        zp = np.where(ls < tol, (tol-ls)/tol*z_rot + ls/tol*zp, zp)
        xp = np.where(ls < 0.0, x_rot, xp)
        zp = np.where(ls < 0.0, z_rot, zp)

        return xp,zp
