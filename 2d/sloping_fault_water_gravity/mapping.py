import numpy as np
from pylab import *
from clawpack.seismic.mappings import Mapping2D
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
from make_topo_and_grid import get_oceanfloor_parameters

def test(mfault):

    from clawpack.clawutil.data import ClawData

    probdata = ClawData()
    probdata.read('setprob.data',force=True)
    sink_depth = probdata.sink_depth

    fault = dtopotools.Fault()
    fault.read('fault.data')
    mapping = Mapping(fault, sink_depth)
    zlower_ocean = mapping.zlower_ocean
    fault_center = mapping.xcenter
    fault_depth = mapping.fault_depth
    fault_width = mapping.fault_width

    # num of cells here determined in an identical fashion to that in setrun.py
    # additional comments can be found there
    xupper_domain = mapping.xlower_shore
    # determine cell number and set computational boundaries
    target_dh = fault_width/mfault
    # x direction
    num_cells_above_fault = np.rint((xupper_domain - (fault_center+0.5*fault_width))/target_dh)
    dx = (xupper_domain - (fault_center+0.5*fault_width))/num_cells_above_fault
    target_num_cells = np.rint(probdata.domain_width/dx)
    num_cells_below_fault = target_num_cells - num_cells_above_fault - mfault
    mx = int(target_num_cells)
    xlower_domain = fault_center-0.5*fault_width - num_cells_below_fault*dx
    xupper_domain = fault_center+0.5*fault_width + num_cells_above_fault*dx
    # z direction
    num_cells_above_fault = np.rint(fault_depth/target_dh)
    dz = fault_depth/num_cells_above_fault
    target_num_cells = np.rint(probdata.domain_depth/dz)
    num_cells_below_fault = target_num_cells - num_cells_above_fault
    mz = int(target_num_cells)
    zlower_domain = -fault_depth - num_cells_below_fault*dz
    zupper_domain = -fault_depth + num_cells_above_fault*dz

    mapping.zclower_ocean = np.floor(mapping.zlower_ocean/dz)*dz
    print mapping.zclower_ocean
    print fault_depth

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

    def __init__(self, fault, sink_depth):
        super(Mapping,self).__init__(fault)
        # Obtain topography parameters to match 1D Geoclaw
        self.xlower_domain, \
        self.xlower_slope, \
        self.xlower_shelf, \
        self.xlower_beach, \
        self.xlower_shore, \
        self.xupper_domain, \
        self.zlower_ocean, \
        self.zlower_shelf, \
        self.zlower_beach, \
        self.zlower_shore = get_oceanfloor_parameters()
        # set ocean computational depth to null
        self.zclower_ocean = 1.0

    def mapc2p(self,xc,zc):
        """
        map computational grid to physical grid that rotates near the fault
        so cell edges match the fault line.  Linear interpolation is used to
        adjust the rotation angle based on distance from fault in computational space.
        The variable tol ensures the physical grid also lines up with a horizontal sea floor
        """

        # define grid that lines up with ocean floor
        x_flr = xc
        z_flr = zc*self.zlower_ocean/self.zclower_ocean
        slope = (self.zlower_shelf - self.zlower_ocean)/(self.xlower_shelf - self.xlower_slope)
        z_flr = where(xc > self.xlower_slope,
            zc*(self.zlower_ocean + (xc-self.xlower_slope)*slope)/self.zclower_ocean, z_flr)
        z_flr = where(xc > self.xlower_shelf, zc*self.zlower_shelf/self.zclower_ocean, z_flr)

        # construct distance funtion in computation domain_width
        ls_flr = self.zclower_ocean - zc
        ls_flr = np.where(ls_flr < 0.0, 0.0, ls_flr)

        # define grid that is rotated to line up with fault
        x_rot = self.xcenter + np.cos(self.theta)*(xc-self.xcenter) + np.sin(self.theta)*(zc-self.zcenter)
        z_rot = self.zcenter - np.sin(self.theta)*(xc-self.xcenter) + np.cos(self.theta)*(zc-self.zcenter)
        # constuct function in computational domain
        ls_rot = np.abs(zc - self.zcenter)
        ls_rot = np.where(xc < self.xcl, np.sqrt((xc-self.xcl)**2 + (zc-self.zcenter)**2), ls_rot)
        ls_rot = np.where(xc > self.xcr, np.sqrt((xc-self.xcr)**2 + (zc-self.zcenter)**2), ls_rot)

        # Interpolate between custom grids and cartesian grid
        xp = xc
        zp = zc
        tol = self.fault_depth + self.zclower_ocean
        # Interpolate between rotated grid and current grid
        xp = np.where(ls_rot < tol, (tol-ls_rot)/tol*x_rot + ls_rot/tol*xp, xp)
        zp = np.where(ls_rot < tol, (tol-ls_rot)/tol*z_rot + ls_rot/tol*zp, zp)
        # Interpolate between floor grid and current grid
        xp = where(ls_flr < tol, (tol-ls_flr)/tol*x_flr + ls_flr/tol*xp, xp)
        zp = where(ls_flr < tol, (tol-ls_flr)/tol*z_flr + ls_flr/tol*zp, zp)

        return xp,zp
