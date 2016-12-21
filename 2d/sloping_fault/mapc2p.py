from clawpack.clawutil.data import ClawData
import numpy
from pylab import *
from math import atan2

probdata = ClawData()
probdata.read('setprob.data',force=True)
data = ClawData()
data.read('claw.data',force=True)

fault_width = probdata.fault_width
theta = probdata.fault_dip
fault_center = [probdata.fault_center,-probdata.fault_depth]
xcl = fault_center[0] - 0.5*fault_width
xcr = fault_center[0] + 0.5*fault_width
xp1 = fault_center[0] - 0.5*fault_width*cos(theta)
xp2 = fault_center[0] + 0.5*fault_width*cos(theta)
yp1 = fault_center[1] + 0.5*fault_width*sin(theta)
yp2 = fault_center[1] - 0.5*fault_width*sin(theta)

basin_center = probdata.basin_center
basin_width = probdata.basin_width
basin_elevation = probdata.basin_elevation

slope_center = probdata.slope_center
slope_width = probdata.slope_width

shelf_center = probdata.shelf_center
shelf_width = probdata.shelf_width
shelf_elevation = probdata.shelf_elevation

beach_center = probdata.beach_center
beach_width = probdata.beach_width

shore_elevation = probdata.shore_elevation

def surf_fault_ls(x,y):
    ls = numpy.abs(y - fault_center[1])
    ls = numpy.where(x < xcl, numpy.sqrt((x-xcl)**2 + (y-fault_center[1])**2), ls)
    ls = numpy.where(x > xcr, numpy.sqrt((x-xcr)**2 + (y-fault_center[1])**2), ls)
    return ls

def topo(x):
    y = basin_elevation
    y = numpy.where(abs(x - slope_center) < 0.5*slope_width, \
        basin_elevation + ((x - slope_center)/slope_width + 0.5)*(shelf_elevation - basin_elevation), y)
    y = numpy.where(abs(x - shelf_center) < 0.5*shelf_width, shelf_elevation, y)
    y = numpy.where(abs(x - beach_center) < 0.5*beach_width, \
        shelf_elevation + ((x - beach_center)/beach_width + 0.5)*(shore_elevation - shelf_elevation), y)
    y = numpy.where(x > beach_center + 0.5*beach_width, shore_elevation, y)
    return y

# Will need to do better than this once topography is entered via topo files
x = linspace(data.lower[0], data.upper[0], 10000)
y = topo(x)
ls = surf_fault_ls(x,y)
tol = min(ls)

def mapc2p(xc,yc):
    """
    map computational grid to physical grid that rotates near the fault
    so cell edges match the fault line.  Linear interpolation is used to
    adjust the rotation angle based on distance from fault in computational space.
    The variable tol ensures the physical grid also lines up with sea floor
    """

    # construct grid that lines up with sea floor
    xsf = xc
    ysf = yc + topo(xsf)

    # constucted signed distance function in computational domain
    ls = numpy.abs(yc - fault_center[1])
    ls = numpy.where(xc < xcl, numpy.sqrt((xc-xcl)**2 + (yc-fault_center[1])**2), ls)
    ls = numpy.where(xc > xcr, numpy.sqrt((xc-xcr)**2 + (yc-fault_center[1])**2), ls)

#    # define grid that is rotated to line up with fault
    xrot = fault_center[0] + numpy.cos(theta)*(xc-fault_center[0]) + numpy.sin(theta)*(yc-fault_center[1])
    yrot = fault_center[1] - numpy.sin(theta)*(xc-fault_center[0]) + numpy.cos(theta)*(yc-fault_center[1])

    # Interpolate between rotated grid and cartesian grid near the fault,
    # using seafloor grid far away from fault.
    xp = numpy.where(ls < tol, (tol-ls)/tol*xrot + ls/tol*xsf, xsf)
    yp = numpy.where(ls < tol, (tol-ls)/tol*yrot + ls/tol*ysf, ysf)

    return xp,yp

def test(mfault):

    fault_depth = probdata.fault_depth

    # num of cells here determined in an identical fashion to that in setrun.py
    # additional comments can be found there
    dx = fault_width/mfault
    num_cells_above = numpy.rint(fault_depth/dx)
    dy = fault_depth/num_cells_above
    mx = int(numpy.ceil((data.upper[0]-data.lower[0])/dx)) # mx
    my = int(numpy.ceil((data.upper[1]-data.lower[1])/dy)) # my
    mr = mx - mfault

    x = linspace(data.lower[0], data.upper[0], mx+1)
    y = linspace(data.lower[1], data.upper[1], my+1)
    xc,yc = meshgrid(x,y)
    xp,yp = mapc2p(xc,yc)
    figure()
    plot(xp,yp,'k-')
    plot(xp.T,yp.T,'k-')
    plot((xp1,xp2),(yp1,yp2),'-g')
    axis('scaled')
