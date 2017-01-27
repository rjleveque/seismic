import numpy as np
from clawpack.geoclaw.data import LAT2METER

class Mapping2D(object):

    def __init__(self, fault):

        xp1 = 1e10
        xp2 = -1e10
        for subfault in fault.subfaults:
            theta = subfault.dip/180.0*np.pi
            tmp = subfault.longitude*LAT2METER
            if (tmp < xp1):
                xp1 = tmp
                yp1 = -subfault.depth
            tmp = subfault.longitude*LAT2METER + np.cos(theta)*subfault.width
            if (xp2 < tmp):
                xp2 = tmp
                yp2 = -subfault.depth - np.sin(theta)*subfault.width

        xcenter = 0.5*(xp1 + xp2)
        ycenter = 0.5*(yp1 + yp2)
        fault_width = np.sqrt((xp2-xp1)**2 + (yp2-yp1)**2)

        xcl = xcenter - 0.5*fault_width
        xcr = xcenter + 0.5*fault_width

        self.fault_width = fault_width
        self.fault_depth = -ycenter
        self.xcenter = xcenter
        self.ycenter = ycenter
        self.theta = theta
        self.xcl = xcl
        self.xcr = xcr
        self.xp1 = xp1
        self.xp2 = xp2
        self.yp1 = yp1
        self.yp2 = yp2

class Mapping3D(object):

    def __init__(self, fault):

        xp1 = 1e10
        xp2 = -1e10
        yp1 = 1e10
        yp2 = -1e10
        for subfault in fault.subfaults:
            theta = subfault.dip/180.0*np.pi
            tmpx = subfault.longitude*LAT2METER
            tmpy = subfault.latitude - 0.5*subfault.length
            if (tmpx <= xp1 and tmpy <= yp1):
                xp1 = tmpx
                yp1 = tmpy
                zp1 = -subfault.depth
            tmpx = subfault.longitude*LAT2METER + np.cos(theta)*subfault.width
            tmpy = subfault.latitude + 0.5*subfault.length
            if (xp2 <= tmpx and yp2 <= tmpy):
                xp2 = tmpx
                yp2 = tmpy
                zp2 = -subfault.depth - np.sin(theta)*subfault.width

        xcenter = 0.5*(xp1 + xp2)
        ycenter = 0.5*(yp1 + yp2)
        zcenter = 0.5*(zp1 + zp2)
        fault_width = np.sqrt((xp2-xp1)**2 + (zp2-zp1)**2)
        fault_length = yp2 - yp1

        xcl = xcenter - 0.5*fault_width
        xcr = xcenter + 0.5*fault_width

        self.fault_width = fault_width
        self.fault_length = fault_length
        self.fault_depth = -zcenter
        self.xcenter = xcenter
        self.ycenter = ycenter
        self.zcenter = zcenter
        self.theta = theta
        self.xcl = xcl
        self.xcr = xcr
        self.xp1 = xp1
        self.xp2 = xp2
        self.zp1 = zp1
        self.zp2 = zp2