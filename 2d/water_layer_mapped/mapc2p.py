from clawpack.clawutil.data import ClawData
import numpy
from pylab import *

#probdata = ClawData()
#probdata.read('setprob.data',force=True)

xp1 = 40e3
yp0 = -40e3
ypi1 = -10e3
ypi2 = -4e3
yci = 0.8

def mapc2p(xc,yc):
    spi = (ypi2 - ypi1)/(2.*xp1)
    xp = xp1 * xc
    ypi = ypi1 + spi*(xp + xp1)
    s1 = (ypi - yp0)/yci
    s2 = -ypi/(1.-yci)
    yp = where(yc<yci, yp0+s1*yc, -s2*(1.-yc))
    return xp,yp


def test(mx,my):

    x = linspace(-1,1,mx)
    y = linspace(0,1,my)
    xc,yc = meshgrid(x,y)
    xp,yp = mapc2p(xc,yc)
    figure(500)
    clf()
    plot(xp,yp,'k-')
    plot(xp.T,yp.T,'k-')
    plot((-xp1,xp1),(ypi1,ypi2),'-g', linewidth=2.0)
    #axis('scaled')
