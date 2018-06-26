
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
from clawpack.clawutil.data import ClawData
from clawpack.visclaw.data import ClawPlotData
from mapping import Mapping
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
from clawpack.geoclaw.data import LAT2METER
reload(dtopotools)

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    fault = dtopotools.Fault()
    fault.read(plotdata.outdir + '/fault.data')
    probdata = ClawData()
    probdata.read(plotdata.outdir + '/setprob.data',force=True)
    clawdata = ClawData()
    clawdata.read(plotdata.outdir + '/claw.data', force=True)

    mapping = Mapping(fault, probdata, clawdata)
    fault_width = mapping.fault_width
    xcenter = mapping.xcenter
    zcenter = mapping.zcenter
    xp1 = mapping.xp1
    xp2 = mapping.xp2
    zp1 = mapping.zp1
    zp2 = mapping.zp2

    xlimits = [xcenter-0.5*probdata.domain_width+probdata.abl_depth,
               xcenter+0.5*probdata.domain_width-probdata.abl_depth]
    zlimits = [-probdata.domain_depth+probdata.abl_depth, 0.0]
    xlimitsW = xlimits
    zlimitsW = [probdata.zlower_ocean, 0.0]
    abl_depth = probdata.abl_depth

    gaugedata = ClawData()
    gaugedata.read(plotdata.outdir + '/gauges.data',force=True)
    ngauges = gaugedata.ngauges/2
    xc = np.zeros(ngauges)
    for j in range(ngauges):
        g = plotdata.getgauge(j)
        xc[j] = g.location[0]
    fault.create_dtopography(xc/LAT2METER,np.array([0.]),[1.0],horiz_disp=True)


    plotdata_avg = ClawPlotData()
    plotdata_avg.outdir = '/home/vogl2/workspace/clawpack/seismic/2d/sloping_fault_water_gravity_averaged/output_topo'

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'
    plotdata.parallel = True

    def plot_seasurface(current_data):
        from pylab import plot,zeros,gca, legend
        t = current_data.t

        zs = zeros(ngauges)
        zs_avg = zeros(ngauges)
        for gaugeno in range(ngauges):
            g = plotdata.getgauge(ngauges+gaugeno)
            for k in range(1,len(g.t)):
                if (g.t[k] > t and g.t[k-1] <= t):
                    dt = g.t[k] - g.t[k-1]
                    zs[gaugeno] = g.q[5,k]*(t-g.t[k-1])/dt +\
                                  g.q[5,k-1]*(g.t[k]-t)/dt
            g = plotdata_avg.getgauge(ngauges+gaugeno)
            for k in range(1,len(g.t)):
                if (g.t[k] > t and g.t[k-1] <= t):
                    dt = g.t[k] - g.t[k-1]
                    zs_avg[gaugeno] = g.q[5,k]*(t-g.t[k-1])/dt +\
                                  g.q[5,k-1]*(g.t[k]-t)/dt

        ax = gca()
        plot(xc[:ngauges],zs,'-b',linewidth=2,label='mapped topography')
        plot(xc[:ngauges],zs_avg,'-r',linewidth=2,label='averaged parameters')
        ax.legend(fontsize=18, loc='lower left')
        ax.set_ylabel('ocean surface',fontsize=18)


    def plot_seafloor(current_data):
        from pylab import plot,zeros,gca, legend, colorbar, gcf
        t = current_data.t

        zf = zeros(ngauges)
        zf_avg = zeros(ngauges)
        for gaugeno in range(ngauges):
            g = plotdata.getgauge(gaugeno)
            for k in range(1,len(g.t)):
                if (g.t[k] > t and g.t[k-1] <= t):
                    dt = g.t[k] - g.t[k-1]
                    zf[gaugeno] = g.q[5,k]*(t-g.t[k-1])/dt +\
                                  g.q[5,k-1]*(g.t[k]-t)/dt
            g = plotdata_avg.getgauge(gaugeno)
            for k in range(1,len(g.t)):
                if (g.t[k] > t and g.t[k-1] <= t):
                    dt = g.t[k] - g.t[k-1]
                    zf_avg[gaugeno] = g.q[5,k]*(t-g.t[k-1])/dt +\
                                  g.q[5,k-1]*(g.t[k]-t)/dt

        ax = gca()
        plot(xc[:ngauges],zf,'-k',linewidth=2,label='mapped topography')
        plot(xc[:ngauges],zf_avg,'-g',linewidth=2,label='averaged parameters')
        ax.legend(fontsize=18,loc='lower left')
        ax.set_ylabel('ocean floor',fontsize=18)
        f = gcf()
        f.tight_layout()

    # Figure
    plotfigure = plotdata.new_plotfigure(name='results', figno=1)
    plotfigure.kwargs = {'figsize':(10,8)}

    # Set axes for seasurface:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimitsW
    plotaxes.ylimits = [-0.3,0.5]
    plotaxes.title = ''
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_seasurface

    # Set axes for seafloor:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimitsW
    plotaxes.ylimits = [-0.3,0.5]
    plotaxes.title = ''
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_seafloor

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
#    plotdata.parallel = True

    return plotdata
