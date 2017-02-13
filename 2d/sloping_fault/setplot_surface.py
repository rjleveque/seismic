
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
from mapping import Mapping
from clawpack.clawutil.data import ClawData
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
    
    updip_gaugeno = 30
    downdip_gaugeno = 305

    fault = dtopotools.Fault()
    fault.read(plotdata.outdir + '/fault.data')

    mapping = Mapping(fault)
    fault_width = mapping.fault_width
    xcenter = mapping.xcenter
    ycenter = mapping.ycenter
    xp1 = mapping.xp1
    xp2 = mapping.xp2
    yp1 = mapping.yp1
    yp2 = mapping.yp2

    probdata = ClawData()
    probdata.read(plotdata.outdir + '/setprob.data',force=True)
    xlimits = [xcenter-0.5*probdata.domain_width,xcenter+0.5*probdata.domain_width]
    ylimits = [-probdata.domain_depth,0.0]
    abl_depth = probdata.abl_depth
    xlimits_trunc = [xlimits[0]+abl_depth,xlimits[1]-abl_depth]
    ylimits_trunc = [ylimits[0]+abl_depth,ylimits[1]]

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'

    gaugedata = ClawData()
    gaugedata.read(plotdata.outdir + '/gauges.data',force=True)
    ngauges = gaugedata.ngauges
    xc = np.zeros(ngauges)
    for j in range(ngauges):
        g = plotdata.getgauge(j)
        xc[j] = g.location[0]
    fault.create_dtopography(xc/LAT2METER,np.array([0.]),[1.0],y_disp=True)

    def plot_vertical_displacement(current_data):
        from pylab import plot,zeros,gca,legend
        t = current_data.t

        ys = zeros(ngauges)
        for gaugeno in range(ngauges):
            g = plotdata.getgauge(gaugeno)
            for k in range(1,len(g.t)):
                if g.t[k] > t:
                    break
                dt = g.t[k] - g.t[k-1]
                v = 0.5*(g.q[4,k]+g.q[4,k-1])
                ys[gaugeno] += dt*v
        ax = gca()
        kwargs ={'linestyle':'-','color':'blue','label':'numerical'}
        plot(xc[:ngauges],ys,**kwargs)
        kwargs = {'linestyle':'--','color':'r','label':'Okada'}
        fault.plot_okada(ax,displacement='vertical',kwargs=kwargs)
        legend()

    def plot_fault(current_data):
        from pylab import linspace, plot
        xl = linspace(xp1,xp2,100)
        yl = linspace(yp1,yp2,100)
        plot(xl,yl,'k')

    def plot_gauge_vertical_displacement(t,gaugeno):
        from pylab import plot,zeros
 
        g = plotdata.getgauge(gaugeno)
        ys = zeros(len(g.t))
        for k in range(1,len(g.t)):
            if (g.t[k] > t):
                break
            ys[k] = ys[k-1] + 0.5*(g.t[k]-g.t[k-1])*(g.q[4,k]+g.q[4,k-1])

        plot(g.t[0:k],ys[0:k],'-b')

    def plot_gauge_horizontal_displacement(t,gaugeno):
        from pylab import plot,zeros
 
        g = plotdata.getgauge(gaugeno)
        xs = zeros(len(g.t))
        for k in range(1,len(g.t)):
            if (g.t[k] > t):
                break
            xs[k] = xs[k-1] + 0.5*(g.t[k]-g.t[k-1])*(g.q[3,k]+g.q[3,k-1])

        plot(g.t[0:k],xs[0:k],'-b')

    def plot_gauge_vertical_displacement_updip(current_data):
        plot_gauge_vertical_displacement(current_data.t,updip_gaugeno)

    def plot_gauge_vertical_displacement_downdip(current_data):
        plot_gauge_vertical_displacement(current_data.t,downdip_gaugeno)

    def plot_gauge_horizontal_displacement_updip(current_data):
        plot_gauge_horizontal_displacement(current_data.t,updip_gaugeno)

    def plot_gauge_horizontal_displacement_downdip(current_data):
        plot_gauge_horizontal_displacement(current_data.t,downdip_gaugeno)

    def plot_gauge_vertical_acceleration(t,gaugeno):
        from pylab import plot,zeros

        g = plotdata.getgauge(gaugeno)
        dvs = zeros(len(g.t))
        for k in range(1,len(g.t)):
            if (g.t[k] > t):
                break
            dvs[k] = (g.q[4,k]-g.q[4,k-1])/(g.t[k]-g.t[k-1])

        plot(g.t[0:k],dvs[0:k],'-r')

    def plot_gauge_horizontal_acceleration(t,gaugeno):
        from pylab import plot,zeros

        g = plotdata.getgauge(gaugeno)
        dus = zeros(len(g.t))
        for k in range(1,len(g.t)):
            if (g.t[k] > t):
                break
            dus[k] = (g.q[3,k]-g.q[3,k-1])/(g.t[k]-g.t[k-1])

        plot(g.t[0:k],dus[0:k],'-r')

    def plot_gauge_vertical_acceleration_updip(current_data):
        plot_gauge_vertical_acceleration(current_data.t,updip_gaugeno)

    def plot_gauge_vertical_acceleration_downdip(current_data):
        plot_gauge_vertical_acceleration(current_data.t,downdip_gaugeno)

    def plot_gauge_horizontal_acceleration_updip(current_data):
        plot_gauge_horizontal_acceleration(current_data.t,updip_gaugeno)

    def plot_gauge_horizontal_acceleration_downdip(current_data):
        plot_gauge_horizontal_acceleration(current_data.t,downdip_gaugeno)

    def dip_direction_vel(current_data):
        # return vel dot tau, where tau is tangent to fault
        tau_x = (xp2 - xp1)/fault_width
        tau_y = (yp2 - yp1)/fault_width
        u = current_data.q[3,:,:]
        v = current_data.q[4,:,:]
        return u*tau_x + v*tau_y

    # Figure for surfaces and s-waves
    plotfigure = plotdata.new_plotfigure(name='surface_and_s_waves', figno=1)
    plotfigure.kwargs = {'figsize':(8,8)}

    # Set axes for vertical displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-0.3, 0.5]
    plotaxes.title = 'vertical displacement'
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_vertical_displacement

    # Set axes for slip-direction velocity:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = 'dip-direction velocity'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_fault

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = dip_direction_vel
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.01
    plotitem.pcolor_cmax = 0.01
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

    # Figure for surfaces and s-waves
    plotfigure = plotdata.new_plotfigure(name='surface_and_s_waves (truncated)', figno=2)
    plotfigure.kwargs = {'figsize':(8,8)}

    # Set axes for vertical displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits_trunc
    plotaxes.ylimits = [-0.3, 0.5]
    plotaxes.title = 'vertical displacement'
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_vertical_displacement

    # Set axes for slip-direction velocity:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits_trunc
    plotaxes.ylimits = ylimits_trunc
    plotaxes.title = 'dip-direction velocity'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_fault

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = dip_direction_vel
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.01
    plotitem.pcolor_cmax = 0.01
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

    # Figure for up-dip ground motion
    plotfigure = plotdata.new_plotfigure(name='up-dip ground motion', figno=3)
    plotfigure.kwargs = {'figsize':(8,8)}

    # Set axes for surface profile:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = xlimits_trunc
    plotaxes.ylimits = [-0.3, 0.5]
    plotaxes.title = 'surface'
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_vertical_displacement

    # Set axes for up-dip vertical displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = [0.0,200.0]
    plotaxes.ylimits = [-0.02, 0.02]
    plotaxes.title = ('vertical displacement at x=%s' \
                       % plotdata.getgauge(updip_gaugeno).location[0])
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_gauge_vertical_displacement_updip

    # Set axes for up-dip horizontal displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = [0.0,200.0]
    plotaxes.ylimits = [-0.01, 0.02]
    plotaxes.title = ('horizontal displacement at x=%s' \
                       % plotdata.getgauge(updip_gaugeno).location[0])
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_gauge_horizontal_displacement_updip

    # Figure for down-dip ground motion
    plotfigure = plotdata.new_plotfigure(name='down-dip ground motion', figno=4)
    plotfigure.kwargs = {'figsize':(8,8)}

    # Set axes for surface profile:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = xlimits_trunc
    plotaxes.ylimits = [-0.3, 0.5]
    plotaxes.title = 'surface'
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_vertical_displacement

    # Set axes for down-dip vertical displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = [0.0,200.0]
    plotaxes.ylimits = [-0.03, 0.02]
    plotaxes.title = ('vertical displacement at x=%s' \
                       % plotdata.getgauge(downdip_gaugeno).location[0])
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_gauge_vertical_displacement_downdip

    # Set axes for down-dip horizontal displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = [0.0,200.0]
    plotaxes.ylimits = [-0.08, 0.01]
    plotaxes.title = ('horizontal displacement at x=%s' \
                       % plotdata.getgauge(downdip_gaugeno).location[0])
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_gauge_horizontal_displacement_downdip

    # Figure for up-dip ground acceleration
    plotfigure = plotdata.new_plotfigure(name='up-dip ground acceleration', figno=5)
    plotfigure.kwargs = {'figsize':(8,8)}

    # Set axes for surface profile:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = xlimits_trunc
    plotaxes.ylimits = [-0.3, 0.5]
    plotaxes.title = 'surface'
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_vertical_displacement

    # Set axes for up-dip vertical acceleration:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = [0.0,200.0]
    plotaxes.ylimits = [-0.015, 0.015]
    plotaxes.title = ('vertical acceleration at x=%s' \
                       % plotdata.getgauge(updip_gaugeno).location[0])
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_gauge_vertical_acceleration_updip

    # Set axes for up-dip horizontal acceleration:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = [0.0,200.0]
    plotaxes.ylimits = [-0.015, 0.015]
    plotaxes.title = ('horizontal acceleration at x=%s' \
                       % plotdata.getgauge(updip_gaugeno).location[0])
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_gauge_horizontal_acceleration_updip

    # Figure for down-dip ground acceleration
    plotfigure = plotdata.new_plotfigure(name='down-dip ground acceleration', figno=6)
    plotfigure.kwargs = {'figsize':(8,8)}

    # Set axes for surface profile:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = xlimits_trunc
    plotaxes.ylimits = [-0.3, 0.5]
    plotaxes.title = 'surface'
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_vertical_displacement

    # Set axes for down-dip vertical acceleration:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = [0.0,200.0]
    plotaxes.ylimits = [-0.015, 0.015]
    plotaxes.title = ('vertical acceleration at x=%s' \
                       % plotdata.getgauge(downdip_gaugeno).location[0])
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_gauge_vertical_acceleration_downdip

    # Set axes for down-dip horizontal acceleration:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = [0.0,200.0]
    plotaxes.ylimits = [-0.015, 0.015]
    plotaxes.title = ('horizontal acceleration at x=%s' \
                       % plotdata.getgauge(downdip_gaugeno).location[0])
    plotaxes.title_with_t = False
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_gauge_horizontal_acceleration_downdip

    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=7)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = 'Level 3 grid patches'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_fault

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#effeee', '#eeffee', '#eeeffe',
                                  '#eeeeff', '#ffffff']
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0,0,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

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
    plotdata.parallel = True

    return plotdata
