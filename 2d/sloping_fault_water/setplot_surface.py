
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
from clawpack.clawutil.data import ClawData
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

    fault = dtopotools.Fault()
    fault.read(plotdata.outdir + '/fault.data')
    probdata = ClawData()
    probdata.read('setprob.data',force=True)

    mapping = Mapping(fault,probdata.water_scaling)
    fault_width = mapping.fault_width
    xcenter = mapping.xcenter
    ycenter = mapping.ycenter
    xp1 = mapping.xp1
    xp2 = mapping.xp2
    yp1 = mapping.yp1
    yp2 = mapping.yp2

    xlimits = [xcenter-0.5*probdata.domain_width,xcenter+0.5*probdata.domain_width]
    ylimits = [-probdata.domain_depth,0.0]
    xlimitsW = [xp1-10.0*probdata.water_depth,xp2+10.0*probdata.water_depth]
    ylimitsW = [-5.0*probdata.water_depth,probdata.water_depth]
    abl_depth = probdata.abl_depth
    xlimits_trunc = [xlimits[0]+abl_depth,xlimits[1]-abl_depth]
    ylimits_trunc = [ylimits[0]+abl_depth,ylimits[1]]



    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'
    plotdata.parallel = True

    gaugedata = ClawData()
    gaugedata.read('gauges.data',force=True)
    ngauges = gaugedata.ngauges/2
    xc = np.zeros(ngauges)
    for j in range(ngauges):
        g = plotdata.getgauge(j)
        xc[j] = g.location[0]
    fault.create_dtopography(xc/LAT2METER,np.array([0.]),[1.0],y_disp=True)

    def plot_vertical_displacement(current_data):
        from pylab import plot,zeros,gca, legend
        t = current_data.t

        yf = zeros(ngauges)
        for gaugeno in range(ngauges):
            g = plotdata.getgauge(gaugeno)
            for k in range(1,len(g.t)):
                if g.t[k] > t:
                    break
                dt = g.t[k] - g.t[k-1]
                v = 0.5*(g.q[4,k]+g.q[4,k-1])
                yf[gaugeno] += dt*v

        ys = zeros(ngauges)
        for gaugeno in range(ngauges):
            g = plotdata.getgauge(ngauges+gaugeno)
            for k in range(1,len(g.t)):
                if g.t[k] > t:
                    break
                dt = g.t[k] - g.t[k-1]
                v = 0.5*(g.q[4,k]+g.q[4,k-1])
                ys[gaugeno] += dt*v

        ax = gca()
        kwargs ={'linestyle':'-','color':'black','label':'sea floor'}
        plot(xc[:ngauges],yf,**kwargs)
        kwargs ={'linestyle':'-','color':'blue','label':'sea surface'}
        plot(xc[:ngauges],ys,**kwargs)
        kwargs = {'linestyle':'--','color':'r','label':'Okada'}
        fault.plot_okada(ax,displacement='vertical',kwargs=kwargs)
        legend()

    def plot_interfaces(current_data):
        from pylab import linspace, plot
        xl = linspace(xp1,xp2,100)
        yl = linspace(yp1,yp2,100)
        plot(xl,yl,'g')
        xl = linspace(xlimits[0],xlimits[1],100)
        plot(xl,0.0*xl,'b')

    def sigmatr(current_data):
        # return -trace(sigma)
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:])

    # Figure for surface and p waves
    plotfigure = plotdata.new_plotfigure(name='surface_and_p_waves', figno=1)
    plotfigure.kwargs = {'figsize':(8,8)}

    # Set axes for vertical displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-0.3,0.5]
    plotaxes.title = 'vertical displacement'
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_vertical_displacement

    # Set axes for vertical displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = '-trace(sigma)'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -1e6
    plotitem.pcolor_cmax = 1e6
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

    # Figure for surface and p waves (water scale)
    plotfigure = plotdata.new_plotfigure(name='surface_and_p_waves (water scale)', figno=2)
    plotfigure.kwargs = {'figsize':(8,8)}

    # Set axes for vertical displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimitsW
    plotaxes.ylimits = [-0.3,0.5]
    plotaxes.title = 'vertical displacement'
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_vertical_displacement

    # Set axes for vertical displacement:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimitsW
    plotaxes.ylimits = ylimitsW
    plotaxes.title = '-trace(sigma)'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -1e6
    plotitem.pcolor_cmax = 1e6
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = 'Level 3 grid patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffffff']
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0,0,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Horizontal velocity'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Vertical velocity'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 4
    plotitem.plotstyle = 'b-'



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
