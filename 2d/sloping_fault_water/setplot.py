
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
from clawpack.clawutil.data import ClawData
from mapping import Mapping
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
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

    mapping = Mapping(fault, probdata.water_scaling)
    fault_width = mapping.fault_width
    xcenter = mapping.xcenter
    ycenter = mapping.ycenter
    xp1 = mapping.xp1
    xp2 = mapping.xp2
    yp1 = mapping.yp1
    yp2 = mapping.yp2

    xlimits = [xcenter-0.5*probdata.domain_width,xcenter+0.5*probdata.domain_width]
    ylimits = [-probdata.domain_depth,probdata.water_depth]
    xlimitsW = [xp1-10.0*probdata.water_depth,xp2+10.0*probdata.water_depth]
    ylimitsW = [-5.0*probdata.water_depth,probdata.water_depth]
    abl_depth = probdata.abl_depth
    xlimits_trunc = [xlimits[0]+abl_depth, xlimits[1]-abl_depth]
    ylimits_trunc = [ylimits[0]+abl_depth, ylimits[1]]

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'
    plotdata.parallel = True

    def plot_interfaces(current_data):
        from pylab import linspace, plot
        xl = linspace(xp1,xp2,100)
        yl = linspace(yp1,yp2,100)
        plot(xl,yl,'g')
        xl = linspace(xlimits[0],xlimits[1],1000)
        plot(xl,0.0*xl,'b')


    def sigmatr(current_data):
        # return -trace(sigma)
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:])

    def slip_direction_vel(current_data):
        # return vel dot tau, where tau is tangent to fault
        tau_x = (xp2 - xp1)/fault_width
        tau_y = (yp2 - yp1)/fault_width
        tau_x = np.where(current_data.y > ycenter, -tau_x, tau_x)
        tau_y = np.where(current_data.y > ycenter, -tau_y, tau_y)
        u = current_data.q[3,:,:]
        v = current_data.q[4,:,:]
        return u*tau_x + v*tau_y

    # Figure for waves
    plotfigure = plotdata.new_plotfigure(name='trace', figno=1)
    plotfigure.kwargs = {'figsize':(10,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
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

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = 'slip-direction-velocity'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = slip_direction_vel
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.1
    plotitem.pcolor_cmax = 0.1
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

    plotfigure = plotdata.new_plotfigure(name='trace (water scale)', figno=2)
    plotfigure.kwargs = {'figsize':(10,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
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


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimitsW
    plotaxes.ylimits = ylimitsW
    plotaxes.title = 'slip-direction-velocity'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = slip_direction_vel
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -1e-1
    plotitem.pcolor_cmax = 1e-1
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
    plotaxes.title = 'Level 4 grid patches'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffffff']
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0,0,0,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

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
#    plotdata.parallel = True

    return plotdata
