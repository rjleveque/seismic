
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np
from mapc2p import mapc2p
import os

csig = 1e7
cv = 0.5 

outdir2 = os.path.abspath('../water_layer_averaged/_output')

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    def fixsubplots(current_data):
        from pylab import tight_layout, figure
        figure(31)
        tight_layout()

    #plotdata.afterframe = fixsubplots
    
    def plot_interfaces(current_data):
        from pylab import linspace, plot
        xl = linspace(-40e3, 40e3, 101)
        #yl = -4000. + 0*xl
        yl = -10000. + 6000*(xl+40e3)/80e3
        plot(xl,yl,'k')
    

    def sigmatr(current_data):
        q = current_data.q
        return q[0,:,:] + q[1,:,:]

    def div(current_data):
        from numpy import array,zeros,hstack,vstack
        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        mx, my = u.shape
        if (mx<3) or (my<3):
            d = zeros(u.shape)
            return d
        dx, dy = current_data.dx, current_data.dy
        I = array(range(1,mx-1))
        J = array(range(1,my-1))
        ux = (u[I+1,:][:,J] - u[I-1,:][:,J]) / (2*dx)
        vy = (v[:,J+1][I,:] - v[:,J-1][I,:]) / (2*dy)
        dint = ux + vy
        
        #zx = zeros((mx-2,1))
        #zy = zeros((1,my))
        #d = vstack((zy, hstack((zx, ux+vy, zx)), zy))
        
        d0 = dint[:,0]
        d1 = dint[:,-1]
        d2 = vstack((d0, dint.T, d1)).T
        d0 = d2[0,:]
        d1 = d2[-1,:]
        d = vstack((d0,d2,d1))      
        return d

    def curl(current_data):
        from numpy import array,zeros,hstack,vstack
        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        mx, my = u.shape
        if (mx<3) or (my<3):
            c = zeros(u.shape)
            return c
        dx, dy = current_data.dx, current_data.dy
        I = array(range(1,mx-1))
        J = array(range(1,my-1))
        vx = (v[I+1,:][:,J] - v[I-1,:][:,J]) / (2*dx)
        uy = (u[:,J+1][I,:] - u[:,J-1][I,:]) / (2*dy)
        cint = vx - uy

        c0 = cint[:,0]
        c1 = cint[:,-1]
        c2 = vstack((c0, cint.T, c1)).T
        c0 = c2[0,:]
        c1 = c2[-1,:]
        c = vstack((c0,c2,c1))      
        return c

    # Figure for trace(sigma), sigma_12, v side by side
    plotfigure = plotdata.new_plotfigure(name='Comparisons', figno=31)
    plotfigure.kwargs = {'figsize':(13,6)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(231)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'trace(sigma)'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -csig
    plotitem.pcolor_cmax = csig
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(234)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'trace(sigma)'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.outdir = outdir2
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -csig
    plotitem.pcolor_cmax = csig
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]

    # sig12

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(232)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'sigma_12'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 2
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -csig
    plotitem.pcolor_cmax = csig
    plotitem.add_colorbar = False
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(235)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'sigma_12'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.outdir = outdir2
    plotitem.plot_var = 2
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -csig
    plotitem.pcolor_cmax = csig
    plotitem.add_colorbar = False
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]

    # v

    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.show = False
    plotaxes.axescmd = 'subplot(233)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'v'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 4
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -cv
    plotitem.pcolor_cmax = cv
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p


    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.show = False
    plotaxes.axescmd = 'subplot(236)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'v'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.outdir = outdir2
    plotitem.plot_var = 4
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -cv
    plotitem.pcolor_cmax = cv
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]



    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'Grid patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffffff']
    plotitem.amr_celledges_show = [1,0,0]
    plotitem.amr_patchedges_show = [1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p


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


    plotfigure = plotdata.new_plotfigure(name='gauge plot sigma', figno=301, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'sigma11'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 1
    plotitem.plotstyle = 'b-'

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'sigma22'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 2
    plotitem.plotstyle = 'b-'


    #-----------------------------------------
    # Figures for gauge comparison
    #-----------------------------------------
    otherfigure = plotdata.new_otherfigure(name='v gauges',
                    fname='v_gauges.png')
    otherfigure = plotdata.new_otherfigure(name='sig11 gauges',
                    fname='sig11_gauges.png')



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

    
