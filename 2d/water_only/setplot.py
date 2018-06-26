
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np

csig = 10
cdivcurl = 1e-5

h0 = 4500.  # ocean depth
cw = 1500.  # speed of sound in water
g = 9.81
rho = 1025. # density of water
rhog = rho*g
period = 4 * h0 / cw # time for signal transit depth 4 times
ybot = -h0

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
    
    def plot_interfaces(current_data):
        return
        from pylab import linspace, plot
        xl = linspace(-40e3, 40e3, 101)
        #yl = -4000. + 0*xl
        yl = -10000. + 6000*(xl+40e3)/80e3
        plot(xl,yl,'k')
    

    def sigmatr(current_data):
        q = current_data.q
        return (q[0,:,:] + q[1,:,:])

    def pressure(current_data):
        """ pressure in units of meters of water"""
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:])/(2.*rhog)

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

    # Figure for trace(sigma) and sigma_12 side by side
    plotfigure = plotdata.new_plotfigure(name='P and S waves', figno=11)
    plotfigure.kwargs = {'figsize':(13,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.6,.35,.35])' # 'subplot(221)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'pressure'
    #plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = pressure
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -csig
    plotitem.pcolor_cmax = csig
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.5,.6,.45,.35])' # 'subplot(222)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'delta' #'sigma_12'
    #plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 5 #2
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -5. #-csig
    plotitem.pcolor_cmax =  5. #csig
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.1,.35,.35])' # 'subplot(223)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'div(u)'
    #plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = div
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -cdivcurl
    plotitem.pcolor_cmax = cdivcurl
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]


    # Figure for curl:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.5,.1,.45,.35])' # 'subplot(224)'
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'curl(u)'
    #plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = curl
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -cdivcurl
    plotitem.pcolor_cmax = cdivcurl
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [0]



    # Figure for surface displacement
    plotfigure = plotdata.new_plotfigure(name='surface displacement', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = [0,2]
    plotaxes.ylimits = [-15,15]
    plotaxes.title = 'surface displacement'

    def fixup(current_data):
        from pylab import grid
        from clawpack.visclaw import legend_tools
        grid(True)
        labels = ['surface displacement','bottom displacement',
                  'bottom pressure diff']
        colors = ['b','g','r']
        linestyles = '-'
        markers = ''
        legend_tools.add_legend(labels,colors,linestyles,markers,
                loc='upper left')
        
    plotaxes.afteraxes = fixup
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    #plotitem.show = False

    def xsec(current_data):
        # Return x value and surface eta at this point, along y=0
        global surf_min, surf_max
        from pylab import find,ravel,nan
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q

        ij = find((y <= 0.) & (y > -dy))
        x_slice = ravel(x)[ij]

        # extrapolate displacement from cell centers to top surface:
        eta_slice_1 = ravel(q[5,:,:])[ij]
        ij = find((y <= -dy) & (y > -2*dy))
        eta_slice_2 = ravel(q[5,:,:])[ij]
        eta_slice = 1.5*eta_slice_1 - 0.5*eta_slice_2

        #if current_data.level < 3:
            #eta_slice = nan*eta_slice
        #if len(eta_slice) > 0:
            #print('min,max = ',eta_slice.min(),eta_slice.max())
        try:
            surf_min = min(surf_min,eta_slice.min())
            surf_max = max(surf_max,eta_slice.max())
        except:
            pass
        if current_data.level != 2:
            eta_slice *= nan
        return x_slice, eta_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = 'b-'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':5}
    plotitem.amr_data_show = [1]  

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    #plotitem.show = False

    def xsec_bottom_delta(current_data):
        from pylab import find,ravel,nan
        
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q

        ij = find((y <= ybot+dy) & (y > ybot))
        x_slice = ravel(x)[ij]

        # bottom displacement:
        delta_slice = ravel(q[5,:,:])[ij]
        if current_data.level != 2:
            delta_slice *= nan
        return x_slice, delta_slice

    plotitem.map_2d_to_1d = xsec_bottom_delta
    plotitem.plotstyle = 'g-'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':5}
    plotitem.amr_data_show = [1]  

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    #plotitem.show = False

    def xsec_bottom_pressure(current_data):
        from pylab import find,ravel,nan
        
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q

        ij = find((y <= ybot+dy) & (y > ybot))
        x_slice = ravel(x)[ij]

        # bottom pressure in units of meters of water, corrected for delta:
        p_slice = -(ravel(q[0,:,:])[ij] + ravel(q[1,:,:])[ij]) / (2*rhog) \
                    - ravel(q[5,:,:])[ij]
        if current_data.level != 2:
            p_slice *= nan
        return x_slice, p_slice

    plotitem.map_2d_to_1d = xsec_bottom_pressure
    plotitem.plotstyle = 'r-'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':5}
    plotitem.amr_data_show = [1]  



    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = [0,2]
    #plotaxes.ylimits = [0,1]
    plotaxes.title = 'Grid patches'
    #plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffffff']
    plotitem.amr_celledges_show = [1,0,0]
    plotitem.amr_patchedges_show = [1]


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    def aa_gauge(current_data):
        from pylab import grid, tight_layout
        grid(True)
        tight_layout()

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Horizontal velocity'
    plotaxes.afteraxes = aa_gauge

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Vertical velocity'
    plotaxes.afteraxes = aa_gauge

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 4
    plotitem.plotstyle = 'b-'


    plotfigure = plotdata.new_plotfigure(name='gauge plot sigma', figno=301, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    def pressure_diff(current_data):
        q = current_data.q
        p = (q[0,:]+q[1,:]) / (2*rhog)
        delta = q[5,:]
        return p - delta

    def filtered_pressure_diff(current_data):
        from pylab import find,len,sum,zeros,hstack
        q = current_data.q
        import pdb; pdb.set_trace()
        p = (q[0,:]+q[1,:]) / (2*rhog)
        delta = q[5,:]
        d = p - delta
        t = q.t
        fpd = []
        for k1 in range(len(t)):
            t1 = t[k1]
            try:
                k2 = find(t > t1+period).min()
            except:
                break
            if k1==0:
                kmid = int((k1+k2)/2.)
            fpd.append(sum(d[k1:k2])/(k2-k1+1))
        fpd = hstack((zeros(kmid),fpd))
        kz = len(t) - len(fpd)
        if kz>0:
            fpd = hstack((fpd,zeros(kz)))
        return fpd

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'pressure difference'
    plotaxes.afteraxes = aa_gauge

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = pressure_diff
    plotitem.plotstyle = 'b-'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = filtered_pressure_diff
    plotitem.plotstyle = 'k-'



    def displacement(current_data):
        q = current_data.q
        delta = q[5,:]
        return delta

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'displacement'
    plotaxes.afteraxes = aa_gauge

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = displacement
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

    
