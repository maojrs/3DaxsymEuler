
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import os
import numpy as np

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
    
    
    # Plot outline of interface
    def aa(current_data):
      from pylab import linspace,plot,annotate,text
      #gcs = 2.0/200.0
      # Wall pwidth
      #x = [-2.0,-2.0,2.0,2.0]
      #y = [-6,6,6,-6]
      x = [-0.85, -0.85, 0.85, 0.85]
      y = [0.0, 0.85, 0.85, 0.0]
      #y[:] = [xx - gcs for xx in y]
      plot(x,y,'k',linewidth=2.0)

    # Figure for Density
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Density'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    #plotitem.pcolor_cmin = 0.8
    #plotitem.pcolor_cmax = 3.0
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = 1.0
    plotitem.pcolor_cmax = 2.0
    plotitem.show = True       # show on plot?
    
    plotaxes.afteraxes = aa
    
        # Figure for momentum x
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Momentum x', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Momentum x'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 1
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 160.0
    plotitem.show = True       # show on plot?
    
    plotaxes.afteraxes = aa
    
            # Figure for momentum y
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Momentum y', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Momentum y'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 2
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 160.0
    plotitem.show = True       # show on plot?
    
    plotaxes.afteraxes = aa
    
      # Figure for Energy
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Energy', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Energy'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 3
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    plotitem.pcolor_cmin = 200000
    plotitem.pcolor_cmax = 400000
    
    plotaxes.afteraxes = aa
    
    # Figure for Pressure
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=4)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-3,3] #[-8.5,16] #'auto' -16
    plotaxes.ylimits = [-5,5]
    plotaxes.title = 'Pressure'
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.pcolor_cmin = 90000
    plotitem.pcolor_cmax = 200000

    def Pressure(current_data):
        q = current_data.q   # solution when this function called
        aux = current_data.aux
        gamma = aux[0,:,:]
        gamma1 = aux[0,:,:] - 1.0
        pinf = aux[1,:,:]
        omega = aux[2,:,:]
        rho = q[0,:,:]           # density
        momx = q[1,:,:]           # momentum
        momy = q[2,:,:]
        ene = q[3,:,:]           # energy
        P = gamma1*(ene - 0.5*(momx*momx + momy*momy)/rho)/(1.0 - omega*rho)
        P = P - gamma*pinf 
        return P

    plotitem.plot_var = Pressure  # defined above
    #plotitem.plotstyle = '-o'
    #plotitem.color = 'r'
    
    plotaxes.afteraxes = aa

    

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

    return plotdata

    
