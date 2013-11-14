
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
      x = [-0.0085, -0.0085, 0.0085, 0.0085]
      x = [y - 0.0 for y in x]
      y = [0.0, 0.0085, 0.0085, 0.0]
      xborder = [-0.05, 0.05, 0.05, -0.05]
      yborder = [0.0, 0.0]
      #y[:] = [xx - gcs for xx in y]
      plot(x,y,'k',linewidth=2.0)
      
    #When only using SGEOS
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
        P = gamma1*(ene - 0.5*(momx*momx + momy*momy)/rho) #/(1.0 - omega*rho)
        P = P - gamma*pinf
        return P
      
    # Mirrored pressure pcolor plor
    def MirrorPressure_pcolor(current_data):
        from pylab import linspace,plot,pcolor,annotate,text,cm
        xx = current_data.x
        yy = current_data.y
        dy = abs(yy[1,1] - yy[1,2])
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
        P = gamma1*(ene - 0.5*(momx*momx + momy*momy)/rho) #/(1.0 - omega*rho)
        P = P - gamma*pinf
        #xborder = [-0.05, 0.05, 0.05, -0.05]
        #yborder = [0.0, 0.0, 0.02, 0.02]
        #plot(xborder,yborder,'k',linewidth=2.0)
        x = [-0.0085, -0.0085, 0.0085, 0.0085]
        x = [zz - 0.0 for zz in x]
        y = [0.0, 0.0085, 0.0085, 0.0]
        y2 = [-zz for zz in y]
        #xborder = [-0.05, 0.05, 0.05, -0.05]
        #yborder = [-0.02, -0.02, 0.]
        plot(x,y,'k',linewidth=2.0)
        plot(x,y2,'k',linewidth=2.0)
        
        pcolor(xx,yy-0.5*dy,P,cmap=cm.jet, vmin=90000, vmax=300000)
        pcolor(xx,-(yy-0.5*dy),P,cmap=cm.jet, vmin=90000, vmax=300000)
        
         # Mirrored pressure pcolor plor
    def MirrorPressure_contour(current_data):
        from pylab import linspace,plot,contour,contourf,annotate,text,cm,colorbar,show
        import matplotlib.ticker as ticker 
        xx = current_data.x
        yy = current_data.y
        dy = abs(yy[1,1] - yy[1,2])
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
        P = gamma1*(ene - 0.5*(momx*momx + momy*momy)/rho) #/(1.0 - omega*rho)
        P = P - gamma*pinf
        #xborder = [-0.05, 0.05, 0.05, -0.05]
        #yborder = [0.0, 0.0, 0.02, 0.02]
        #plot(xborder,yborder,'k',linewidth=2.0)
        x = [-0.0085, -0.0085, 0.0085, 0.0085]
        x = [zz - 0.0 for zz in x]
        y = [0.0, 0.0085, 0.0085, 0.0]
        y2 = [-zz for zz in y]
        #xborder = [-0.05, 0.05, 0.05, -0.05]
        #yborder = [-0.02, -0.02, 0.]
        plot(x,y,'k',linewidth=2.0)
        plot(x,y2,'k',linewidth=2.0)
        
        locator = ticker.MaxNLocator(20) # if you want no more than 10 contours 
        locator.create_dummy_axis()
        locator.set_bounds(100000, 300000) 
        levs = locator() 
        
        contourf(xx,yy-0.5*dy,P, levs, alpha=.75, cmap=cm.Blues)
        contourf(xx,-(yy-0.5*dy),P, levs,alpha=.75, cmap=cm.Blues)
        colorbar()
        contour(xx,yy-0.5*dy,P,levs, colors='black', linewidth=0.5)
        contour(xx,-(yy-0.5*dy),P,levs, colors='black', linewidth=0.5)

    # Figure for Density
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-0.03,0.03] #'auto'
    plotaxes.ylimits = [-0.05,0.05]#'auto'
    plotaxes.title = 'Density'
    #plotaxes.scaled = True      # so aspect ratio is 1

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
    plotaxes.xlimits = [-0.03,0.03] #'auto'
    plotaxes.ylimits = [-0.05,0.05] #'auto'
    plotaxes.title = 'Momentum x'
    #plotaxes.scaled = True      # so aspect ratio is 1

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
    plotaxes.xlimits = [-0.03,0.03]#'auto'
    plotaxes.ylimits = [-0.05,0.05]#'auto'
    plotaxes.title = 'Momentum y'
    #plotaxes.scaled = True      # so aspect ratio is 1

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
    plotaxes.xlimits = [-0.03,0.03]#'auto'
    plotaxes.ylimits = [-0.05,0.05]#'auto'
    plotaxes.title = 'Energy'
    #plotaxes.scaled = True      # so aspect ratio is 1

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
    plotaxes.xlimits = [-0.03,0.03] #[-3,3] #[-8.5,16] #'auto' -16
    plotaxes.ylimits = [-0.02,0.04]#[-5,5]
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True      # so aspect ratio is 1
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.pcolor_cmin = 90000
    plotitem.pcolor_cmax = 300000
    plotitem.plot_var = Pressure  # defined above
    #plotitem.plotstyle = '-o'
    #plotitem.color = 'r'

    plotaxes.afteraxes = aa
    
     ## Figure for Pressure Rotated
    ## -------------------

    #plotfigure = plotdata.new_plotfigure(name='Pressure', figno=4)

    ##yy = plotdata.current_data.y
    ##dy = abs(yy[1,1] - yy[1,2])

    ## Set up for axes in this figure:
    #plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = [-0.03,0.03] #[-3,3] #[-8.5,16] #'auto' -16
    #plotaxes.ylimits = [-0.02,0.02]#[-5,5]
    #plotaxes.title = 'Pressure'
    #plotaxes.scaled = True      # so aspect ratio is 1
    
    #plotaxes.afteraxes = MirrorPressure_pcolor
    
    ##plotaxes.afteraxes = aa
    
        # Figure for Pressure Contour
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Pressure Contour', figno=5)

        # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-0.03,0.03] #[-3,3] #[-8.5,16] #'auto' -16
    plotaxes.ylimits = [-0.02,0.02]#[-5,5]
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True      # so aspect ratio is 1
    
    plotaxes.afteraxes = MirrorPressure_contour
    
        # Figure for Pressure slice
    # -------------------
    
    plotfigure = plotdata.new_plotfigure(name='Pressure slice', figno=6)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-0.03,0.03] #[-3,3] #[-8.5,16] #'auto' -16
    plotaxes.ylimits = [90000,300000]
    plotaxes.title = 'Pressure slice'
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def xsec(current_data):
        # Return x value and surface eta at this point, along y=0
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        aux = current_data.aux

        ij = find((y <= dy/2.) & (y > -dy/2.))
        x_slice = ravel(x)[ij]
        gamma_slice = ravel(aux[0,:,:])[ij]
        pinf_slice = ravel(aux[1,:,:])[ij]
        rho_slice = ravel(q[0,:,:])[ij]
        momx_slice = ravel(q[1,:,:])[ij]
        momy_slice = ravel(q[2,:,:])[ij]
        ene_slice = ravel(q[3,:,:])[ij]
        P_slice = (gamma_slice - 1.0)*(ene_slice - 0.5*(momx_slice**2 + momy_slice**2)/rho_slice)
        P_slice = P_slice - gamma_slice*pinf_slice
        return x_slice, P_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = '-kx'
    plotitem.kwargs = {'markersize':3}


      ## Figure for Something slice
    ## -------------------
    
    #plotfigure = plotdata.new_plotfigure(name='Something slice', figno=6)
    ## Set up for axes in this figure:
    #plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = [-0.03,0.03] #[-3,3] #[-8.5,16] #'auto' -16
    ##plotaxes.ylimits = [90000,300000]
    #plotaxes.title = 'Something slice'
    #plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    #def xsec(current_data):
        ## Return x value and surface eta at this point, along y=0
        #from pylab import find,ravel
        #x = current_data.x
        #y = current_data.y
        #dy = current_data.dy
        #q = current_data.q
        #aux = current_data.aux

        #ij = find((y <= dy/2.) & (y > -dy/2.))
        #x_slice = ravel(x)[ij]
        #ene_slice = ravel(q[3,:,:])[ij]
        #return x_slice, ene_slice

    #plotitem.map_2d_to_1d = xsec
    #plotitem.plotstyle = '-kx'
    #plotitem.kwargs = {'markersize':3}
    

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

    
