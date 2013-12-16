
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

# Note: To change plotted time scale, edit frametools.py in visclaw

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
      
          # Plot outline of interface
    def aa1D(current_data):
      from pylab import linspace,plot,annotate,text
      #gcs = 2.0/200.0
      x = [-0.0085,-0.0085,0.0085,0.0085]
      y = [-1000000,10000000,1000000,-1000000]
      #y[:] = [xx - gcs for xx in y]
      plot(x,y,'k',linewidth=2.0)
      #plot(-8.0, 180000, 'vk', markersize=10) 
      #plot(-2.0, 180000, 'vk', markersize=10) 
      #plot(0.0, 180000, 'vk', markersize=10) 
      #plot(2.0, 180000, 'vk', markersize=10)
      text(-0.0075,285000,'Water',fontweight='bold',fontsize=20)
      #text(-0.8,285000,'PS',fontweight='bold',fontsize=20)
      text(-0.029,285000,'Air',fontweight='bold',fontsize=20)
      text(0.0095,285000,'Air',fontweight='bold',fontsize=20)
      
                # Plot outline of interface
    def aa1DPSIcm(current_data):
      from pylab import linspace,plot,annotate,text,xlabel,ylabel
      #gcs = 2.0/200.0
      x = [-0.85,-0.85,0.85,0.85]
      y = [-100,100,100,-100]
      #y[:] = [xx - gcs for xx in y]
      plot(x,y,'k',linewidth=2.0)
      xlabel('cm',fontsize='16')
      ylabel('psi',fontsize='16')
      xcav = [-3.0,3.0]
      ycav = [-14.334351113,-14.334351113] #Water vapour pressure for cavitation at room temp in 1atm=0 ref system
      plot(xcav,ycav,'b--')
      #plot(-8.0, 180000, 'vk', markersize=10) 
      #plot(-2.0, 180000, 'vk', markersize=10) 
      #plot(0.0, 180000, 'vk', markersize=10) 
      #plot(2.0, 180000, 'vk', markersize=10)
      text(-0.75,27,'Water',fontweight='bold',fontsize=20)
      #text(-0.8,285000,'PS',fontweight='bold',fontsize=20)
      text(-2.9,27,'Air',fontweight='bold',fontsize=20)
      text(0.95,27,'Air',fontweight='bold',fontsize=20)
      text(1.0,-13,'Vapor pressure',fontsize=15,color='blue')
      
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
        from pylab import linspace,plot,pcolor,contour,contourf,annotate,text,cm,colorbar,show,xlabel,ylabel,xticks,yticks
        import matplotlib.ticker as ticker
        from clawpack.visclaw import colormaps as cmaps
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
        # Convert to PSI
        P = P*0.000145038 - 14.6959488
        x = [-0.0085, -0.0085, 0.0085, 0.0085]
        x = [zz - 0.0 for zz in x]
        y = [0.0, 0.0085, 0.0085, 0.0]
        y2 = [-zz for zz in y]
        plot(x,y,'k',linewidth=2.0)
        plot(x,y2,'k',linewidth=2.0)
        locator = ticker.MaxNLocator(20) # if you want no more than 10 contours 
        locator.create_dummy_axis()
        #For Pa
        #locator.set_bounds(90000, 300000) 
        # For PSI
        locator.set_bounds(-20, 30) 
        levs = locator() 
        #OTHER colormap: cmap = camps.schlieren_grays
        contourf(xx,yy-0.5*dy,P, levs, alpha=.75, cmap=cm.Blues)
        contourf(xx,-(yy-0.5*dy),P, levs,alpha=.75, cmap=cm.Blues)
        cbar = colorbar(shrink=0.65)
        cbar.ax.set_xlabel('psi', fontsize='16', rotation='horizontal')
        # for PSI or Pascals
        pcolor(xx,yy-0.5*dy,P,cmap=cm.Blues, vmin=-20, vmax=30)
        pcolor(xx,-(yy-0.5*dy),P,cmap=cm.Blues, vmin=-20, vmax=30)
        # Convert axis to cm
        xxticks = np.arange(xx.min(), xx.max(), 0.01)
        labelsx = range(xxticks.size) 
        labelsx[:] = [x - 3 for x in labelsx]
        xticks(xxticks, labelsx)
        yyticks = np.arange(-yy.max(), yy.max(), 0.01)
        labelsy = range(yyticks.size)
        labelsy[:] = [x - 2 for x in labelsy]
        yticks(yyticks, labelsy)
        xlabel('cm',fontsize='16')
        ylabel('cm',fontsize='16')
        #pcolor(xx,yy-0.5*dy,P,cmap=cm.Blues, vmin=90000, vmax=300000)
        #pcolor(xx,-(yy-0.5*dy),P,cmap=cm.Blues, vmin=90000, vmax=300000)
        #colorbar()
        contour(xx,yy-0.5*dy,P,levs, colors='black', linewidth=0.5)
        contour(xx,-(yy-0.5*dy),P,levs, colors='black', linewidth=0.5)
        
         # Mirrored pressure and pressure slice
    def MirrorPressurecontour_N_Pressureslice(current_data):
        from pylab import linspace,plot,subplot,pcolor,contour,contourf,annotate,text,cm,colorbar,show,xlabel,ylabel,xticks,yticks
        from pylab import subplots,ylim, xlim, setp,  annotate,text, get, subplot2grid, axes,gca,title
        import matplotlib.ticker as ticker
        from clawpack.visclaw import colormaps as cmaps
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
        # Convert to PSI
        P = P*0.000145038 - 14.6959488
        x = [-0.0085, -0.0085, 0.0085, 0.0085]
        x = [zz - 0.0 for zz in x]
        y = [0.0, 0.0085, 0.0085, 0.0]
        y2 = [-zz for zz in y]
        
        s1 = subplot2grid((5,16), (0,0), colspan=14, rowspan=3) # subplot(211)
        plot(x,y,'k',linewidth=2.0)
        plot(x,y2,'k',linewidth=2.0)
        locator = ticker.MaxNLocator(20) # if you want no more than 10 contours 
        locator.create_dummy_axis()
        #For Pa
        #locator.set_bounds(90000, 300000) 
        # For PSI
        locator.set_bounds(-20, 30) 
        levs = locator()
        #f, (ax1, ax2) = subplots(2, sharex=True,sharey=True,subplot_kw=dict(adjustable='datalim', aspect='equal'))
        #OTHER colormap: cmap = camps.schlieren_grays
        cont1 = contourf(xx,yy-0.5*dy,P, levs, alpha=.75, cmap=cm.Blues)
        cont2 = contourf(xx,-(yy-0.5*dy),P, levs,alpha=.75, cmap=cm.Blues)
        # for PSI or Pascals
        pcolor(xx,yy-0.5*dy,P,cmap=cm.Blues, vmin=-20, vmax=30)
        pcolor(xx,-(yy-0.5*dy),P,cmap=cm.Blues, vmin=-20, vmax=30)
        s1.set_xlim([-0.03,0.03])
        s1.set_ylim([-0.02,0.02])
        # Convert axis to cm
        xxticks = np.arange(xx.min(), xx.max(), 0.01)
        labelsx = range(0)#range(xxticks.size) 
        labelsx[:] = [x - 3 for x in labelsx]
        xticks(xxticks, labelsx)
        yyticks = np.arange(-yy.max(), yy.max(), 0.01)
        labelsy = range(yyticks.size)
        labelsy[:] = [x - 2 for x in labelsy]
        yticks(yyticks, labelsy)
        #xlabel('cm',fontsize='13',fontweight='bold')
        ylabel('cm',fontsize='13',fontweight='bold')
        title('Pressure contours (2D) & pressure slice (1D)', fontweight='bold')
        #pcolor(xx,yy-0.5*dy,P,cmap=cm.Blues, vmin=90000, vmax=300000)
        #pcolor(xx,-(yy-0.5*dy),P,cmap=cm.Blues, vmin=90000, vmax=300000)
        #colorbar()
        contour(xx,yy-0.5*dy,P,levs, colors='black', linewidth=0.5)
        contour(xx,-(yy-0.5*dy),P,levs, colors='black', linewidth=0.5)
        
        s2 = subplot2grid((5,16), (3,0), colspan=14,rowspan=2) #subplot(212)
        x_slice, P_slice = xsec(current_data)
        plot(x_slice,P_slice, 'k-')
        s2.set_xlim([-3,3])
        s2.set_ylim([-20,30])
        #ax1.ylim(-20,30)
        #s2.set_xlim([-3.0,3.0])
        #s2.set_ylim([-20,30])
        #gcs = 2.0/200.0
        x = [-0.85,-0.85,0.85,0.85]
        y = [-100,100,100,-100]
        #y[:] = [xx - gcs for xx in y]
        plot(x,y,'k',linewidth=2.0)
        xlabel('cm',fontsize='13',fontweight='bold')
        ylabel('psi',fontsize='13',fontweight='bold')
        #title('Pressure slice')
        xcav = [-3.0,3.0]
        ycav = [-14.334351113,-14.334351113] #Water vapour pressure for cavitation at room temp in 1atm=0 ref system
        plot(xcav,ycav,'b--')
        #plot(-8.0, 180000, 'vk', markersize=10) 
        #plot(-2.0, 180000, 'vk', markersize=10) 
        #plot(0.0, 180000, 'vk', markersize=10) 
        #plot(2.0, 180000, 'vk', markersize=10)
        text(-0.75,25,'Water',fontweight='bold',fontsize=13)
        #text(-0.8,285000,'PS',fontweight='bold',fontsize=20)
        text(-2.9,25,'Air',fontweight='bold',fontsize=13)
        text(0.95,25,'Air',fontweight='bold',fontsize=13)
        text(1.2,-13,'Vapor pressure',fontweight='bold',fontsize=13,color='blue')
        
        s3 = subplot2grid((5,16), (0,15), rowspan=5)
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(gca())
        cax = divider.append_axes("right", "5%", pad="3%")
        cax.axis('off')
        cbar = colorbar(cont1,cax=s3) #,orientation='horizontal') #, shrink=0.99) #,location='eastoutside') #orientation='horizontal')
        cbar.ax.set_xlabel('psi', fontsize='13', fontweight='bold', rotation='horizontal')
        

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
    #plotitem.pcolor_cmap = 'schlieren'
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
    # Axes for m vs Pa or cm vs PSI
    #plotaxes.xlimits = [-0.03,0.03] #[-3,3] #[-8.5,16] #'auto' -16
    #plotaxes.ylimits = [0.00000,300000]
    plotaxes.xlimits = [-3.0,3.0]
    plotaxes.ylimits = [-20,30]
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
        # Convert to Psi and centimeters
        P_slice = P_slice*0.000145038 - 14.6959488
        x_slice = 100*x_slice
        return x_slice, P_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = '-kx'
    plotitem.kwargs = {'markersize':3}
    
    plotaxes.afteraxes = aa1DPSIcm
    
    
    # Pressure contour(2D) and pressure slice(1D) in one figure
    plotfigure = plotdata.new_plotfigure(name='Contour & Slice', figno=7)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-3,3] #[-3,3] #[-8.5,16] #'auto' -16
    plotaxes.ylimits = [-20,30]#[-5,5]
    plotaxes.title = 'Pressure'
    #plotaxes.scaled = True      # so aspect ratio is 1
    
    plotaxes.afteraxes = MirrorPressurecontour_N_Pressureslice    


   ### Figure for Pressure contour and pressure slice in one figure
    ### -------------------
    
    #plotfigure = plotdata.new_plotfigure(name='Pressure contour and slice', figno=7)

    ## Set up for axes in this sub-figure:
    #plotaxes = plotfigure.new_plotaxes()
    #plotaxes.axescmd = 'subplot(211)'
    #plotaxes.xlimits = 'auto' #[-3,3] #[-8.5,16] #'auto' -16
    #plotaxes.ylimits = [-0.02,0.02]#[-5,5]
    #plotaxes.title = 'Pressure'
    ##plotaxes.scaled = True      # so aspect ratio is 1
    
    #plotaxes.afteraxes = MirrorPressure_contour2

    ## Set up for axes in this subfigure:
    #plotaxes = plotfigure.new_plotaxes()
    #plotaxes.axescmd = 'subplot(212)'
    #plotaxes.xlimits = 'auto' #[-3,3] #[-8.5,16] #'auto' -16
    #plotaxes.ylimits = [-20,30]
    #plotaxes.title = 'Pressure slice'
    #plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    #plotitem.map_2d_to_1d = xsec
    #plotitem.plotstyle = '-kx'
    #plotitem.kwargs = {'markersize':3}
    
    #plotaxes.afteraxes = aa1DPSIcm


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
    plotdata.print_framenos = 'all'          # list of frames to print # [54,90,95,110,126,171,186,199,290]
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
