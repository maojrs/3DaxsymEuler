""" 
Module to set up run time parameters for Clawpack -- classic code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
import numpy as np

#------------------------------
def setrun(claw_pkg='classic'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "classic" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    from clawpack.clawutil import data 
    
    
    assert claw_pkg.lower() == 'classic',  "Expected claw_pkg = 'classic'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('gammagas',     1.4,  'gamma for ideal gas')
    probdata.add_param('gammaplas',    1.1, 'gamma est. for polystirene')
    #probdata.add_param('gammawat',     2.5,  'gamma for water') #7.15 before 
    probdata.add_param('gammawat',     7.15,  'gamma for water')
    #probdata.add_param('gammawat',     1.9,  'gamma for water')



    probdata.add_param('pinfgas',     0.0,  'pinf for stiffend gas/plastic')
    # pinfplas Calculated with c^2=gamma*(p+pinf)/rho to make c =2240m/s (polyestirene speed of sound), 
    # (c= 1484 in water). Values from water obtained fron kirsten's paper
    probdata.add_param('pinfplas',    4789425947.72,  'pinf for stiffend gas/plastic') 
    #probdata.add_param('pinfplas',    1000000000.0,  'pinf for stiffend gas/plastic') 
    #probdata.add_param('pinfwat',       880801075.0,  'pinf for stiffend water') #Previously 300000000.0
    probdata.add_param('pinfwat',        300000000.0,  'pinf for stiffend water') #Previously 300000000.0
    #probdata.add_param('pinfwat',           0.0,  'pinf for stiffend water') #Previously 300000000.0



    probdata.add_param('omegas',     0.0,  'omega (specific excluded volume) for stiffend gas/plastic')
    probdata.add_param('omeplas',     0.0,  'omega (specific excluded volume) for stiffend gas/plastic')
    probdata.add_param('omewat',     0.0,  'omega (specific excluded volume) for stiffend water')
    probdata.add_param('rhog',     1.0,  'air density in kg/m^3')
    probdata.add_param('rhop',     1050.0,  'polystirene density in kg/m^3')
    probdata.add_param('rhow',     1000.0,  'water density in kg/m^3')
    #probdata.add_param('rhow',    1.0,  'water density in kg/m^3')
    
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim
    
    # Lower and upper edge of computational domain:
    clawdata.lower[0] = -0.03 #-7.0           # xlower
    clawdata.upper[0] = 0.03                  # xupper
    clawdata.lower[1] = 0.000000e+00          # ylower
    clawdata.upper[1] = 0.020000e+00          # yupper
    
    # Number of grid cells:
    clawdata.num_cells[0] = 400 #400 #800 #1400 #1400 #3000 #700 #1500      # mx
    clawdata.num_cells[1] = 40 #40 #80 #40 #40 # 40  #20      # my
    

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 4

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 5
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 0
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.000000
    

    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.qNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.q0006'   # File to use for restart data
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
 
    clawdata.output_style = 1
 
    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 300 #500 #400 #300 #200 #50
        clawdata.tfinal = 0.0002 #0.00020 #0.0002 #0.0015 #0.0002 #0.0003 #0.0004 #0.03 #0.01 #0.050000
        clawdata.output_t0 = True  # output at initial (or restart) time?
        
    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  [0., 0.1]
 
    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 2
        clawdata.total_steps = 4
        clawdata.output_t0 = True  # output at initial (or restart) time?
        

    clawdata.output_format == 'ascii'      # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'all'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0
    

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==False: fixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True
    
    # Initial time step for variable dt.  
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 1.000000e-01
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1.000000e+99
    
    # Desired Courant number if variable dt used 
    clawdata.cfl_desired = 0.4000 #0.200000 #0.900000
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 0.490000 #0.350000 #1.000000
    
    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 1000 #500


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR) !godunov or strang(dim splittin), unsplit (for transverse)
    clawdata.dimensional_split = 'godunov'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 0
    
    
    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = [1, 1, 1]
    
    clawdata.use_fwaves = False    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 1
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2
    
    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => soli                print*, i,ixy,mx, pstard wall for systems where q(2) is normal velocity
    
    clawdata.bc_lower[0] = 'user' #''user'  # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'extrap'   # at yupper
                  
    return rundata

    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
    
