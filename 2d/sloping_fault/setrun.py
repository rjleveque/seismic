"""
Module to set up run time parameters for Clawpack -- amrclaw code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np

#------------------------------
def setrun(claw_pkg='amrclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "amrclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data


    assert claw_pkg.lower() == 'amrclaw',  "Expected claw_pkg = 'amrclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    # Specify which topography to use (set topography to False for flat surface)
    topography = True
    if (topography):
        import imp
        make_topo_and_grid = imp.load_source('make_topo_and_grid',\
        '../../tsunami/1d/pwlinear1/make_topo_and_grid.py')
        from make_topo_and_grid import get_seafloor_parameters

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    # Sample setup to write one line to setprob.data ...
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    if (topography):
        xlower, xslope_lower, xshelf_lower, xbeach_lower, xshore_lower, xupper, \
            ybasin_lower, yshelf_lower, ybeach_lower, yshore_lower \
             = get_seafloor_parameters()
        xlower = -250.0e3
        xupper = 150.0e3
        fault_offset = xshelf_lower
        domain_depth = 200e3
    else:
        xlower = -250.0e3
        xslope_lower = -65.0e3
        xshelf_lower = -45.0e3
        xbeach_lower = -5.0e3
        xshore_lower = 0.0
        xupper = 150.0e3
        ybasin_lower = 0.0
        yshelf_lower = 0.0
        ybeach_lower = 0.0
        yshore_lower = 0.0
        fault_offset = -45.0e3
        domain_depth = 200.0e3

    probdata.add_param('fault_center', fault_offset + 25.0e3, 'center of fault')
    probdata.add_param('fault_width', 50735, 'width of fault')
    probdata.add_param('fault_dip', 0.17, 'angle of fault dip')
    probdata.add_param('fault_depth', 19.0e3, 'depth of fault')
    probdata.add_param('basin_center', 0.5*(xlower + xslope_lower), 'center of basin')
    probdata.add_param('basin_width', xslope_lower-xlower, 'width of basin')
    probdata.add_param('basin_elevation', ybasin_lower, 'elevation of basin')
    probdata.add_param('slope_center', 0.5*(xslope_lower + xshelf_lower), 'center of slope')
    probdata.add_param('slope_width', xshelf_lower - xslope_lower, 'width of slope')
    probdata.add_param('shelf_center', 0.5*(xshelf_lower + xbeach_lower), 'center pf shelf')
    probdata.add_param('shelf_width', xbeach_lower - xshelf_lower, 'width of shelf')
    probdata.add_param('shelf_elevation', yshelf_lower, 'elevation of shelf')
    probdata.add_param('beach_center',0.5*(xbeach_lower + xshore_lower), 'center of beach')
    probdata.add_param('beach_width', xshore_lower - xbeach_lower, 'width of beach')
    probdata.add_param('shore_elevation', yshore_lower, 'elevation of shore')



    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Number of grid cells:
    num_cells_fault = 10
    dx = probdata.fault_width/num_cells_fault
    ## specify dy using dx
    num_cells_above = np.rint(probdata.fault_depth/dx)
    dy = probdata.fault_depth/num_cells_above
    num_cells_left = int(np.ceil((probdata.fault_center - 0.5*probdata.fault_width - xlower)/dx))
    num_cells_right = int(np.ceil((xupper - probdata.fault_center - 0.5*probdata.fault_width)/dx))
    clawdata.num_cells[0] = num_cells_left + num_cells_fault + num_cells_right
    clawdata.num_cells[1] = int(np.ceil(domain_depth/dy)) # my

    # Lower and upper edge of computational domain:
    ## note the size of domain is likely expanded here
    clawdata.lower[0] = probdata.fault_center - 0.5*probdata.fault_width - num_cells_left*dx   # xlower
    clawdata.upper[0] = probdata.fault_center + 0.5*probdata.fault_width + num_cells_right*dx     # xupper
    clawdata.lower[1] = -clawdata.num_cells[1]*dy       # ylower
    clawdata.upper[1] = 0.0          # yupper
    domain_depth = clawdata.upper[1] - clawdata.lower[1]

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 5

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 13

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 12


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
        clawdata.num_output_times = 50
        clawdata.tfinal = 100.0
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  list(np.linspace(0,5,11)) + \
            range(6,61)

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 40
        clawdata.output_t0 = True  # output at initial (or restart) time?


    clawdata.output_format = 'ascii'      # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
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
    clawdata.dt_initial = 1e-6

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1.000000e+99

    # Desired Courant number if variable dt used
    clawdata.cfl_desired = 0.900000
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.000000

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 1000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2


    # Number of waves in the Riemann solution:
    clawdata.num_waves = 4

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = ['mc', 'mc', 'mc', 'mc']

    clawdata.use_fwaves = False    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 0


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'user'   # at yupper



    # ---------------
    # Gauges:
    # ---------------
    gauges = rundata.gaugedata.gauges
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    # top edge:
    xgauges = np.linspace(clawdata.lower[0]+1, clawdata.upper[0]-1, 100)
    for gaugeno,x in enumerate(xgauges):
        gauges.append([gaugeno,x,clawdata.upper[1]-1,0,1e10])
    ## bottom edge:
    #for gaugeno,x in enumerate(xgauges):
    #    gauges.append([200+gaugeno,x,clawdata.lower[1]+1000,0,1e10])

    ## above fault plane:
    #xgauges = np.linspace(probdata.fault_center-0.5*probdata.fault_width,
    #                        probdata.fault_center+0.5*probdata.fault_width)
    #for gaugeno,x in enumerate(xgauges):
    #    gauges.append([300+gaugeno,x,-probdata.fault_depth+1,0,1e10])

    # below fault plane:
    #for gaugeno,x in enumerate(xgauges):
    #    gauges.append([400+gaugeno,x,-probdata.fault_depth-1,0,1e10])


    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 1

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.
      clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5



    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 4

    # List of refinement ratios at each level (length at least
    # amr_level_max-1)
    amrdata.refinement_ratios_x = [4,4,2]
    amrdata.refinement_ratios_y = [4,4,2]
    amrdata.refinement_ratios_t = [4,4,2]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one
    # of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center', 'center', 'center', 'center', 'center', \
        'center', 'center', 'center', 'center', 'center', 'center', \
        'capacity','yleft']



    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.  # Richardson tolerance

    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 0.0001  # tolerance used in this routine
    # User can modify flag2refine to change the criterion for flagging.
    # Default: check maximum absolute difference of first component of q
    # between a cell and each of its neighbors.

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 2

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells
    # refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged
    # cells)
    amrdata.clustering_cutoff =0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0


    # ---------------
    # Regions:
    # ---------------
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    xbuffer = 0.5*probdata.fault_width
    ybuffer = xbuffer
    xcb = [probdata.fault_center-xbuffer,probdata.fault_center+xbuffer]

    # high-resolution region to surround the fault during slip
    regions.append([amrdata.amr_levels_max,amrdata.amr_levels_max, 0,1,
                    xcb[0],xcb[1],
                    -probdata.fault_depth-dx, -probdata.fault_depth+dx])

    for j in range(amrdata.amr_levels_max-1):
        regions.append([1,amrdata.amr_levels_max-j, 0,1e9,
                    xcb[0]-(j+1)*xbuffer,xcb[1]+(j+1)*xbuffer,
                    -probdata.fault_depth-(j+1)*ybuffer,0.0])

    regions.append([1,1, 0,1e9, -1.e9,1.e9, -1.e9,0.0])
    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting

    return rundata

    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
