"""
Module to set up run time parameters for Clawpack -- amrclaw code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
reload(dtopotools)
from clawpack.seismic.mappings import Mapping2D
from make_topo_and_grid import get_oceanfloor_parameters

USE_TOPO = True

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

    if (USE_TOPO):
        # Obtain topography parameters to match 1D Geoclaw
        xlower_domain, xlower_slope, xlower_shelf, xlower_beach, xlower_shore, xupper_domain, \
          zlower_ocean, zlower_shelf, zlower_beach, zlower_shore = get_oceanfloor_parameters()
    else:
        xlower_domain = -150e3
        xlower_slope = -65e3
        xlower_shelf = -45e3
        xlower_beach = -5e3
        xlower_shore = 0.0
        xupper_domain = 0.0
        zlower_ocean = -4500.0
        zlower_shelf = -4500.0
        zlower_shore = -4500.0


    # Adjust topo parameters to remove shore
    zlower_shore = min(zlower_shore, -100.0)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('zlower_ocean', zlower_ocean, 'z-coord of ocean floor')
    probdata.add_param('xlower_slope', xlower_slope, 'x-coord of beginning of slope')
    probdata.add_param('xlower_shelf', xlower_shelf, 'x-coord of beginning of shelf')
    probdata.add_param('zlower_shelf', zlower_shelf, 'z-coord of shelf')
    probdata.add_param('xlower_beach', xlower_beach, 'x-coord of beginning of beach')
    probdata.add_param('xlower_shore', xlower_shore, 'x-coord of beginning of shore')
    probdata.add_param('zlower_shore', zlower_shore, 'z-coord of shore')
    probdata.add_param('abl_depth', 30e3, 'depth of absorbing layer')
    probdata.add_param('fault_zshift', 0.0, 'vertical shift of comp domain to match fault depth')
    probdata.add_param('domain_depth', 50e3, 'depth of domain')
    probdata.add_param('domain_width', xupper_domain-xlower_domain, 'width of domain')

    #------------------------------------------------------------------
    # Read in fault information
    #------------------------------------------------------------------
    fault = dtopotools.Fault()
    fault.read('fault.data')

    mapping = Mapping2D(fault)
    fault_width = mapping.fault_width
    fault_depth = mapping.fault_depth
    fault_center = mapping.xcenter

    rupture_rise_time = 0.0
    for subfault in fault.subfaults:
        rupture_rise_time = max(rupture_rise_time,subfault.rupture_time
                                    + subfault.rise_time)

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
    num_cells_fault = 20

    # determine cell number and set computational boundaries
    dx = fault_width/num_cells_fault
    # x direction
    num_cells_above_fault = np.ceil((xupper_domain - (fault_center+0.5*fault_width))/dx)
    target_num_cells = np.rint(probdata.domain_width/dx)
    num_cells_below_fault = target_num_cells - num_cells_above_fault - num_cells_fault
    clawdata.num_cells[0] = int(target_num_cells)
    clawdata.lower[0] = fault_center-0.5*fault_width - num_cells_below_fault*dx
    clawdata.upper[0] = fault_center+0.5*fault_width + num_cells_above_fault*dx
    # z direction
    num_cells_across_ocean = np.ceil(-zlower_ocean/dx)
    dz = -zlower_ocean/num_cells_across_ocean
    clawdata.num_cells[1] = int(probdata.domain_depth/dz)
    clawdata.lower[1] = -clawdata.num_cells[1]*dz
    clawdata.upper[1] = 0.0

    probdata.fault_zshift = -fault_depth + np.ceil(fault_depth/dz)*dz

    # add absorbing layer
    target_num_cells = np.rint(probdata.abl_depth/dx)
    clawdata.lower[0] -= target_num_cells*dx
    clawdata.upper[0] += target_num_cells*dx
    clawdata.num_cells[0] += 2*int(target_num_cells)
    target_num_cells = np.rint(probdata.abl_depth/dz)
    clawdata.lower[1] -= target_num_cells*dz
    clawdata.num_cells[1] += int(target_num_cells)

    # Note adjustments in computational domain size
    probdata.domain_width = clawdata.upper[0] - clawdata.lower[0]
    probdata.domain_depth = clawdata.upper[1] - clawdata.lower[1]

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 6

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 15

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
        clawdata.num_output_times = 300
        clawdata.tfinal = 600.0
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  list(np.linspace(0,10,11))

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 10
        clawdata.output_t0 = True  # output at initial (or restart) time?


    clawdata.output_format = 'binary'      # 'ascii', 'binary', 'netcdf'

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
    clawdata.dt_initial = 0.001

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1.000000e+99

    # Desired Courant number if variable dt used
    clawdata.cfl_desired = 0.700000
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
    clawdata.source_split = 'godunov'


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

    clawdata.bc_lower[1] = 'extrap'   # at zlower
    clawdata.bc_upper[1] = 'user'   # at yupper

    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.
      clawdata.checkpt_times = []

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5



    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 3

    # List of refinement ratios at each level (length at least
    # amr_level_max-1)
    amrdata.refinement_ratios_x = [8,2]
    amrdata.refinement_ratios_y = [8,2]
    amrdata.refinement_ratios_t = [8,2]

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one
    # of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center', 'center', 'center', 'center', 'center', \
        'center', 'center', 'center', 'center', 'center', 'center', \
        'capacity','yleft','xleft','yleft']



    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.  # Richardson tolerance

    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 1.0e-5 # tolerance used in this routine
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
    amrdata.clustering_cutoff = 0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0


    # ---------------
    # Gauges:
    # ---------------
    gauges = rundata.gaugedata.gauges
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]

    xgauges = np.linspace(clawdata.lower[0]+1, clawdata.upper[0]-1,
                            np.rint(probdata.domain_width/1e3))
    ngauges = len(xgauges)

    # ocean floor:
    dz_level2 = dz/amrdata.refinement_ratios_y[0]

    for gaugeno,x in enumerate(xgauges):
        gauges.append([gaugeno,x,zlower_ocean - 0.5*dz_level2,0,1e10])

    # ocean surface:
    dz_level3 = dz_level3/amrdata.refinement_ratios_y[1]
    for gaugeno,x in enumerate(xgauges):
        gauges.append([ngauges+gaugeno,x,-0.5*dz_level3,0,1e10])

    # set gauge output increment to match rest of domain
    if clawdata.output_style==1:
        rundata.gaugedata.min_time_increment = clawdata.tfinal/clawdata.num_output_times

    elif clawdata.output_style == 2:
        rundata.gauagedata.min_time_increment = min(clawdata.output_times)


    # ---------------
    # Regions:
    # ---------------
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    # Region for the fault
    regions.append([amrdata.amr_levels_max-1, amrdata.amr_levels_max-1,
                    0,clawdata.dt_initial, #rupture_rise_time,
                    fault_center-0.5*fault_width,fault_center+0.5*fault_width,
                    -fault_depth-dx, -fault_depth+dx])

    # Debug
    # regions.append([1,amrdata.amr_levels_max-1,
    #                 0,1e9,
    #                 -1e9,1e9,
    #                 zlower_ocean,1e9])


    # Region for shelf (if exists)
    if (zlower_shelf > zlower_ocean):
        regions.append([1, amrdata.amr_levels_max-1,
                        0,1e9,
                        -1e9, 1e9,
                        -1e9, zlower_ocean])


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
