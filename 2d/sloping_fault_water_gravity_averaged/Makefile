
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/seismic/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= -fopenmp -ffpe-trap='overflow','invalid','zero','underflow','inexact','denormal'

# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
include $(AMRLIB)/Makefile.amr_2d
SEISLIB = $(CLAW)/seismic/src/2d


# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \
    $(SEISLIB)/fault_module.f90

SOURCES = \
  cellave.f \
  qinit.f \
  setprob.f90 \
  fdisc.f \
  bc2amr.f \
  setaux.f90 \
  src2.f90 \
  b4step2.f90 \
  mapc2p.f90 \
  $(SEISLIB)/tick.f \
  $(CLAW)/seismic/src/rp/rpn2_vc_elasticity_mapped_fault_abl.f90 \
	$(CLAW)/seismic/src/rp/rpt2_vc_elasticity_mapped_fault_abl.f90 \
  $(AMRLIB)/fss.f \
  $(AMRLIB)/zeroin.f
#	$(AMRLIB)/cellave.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)
