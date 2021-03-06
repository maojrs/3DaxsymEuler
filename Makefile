
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = classic                  # Clawpack package to use
EXE = xclaw                         # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= -I/usr/include -L/usr/lib -lsilo -lm -lstdc++

# ---------------------------------
# List of sources for this program:
# ---------------------------------

MODULES = \

SOURCES = \
  auxmodule.f90 \
  qinit.f \
  rpn2HLLC-exact.f \
  rpt2ac.f \
  setprob.f \
  setaux.f \
  src2.f \
  b4step2.f \
  out2_2D.f \
  out2silo_2D.f90 \
  bc2.f \
  phi_exact_tamman.f90 \
  $(CLAW)/classic/src/2d/driver.f90 \
  $(CLAW)/classic/src/2d/claw2ez.f \
  $(CLAW)/classic/src/2d/claw2.f \
  $(CLAW)/classic/src/2d/step2.f90 \
  $(CLAW)/classic/src/2d/step2ds.f90 \
  $(CLAW)/classic/src/2d/dimsp2.f \
  $(CLAW)/classic/src/2d/flux2.f90 \
  $(CLAW)/classic/src/2d/copyq2.f \
  $(CLAW)/classic/src/2d/inlinelimiter.f90 \
  $(CLAW)/classic/src/2d/restart2.f \
  $(CLAW)/classic/src/2d/opendatafile.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

