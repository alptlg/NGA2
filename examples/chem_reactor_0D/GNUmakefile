# 0D homogeneous reactor (NGA2 example structure)
# chem_data_fc.f90 in src/ - generate with: python ../../tools/scripts/chemistry/yaml2nga.py <mechanism.yaml> src/chem_data_fc.f90
# Run: ./chem_reactor_0D.dp.gnu.opt.exe -i input
#      (or with MPI: mpirun -np 1 ./chem_reactor_0D.dp.gnu.opt.mpi.exe -i input)

NGA_HOME ?= ../..
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_CVODE = TRUE
COMP      = gnu
DEBUG     = FALSE
EXEBASE   = chem_reactor_0D

# User code in src/
Udirs   := src
Upack   += $(foreach dir, $(Udirs), $(wildcard $(dir)/Make.package))
Ulocs   += $(foreach dir, $(Udirs), $(wildcard $(dir)))
include $(Upack)
INCLUDE_LOCATIONS += $(Ulocs)
VPATH_LOCATIONS   += $(Ulocs)

include $(NGA_HOME)/tools/GNUMake/Make.defs

# Only libraries (param, parallel, messager, precision, string) - no core/grid/config
Bdirs   := libraries
Bpack   += $(foreach dir, $(Bdirs), $(NGA_HOME)/src/$(dir)/Make.package)
include $(Bpack)

# Restrict to minimal deps + our sources (no main.f90, no full NGA stack)
f90EXE_sources := $(filter precision_$(PREC).f90 string.f90 parallel.f90 messager.f90 param.f90, $(f90EXE_sources))
f90EXE_sources += simulation.f90 chem_data_fc.f90

all: $(executable)
	@echo "Build complete. Run: ./$(executable) -i input"

include $(NGA_HOME)/tools/GNUMake/Make.rules
