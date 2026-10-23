# chem_reactor_0D template: 0D reactor + ignition delay (see README.md)
# One executable per mechanism: mechs/<MECH>/ holds chem_data_fc.f90 (+ customkinetics.f90 for ARCANE QSS),
# its build tree and chem_reactor_0D.<suffix>.exe (one executable per mechanism, hence EXEBASE inside mechs/).
# Normally driven by idt.py ("python idt.py build [mech] [-i input]"); by hand:
#   make MECH=<name>        (requires CVODE_DIR)

# NGA2 root: found automatically inside NGA2/examples/ (or in a sibling folder of Repositories/NGA2); else set NGA_HOME
NGA_HOME ?= $(if $(wildcard ../../tools/GNUMake/Make.defs),../..,../../Repositories/NGA2)
MECH     ?= detailed
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_CVODE = TRUE
COMP      = gnu
DEBUG     = FALSE
EXEBASE   = mechs/$(MECH)/chem_reactor_0D

ifneq ($(MAKECMDGOALS),clean)
ifeq ($(wildcard mechs/$(MECH)/chem_data_fc.f90),)
  $(error mechs/$(MECH)/chem_data_fc.f90 not found: run "python idt.py build $(MECH)" first)
endif
endif

Udirs := src mechs/$(MECH)
INCLUDE_LOCATIONS += $(Udirs)
VPATH_LOCATIONS   += $(Udirs)

include $(NGA_HOME)/tools/GNUMake/Make.defs
BuildDir = mechs/$(MECH)/build

# NGA2 libraries needed here only (no grid/config stack; sources listed here, no src/Make.package)
Bdirs := libraries
include $(foreach dir, $(Bdirs), $(NGA_HOME)/src/$(dir)/Make.package)
f90EXE_sources := $(filter precision_$(PREC).f90 string.f90 parallel.f90 messager.f90 param.f90, $(f90EXE_sources))
f90EXE_sources += simulation.f90 chem_data_fc.f90 $(notdir $(wildcard mechs/$(MECH)/customkinetics.f90))

# keep IEEE NaN handling under -ffast-math (ARCANE QSS code tests cqss /= cqss)
ifeq ($(COMP),gnu)
  F90FLAGS += -fno-finite-math-only
endif

all: $(executable)
	@echo "Build complete: $(executable)"

include $(NGA_HOME)/tools/GNUMake/Make.rules
