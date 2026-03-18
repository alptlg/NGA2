# chem_reactor_0D — 0D Homogeneous Reactor

0D adiabatic constant-pressure chemistry solver. Integrates mass fraction and temperature evolution using CVODE (SUNDIALS) with adaptive stepping. Chemical source terms come from the `fcmech` kinetics module.

## Naming convention (for future variants)

Internal prefixes for homogeneous reactor variants:

| Prefix | Variant | Description |
|--------|---------|-------------|
| `hr_ib_` | Isobaric | Constant pressure (current) |
| `hr_ic_` | Isochoric | Constant volume (future) |

## Features

- **Adiabatic, constant pressure**: ideal gas, NASA7 thermo
- **Initial composition**: Fuel + equivalence ratio (standard air) or per-species mass fractions
- **Kinetics**: Arrhenius, three-body, falloff (Troe), plog
- **Output**: Time series of species mass fractions, T, P, ρ

## Requirements

- **Fortran compiler**: gfortran (or Intel, etc.)
- **SUNDIALS/CVODE**: For stiff ODE integration
  - macOS: `brew install sundials`
  - Linux: install via package manager or build from source

## Directory Structure

```
chem_reactor_0D/
├── GNUmakefile
├── README.md
├── input              # Input parameters
├── plot_results.py    # Plot results
└── src/
    ├── Make.package
    ├── simulation.f90      # Main program
    ├── chem_data_fc.f90      # Kinetics (generated)
    └── chem_data_reactions.txt # Reaction index (generated)
```

## Quick Start

### 1. Generate the kinetics module

`chem_data_fc.f90` must be generated from a Cantera YAML mechanism:

```bash
cd examples/chem_reactor_0D

# Generate from a mechanism (e.g. Nakamura ammonia, gri30)
python ../../tools/scripts/chemistry/yaml2nga.py \
    ../../tools/scripts/chemistry/kinetics/Nakamura.yaml \
    src/chem_data_fc.f90
```

Or with gri30:
```bash
python ../../tools/scripts/chemistry/yaml2nga.py \
    ../../tools/scripts/chemistry/kinetics/gri30.yaml \
    src/chem_data_fc.f90
```

This creates `src/chem_data_fc.f90` and `src/chem_data_reactions.txt`.

### 2. Build

```bash
make CVODE_DIR=$(brew --prefix sundials)
```

If SUNDIALS is elsewhere, set `CVODE_DIR` to the install root (with `include/`, `lib/`, `fortran/`).

### 3. Run

```bash
./chem_reactor_0D.dp.gnu.opt.exe -i input
```

With MPI (single rank):

```bash
mpirun -np 1 ./chem_reactor_0D.dp.gnu.opt.mpi.exe -i input
```

### 4. Plot results

```bash
python plot_results.py results_hr.out
```

## Input File

The `input` file uses NGA2 param format (`key : value`). Lines starting with `#` are comments.

### Required parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `Temperature` | Initial temperature (K) | `1000.0` |
| `Pressure` | Pressure (Pa) | `100000.0` |
| `End time` | Integration end time (s) | `0.1` |
| `Time step` | Output interval (s) | `1e-6` |
| `Output file` | Output filename | `results_hr.out` |

### Reactor type

| Parameter | Description | Example |
|-----------|-------------|---------|
| `Reactor` | `isobar` (constant P) or `isochor` (constant V) | `isobar` |

Default: `isobar`. Case independent.

### Initial composition (choose one)

**Option A — Fuel + equivalence ratio**

```text
Fuel : NH3
Equivalence ratio : 1.0
```

Uses standard air (N2=0.79, O2=0.21 mole frac). Works for fuels without carbon (e.g. NH3) or with carbon (e.g. CH4).

**Option B — Per-species mass fractions**

```text
Initial Y CH4 : 0.05
Initial Y O2 : 0.2
Initial Y N2 : 0.75
```

`Initial Y <species>` for each species. Sum must be 1.0.

### Command-line input

```bash
./chem_reactor_0D.dp.gnu.opt.exe -i input
```

## Output Format

`results_hr.out` is a space-separated file:

- **Header**: `time`, `species1-Y1`, `species2-Y2`, ..., `T`, `P`, `rho`
- **Columns**: Fixed-width columns (18 chars for time, species, T, P, rho)

## Plotting

`plot_results.py` reads the output and plots temperature and selected species:

```bash
python plot_results.py [results_hr.out] [-o output_prefix]
```

Requires: `numpy`, `matplotlib`

## Available Mechanisms

YAML mechanisms are in `tools/scripts/chemistry/kinetics/`:

- `gri30.yaml` — GRI-Mech 3.0 (CH4/hydrocarbon)
- `Nakamura.yaml` — Ammonia (NH3)
- `h2o2.yaml` — H2/O2
- `ONE.yaml` — Minimal single reaction

## Build Options

| Variable | Default | Description |
|----------|---------|-------------|
| `PRECISION` | `DOUBLE` | `DOUBLE` or `SINGLE` |
| `USE_MPI` | `TRUE` | MPI support |
| `COMP` | `gnu` | Compiler: `gnu`, `intel`, etc. |
| `DEBUG` | `FALSE` | Debug build |

Example:

```bash
make PRECISION=DOUBLE COMP=intel CVODE_DIR=/path/to/sundials
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `CVODE_DIR not set` | Install SUNDIALS and pass `CVODE_DIR=$(brew --prefix sundials)` |
| `Fuel, O2, or N2 not found` | Mechanism must define O2 and N2; fuel name must match species name |
| `No Initial Y specified` | Provide either Fuel+Equivalence ratio or Initial Y per species |
| `chem_data_fc.f90` missing | Run `yaml2nga.py` to generate from a YAML mechanism |

## See Also

- `tools/scripts/chemistry/yaml2nga.py` — YAML → chem_data_fc converter
- `tools/scripts/finite_chemistry/V2/` — Cantera comparison utilities
