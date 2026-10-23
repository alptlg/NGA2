# 1D laminar flame speed template (NGA2)

Computes laminar flame speeds S_L of several mechanisms (detailed, skeletal, ARCANE QSS-reduced) and compares them.
The solver is a Fortran port of Cantera 3.0 `FreeFlame` on NGA2's chemistry module: with the default
`Convection scheme : upwind` it reproduces Cantera's flame speeds (to < 0.1 %; `python lfs.py compare` checks it);
`second-order` gives flame speeds that depend much less on the grid.

## Requirements
* Python with numpy, pandas, matplotlib and PyYAML (or ruamel.yaml). No Cantera is needed (only to compare with
  Cantera results).
* gfortran, MPI and LAPACK/BLAS (from `/usr`; otherwise `export LAPACK_DIR=...`).
* NGA2 with `tools/scripts/chemistry/yaml2nga_lib_cf.py` (writes the Fortran chemistry module without Cantera, with
  the same rates, thermo and transport as Cantera).

## Use it
1. Copy this folder into `NGA2/examples/` (NGA2 is then found automatically); elsewhere set `NGA_HOME` to the NGA2
   repository.
2. Put each mechanism in its own folder under `mechs/` and list it in the study file `input` (NGA2 input format):
   ```
   Mechanism detailed : mechs/detailed/detailed.yaml              # detailed / skeletal: Cantera YAML
   Mechanism reduced : mechs/reduced/reduced.yaml mechs/reduced/reduced.f90   # ARCANE QSS: YAML + Fortran
   ```
3. Set fuels, pressures and equivalence ratios in `input`, then:
   ```bash
   python lfs.py all
   ```
4. Look at `results/` (below).

Folder layout: `input` (study settings, read by `lfs.py`), `input_single` (one flame, read by the executable),
`lfs.py` (driver), `src/simulation.f90` (flame solver), `GNUmakefile` (`make MECH=<name>`), `mechs/<name>/` (your
mechanism files + generated chemistry and executable), `results/`.
The example study uses Li-ion battery electrolyte solvent blends (EC, DMC, EMC, DEC in air) with the LIB electrolyte
surrogate mechanism of Kanayama et al. (Combust. Flame), which NGA2 ships in `tools/scripts/chemistry/kinetics/`
(the example `input` points there; generated files still go to `mechs/detailed/`). Replace both with yours.

### What happens to a mechanism
1. You provide the Cantera YAML (plus ARCANE's Fortran file for a QSS mechanism) and name it in `input`. The study
   file points to these source files, not to generated ones, so `build` can regenerate whenever they change;
   `lfs.py` also reads the species list from the YAML.
2. `python lfs.py build` runs NGA2's `tools/scripts/chemistry/yaml2nga_lib_cf.py`, which writes
   `mechs/<name>/chem_data_fc.f90` (species, thermo, transport properties and
   reaction rates; for a QSS mechanism the rates come from a patched
   copy of ARCANE's file, `mechs/<name>/customkinetics.f90`).
3. `make MECH=<name>` (called by `build`) compiles it with `src/simulation.f90` into
   `mechs/<name>/chem_reactor_1D_lfs.*.exe`: the chemistry is compiled in, so there is one executable per mechanism.
4. `python lfs.py run` starts that executable for every point, each with a one-run input file like `input_single`.

NGA2, and with it `yaml2nga_lib_cf.py`, is taken from `$NGA_HOME` if set, otherwise from `../..` (this folder inside
`NGA2/examples/`); the GNUmakefile uses the same rule.

## Commands
| command | what it does |
|---|---|
| `python lfs.py all` | build + run + plot |
| `python lfs.py build [mech]` | write the Fortran chemistry (with transport) and compile (skips up-to-date ones; `--force`) |
| `python lfs.py run [mech]` | flame-speed sweep in parallel (resumes; retries failed points) |
| `python lfs.py plot` | figures and `summary.md` |
| `python lfs.py compare <lfs.csv>` | compare with a Cantera `lfs.csv` (columns mech, fuel, P_atm, phi, sl_cm_s, cpu_s, npts) |
| `mechs/<mech>/chem_reactor_1D_lfs.*.exe -i input_single` | one flame by hand (`input_single` lists the keywords) |

Every Python command reads the study file `input`; `-i <file>` uses another one (e.g. `python lfs.py all -i input_5atm`).
`run` keeps the flames already in `results/lfs.csv` and only adds missing ones. It stops if the unburnt temperature,
refinement, convection scheme or mixtures differ from those of the existing results (`results/run_settings.json`):
use another `Output directory` then. After rebuilding a mechanism with other kinetics, delete its rows in `lfs.csv`.

## Results
| file | content |
|---|---|
| `results/lfs_P<P>atm.png` | S_L vs phi for every fuel, error vs the reference mechanism |
| `results/summary.md` | per mechanism: S_L error (mean, max), failures, CPU, speed-up, grid points |
| `results/lfs.csv` | every flame: S_L, CPU, grid points, how it was started |
| `results/solutions/<mech>_<fuel>_P<P>_phi<phi>.dat` | converged flames: z, u, T, Y_k |
| `results/compare_cantera.md` | written by `compare` |

## How the flame speed is calculated
**The problem.** A flat flame moves into the unburnt mixture at the laminar flame speed S_L. In the frame of the
flame the problem is steady: unburnt gas enters on the left at T0 with velocity S_L, burns, and leaves on the right
at the adiabatic flame temperature. S_L is the inlet velocity for which this steady flame exists; it is part of
the solution.

**Assumptions.** One-dimensional, planar, steady, adiabatic; ideal gas; low Mach number (constant pressure, no
viscous heating or kinetic energy); mixture-averaged diffusion with a correction velocity; no thermal diffusion
(Soret), no radiation, no gravity.

**Equations** (z = position, u = velocity, rho = density, Y_k = mass fractions, W_k = molar masses,
w_k = molar production rates, h_k = molar enthalpies):
* mass: d(rho u)/dz = 0, so rho u = rho_0 S_L everywhere
* species: rho u dY_k/dz = - dj_k/dz + W_k w_k
* energy: rho c_p u dT/dz = d/dz(lambda dT/dz) - sum_k h_k w_k - sum_k (j_k / W_k) dh_k/dz
* diffusion: j_k = -rho (W_k / W) D_km dX_k/dz, minus Y_k sum_j j_j so that the fluxes sum to zero;
  D_km = (1 - Y_k) / sum_{j != k} X_j / D_kj; mixture conductivity lambda = (sum X_k lambda_k + 1 / sum X_k / lambda_k) / 2
* ideal gas: rho = P W / (R T) with P fixed

**Boundary conditions.** Inlet (left): T = T0, and the inflow of each species (convection + diffusion) equals that
of the unburnt mixture. Outlet (right): zero gradients of T and Y_k.

**Holding the flame in place.** A steady flame can sit anywhere, so the temperature is fixed at one grid point
(about a quarter of the way from T0 to the flame temperature, snapped to a grid point). This extra condition is
what lets the solver find S_L.

**Discretisation.** Grid points that are refined around the flame; convection by first-order upwind differences
(as Cantera) or second order (below); diffusion by central differences with transport properties evaluated midway
between points.

**Solution method.** The equations at all points are solved together with a damped Newton method (Jacobian by
finite differences, reused for up to 20 iterations). If Newton fails, a few small implicit time steps are taken to
bring the solution closer and Newton is tried again. When converged, the grid is refined where the solution changes
too much between neighbouring points, and the problem is solved again, until no new points are needed.

**Convection scheme.** `upwind` (default, as Cantera) takes the slope at a point from the point upstream. It is
very robust but acts like extra, unphysical diffusion (about u h / 2, h = grid spacing), so S_L depends on the grid:
several % with the default refinement, more for slow lean or high-pressure flames, and ~1 % even with 1000 points.
`second-order` uses three-point central slopes (error ~ h^2) and removes most of this:
* deferred correction: Newton always solves the robust upwind equations plus a source equal to the difference between
  upwind and second-order convection of the previous solution, repeated until the solution stops changing (each pass
  is a quick Newton solve); the converged flame satisfies the second-order equations;
* the second-order slope falls back to upwind where a variable is not resolved (coarse upstream cells, local extrema,
  tiny upstream tails), so no wiggles or negative mass fractions; cold starts are solved upwind first;
* with the default refinement S_L is then typically within ~1 % of its grid-converged value; a sweep costs about
  twice the CPU time of `upwind`;
* the error of a reduced mechanism vs detailed changes little between the two schemes, so `upwind` is fine for
  comparing mechanisms; use `second-order` when absolute flame speeds matter. Results then no longer match Cantera
  point for point (they approach its grid-converged limit).

**Starting the solution.**
* Cold start (no previous flame): a coarse grid over the domain width with straight-line profiles from the unburnt
  gas to its adiabatic equilibrium between 30 % and 50 % of the domain; if that fails, the grid is made finer
  (12, 24, 48 points) and/or the temperature is first held fixed. If the converged flame touches either end of the
  domain, the domain is doubled and the solve repeated. The equilibrium (as Cantera `equilibrate('HP')`) is computed
  by the solver itself: element-potential method, Newton with an analytic Jacobian, started from complete-combustion
  products.
* Continuation (sweep): the reference mechanism starts cold at the first phi >= 1, then each next phi (up to the
  richest, then down to the leanest) starts from its neighbour. Every other mechanism starts from the reference
  flame at the same phi. If a start from a profile fails, cold starts on the `Cold start widths` are tried.

## Errors and comparison
Every mechanism is solved at exactly the same points (fuel, P, phi) and compared with the reference mechanism
(`Reference mechanism`) at the same point: S_L error = 100 x (S_L,mech - S_L,ref) / S_L,ref [%]. Plots: signed error
vs phi; `summary.md`: mean and max of |error| and the failures.

## CPU time
* What is measured: the process CPU time of each flame solve (`cpu_s`), written by the solver.
* Not the same basis for every flame: in a sweep the reference flame starts from its neighbour phi while the other
  mechanisms start from the converged reference flame at the same phi, and grids differ. CPU times and speed-ups in
  `summary.md` are therefore indicative only; runs in parallel also disturb each other (set `Workers : 1` for timings
  to quote).

## Defaults
| setting | default | where (`input` keyword) |
|---|---|---|
| unburnt temperature T0 | 300 K | `Unburnt temperature` |
| pressures | 1 and 5 atm | `Pressures` |
| equivalence ratios | 0.5 to 2.0, step 0.1 | `Equivalence ratio range` (or `Equivalence ratios`) |
| domain width for cold starts | 4 cm (then 8 cm, 12 cm if needed; doubled automatically if the flame touches an edge) | `Cold start widths` |
| initial grid | 8 points at 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1 x width | solver |
| grid refinement | ratio 3, slope 0.08, curve 0.08, prune 0.05 | `Refine ratio / slope / curve / prune` |
| maximum grid points | 400 (a flame whose grid reaches it fails, as in Cantera) | `Max grid points` |
| convection | first-order upwind (as Cantera) | `Convection scheme` (`upwind` or `second-order`) |
| Newton tolerances | relative 1e-4, absolute 1e-9 (time steps: 1e-4, 1e-11) | solver |
| time steps when Newton fails | 10 steps starting at 1e-5 s | solver |
| bounds during iterations | T between 200 K and 2 x the lowest species T_max (by hand: `T upper bound`, default 4000 K); Y_k >= -1e-7 | solver |
| accepted flame speeds | 0 < S_L <= 200 cm/s | `Max flame speed` |
| time limit per flame | 5400 s | `Run timeout` |

What the refinement criteria mean: **slope** = largest change of any variable between two neighbouring points, as a
fraction of its total change across the flame; **curve** = same for the change of slope; **ratio** = largest size
ratio of neighbouring grid intervals; **prune** = points where all changes are below this fraction are removed.

## Good to know
* Mechanisms: detailed/skeletal YAMLs are converted by `yaml2nga_lib_cf.py` (exact reverse rates and PLOG, Cantera's
  mixture-averaged transport fits). ARCANE QSS Fortran is patched at build time (QSS concentration floor and 1e3
  ceiling removed, inclusive PLOG pressure brackets). Mechanisms whose YAML holds placeholder reactions (rates from a
  separate hand-written routine) are not supported.
* Accuracy: with `upwind`, two converged runs that end on slightly different grids can differ by ~1 %; see
  "Convection scheme".
* Memory: each running flame needs about 5 x (species + 2)^2 x points x 8 bytes for the Jacobian (~0.4 GB for 160
  species and 400 points); choose `Workers` so that Workers x that fits in memory.
* Species names cannot contain `! # % : =` (NGA2 input format).
* Unlike most NGA2 examples, `src/simulation.f90` is a standalone program (NGA2 libraries + LAPACK, no grid/config
  stack) and there is one executable per mechanism, in `mechs/<name>/`; each flame runs serially (no `mpiexec`).
