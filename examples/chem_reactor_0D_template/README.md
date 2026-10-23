# 0D ignition-delay template (NGA2)

Compares ignition delay times (IDT), end states and time histories of several mechanisms (detailed, skeletal,
ARCANE QSS-reduced) in an adiabatic 0D reactor. The IDT definition is the one of a Cantera reactor sweep, and
`python idt.py check` recomputes the same points with Cantera to confirm it.

## Requirements
* Python with numpy, pandas, matplotlib and PyYAML (or ruamel.yaml); activate it first. Cantera (>= 3.0) is needed
  only for `check`.
* gfortran, MPI and SUNDIALS/CVODE, as for the other NGA2 chemistry examples; `CVODE_DIR` = SUNDIALS install root
  (the same variable NGA2's `Make.defs` uses).
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
3. Set fuels, pressures, equivalence ratios and temperatures in `input`, then:
   ```bash
   export CVODE_DIR=/path/to/sundials
   python idt.py all
   ```
4. Look at `results/` (below).

Folder layout: `input` (study settings, read by `idt.py`), `input_single` (one run, read by the executable),
`idt.py` and `plot_results.py` (drivers), `src/simulation.f90` (reactor), `GNUmakefile` (`make MECH=<name>`),
`mechs/<name>/` (your mechanism files + generated chemistry and executable), `results/`.
The example study uses Li-ion battery electrolyte solvent blends (EC, DMC, EMC, DEC in air) with the LIB electrolyte
surrogate mechanism of Kanayama et al. (Combust. Flame), which NGA2 ships in `tools/scripts/chemistry/kinetics/`
(the example `input` points there; generated files still go to `mechs/detailed/`). Replace both with yours.

### What happens to a mechanism
1. You provide the Cantera YAML (plus ARCANE's Fortran file for a QSS mechanism) and name it in `input`. The study
   file points to these source files, not to generated ones, so `build` can regenerate whenever they change;
   `check` also runs Cantera on the YAML.
2. `python idt.py build` runs NGA2's `tools/scripts/chemistry/yaml2nga_lib_cf.py`, which writes
   `mechs/<name>/chem_data_fc.f90` (species, thermo and reaction rates; for a QSS mechanism the rates come from a patched
   copy of ARCANE's file, `mechs/<name>/customkinetics.f90`).
3. `make MECH=<name>` (called by `build`) compiles it with `src/simulation.f90` into
   `mechs/<name>/chem_reactor_0D.*.exe`: the chemistry is compiled in, so there is one executable per mechanism.
4. `python idt.py run` starts that executable for every point, each with a one-run input file like `input_single`.

NGA2, and with it `yaml2nga_lib_cf.py`, is taken from `$NGA_HOME` if set, otherwise from `../..` (this folder inside
`NGA2/examples/`); the GNUmakefile uses the same rule.

## Commands
| command | what it does |
|---|---|
| `python idt.py all` | build + run + plot |
| `python idt.py build [mech]` | write the Fortran chemistry and compile (skips up-to-date ones; `--force`) |
| `python idt.py run [mech]` | run the sweep in parallel (resumes; retries failed points) |
| `python idt.py plot` | figures and `summary.md` |
| `python idt.py check [--every N]` | every N-th point (default 5) again with Cantera (not for QSS mechanisms) |
| `python plot_results.py LP40 1 1 650` | time history of one case: fuel, P [atm], phi, T0 [K] |
| `mechs/<mech>/chem_reactor_0D.*.exe -i input_single` | a single run by hand (`input_single` is an example) |

Every Python command reads the study file `input`; `-i <file>` uses another one (e.g. `python idt.py all -i input_rich`).
`run` keeps the points already in `results/idt.csv` and only adds missing ones. It stops if the reactor, end mode,
tolerances, IDT criterion or mixtures differ from those of the existing results (`results/run_settings.json`): use
another `Output directory` then. After rebuilding a mechanism with other kinetics, delete its rows in `idt.csv`.

## Results
| file | content |
|---|---|
| `results/summary.md` | per mechanism: IDT error (mean, max, points > 20 %), end T / CO / CO2 / H2O deviation, CPU, speed-up |
| `results/idt_P<P>atm.png` | IDT vs 1000/T0 for every fuel and phi, error vs the reference, CPU per run |
| `results/endstate/endstate_<fuel>_P<P>atm_phi<phi>[_<n>].png` | end state vs T0: T, P, rho and species, all mechanisms |
| `results/history/history_<fuel>_P<P>atm_phi<phi>_T<T0>K[_<n>].png` | time evolution of one case, all mechanisms |
| `results/idt.csv` | every run: IDT, end state, CPU |
| `results/check_cantera.csv` | written by `check` |

`_<n>` numbers the pages when a figure set needs more than one png (28 panels each).

## What is computed
* **Mixture**: fuel and oxidizer mole fractions + equivalence ratio, exactly as Cantera `set_equivalence_ratio`.
* **Reactor**: adiabatic, solved either isochoric (constant volume, the default) or isobaric (constant pressure).
  Choose it with `Reactor : isochor` or `Reactor : isobar` (also `cv` / `cp`) in `input` (whole study) or
  `input_single` (one run). Isochoric: density fixed, pressure rises during ignition. Isobaric: pressure fixed, the
  gas expands. Constant pressure gives longer IDTs than constant volume, so only compare like with like.
* **End of run** (`End mode`): `flamemaster` (default; the stop rule of FlameMaster's 0D reactor) stops once T has
  stopped changing and t > 3 IDT, i.e. at equilibrium within ~0.3 K; `fixed` stops at max(4 IDT, 0.5 s);
  `maxtime` integrates to `Max time`.
* **Kinetics**: detailed/skeletal mechanisms are converted by `yaml2nga_lib_cf.py` (exact reverse rates and PLOG, same
  rates as Cantera). QSS mechanisms use ARCANE's Fortran; `build` patches two generator artefacts (the floor and the
  1e3 ceiling on QSS concentrations are removed, PLOG pressure brackets made inclusive). Mechanisms whose YAML holds
  placeholder reactions (all rate constants zero, rates from a separate hand-written routine) are not supported:
  `build` stops with a message.
* **Species shown**: species of the reference mechanism that another mechanism also has. Very small values are not
  meaningful and are left out: end states below `End state x_min` (1e-10); in time histories, species whose peak
  stays below 10 x `--x-min` (default 1e-15).

### IDT criterion
The IDT is the time of the steepest temperature rise, taken on the solver's own steps inside the Fortran (not from
the saved trajectory, which is thinned):
1. after every CVODE step i: slope dT/dt = (T_i - T_i-1)/(t_i - t_i-1), at the mid time (t_i + t_i-1)/2;
2. only steps with mid time >= `IDT min time` (1 us) and T_i > T0 + `IDT min rise` (55 K) count, which skips the
   initial transient and weak low-temperature heat release;
3. the search stops at the first step with T > T0 + `IDT window` (400 K): ignition has happened; the largest slope
   found so far is the IDT;
4. if T never exceeds T0 + 400 K before `Max time` (1e4 s), the point is "not ignited" (an IDT is still reported if
   T rose by more than 55 K).

At high T0 the temperature trace can have two nearly equal maxima of dT/dt, so IDTs there are less robust.

## Errors and comparison
Every mechanism is run at exactly the same points (fuel, P, phi, T0), and each point is compared with the reference
mechanism (`Reference mechanism`) at the same point:
* **IDT error** per point: 100 x (IDT_mech - IDT_ref) / IDT_ref [%] (positive = the mechanism ignites later).
  Plots: signed error vs 1000/T0 (green band = +-20 %). `summary.md`: mean and max of |error| over all points,
  number of points above 20 %.
* **Ignition mismatch**: points where one mechanism ignites and the other does not. **Failed**: solver failures
  (retried by `run`).
* **End state** (points where both ignited): temperature as a difference in K, CO / CO2 / H2O (and every species in
  the end-state figures) as 100 x (X_mech - X_ref) / X_ref [%]. `summary.md` gives the maximum over all points; the
  figures compare species only where the reference mole fraction is above `End state x_min`.
* **Time histories**: the peak mole fraction of each species, 100 x (peak_mech - peak_ref) / peak_ref [%].

## CPU time
* What is measured: the process CPU time of the integration loop only (first to last CVODE step), written as
  `cpu_s` for every run; setup and the result file are not included (the optional trajectory of `plot_results.py`
  is).
* Same basis for every mechanism: same solver (CVODE BDF, dense finite-difference Jacobian), same tolerances, same
  start state and the same stopping rule. The difference is the cost of the chemistry (number of species and
  reactions, QSS algebra).
* `summary.md` gives the median CPU per mechanism and the **speed-up** = median, over the points ignited in both, of
  CPU(reference) / CPU(mech) (2 = twice as fast). The IDT figures show CPU per run vs 1000/T0.
* Caveat: the sweep runs `Workers` simulations at the same time, so CPU times are disturbed by the other runs
  (shared caches and memory bandwidth). For numbers to quote, set `Workers : 1` on an otherwise idle machine and
  rerun the mechanisms you compare (delete their rows in `results/idt.csv`).

## Good to know
* Solver tolerances `CVODE rtol` 1e-10, `CVODE atol` 1e-20. Keep `atol` ~1e-20: larger values can give wrong IDTs
  and solver failures at low T0, where radical mass fractions are tiny.
* Pressures are in atm in the study file (`Pressures`) but in Pa in one-run files (`Pressure`).
* Time histories: linear time axis by default (ignition at ~1/3, then the approach to equilibrium);
  `--time-axis log` shows the early low-temperature chemistry; `--keep` keeps the trajectory files; `--rows` should
  exceed the number of CVODE steps (printed) so that short ignition peaks are not dropped.
* Species names cannot contain `! # % : =` (NGA2 input format).
* Unlike most NGA2 examples, `src/simulation.f90` is a standalone program (NGA2 libraries + CVODE, no grid/config
  stack) and there is one executable per mechanism, in `mechs/<name>/`; runs are serial (no `mpiexec`).
