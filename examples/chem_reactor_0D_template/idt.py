#!/usr/bin/env python3
"""Ignition-delay study with the NGA2 0D reactor (chem_reactor_0D template).

    python idt.py build [mech ...]   generate mechs/<mech>/chem_data_fc.f90 (yaml2nga_lib_cf.py) and compile
                                     (skips mechanisms whose executable is up to date; --force)
    python idt.py run   [mech ...]   sweep fuels x P x phi x T0 -> <Output directory>/idt.csv (resumes)
    python idt.py plot               figures + summary.md from idt.csv
    python idt.py check [--every N]  same points with Cantera (same IDT definition), non-QSS mechanisms
    python idt.py all                build + run + plot
    -i <file> selects the study input file (default: input, NGA2 input format; see that file).

Mechanisms ('Mechanism <name> : <yaml> [<f90>]'): a Cantera YAML for detailed/skeletal mechanisms, YAML + ARCANE Fortran
for QSS ones. QSS files are copied to mechs/<mech>/customkinetics.f90 with two ARCANE generator artefacts patched: the
floor (TINY) and the 1e3 ceiling on QSS concentrations are removed (they bias the rates), and the strict PLOG bracket
'P > P_i' is made inclusive (otherwise a pressure exactly on a PLOG node matches no bracket). The source terms are
ARCANE's; thermo and the reactor are NGA2's. Detailed/skeletal kinetics are exact (k_b = k_f/K_c, Cantera PLOG).
Building needs CVODE_DIR, PyYAML and NGA2's tools/scripts/chemistry/yaml2nga_lib_cf.py (no Cantera); running/plotting
needs numpy, pandas, matplotlib; only 'check' needs Cantera.
"""
import argparse, json, os, re, shutil, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ATM = 101325.0
NGA_HOME = Path(os.environ.get("NGA_HOME") or next(
    (p for p in (HERE / "../..", HERE / "../../Repositories/NGA2") if (p / "tools/GNUMake/Make.defs").exists()),
    HERE / "../..")).resolve()                     # NGA2/examples/<this>, or a sibling of Repositories/NGA2; else NGA_HOME
CEILING_RE = re.compile(r"^(\s*)(cqss\s*=\s*MIN\(\s*cqss\s*,\s*1e03_pr\s*\))\s*$", re.M | re.I)
FLOOR_RE = re.compile(r"^(\s*)(cqss\s*=\s*MAX\(\s*cqss\s*,\s*TINY\(1\.0_pr\)\s*\))\s*$", re.M | re.I)
PLOG_RE = re.compile(r"IF \( Ploc>P_pdep\(i\) \) THEN")


# ------------------------------------------------------------------------------------------------ config
def read_input(path):
    """NGA2 input file, read as src/libraries/param.f90 does: 'Tag : value' (or 'Tag = value'); text after ! # % is a
    comment; a line without ':' or '=' continues the previous value; commas, semicolons and tabs count as blanks;
    the first occurrence of a tag wins. Returns {tag: [tokens]} in file order."""
    params, tag = {}, None
    for ln in open(path):
        cut = [ln.find(c) for c in "!#%" if c in ln]
        ln = (ln[:min(cut)] if cut else ln).strip()
        if not ln:
            continue
        k = [ln.find(c) for c in ":=" if c in ln]
        if k:
            tag, val = ln[:min(k)].strip(), ln[min(k) + 1:]
            if tag in params:
                print(f"[input] '{tag}' appears more than once: only the first is used")
                tag = None
                continue
            params[tag] = val
        elif tag:
            params[tag] += " " + ln
    return {t: re.sub(r"[\t,;]", " ", v).split() for t, v in params.items()}


def floats(p, tag, default=None):
    if tag in p:
        return [float(v) for v in p[tag]]
    if default is None:
        sys.exit(f"[input] missing '{tag}'")
    return default


def pairs(tokens):
    """'A 0.5 B 0.5' -> {'A': 0.5, 'B': 0.5}"""
    if len(tokens) % 2:
        sys.exit(f"[input] expected '<species> <value>' pairs, got {' '.join(tokens)}")
    return {tokens[i]: float(tokens[i + 1]) for i in range(0, len(tokens), 2)}


def load_config(path):
    """Study settings from the NGA2 input file -> internal dict."""
    p = read_input(path)
    base = Path(path).resolve().parent
    absp = lambda f: str((base / os.path.expanduser(f)).resolve())
    cfg = {"_base": base, "_path": str(path), "name": " ".join(p.get("Name", [])),
           "output_dir": (base / p.get("Output directory", ["results"])[0]).resolve(),
           "workers": int(floats(p, "Workers", [6])[0]), "mechanisms": {}, "fuels": {}, "species_labels": {}}
    for t, v in p.items():
        if t.startswith("Mechanism "):
            m = {"yaml": absp(v[0])}
            if len(v) > 1:
                m["f90"] = absp(v[1])
            cfg["mechanisms"][t.split(None, 1)[1]] = m
        elif t.startswith("Fuel "):
            cfg["fuels"][t.split(None, 1)[1]] = pairs(v)
        elif t.startswith("Label "):
            cfg["species_labels"][t.split(None, 1)[1]] = " ".join(v)
    if not cfg["mechanisms"] or not cfg["fuels"]:
        sys.exit("[input] needs at least one 'Mechanism <name> :' and one 'Fuel <name> :' line")
    cfg["reference"] = p.get("Reference mechanism", [next(iter(cfg["mechanisms"]))])[0]
    cfg["oxidizer"] = pairs(p.get("Oxidizer", ["O2", "1.0", "N2", "3.76"]))
    for scale in (1000.0, 10000.0):                 # '1000/T0 range' (as the plots) or the older '10000/T0 range'
        key = f"{scale:g}/T0"
        if key in p or key + " range" in p:
            break
    else:
        sys.exit("[input] missing '1000/T0 range' (or '1000/T0 : <list>')")
    if key in p:
        inv = floats(p, key)
    else:
        a, b, st = floats(p, key + " range")
        inv = list(np.round(np.arange(a, b + 0.5 * st, st), 10))
    reactor = {"isochor": "isochor", "cv": "isochor", "isobar": "isobar", "cp": "isobar"}.get(
        " ".join(p.get("Reactor", ["isochor"])).lower())
    end_mode = p.get("End mode", ["flamemaster"])[0].lower()
    if reactor is None:
        sys.exit("[input] Reactor must be isochor (cv) or isobar (cp)")
    if end_mode not in ("flamemaster", "fixed", "maxtime"):
        sys.exit("[input] End mode must be flamemaster, fixed or maxtime")
    i = {"reactor": reactor, "pressures_atm": floats(p, "Pressures"),
         "phis": floats(p, "Equivalence ratios"), "T0_K": [scale / v for v in inv],
         "max_time_s": floats(p, "Max time", [1.0e4])[0], "min_time_s": floats(p, "IDT min time", [1.0e-6])[0],
         "min_rise_K": floats(p, "IDT min rise", [55.0])[0], "window_K": floats(p, "IDT window", [400.0])[0],
         "end_mode": end_mode, "rtol": floats(p, "CVODE rtol", [1.0e-10])[0],
         "atol": floats(p, "CVODE atol", [1.0e-20])[0], "timeout_s": floats(p, "Run timeout", [7200.0])[0]}
    cfg["idt"] = i
    e = {}
    if "End state fuels" in p:
        e["fuels"] = p["End state fuels"]
        if set(e["fuels"]) - set(cfg["fuels"]):
            sys.exit(f"[input] End state fuels not defined as 'Fuel <name>': {sorted(set(e['fuels']) - set(cfg['fuels']))}")
    if "End state pressures" in p:
        e["pressures_atm"] = floats(p, "End state pressures")
    if "End state equivalence ratios" in p:
        e["phis"] = floats(p, "End state equivalence ratios")
    if "End state x_min" in p:
        e["x_min"] = floats(p, "End state x_min")[0]
    cfg["endstate"] = e
    return cfg


def mech_dir(name):
    return HERE / "mechs" / name


def executable(name):
    exes = list(mech_dir(name).glob("chem_reactor_0D*.exe"))
    return max(exes, key=lambda p: p.stat().st_mtime) if exes else None


def run_env():
    env = os.environ.copy()
    cv = env.get("CVODE_DIR")
    if cv:
        libs = [str(Path(cv) / d) for d in ("lib", "lib64") if (Path(cv) / d).is_dir()]
        env["LD_LIBRARY_PATH"] = ":".join(libs + [env.get("LD_LIBRARY_PATH", "")]).strip(":")
    env["OMP_NUM_THREADS"] = "1"
    return env


# ------------------------------------------------------------------------------------------------ build
def patch_qss(src):
    n_c = len(CEILING_RE.findall(src)); src = CEILING_RE.sub(r"\1! idt.py: QSS 1e3 ceiling removed\n\1! \2", src)
    n_f = len(FLOOR_RE.findall(src)); src = FLOOR_RE.sub(r"\1! idt.py: QSS floor removed\n\1! \2", src)
    n_p = len(PLOG_RE.findall(src)); src = PLOG_RE.sub("IF ( Ploc>=P_pdep(i) ) THEN", src)
    return src, f"ceiling removed {n_c}, floor removed {n_f}, PLOG brackets fixed {n_p}"


GENERATED = ["build", "chem_data_fc.f90", "chem_data_reactions.txt", "customkinetics.f90", "chem_reactor_0D*.exe",
             "source.txt"]


def up_to_date(name, m):
    """Executable newer than the mechanism files, the reactor source, the GNUmakefile and the generator."""
    exe = executable(name)
    srcs = [Path(m["yaml"]), HERE / "src/simulation.f90", HERE / "GNUmakefile",
            NGA_HOME / "tools/scripts/chemistry/yaml2nga_lib_cf.py"] + ([Path(m["f90"])] if m.get("f90") else [])
    return exe is not None and all(exe.stat().st_mtime > p.stat().st_mtime for p in srcs if p.exists())


def build(cfg, names, jobs=8, force=False):
    sys.path.insert(0, str(NGA_HOME / "tools/scripts/chemistry"))
    from yaml2nga_lib_cf import read_yaml, check_custom_f90, all_rates_zero, write_fcmech
    if not os.environ.get("CVODE_DIR"):
        sys.exit("CVODE_DIR is not set (SUNDIALS install root)")
    for name in names:
        m = cfg["mechanisms"][name]
        if not force and up_to_date(name, m):
            print(f"[{name}] up to date ({executable(name).relative_to(HERE)}); --force to rebuild")
            continue
        d = mech_dir(name)
        d.mkdir(parents=True, exist_ok=True)
        if m.get("f90") and Path(m["f90"]).resolve() == (d / "customkinetics.f90").resolve():
            sys.exit(f"[{name}] rename the ARCANE file: customkinetics.f90 is the name of the generated copy")
        for g in GENERATED:                         # remove what a previous build made; keep the user's inputs
            for f in d.glob(g):
                shutil.rmtree(f) if f.is_dir() else f.unlink()
        custom = None
        if m.get("f90"):
            src, msg = patch_qss(Path(m["f90"]).read_text())
            custom = d / "customkinetics.f90"
            custom.write_text(src)
            print(f"[{name}] ARCANE QSS kinetics {Path(m['f90']).name}: {msg}")
        mech = read_yaml(m["yaml"], transport=False, custom=bool(custom))
        if custom:
            check_custom_f90(custom, mech)
        elif all_rates_zero(mech):                  # placeholder mechanisms need their own rate routine
            sys.exit(f"[{name}] {Path(m['yaml']).name}: all rate constants are zero (placeholder reactions). Its "
                     "rates come from a separate routine; only Cantera YAML kinetics or ARCANE QSS "
                     "(yaml + f90 with customkinetics) are supported")
        info = write_fcmech(mech, d / "chem_data_fc.f90", custom)
        print(f"[{name}] {info['nS']} transported species, " +
              ("ARCANE QSS source terms" if custom else f"{info['nRf']} reactions, exact kinetics"))
        r = subprocess.run(["make", f"MECH={name}", f"-j{jobs}"], cwd=HERE, capture_output=True, text=True)
        if r.returncode != 0 or executable(name) is None:
            print(r.stdout[-3000:], r.stderr[-3000:])
            sys.exit(f"[{name}] build failed")
        print(f"[{name}] built {executable(name).relative_to(HERE)}")


# ------------------------------------------------------------------------------------------------ run
def stream_lines(tag, comp):
    names = list(comp)
    return [f"{tag} : {' '.join(names)}", f"{tag} mole fractions : {' '.join(repr(float(comp[s])) for s in names)}"]


def input_text(cfg, fuel, P_atm, phi, T0):
    i = cfg["idt"]
    L = stream_lines("Fuel", cfg["fuels"][fuel]) + stream_lines("Oxidizer", cfg["oxidizer"])
    L += [f"Equivalence ratio : {float(phi)!r}", f"Temperature : {float(T0)!r}", f"Pressure : {float(P_atm) * ATM!r}",
          f"Reactor : {i.get('reactor', 'isochor')}", f"Max time : {float(i.get('max_time_s', 1e4))!r}",
          f"IDT min time : {float(i.get('min_time_s', 1e-6))!r}", f"IDT min rise : {float(i.get('min_rise_K', 55))!r}",
          f"IDT window : {float(i.get('window_K', 400))!r}", f"End mode : {i.get('end_mode', 'flamemaster')}",
          f"CVODE rtol : {float(i.get('rtol', 1e-10))!r}", f"CVODE atol : {float(i.get('atol', 1e-20))!r}",
          "Result file : idt.res"]
    return "\n".join(L) + "\n"


def ffloat(v):
    """Fortran real, also when the exponent letter was dropped (e.g. 1.5-100)."""
    v = v.strip()
    return float(re.sub(r"(\d)([+-]\d{3})$", r"\1e\2", v))


def read_result(path):
    out, X = {}, {}
    for ln in open(path):
        k, v = ln.split(None, 1)
        if k.startswith("X_"):
            X[k] = ffloat(v)
        elif k == "reactor":
            out[k] = v.strip()
        else:
            out[k] = ffloat(v)
    return out, X


def run_point(cfg, exe, fuel, P_atm, phi, T0, env):
    with tempfile.TemporaryDirectory(dir=cfg["output_dir"]) as d:
        Path(d, "input").write_text(input_text(cfg, fuel, P_atm, phi, T0))
        try:
            r = subprocess.run([str(exe), "-i", "input"], cwd=d, env=env, capture_output=True, text=True,
                               timeout=cfg["idt"].get("timeout_s", 7200))
        except subprocess.TimeoutExpired:
            return {"status": -999, "note": "timeout"}, {}
        log = (r.stdout + r.stderr).strip()
        try:
            return read_result(Path(d, "idt.res"))
        except Exception as e:
            return {"status": -998, "note": (log.splitlines()[-1][:120] if log else f"no result ({e})")}, {}


def check_settings(out, settings, have_results):
    """Results are reused point by point: refuse to mix them with results of different run settings."""
    f = out / "run_settings.json"
    settings = json.loads(json.dumps(settings))
    if have_results and f.exists():
        old = json.loads(f.read_text())
        diff = [k for k in sorted(set(old) | set(settings)) if old.get(k) != settings.get(k)]
        if diff:
            sys.exit(f"[run] {out} holds results made with other settings ({', '.join(diff)}): use another "
                     "'Output directory' or delete it")
    f.write_text(json.dumps(settings, indent=1) + "\n")


def run(cfg, names):
    out = cfg["output_dir"]; out.mkdir(parents=True, exist_ok=True)
    csv = out / "idt.csv"
    old = pd.read_csv(csv) if csv.exists() else pd.DataFrame()
    check_settings(out, {k: v for k, v in cfg["idt"].items() if k not in ("pressures_atm", "phis", "T0_K", "timeout_s")}
                   | {"oxidizer": cfg["oxidizer"], "fuels": cfg["fuels"]}, len(old) > 0)
    done = set()
    if len(old):
        old = old[(old.status == 0) | ~old.mech.isin(names)]           # failed points are run again
        done = set(zip(old.mech, old.fuel, old.P_atm, old.phi, old.T0_K.round(6)))
    i = cfg["idt"]
    jobs = [(m, f, P, phi, T0) for m in names for f in cfg["fuels"] for P in i["pressures_atm"] for phi in i["phis"]
            for T0 in i["T0_K"] if (m, f, P, phi, round(T0, 6)) not in done]
    for m in names:
        if executable(m) is None:
            sys.exit(f"[{m}] no executable: run 'python idt.py build {m} -i {cfg['_path']}'")
    print(f"{len(jobs)} points to run ({len(done)} already in {csv.name})")
    env = run_env()
    rows = []

    def job(j):
        m, f, P, phi, T0 = j
        res, X = run_point(cfg, executable(m), f, P, phi, T0, env)
        idt = res.get("idt_s", -1.0)
        row = dict(mech=m, fuel=f, P_atm=P, phi=phi, T0_K=T0, idt_s=idt if idt > 0 else np.nan,
                   ignited=int(res.get("ignited", 0)), status=int(res.get("status", 0)), note=res.get("note", ""),
                   cpu_s=res.get("cpu_s", np.nan), nsteps=res.get("nsteps", np.nan), t_end=res.get("t_end", np.nan),
                   T_end=res.get("T_end", np.nan), P_end=res.get("P_end", np.nan), rho_end=res.get("rho_end", np.nan))
        row.update(X)
        return row

    with ThreadPoolExecutor(cfg.get("workers", 6)) as ex:
        for n, row in enumerate(ex.map(job, jobs), 1):
            rows.append(row)
            print(f"[{n}/{len(jobs)}] {row['mech']} {row['fuel']} {row['P_atm']:g} atm phi {row['phi']:g} "
                  f"T0 {row['T0_K']:.0f} K: IDT {row['idt_s']:.4g} s, T_end {row['T_end']:.1f} K", flush=True)
            if n % 50 == 0 or n == len(jobs):
                pd.concat([old, pd.DataFrame(rows)], ignore_index=True).to_csv(csv, index=False)
    print(f"wrote {csv}")


# ------------------------------------------------------------------------------------------------ check (Cantera)
def peak_dtdt(t, T, T0, min_time, min_rise):
    t, T = np.asarray(t), np.asarray(T)
    d = np.diff(T) / np.diff(t); tm = 0.5 * (t[1:] + t[:-1])
    q = np.where((tm >= min_time) & (T[1:] > T0 + min_rise))[0]
    return float(tm[q[np.argmax(d[q])]]) if q.size else np.nan


def check(cfg, every):
    import cantera as ct
    ct.suppress_thermo_warnings()
    csv = cfg["output_dir"] / "idt.csv"
    d = pd.read_csv(csv)
    i = cfg["idt"]
    rows = []
    for name, m in cfg["mechanisms"].items():
        if m.get("f90"):
            print(f"[{name}] skipped: ARCANE QSS mechanisms need a Cantera build with ARCANE's custom kinetics")
            continue
        g = ct.Solution(m["yaml"])
        sub = d[d.mech == name].sort_values(["fuel", "P_atm", "phi", "T0_K"]).iloc[::every]
        for _, r in sub.iterrows():
            g.set_equivalence_ratio(r.phi, cfg["fuels"][r.fuel], cfg["oxidizer"]); g.TP = r.T0_K, r.P_atm * ATM
            reac = ct.IdealGasReactor(g) if i.get("reactor", "isochor") == "isochor" else ct.IdealGasConstPressureReactor(g)
            net = ct.ReactorNet([reac]); ts, Ts = [0.0], [r.T0_K]
            while net.time < i.get("max_time_s", 1e4) and reac.T <= r.T0_K + i.get("window_K", 400) and len(ts) < 400000:
                ts.append(net.step()); Ts.append(reac.T)
            idt = peak_dtdt(ts, Ts, r.T0_K, i.get("min_time_s", 1e-6), i.get("min_rise_K", 55))
            net.advance(max(r.t_end, net.time)) if np.isfinite(r.t_end) else None
            rows.append(dict(mech=name, fuel=r.fuel, P_atm=r.P_atm, phi=r.phi, T0_K=r.T0_K, idt_nga=r.idt_s, idt_ct=idt,
                             T_end_nga=r.T_end, T_end_ct=reac.T))
        c = pd.DataFrame([x for x in rows if x["mech"] == name])
        e = 100 * (c.idt_nga / c.idt_ct - 1)
        print(f"[{name}] {len(c)} points: IDT NGA2 vs Cantera mean|e| {e.abs().mean():.3g} %, max|e| {e.abs().max():.3g} %;"
              f" end T max|dT| {(c.T_end_nga - c.T_end_ct).abs().max():.3g} K")
    pd.DataFrame(rows).to_csv(cfg["output_dir"] / "check_cantera.csv", index=False)


# ------------------------------------------------------------------------------------------------ plot
def plot(cfg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter
    out = cfg["output_dir"]
    d = pd.read_csv(out / "idt.csv")
    i = cfg["idt"]                                 # only the points of the current sweep (idt.csv may hold older ones)
    d = d[d.fuel.isin(cfg["fuels"]) & d.P_atm.isin(i["pressures_atm"]) & d.phi.isin(i["phis"]) &
          d.T0_K.round(6).isin([round(t, 6) for t in i["T0_K"]])].copy()
    bad = d.status != 0                            # CVODE failure / timeout: no IDT or end state
    d.loc[bad, [c for c in d.columns if c in ("idt_s", "t_end", "T_end", "P_end", "rho_end") or c.startswith("X_")]] = np.nan
    d.loc[bad, "ignited"] = -1
    mechs = [m for m in cfg["mechanisms"] if m in set(d.mech)]
    ref = cfg.get("reference", mechs[0] if mechs else "")
    if ref not in mechs:
        sys.exit(f"reference '{ref}' has no results for this sweep: python idt.py run {ref}")
    fuels, phis = list(cfg["fuels"]), list(cfg["idt"]["phis"])
    labels = cfg.get("species_labels", {})
    comp = {f: " / ".join(f"{labels.get(k, k)} {v:.3f}" for k, v in sorted(c.items(), key=lambda kv: -kv[1]))
            for f, c in cfg["fuels"].items()}
    styles = [("k", "o", "-"), ("tab:purple", "s", "--"), ("tab:blue", "^", "-."), ("tab:green", "v", ":"),
              ("tab:red", "D", "--"), ("tab:orange", "x", "-."), ("tab:olive", "P", ":"), ("tab:cyan", "*", "--")]
    sty = {m: styles[k % len(styles)] for k, m in enumerate([ref] + [m for m in mechs if m != ref])}
    key = ["fuel", "P_atm", "phi", "T0_K"]
    R = d[d.mech == ref].set_index(key)

    def decades(s):
        s = s[(s > 0) & np.isfinite(s)]
        return 10 ** np.floor(np.log10(s.min())), 10 ** np.ceil(np.log10(s.max()))

    def logaxis(ax, lo, hi):
        ax.set_ylim(lo, hi); ax.yaxis.set_major_locator(LogLocator(base=10, numticks=30))
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=300))
        ax.yaxis.set_minor_formatter(NullFormatter()); ax.grid(True, which="both", alpha=0.3)

    ylo, yhi = decades(d.idt_s); clo, chi = decades(d.cpu_s)
    err = {m: 100 * (d[d.mech == m].set_index(key).idt_s / R.idt_s - 1) for m in mechs if m != ref}
    emax = max([e.abs().max() for e in err.values()] + [0.0])
    elim = min(100, max(20, 10 * np.ceil(emax / 10))) if np.isfinite(emax) else 100
    # IDT grid, one figure per pressure: per phi an IDT, error and CPU row; columns = fuels
    for P in cfg["idt"]["pressures_atm"]:
        fig, axes = plt.subplots(3 * len(phis), len(fuels), figsize=(max(4.3 * len(fuels), 9.0), 7.3 * len(phis)), sharex=True,
                                 squeeze=False, gridspec_kw=dict(height_ratios=[3, 1.3, 1.3] * len(phis)))
        for i, phi in enumerate(phis):
            for j, f in enumerate(fuels):
                ax, ae, ac = axes[3 * i, j], axes[3 * i + 1, j], axes[3 * i + 2, j]
                for m in [ref] + [m for m in mechs if m != ref]:
                    c, mk, ls = sty[m]
                    s = d[(d.mech == m) & (d.fuel == f) & (d.P_atm == P) & (d.phi == phi)].sort_values("T0_K")
                    mfc = "none" if mk in "osD" else None
                    ax.semilogy(1000 / s.T0_K, s.idt_s, ls, marker=mk, color=c, ms=3.2, lw=1.1, mfc=mfc, label=m)
                    ac.semilogy(1000 / s.T0_K, s.cpu_s, ls, marker=mk, color=c, ms=2.6, lw=0.9, mfc=mfc)
                    if m != ref:
                        e = 100 * (s.set_index("T0_K").idt_s / R.loc[(f, P, phi)].idt_s.reindex(s.T0_K) - 1)
                        ae.plot(1000 / e.index, e.clip(-elim, elim), ls, marker=mk, color=c, ms=3, lw=1.0, mfc=mfc)
                logaxis(ax, ylo, yhi); logaxis(ac, clo, chi)
                ae.grid(True, alpha=0.3); ae.axhline(0, color="k", lw=0.7); ae.axhspan(-20, 20, color="tab:green", alpha=0.07)
                ae.set_ylim(-elim, elim)
                if i == 0:
                    ax.set_title(f, fontsize=14, fontweight="bold", pad=34)
                    ax.annotate(comp[f], xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 14),
                                textcoords="offset points", ha="center", va="bottom", fontsize=9)
                if j == 0:
                    ax.set_ylabel(f"IDT [s]  ($\\phi$={phi:g})"); ae.set_ylabel(f"vs {ref} [%]"); ac.set_ylabel("CPU [s]")
                if i == len(phis) - 1:
                    ac.set_xlabel("1000 / T$_0$ [1/K]")
        axes[0, 0].legend(fontsize=8, loc="lower right")
        fig.suptitle(f"Ignition delay (NGA2 0D, {cfg['idt'].get('reactor', 'isochor')}), P = {P:g} atm; IDT = peak dT/dt\n"
                     f"green band $\\pm$20 %, errors clipped at $\\pm${elim:.0f} %; CPU per run (parallel sweep)", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        fig.savefig(out / f"idt_P{P:g}atm.png", dpi=110); plt.close(fig)
    plot_endstate(cfg, d, mechs, ref, sty, comp)
    # summary
    L = [f"# IDT study: {cfg.get('name', '')}", "", f"Reference: {ref}; reactor {cfg['idt'].get('reactor', 'isochor')}, "
         f"IDT = peak dT/dt (t >= {cfg['idt'].get('min_time_s', 1e-6):g} s, T > T0 + {cfg['idt'].get('min_rise_K', 55):g} K, "
         f"until T > T0 + {cfg['idt'].get('window_K', 400):g} K), end mode {cfg['idt'].get('end_mode', 'flamemaster')}.", "",
         "| mech | points | IDT mean abs err % | IDT max abs err % | > 20 % | ignition mismatch | failed | "
         "max abs dT_end K | max abs dCO % | max abs dCO2 % | max abs dH2O % | median CPU s | speed-up vs ref |",
         "|" + "---|" * 13]
    for m in mechs:
        s = d[d.mech == m].set_index(key)
        r = R.reindex(s.index)
        both = (s.ignited == 1) & (r.ignited == 1)
        e = (100 * (s.idt_s / r.idt_s - 1)).abs()
        dX = lambda c: (100 * (s[c] / r[c] - 1))[both].abs().max() if c in s and c in r else np.nan
        valid = (s.ignited >= 0) & (r.ignited >= 0)
        L.append(f"| {m} | {len(s)} | {e.mean():.2f} | {e.max():.1f} | {(e > 20).sum()} | {((s.ignited != r.ignited) & valid).sum()} | "
                 f"{(s.status != 0).sum()} | {(s.T_end - r.T_end)[both].abs().max():.2f} | {dX('X_CO'):.2f} | {dX('X_CO2'):.2f} | "
                 f"{dX('X_H2O'):.2f} | {s.cpu_s.median():.3g} | {(r.cpu_s / s.cpu_s)[both].median():.2f} |")
    L += ["", f"CPU: process CPU time of the integration of each run (median over points); speed-up = median over the "
          f"points ignited in both of CPU({ref}) / CPU(mech). Runs were {cfg.get('workers', 6)} in parallel: "
          f"{'serial, comparable' if cfg.get('workers', 6) == 1 else 'indicative only; set Workers : 1 for a fair timing'}."]
    (out / "summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L)); print(f"figures in {out}")


X_MIN = 1e-10        # end-state species floor (input: End state x_min): smaller mole fractions are not meaningful
                     # (trace species not relaxed to equilibrium, or close to the solver's absolute tolerance).


def plot_endstate(cfg, d, mechs, ref, sty, comp, per_png=28, ncol=4):
    """End state vs T0, one png set per (fuel, P, phi): T, P, rho and the species of the reference mechanism that at
    least one other mechanism also carries (all species if only the reference is plotted) and that some mechanism has
    above x_min in that case; ordered by the reference end mole fraction; species axes from x_min up; deviations only
    where the reference is above x_min. 28 panels per png -> results/endstate/endstate_<fuel>_P<P>atm_phi<phi>[_<n>].png.
    Cases and floor: 'End state ...' lines of the input file; default = whole sweep grid, 1e-10."""
    from concurrent.futures import ProcessPoolExecutor
    out = cfg["output_dir"] / "endstate"; out.mkdir(exist_ok=True)
    for f in out.glob("endstate_*.png"):
        f.unlink()
    sel = cfg.get("endstate", {})
    xcols = [c for c in d.columns if c.startswith("X_")]
    has = {m: {c for c in xcols if d.loc[d.mech == m, c].notna().any()} for m in mechs}
    others = set().union(*[has[m] for m in mechs if m != ref]) if len(mechs) > 1 else has[ref]
    dr = d[(d.mech == ref) & (d.ignited == 1)]
    spc = sorted([c for c in xcols if c in has[ref] and c in others], key=lambda c: -np.nanmedian(dr[c].clip(lower=0)))
    Q = [("T_end", "T$_{end}$ [K]", False, 1.0), ("P_end", "P$_{end}$ [atm]", False, 1 / ATM),
         ("rho_end", "$\\rho_{end}$ [kg/m$^3$]", False, 1.0)] + [(c, f"X$_{{end}}$ {c[2:]}", True, 1.0) for c in spc]
    order = [ref] + [m for m in mechs if m != ref]
    head = dict(reactor=cfg["idt"].get("reactor", "isochor"), end=cfg["idt"].get("end_mode", "flamemaster"),
                x_min=float(sel.get("x_min", X_MIN)))
    cases = [(f, P, phi) for f in sel.get("fuels", list(cfg["fuels"])) for P in sel.get("pressures_atm", cfg["idt"]["pressures_atm"])
             for phi in sel.get("phis", cfg["idt"]["phis"])]
    jobs = [(d[(d.fuel == f) & (d.P_atm == P) & (d.phi == phi)][["mech", "T0_K", "ignited", "T_end", "P_end", "rho_end"] + spc],
             f, P, phi, comp[f], Q, order, ref, sty, {m: sorted(has[m]) for m in mechs}, head, out, per_png, ncol)
            for f, P, phi in cases]
    with ProcessPoolExecutor(cfg.get("workers", 6)) as ex:
        res = list(ex.map(_endstate_case, jobs))
    nq = [r[1] for r in res if r[0]] or [3]
    print(f"end-state figures: {len(spc)} shared species; per case T, P, rho + {min(nq) - 3}-{max(nq) - 3} species above "
          f"x_min = {head['x_min']:g}; {sum(r[0] for r in res)} png in {out}")


def _endstate_case(job):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter
    sub, f, P, phi, comp, Q, order, ref, sty, has, head, out, per_png, ncol = job
    if sub.empty:
        return 0, 0
    has = {m: set(v) for m, v in has.items()}
    xmin = head["x_min"]
    ign_all = sub[sub.ignited == 1]
    keep = [q for q in Q if not q[2] or (ign_all[q[0]] > xmin).any()]
    nskip = len(Q) - len(keep); Q = keep
    R = sub[sub.mech == ref].set_index("T0_K").sort_index()
    npng = int(np.ceil(len(Q) / per_png)); nrow = int(np.ceil(per_png / ncol))
    for ip in range(npng):
        qs = Q[ip * per_png:(ip + 1) * per_png]
        nrow = int(np.ceil(len(qs) / ncol))                    # only as many rows as panels on this page
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 2.75 * nrow + 0.6), squeeze=False)
        axes = axes.ravel()
        for iq, (ax, (col, lab, logy, sc)) in enumerate(zip(axes, qs)):
            ms = [m for m in order if col in ("T_end", "P_end", "rho_end") or col in has[m]]
            ign = sub[(sub.ignited == 1) & sub.mech.isin(ms)][col] * sc
            vals = ign[ign > xmin] if logy else ign.dropna()
            lo, hi = vals.min(), vals.max()
            if logy:
                lo = max(xmin, 10 ** np.floor(np.log10(lo))); hi = max(10 ** np.ceil(np.log10(hi)), 10 * lo)
            else:
                pad = 0.05 * (hi - lo if hi > lo else abs(hi) + 1e-12); lo, hi = lo - pad, hi + pad
            tit = []
            for m in ms:
                c, mk, ls = sty[m]
                s = sub[sub.mech == m].set_index("T0_K").sort_index(); ig = s.ignited == 1
                y = s[col] * sc
                (ax.semilogy if logy else ax.plot)(s.index[ig], y[ig].clip(lower=lo) if logy else y[ig], ls, color=c,
                                                   lw=1.1, marker=mk, ms=2.8, mfc="none" if mk in "osD" else None, label=m)
                ax.plot(s.index[~ig], y[~ig].clip(lo, hi), mk, color=c, ms=6, mfc="none", mew=1.2, clip_on=False)
                if m != ref:
                    r = R[col].reindex(s.index)
                    ok = ig & (R.ignited.reindex(s.index) == 1) & (r * sc > (xmin if logy else -np.inf))
                    rel = (100 * (s[col][ok] / r[ok] - 1)).abs()
                    tit.append(f"{m} {rel.max():.2g}%" if rel.notna().any() else f"{m} -")
            ax.set_ylim(lo, hi)
            if logy:
                ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
                ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=60))
                ax.yaxis.set_minor_formatter(NullFormatter())
            ax.grid(True, which="both", alpha=0.25)
            ax.set_title(lab + ("\nmax|$\\Delta$|: " + ", ".join(tit) if tit else ""), fontsize=8.5)
            ax.tick_params(labelsize=7.5)
            if iq >= len(qs) - ncol:
                ax.set_xlabel("T$_0$ [K]", fontsize=8.5)
        for ax in axes[len(qs):]:
            ax.axis("off")
        h, l = axes[0].get_legend_handles_labels()
        fig.legend(h, l, loc="upper right", bbox_to_anchor=(0.995, 0.997), ncol=len(order), fontsize=10, frameon=False)
        fig.suptitle(f"{f} ({comp}), P$_0$ = {P:g} atm, $\\phi$ = {phi:g}: end state vs T$_0$ ({head['reactor']}, "
                     f"{head['end']} end)   [{ip + 1}/{npng}]\nOpen markers: no ignition. Species with X > {xmin:g} in some mechanism "
                     f"({nskip} below not shown); max|$\\Delta$| vs {ref} where X$_{{{ref}}}$ > {xmin:g}", fontsize=10.5, x=0.01, y=0.997, ha="left")
        fig.subplots_adjust(left=0.05, right=0.99, top=1 - 0.95 / (2.75 * nrow + 0.6), bottom=0.6 / (2.75 * nrow + 0.6),
                            wspace=0.3, hspace=0.62)
        name = f"endstate_{f}_P{P:g}atm_phi{phi:g}" + (f"_{ip + 1}" if npng > 1 else "") + ".png"
        fig.savefig(out / name, dpi=100); plt.close(fig)
    return npng, len(Q)

# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["build", "run", "plot", "check", "all"])
    ap.add_argument("mechs", nargs="*", help="mechanisms (default: all in the input file)")
    ap.add_argument("-i", "--input", default="input", help="study input file (default: input)")
    ap.add_argument("--every", type=int, default=5, help="check: use every N-th point")
    ap.add_argument("-j", type=int, default=8, help="build: make -j")
    ap.add_argument("--force", action="store_true", help="build: rebuild even if up to date")
    a = ap.parse_args()
    cfg = load_config(a.input)
    names = a.mechs or list(cfg["mechanisms"])
    unknown = [n for n in names if n not in cfg["mechanisms"]]
    if unknown:
        sys.exit(f"unknown mechanisms {unknown}")
    if a.command in ("build", "all"):
        build(cfg, names, a.j, a.force)
    if a.command in ("run", "all"):
        run(cfg, names)
    if a.command in ("plot", "all"):
        plot(cfg)
    if a.command == "check":
        check(cfg, a.every)


if __name__ == "__main__":
    main()
