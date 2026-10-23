#!/usr/bin/env python3
"""Laminar flame speed study with the NGA2 1D flame solver (chem_reactor_1D_lfs).

    python lfs.py build   [mech ...]          generate mechs/<mech>/chem_data_fc.f90 and compile
    python lfs.py run     [mech ...]          flame-speed sweep -> results/lfs.csv (resumes, retries failures)
    python lfs.py plot                        results/lfs_P<P>atm.png + results/summary.md
    python lfs.py compare <cantera lfs.csv>   NGA2 vs a Cantera lfs.csv (columns mech, fuel, P_atm, phi, sl_cm_s,
                                              cpu_s, npts; same mechanism names)
    python lfs.py all                         build + run + plot
    -i <file> selects the study input file (default: input, NGA2 input format; see that file).

Sweep: per (fuel, P) the reference mechanism starts cold at the first phi >= 1 and continues up to the richest and
then down to the leanest phi, each flame seeded by its neighbour; every other
mechanism is seeded by the reference flame at the same phi. A failed seeded solve falls back to cold starts on the
domain widths of the input file ('Cold start widths'). Converged flames: results/solutions/<mech>_<fuel>_P<P>_phi<phi>.dat (z u T Y_k).
"""
import argparse, json, os, re, shutil, subprocess, sys, tempfile, threading
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
GENERATED = ["build", "chem_data_fc.f90", "chem_data_reactions.txt", "customkinetics.f90", "chem_reactor_1D_lfs*.exe"]
HEADER = ["mech", "fuel", "P_atm", "phi", "sl_cm_s", "cpu_s", "npts", "seed", "note"]


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
    cfg = {"_path": str(path), "name": " ".join(p.get("Name", [])),
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
    if "Equivalence ratios" in p:
        phis = floats(p, "Equivalence ratios")
    else:
        a, b, st = floats(p, "Equivalence ratio range")
        phis = [float(x) for x in np.round(np.arange(a, b + 0.5 * st, st), 10)]
    cfg["lfs"] = {"pressures_atm": floats(p, "Pressures"), "phis": phis,
                  "T0_K": floats(p, "Unburnt temperature", [300.0])[0],
                  "refine": {"ratio": floats(p, "Refine ratio", [3.0])[0], "slope": floats(p, "Refine slope", [0.08])[0],
                             "curve": floats(p, "Refine curve", [0.08])[0], "prune": floats(p, "Refine prune", [0.05])[0]},
                  "max_grid_points": int(floats(p, "Max grid points", [400])[0]),
                  "convection": p.get("Convection scheme", ["upwind"])[0],
                  "widths_m": floats(p, "Cold start widths", [0.04, 0.08, 0.12]),
                  "sl_max_cm_s": floats(p, "Max flame speed", [200.0])[0],
                  "timeout_s": floats(p, "Run timeout", [5400.0])[0]}
    return cfg


def mech_dir(name):
    return HERE / "mechs" / name


def executable(name):
    exes = sorted(mech_dir(name).glob("chem_reactor_1D_lfs*.opt*.exe"))
    return exes[-1] if exes else None


def run_env():
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    return env


_SPECIES = {}


def mech_species(m):
    """Species of the mechanism's YAML phase, in mechanism order, and the lowest species T_max (no Cantera needed)."""
    key = m["yaml"]
    if key not in _SPECIES:
        sys.path.insert(0, str(NGA_HOME / "tools/scripts/chemistry"))
        from yaml2nga_lib_cf import read_yaml
        mech = read_yaml(key, transport=False, custom=True)
        _SPECIES[key] = (list(mech["species"]), float(np.min(mech["Tlim"][1])))
    return _SPECIES[key]


def stream_lines(tag, comp):
    names = list(comp)
    return [f"{tag} : {' '.join(names)}", f"{tag} mole fractions : {' '.join(repr(float(comp[s])) for s in names)}"]


# ------------------------------------------------------------------------------------------------ build
def patch_qss(src):
    n_c = len(CEILING_RE.findall(src)); src = CEILING_RE.sub(r"\1! lfs.py: QSS 1e3 ceiling removed\n\1! \2", src)
    n_f = len(FLOOR_RE.findall(src)); src = FLOOR_RE.sub(r"\1! lfs.py: QSS floor removed\n\1! \2", src)
    n_p = len(PLOG_RE.findall(src)); src = PLOG_RE.sub("IF ( Ploc>=P_pdep(i) ) THEN", src)
    return src, f"ceiling removed {n_c}, floor removed {n_f}, PLOG brackets fixed {n_p}"


def up_to_date(name, m):
    exe = executable(name)
    srcs = [Path(m["yaml"]), HERE / "src/simulation.f90", HERE / "GNUmakefile",
            NGA_HOME / "tools/scripts/chemistry/yaml2nga_lib_cf.py"] + ([Path(m["f90"])] if m.get("f90") else [])
    return exe is not None and all(exe.stat().st_mtime > p.stat().st_mtime for p in srcs if p.exists())


def build(cfg, names, jobs=8, force=False):
    sys.path.insert(0, str(NGA_HOME / "tools/scripts/chemistry"))
    from yaml2nga_lib_cf import read_yaml, check_custom_f90, all_rates_zero, write_fcmech
    for name in names:
        m = cfg["mechanisms"][name]
        if not force and up_to_date(name, m):
            print(f"[{name}] up to date ({executable(name).relative_to(HERE)}); --force to rebuild")
            continue
        d = mech_dir(name)
        d.mkdir(parents=True, exist_ok=True)
        if m.get("f90") and Path(m["f90"]).resolve() == (d / "customkinetics.f90").resolve():
            sys.exit(f"[{name}] rename the ARCANE file: customkinetics.f90 is the name of the generated copy")
        for g in GENERATED:
            for f in d.glob(g):
                shutil.rmtree(f) if f.is_dir() else f.unlink()
        custom = None
        if m.get("f90"):
            src, msg = patch_qss(Path(m["f90"]).read_text())
            custom = d / "customkinetics.f90"
            custom.write_text(src)
            print(f"[{name}] ARCANE QSS kinetics {Path(m['f90']).name}: {msg}")
        mech = read_yaml(m["yaml"], transport=True, custom=bool(custom))
        if custom:
            check_custom_f90(custom, mech)
        elif all_rates_zero(mech):
            sys.exit(f"[{name}] all rate constants are zero (placeholder reactions): not supported")
        info = write_fcmech(mech, d / "chem_data_fc.f90", custom)
        if not info["transport"]:
            sys.exit(f"[{name}] the mechanism has no transport data")
        print(f"[{name}] {info['nS']} transported species, " +
              ("ARCANE QSS source terms" if custom else f"{info['nRf']} reactions, exact kinetics"))
        r = subprocess.run(["make", f"MECH={name}", f"-j{jobs}"], cwd=HERE, capture_output=True, text=True)
        if r.returncode != 0 or executable(name) is None:
            print(r.stdout[-3000:], r.stderr[-3000:])
            sys.exit(f"[{name}] build failed")
        print(f"[{name}] built {executable(name).relative_to(HERE)}")


# ------------------------------------------------------------------------------------------------ one flame
def sol_path(cfg, mech, fuel, P, phi):
    return cfg["output_dir"] / "solutions" / f"{mech}_{fuel}_P{P:g}_phi{phi:g}.dat"


def read_solution(path):
    names = open(path).readline().split()
    return names[3:], np.loadtxt(path, skiprows=1, ndmin=2)


def flame(cfg, mech, fuel, P, phi, seed=None, width=None):
    """One flame. Returns dict(ok, sl, cpu, npts, note); writes the solution file on success.
    The Fortran builds the unburnt mixture (Fuel / Oxidizer / Equivalence ratio, as the 0D reactor) and, for cold
    starts, its adiabatic equilibrium; no Cantera at run time."""
    c = cfg["lfs"]; m = cfg["mechanisms"][mech]
    species, tmax_sp = mech_species(m)
    out = sol_path(cfg, mech, fuel, P, phi); out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=cfg["output_dir"]) as d:
        d = Path(d)
        r = c.get("refine", {})
        L = stream_lines("Fuel", cfg["fuels"][fuel]) + stream_lines("Oxidizer", cfg["oxidizer"])
        L += [f"Equivalence ratio : {float(phi)!r}", f"Temperature : {float(c['T0_K'])!r}", f"Pressure : {P * ATM!r}",
              f"T upper bound : {2.0 * tmax_sp!r}",
              f"Refine ratio : {float(r.get('ratio', 3.0))!r}", f"Refine slope : {float(r.get('slope', 0.08))!r}",
              f"Refine curve : {float(r.get('curve', 0.08))!r}", f"Refine prune : {float(r.get('prune', 0.05))!r}",
              f"Max grid points : {int(c.get('max_grid_points', 400))}", f"Convection scheme : {c.get('convection', 'upwind')}",
              "Result file : lfs.res", "Solution file : flame.dat"]
        if seed:                                       # map the seed's species onto this mechanism (missing -> 0)
            names, a = read_solution(seed)
            idx = {s: i for i, s in enumerate(names)}
            cols = [a[:, 3 + idx[s]] if s in idx else np.zeros(len(a)) for s in species]
            prof = np.column_stack([a[:, :3]] + cols)
            with open(d / "profile.dat", "w") as f:
                f.write(f"{len(prof)}\n")
                np.savetxt(f, prof, fmt="%.17e")
            L.append("Initial profile : profile.dat")
        else:                                          # cold start (equilibrium computed by the Fortran)
            L.append(f"Width : {float(width)!r}")
        (d / "input").write_text("\n".join(L) + "\n")
        try:
            p = subprocess.run([str(executable(mech)), "-i", "input"], cwd=d, env=run_env(), capture_output=True,
                               text=True, timeout=c.get("timeout_s", 5400))
        except subprocess.TimeoutExpired:
            return dict(ok=False, note="timeout")
        if not (d / "lfs.res").exists():
            log = (p.stdout + p.stderr).strip().splitlines()
            return dict(ok=False, note=(log[-1] if log else "crashed")[:150])
        res = {k: v for k, v in (ln.split(None, 1) for ln in open(d / "lfs.res"))}
        sl = float(res["sl_m_s"]) * 100.0 if res["status"].strip() == "0" else float("nan")
        info = dict(cpu=float(res["cpu_s"]), npts=int(res["npts"]))
        if np.isfinite(sl) and 0.0 < sl <= c.get("sl_max_cm_s", 200.0):
            shutil.copy(d / "flame.dat", out)
            return dict(ok=True, sl=sl, **info, note="")
        return dict(ok=False, **info, note=f"status {res['status'].strip()}, S_L {sl:.4g} cm/s")


def flame_with_fallback(cfg, mech, fuel, P, phi, seed):
    notes = []
    if seed and Path(seed).is_file():
        r = flame(cfg, mech, fuel, P, phi, seed=seed)
        if r["ok"]:
            return r, Path(seed).name
        notes.append(f"seeded: {r['note']}")
    for w in cfg["lfs"]["widths_m"]:
        r = flame(cfg, mech, fuel, P, phi, width=w)
        if r["ok"]:
            r["note"] = " | ".join(notes); return r, f"fresh (w={w:g})"
        notes.append(f"fresh w={w:g}: {r['note']}")
    return dict(ok=False, note=" | ".join(notes)[:300]), "fresh"


# ------------------------------------------------------------------------------------------------ run
def check_settings(out, settings, have_results):
    """Results are reused point by point: refuse to mix them with results of different solver settings."""
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
    csv = out / "lfs.csv"
    lock = threading.Lock()
    old = pd.read_csv(csv) if csv.exists() else pd.DataFrame(columns=HEADER)
    old = old[old.sl_cm_s.notna() | ~old.mech.isin(names)]           # failed points are run again
    check_settings(out, {k: cfg["lfs"][k] for k in ("T0_K", "refine", "max_grid_points", "convection")}
                   | {"oxidizer": cfg["oxidizer"], "fuels": cfg["fuels"]}, len(old) > 0)
    old.to_csv(csv, index=False)
    done = set(zip(old.mech, old.fuel, old.P_atm.astype(float), old.phi.astype(float).round(6)))
    ref = cfg.get("reference", list(cfg["mechanisms"])[0])
    c = cfg["lfs"]; phis = sorted(c["phis"])
    order = [p for p in phis if p >= 1.0 - 1e-9] + [p for p in phis[::-1] if p < 1.0 - 1e-9]

    def record(mech, fuel, P, phi, r, seed):
        with lock:
            row = [mech, fuel, P, phi, r.get("sl", np.nan) if r["ok"] else np.nan, r.get("cpu", np.nan),
                   r.get("npts", 0), seed, r.get("note", "")]
            pd.DataFrame([row], columns=HEADER).to_csv(csv, mode="a", header=False, index=False)
        print(f"[lfs] {mech:10s} {fuel:6s} P={P:g} phi={phi:4.2f}: " +
              (f"{r['sl']:8.3f} cm/s  {r['npts']} pts  cpu {r['cpu']:7.1f} s  ({seed})" if r["ok"] else f"FAILED {r['note']}"),
              flush=True)

    def chain(fuel, P):                                  # reference mechanism: phi = 1 up, then down
        prev = None
        for phi in order:
            if phi < 1.0 - 1e-9 and (prev is None or prev >= 1.0 - 1e-9):
                prev = 1.0
            seed = sol_path(cfg, ref, fuel, P, prev) if prev is not None else None
            if (ref, fuel, float(P), round(phi, 6)) not in done:
                r, s = flame_with_fallback(cfg, ref, fuel, P, phi, seed if seed and seed.is_file() else None)
                record(ref, fuel, P, phi, r, s)
            if sol_path(cfg, ref, fuel, P, phi).is_file():
                prev = phi

    def point(args):
        mech, fuel, P, phi = args
        r, s = flame_with_fallback(cfg, mech, fuel, P, phi, sol_path(cfg, ref, fuel, P, phi))
        record(mech, fuel, P, phi, r, s)

    for m in names:
        if executable(m) is None:
            sys.exit(f"[{m}] no executable: run 'python lfs.py build {m} -i {cfg['_path']}'")
    workers = cfg.get("workers", 6)
    if ref in names:
        with ThreadPoolExecutor(workers) as ex:
            list(ex.map(lambda fp: chain(*fp), [(f, P) for f in cfg["fuels"] for P in c["pressures_atm"]]))
    others = [(m, f, P, phi) for m in names if m != ref for f in cfg["fuels"] for P in c["pressures_atm"]
              for phi in order if (m, f, float(P), round(phi, 6)) not in done]
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(point, others))
    print(f"wrote {csv}")


# ------------------------------------------------------------------------------------------------ plot / compare
def plot(cfg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = cfg["output_dir"]
    d = pd.read_csv(out / "lfs.csv")
    ref = cfg.get("reference", list(cfg["mechanisms"])[0])
    mechs = [m for m in cfg["mechanisms"] if m in set(d.mech)]
    fuels = list(cfg["fuels"])
    labels = cfg.get("species_labels", {})
    comp = {f: " / ".join(f"{labels.get(k, k)} {v:.3f}" for k, v in sorted(cf.items(), key=lambda kv: -kv[1]))
            for f, cf in cfg["fuels"].items()}
    colors = ["k", "tab:purple", "tab:blue", "tab:green", "tab:red", "tab:orange", "tab:olive", "tab:cyan"]
    marks = ["o", "s", "^", "v", "D", "x", "P", "*"]
    sty = {m: (colors[i % 8], marks[i % 8], "-" if m == ref else "--") for i, m in enumerate([ref] + [x for x in mechs if x != ref])}
    key = ["fuel", "P_atm", "phi"]
    R = d[d.mech == ref].set_index(key)
    for P in cfg["lfs"]["pressures_atm"]:
        fig, axes = plt.subplots(2, len(fuels), figsize=(max(4.3 * len(fuels), 9.0), 7.5), sharex=True, squeeze=False,
                                 gridspec_kw=dict(height_ratios=[3, 1.4]))
        for j, f in enumerate(fuels):
            ax, ae = axes[0, j], axes[1, j]
            for m in [ref] + [x for x in mechs if x != ref]:
                c, mk, ls = sty[m]
                s = d[(d.mech == m) & (d.fuel == f) & (d.P_atm == P)].sort_values("phi")
                ax.plot(s.phi, s.sl_cm_s, ls, marker=mk, color=c, ms=3.5, lw=1.2, mfc="none" if m == ref else None, label=m)
                if m != ref:
                    r = R.reindex(list(zip(s.fuel, s.P_atm, s.phi))).sl_cm_s.values
                    ae.plot(s.phi, 100 * (s.sl_cm_s.values / r - 1), ls, marker=mk, color=c, ms=3, lw=1.0)
            ax.set_title(f"{f}\n{comp[f]}", fontsize=10); ax.grid(True, alpha=0.3)
            ae.axhline(0, color="k", lw=0.7); ae.axhspan(-5, 5, color="tab:green", alpha=0.07); ae.grid(True, alpha=0.3)
            ae.set_xlabel("$\\phi$")
            if j == 0:
                ax.set_ylabel("S$_L$ [cm/s]"); ae.set_ylabel(f"vs {ref} [%]")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(f"Laminar flame speed (NGA2 1D), P = {P:g} atm, T$_0$ = {cfg['lfs']['T0_K']:g} K "
                     f"(green band $\\pm$5 %)", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(out / f"lfs_P{P:g}atm.png", dpi=110); plt.close(fig)
    L = [f"# LFS study: {cfg.get('name', '')}", "", f"Reference: {ref}. Error = 100 x (S_L,mech - S_L,ref) / S_L,ref. "
         f"Convection: {cfg['lfs'].get('convection', 'upwind')}.", "",
         "| mech | points | failed | S_L mean abs err % | S_L max abs err % | median CPU s | speed-up vs ref | median points |",
         "|" + "---|" * 8]
    for m in mechs:
        s = d[d.mech == m].set_index(key); r = R.reindex(s.index)
        e = (100 * (s.sl_cm_s / r.sl_cm_s - 1)).abs()
        sp = (r.cpu_s / s.cpu_s)[s.sl_cm_s.notna() & r.sl_cm_s.notna()].median()
        L.append(f"| {m} | {len(s)} | {s.sl_cm_s.isna().sum()} | {e.mean():.2f} | {e.max():.2f} | {s.cpu_s.median():.3g} | "
                 f"{sp:.2f} | {s.npts.median():.0f} |")
    L += ["", "CPU: process CPU time of each solve; seeded solves start from different flames (reference: neighbour phi, "
          "others: same-phi reference), so CPU times are indicative only."]
    (out / "summary.md").write_text("\n".join(L) + "\n"); print("\n".join(L))


def compare(cfg, ct_csv):
    d = pd.read_csv(cfg["output_dir"] / "lfs.csv"); c = pd.read_csv(ct_csv)
    key = ["mech", "fuel", "P_atm", "phi"]
    for x in (d, c):
        x["phi"] = x.phi.astype(float).round(6); x["P_atm"] = x.P_atm.astype(float)
    m = d.merge(c[key + ["sl_cm_s", "cpu_s", "npts"]], on=key, suffixes=("_nga", "_ct"))
    L = ["| mech | points | both converged | mean abs diff % | max abs diff % | median npts NGA2 / Cantera | "
         "median CPU NGA2 / Cantera [s] |", "|" + "---|" * 7]
    for mech, s in m.groupby("mech", sort=False):
        ok = s.sl_cm_s_nga.notna() & s.sl_cm_s_ct.notna()
        e = (100 * (s.sl_cm_s_nga / s.sl_cm_s_ct - 1))[ok].abs()
        L.append(f"| {mech} | {len(s)} | {ok.sum()} | {e.mean():.4f} | {e.max():.4f} | {s.npts_nga.median():.0f} / "
                 f"{s.npts_ct.median():.0f} | {s.cpu_s_nga.median():.3g} / {s.cpu_s_ct.median():.3g} |")
    worst = m.assign(e=100 * (m.sl_cm_s_nga / m.sl_cm_s_ct - 1)).dropna(subset=["e"]).sort_values("e", key=abs).tail(5)
    L += ["", "Largest differences (e = NGA2 vs Cantera, %):", "```",
          worst[key + ["sl_cm_s_nga", "sl_cm_s_ct", "e", "npts_nga", "npts_ct"]].to_string(index=False), "```"]
    (cfg["output_dir"] / "compare_cantera.md").write_text("\n".join(L) + "\n"); print("\n".join(L))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["build", "run", "plot", "compare", "all"])
    ap.add_argument("args", nargs="*", help="mechanisms (build/run) or the Cantera lfs.csv (compare)")
    ap.add_argument("-i", "--input", default="input", help="study input file (default: input)")
    ap.add_argument("-j", type=int, default=8, help="build: parallel make jobs")
    ap.add_argument("--force", action="store_true", help="build: rebuild even if up to date")
    a = ap.parse_args()
    cfg = load_config(a.input)
    if a.command == "compare":
        if len(a.args) != 1:
            sys.exit("compare needs one argument: a Cantera lfs.csv")
        compare(cfg, a.args[0]); return
    names = a.args or list(cfg["mechanisms"])
    if a.command in ("build", "all"):
        build(cfg, names, a.j, a.force)
    if a.command in ("run", "all"):
        run(cfg, names)
    if a.command in ("plot", "all"):
        plot(cfg)


if __name__ == "__main__":
    main()
