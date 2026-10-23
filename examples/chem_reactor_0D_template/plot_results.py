#!/usr/bin/env python3
"""Time evolution of one 0D case for every mechanism of the study input (detailed / skeletal / reduced in the same panels).

    python plot_results.py LP40 1 1 650 [-i input] [--mechs detailed reduced] [--x-min 1e-15] [--rows 20000] [--keep]
                           <fuel> <P [atm]> <phi> <T0 [K]>

Runs each mechanism's executable (mechs/<mech>/, built by "python idt.py build") with the reactor settings of the input
and a trajectory output, then plots T, P, rho and the shared species vs time:
results/history/history_<fuel>_P<P>atm_phi<phi>_T<T0>K[_<n>].png, 28 panels per png (--keep: trajectories and idt.res
kept in results/history/<case>/<mech>/). Dotted lines: IDT of each mechanism. --rows should exceed the number of CVODE
steps (printed) so that no step is dropped from the trajectory (short ignition peaks).
Species: those of the reference that another mechanism also carries and whose mole fraction exceeds x_min at some
time in some mechanism by at least a factor 10 (default x_min 1e-15: with CVODE atol 1e-20 on mass fractions, values
above ~1e-15 are resolved to better than 1e-5 relative); species axes start at x_min; panel title: difference of the
peak mole fraction vs the reference.
"""
import argparse, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

import idt

X_MIN = 1e-15


def read_traj(path):
    """chem_reactor_0D trajectory: time, <sp>-X (or -Y), T, P, rho."""
    with open(path) as f:
        head = f.readline().split()
        rows = [re.sub(r"(\d)([+-]\d{3})(?=\s|$)", r"\1E\2", ln) for ln in f]
    a = np.loadtxt(rows, ndmin=2)
    sp = [h[:-2] for h in head[1:-3]]
    return dict(t=a[:, 0], X=dict(zip(sp, a[:, 1:-3].T)), T=a[:, -3], P=a[:, -2], rho=a[:, -1])


def run_case(cfg, mech, fuel, P, phi, T0, rows, case_dir, env):
    d = case_dir / mech
    d.mkdir(parents=True, exist_ok=True)
    txt = idt.input_text(cfg, fuel, P, phi, T0) + \
        f"Output file : traj.out\nNumber of outputs : {rows}\nOutput mole fractions : .true.\n"
    (d / "input").write_text(txt)
    r = subprocess.run([str(idt.executable(mech)), "-i", "input"], cwd=d, env=env, capture_output=True, text=True,
                       timeout=cfg["idt"].get("timeout_s", 7200))
    if not (d / "traj.out").exists():
        sys.exit(f"[{mech}] run failed:\n{(r.stdout + r.stderr)[-1500:]}")
    res, _ = idt.read_result(d / "idt.res")
    return mech, read_traj(d / "traj.out"), res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fuel"); ap.add_argument("P_atm", type=float)
    ap.add_argument("phi", type=float); ap.add_argument("T0", type=float)
    ap.add_argument("-i", "--input", default="input", help="study input file (default: input)")
    ap.add_argument("--mechs", nargs="*", help="mechanisms (default: all in the input file)")
    ap.add_argument("--x-min", type=float, default=X_MIN, help="species floor (default 1e-15)")
    ap.add_argument("--rows", type=int, default=20000, help="trajectory rows kept per run (default 20000)")
    ap.add_argument("--keep", action="store_true", help="keep trajectories in results/history/<case>/")
    ap.add_argument("--per-png", type=int, default=28)
    ap.add_argument("--time-axis", choices=["linear", "log"], default="linear",
                    help="linear (default, 0..t_end: ignition and the approach to equilibrium, FlameMaster style) or "
                         "log (early low-temperature chemistry)")
    a = ap.parse_args()
    cfg = idt.load_config(a.input)
    mechs = a.mechs or list(cfg["mechanisms"])
    ref = cfg.get("reference", mechs[0]) if cfg.get("reference", mechs[0]) in mechs else mechs[0]
    for m in mechs:
        if idt.executable(m) is None:
            sys.exit(f"[{m}] no executable: run 'python idt.py build {m} -i {a.input}'")
    tag = f"{a.fuel}_P{a.P_atm:g}atm_phi{a.phi:g}_T{a.T0:g}K"
    out = cfg["output_dir"] / "history"; case_dir = out / tag
    out.mkdir(parents=True, exist_ok=True)
    env = idt.run_env()
    try:
        with ThreadPoolExecutor(len(mechs)) as ex:
            runs = {m: (tr, res) for m, tr, res in
                    ex.map(lambda m: run_case(cfg, m, a.fuel, a.P_atm, a.phi, a.T0, a.rows, case_dir, env), mechs)}
    finally:
        if not a.keep:
            shutil.rmtree(case_dir, ignore_errors=True)
    for m, (tr, res) in runs.items():
        thin = "" if len(tr["t"]) > res["nsteps"] else f"  (thinned: {res['nsteps']:.0f} steps > --rows)"
        print(f"[{m}] {len(tr['t'])} rows, IDT {res['idt_s']:.5g} s, end t {res['t_end']:.4g} s, T_end {res['T_end']:.2f} K{thin}")

    # quantities: T, P, rho + shared species above x_min at some time, ordered by the reference peak mole fraction
    order = [ref] + [m for m in mechs if m != ref]
    has = {m: set(runs[m][0]["X"]) for m in mechs}
    others = set().union(*[has[m] for m in mechs if m != ref]) if len(mechs) > 1 else has[ref]
    peak = lambda m, s: np.nanmax(runs[m][0]["X"][s]) if s in has[m] else np.nan
    spc = [s for s in runs[ref][0]["X"] if s in others and max(np.nan_to_num(peak(m, s)) for m in mechs) > 10 * a.x_min]
    spc.sort(key=lambda s: -np.nan_to_num(peak(ref, s)))
    nskip = len([s for s in runs[ref][0]["X"] if s in others]) - len(spc)
    Q = [("T", "T [K]", False, 1.0), ("P", "P [atm]", False, 1 / idt.ATM), ("rho", "$\\rho$ [kg/m$^3$]", False, 1.0)] + \
        [(s, f"X {s}", True, 1.0) for s in spc]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter
    styles = [("k", "-"), ("tab:purple", "--"), ("tab:blue", "-."), ("tab:green", ":"), ("tab:red", "--"), ("tab:orange", "-.")]
    sty = {m: styles[k % len(styles)] for k, m in enumerate(order)}
    idts = {m: runs[m][1]["idt_s"] for m in order}
    t_lo = min(max(tr["t"][tr["t"] > 0].min(), 1e-9) for tr, _ in runs.values())
    t_hi = max(tr["t"][-1] for tr, _ in runs.values())
    ncol = 4
    npng = int(np.ceil(len(Q) / a.per_png))
    for old in out.glob(f"history_{tag}*.png"):
        old.unlink()
    comp = " / ".join(f"{cfg.get('species_labels', {}).get(k, k)} {v:.3f}"
                      for k, v in sorted(cfg["fuels"][a.fuel].items(), key=lambda kv: -kv[1]))
    for ip in range(npng):
        qs = Q[ip * a.per_png:(ip + 1) * a.per_png]
        nrow = int(np.ceil(len(qs) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 2.75 * nrow + 0.6), squeeze=False)
        axes = axes.ravel()
        for iq, (ax, (key, lab, logy, sc)) in enumerate(zip(axes, qs)):
            ms = [m for m in order if not logy or key in has[m]]
            if logy:
                ax.set_yscale("log")
            ys = {m: (runs[m][0]["X"][key] if logy else runs[m][0][key]) * sc for m in ms}
            for m in ms:
                tr = runs[m][0]; k = tr["t"] > 0; c, ls = sty[m]
                y = ys[m][k]
                ax.plot(tr["t"][k], np.clip(y, a.x_min, None) if logy else y, ls, color=c, lw=1.3, label=m)
                if np.isfinite(idts[m]) and idts[m] > 0:
                    ax.axvline(idts[m], color=c, ls=":", lw=0.8)
            if logy:
                top = max(np.nanmax(v) for v in ys.values())
                lo = max(a.x_min, 10 ** np.floor(np.log10(max(min(np.nanmin(v[v > a.x_min]) if (v > a.x_min).any() else top
                                                                    for v in ys.values()), a.x_min))))
                ax.set_ylim(lo, max(10 ** np.ceil(np.log10(top)), 10 * lo))
                ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
                ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=60))
                ax.yaxis.set_minor_formatter(NullFormatter())
                tit = [f"{m} {100 * (np.nanmax(ys[m]) / np.nanmax(ys[ref]) - 1):+.2g}%" for m in ms if m != ref and ref in ms]
                lab = lab + ("\npeak vs " + ref + ": " + ", ".join(tit) if tit else "")
            if a.time_axis == "log":
                ax.set_xscale("log"); ax.set_xlim(t_lo, t_hi)
            else:
                ax.set_xlim(0.0, t_hi); ax.ticklabel_format(axis="x", style="sci", scilimits=(-2, 3))
            ax.grid(True, which="both", alpha=0.25)
            ax.set_title(lab, fontsize=8.5); ax.tick_params(labelsize=7.5)
            if iq >= len(qs) - ncol:
                ax.set_xlabel("t [s]", fontsize=8.5)
        for ax in axes[len(qs):]:
            ax.axis("off")
        h, l = axes[0].get_legend_handles_labels()
        fig.legend(h, l, loc="upper right", bbox_to_anchor=(0.995, 0.997), ncol=len(order), fontsize=10, frameon=False)
        fig.suptitle(f"{a.fuel} ({comp}), P$_0$ = {a.P_atm:g} atm, $\\phi$ = {a.phi:g}, T$_0$ = {a.T0:g} K: time evolution "
                     f"({cfg['idt'].get('reactor', 'isochor')}, {a.time_axis} time axis to the end of the run)   [{ip + 1}/{npng}]\nDotted: IDT "
                     + ", ".join(f"{m} {v:.4g} s" for m, v in idts.items()) +
                     f". Species whose peak X exceeds 10 x {a.x_min:g} ({nskip} below not shown); axes from {a.x_min:g}", fontsize=10.5,
                     x=0.01, y=0.997, ha="left")
        fig.subplots_adjust(left=0.05, right=0.99, top=1 - 0.95 / (2.75 * nrow + 0.6), bottom=0.6 / (2.75 * nrow + 0.6),
                            wspace=0.3, hspace=0.62)
        name = f"history_{tag}" + (f"_{ip + 1}" if npng > 1 else "") + ".png"
        fig.savefig(out / name, dpi=100); plt.close(fig)
        print(f"saved {out / name}")
    print(f"{len(spc)} shared species with peak X > 10 x_min = {10 * a.x_min:g} ({nskip} below), {npng} png")


if __name__ == "__main__":
    main()
