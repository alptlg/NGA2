#!/usr/bin/env python3
"""
Plot temperature and species mass fraction evolution from chem_reactor_0D output.

Output format: time, species-Y1..Yn, T, P, rho.

Requires: numpy, matplotlib
Usage:
  python plot_results.py [results.out] [-o results]
"""
import argparse
import glob
import re
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


def _fix_scientific_notation(line):
    """Fix malformed scientific notation: 7.6781603489-114 -> 7.6781603489E-114."""
    return re.sub(r"(\d)(-\d+)(?=\s|$)", r"\1E\2", line)


def parse_results(path):
    """Parse chem_reactor_0D output (fixed-width columns)."""
    with open(path) as f:
        header = f.readline().split()
    # Header: time, species-Y1, ..., T [, P, rho]
    has_tpr = len(header) >= 5 and header[-2] == "P" and header[-1] == "rho"
    if has_tpr:
        species_cols = header[1:-3]
    else:
        species_cols = header[1:-1]  # skip time and T
    species_names = []
    for col in species_cols:
        m = re.match(r"(.+)-Y\d+", col)
        if m:
            species_names.append(m.group(1))
    n_species = len(species_names)
    # Preprocess to fix malformed scientific notation (missing E before exponent)
    with open(path) as f:
        next(f)  # skip header
        lines = [_fix_scientific_notation(line) for line in f]
    data = np.loadtxt(lines, dtype=float)
    time = np.asarray(data[:, 0], dtype=float)
    Y = data[:, 1 : 1 + n_species]
    T = data[:, -3] if has_tpr else data[:, -1]
    P = data[:, -2] if has_tpr else np.full_like(time, np.nan)
    rho = data[:, -1] if has_tpr else np.full_like(time, np.nan)
    return {"time": time, "Y": Y, "species": species_names, "T": T, "P": P, "rho": rho}


def plot_results_file(input_file, output_basename="results"):
    """Plot chem_reactor_0D results from file. Returns 0 on success, 1 on error."""
    if not os.path.isfile(input_file):
        print(f"Error: file not found: {input_file}")
        return 1

    data = parse_results(input_file)
    print(f"Parsed {len(data['species'])} species from {input_file}")

    # Remove old output files before creating new ones
    for old in glob.glob(output_basename + "*.png"):
        os.remove(old)
        print(f"Removed {old}")

    # Build list of items to plot: (name, vals)
    plot_items = []
    for i, name in enumerate(data["species"]):
        plot_items.append((name, data["Y"][:, i]))
    plot_items.append(("T", data["T"]))
    if not np.all(np.isnan(data["P"])):
        plot_items.append(("P", data["P"]))
    if not np.all(np.isnan(data["rho"])):
        plot_items.append(("rho", data["rho"]))

    n_items = len(plot_items)
    n_rows, n_cols = 8, 6
    per_page = n_rows * n_cols

    for page_start in range(0, n_items, per_page):
        page_items = plot_items[page_start : page_start + per_page]
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.2 * n_cols, 2.6 * n_rows))
        axes = axes.flatten()

        for k, (name, vals) in enumerate(page_items):
            ax = axes[k]
            ax.plot(data["time"], vals, "b-", lw=1.5)
            ax.set_xlim(data["time"][0], data["time"][-1])
            ax.set_xlabel("time (s)")
            ax.set_ylabel("Y" if name not in ("T", "P", "rho") else name)
            ax.set_title(name)
            ax.grid(True, alpha=0.3)
            sf = ScalarFormatter(useMathText=False)
            sf.set_useOffset(False)
            ax.xaxis.set_major_formatter(sf)
            ax.yaxis.set_major_formatter(sf)
            ax.ticklabel_format(style="scientific", axis="x", scilimits=(0, 0), useOffset=False)

        for k in range(len(page_items), len(axes)):
            axes[k].set_visible(False)

        plt.tight_layout()
        suffix = f"_{page_start // per_page}" if n_items > per_page else ""
        outpath = output_basename + suffix + ".png"
        plt.savefig(outpath, dpi=150, bbox_inches="tight")
        print(f"Saved {outpath}")
        plt.close()
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Plot chem_reactor_0D temperature and mass fractions"
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default="results.out",
        help="chem_reactor_0D output file (default: results.out)",
    )
    parser.add_argument(
        "-o", "--output",
        default="results",
        help="Output plot basename (default: results)",
    )
    args = parser.parse_args()
    return plot_results_file(args.input_file, args.output)


if __name__ == "__main__":
    exit(main())
