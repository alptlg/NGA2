#!/usr/bin/env python3
"""
Generate chem_data_fc.f90 from a YAML chemical mechanism.

Usage: python yaml2nga.py <mechanism.yaml> [output.f90]
       Default output: chem_data_fc.f90 (if not specified)

Efficiency features:
- Vectorized rate coefficient evaluation (group by type)
- Pre-fitted backward Arrhenius (no equilibrium constants at runtime)
"""

import argparse
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    print("Error: PyYAML required. Install with: pip install pyyaml", file=sys.stderr)
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("Error: NumPy required. Install with: pip install numpy", file=sys.stderr)
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# Constants
R_UNIV = 8.314462618  # J/(mol·K)
CAL_TO_J = 4.184  # cal to J
ATM_TO_PA = 101325.0


def parse_pressure(p_val):
    """Parse pressure value (float in Pa or string like '0.1 atm') to Pa."""
    if p_val is None:
        return ATM_TO_PA
    if isinstance(p_val, (int, float)):
        return float(p_val)
    s = str(p_val).strip().lower()
    m = re.match(r"^\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*(\w+)?\s*$", s)
    if m:
        num = float(m.group(1))
        unit = (m.group(2) or "").strip()
        if unit in ("atm", "atmosphere"):
            return num * ATM_TO_PA
        if unit in ("bar",):
            return num * 1e5
        if unit in ("pa", "pascal"):
            return num
        if not unit:
            return num
    try:
        return float(p_val)
    except (ValueError, TypeError):
        return ATM_TO_PA
# cm³ -> m³: 1 m³ = 1e6 cm³; k in (conc)^(1-n)/s => A_SI = A_cgs * 1e-6^(n-1)
CM_TO_M3_FACTOR = 1e-6  # per power of concentration

# Fortran real format (aligned with ARCANE tools.format_to_significant)
def format_to_significant(x, n=15):
    """Format number for Fortran, coherent with double precision."""
    return "0.0" if x == 0 else f"{x:e}"


def f90_real(x):
    """Format scalar for Fortran real(WP) literal."""
    return f"{format_to_significant(x).replace('E', 'e')}_WP"


def f90_array_parameter(name, values, max_per_line=4):
    """Emit Fortran parameter array with continuation lines. Comma before & for continuation."""
    lines = []
    lines.append("  real(WP), parameter, dimension(" + str(len(values)) + ") :: " + name + " = (/ &")
    for i in range(0, len(values), max_per_line):
        chunk = values[i : i + max_per_line]
        vals = ", ".join(f90_real(v) for v in chunk)
        suffix = ", &" if i + max_per_line < len(values) else " /)"
        lines.append("       " + vals + suffix)
    return lines


def f90_int_array_parameter(name, values, max_per_line=8):
    """Emit Fortran integer parameter array. Comma before & for continuation."""
    lines = []
    lines.append("  integer, parameter, dimension(" + str(len(values)) + ") :: " + name + " = (/ &")
    for i in range(0, len(values), max_per_line):
        chunk = values[i : i + max_per_line]
        vals = ", ".join(str(int(v)) for v in chunk)
        suffix = ", &" if i + max_per_line < len(values) else " /)"
        lines.append("       " + vals + suffix)
    return lines


def f90_int_2d_array_parameter(name, arr, nrows, ncols, max_per_line=8):
    """Emit Fortran integer parameter matrix (column-major)."""
    flat = arr.flatten(order="F")
    if len(flat) <= _MAX_CONTINUATION_VALS:
        lines = []
        lines.append("  integer, parameter, dimension(" + str(nrows) + "," + str(ncols) + ") :: " + name + " = reshape((/ &")
        for i in range(0, len(flat), max_per_line):
            chunk = flat[i : i + max_per_line]
            vals = ", ".join(str(int(v)) for v in chunk)
            if i + max_per_line < len(flat):
                lines.append("       " + vals + ", &")
            else:
                lines.append("       " + vals + " /), (/ " + str(nrows) + ", " + str(ncols) + " /))")
        return lines
    # Fallback for very large arrays: chunk
    lines = []
    chunk_size = _MAX_CONTINUATION_VALS
    n_chunks = (len(flat) + chunk_size - 1) // chunk_size
    chunk_names = []
    for c in range(n_chunks):
        start = c * chunk_size
        end = min(start + chunk_size, len(flat))
        chunk_flat = flat[start:end]
        chunk_name = name + "_c" + str(c + 1)
        chunk_names.append(chunk_name)
        lines.append("  integer, parameter, dimension(" + str(len(chunk_flat)) + ") :: " + chunk_name + " = (/ &")
        for i in range(0, len(chunk_flat), max_per_line):
            sub = chunk_flat[i : i + max_per_line]
            vals = ", ".join(str(int(v)) for v in sub)
            suffix = ", &" if i + max_per_line < len(chunk_flat) else " /)"
            lines.append("       " + vals + suffix)
        lines.append("")
    lines.append("  integer, parameter, dimension(" + str(nrows) + "," + str(ncols) + ") :: " + name + " = reshape((/ " + ", ".join(chunk_names) + " /), (/ " + str(nrows) + ", " + str(ncols) + " /))")
    return lines


# Fortran limits continuation lines (F90: 39); use 38*4=152 to stay under
_MAX_CONTINUATION_VALS = 152


def _emit_array_constructor(values, max_per_line=4, indent="       ", bracket_style=True):
    """Emit lines for array constructor. bracket_style: [ ] (F2003) else (/ /) (F90)."""
    lines = []
    close_b = " ]" if bracket_style else " /)"
    for i in range(0, len(values), max_per_line):
        chunk = values[i : i + max_per_line]
        vals = ", ".join(f90_real(v) for v in chunk)
        suffix = ", &" if i + max_per_line < len(values) else " " + close_b
        lines.append(indent + vals + suffix)
    return lines


def f90_2d_array_parameter(name, arr, nrows, ncols, max_per_line=4):
    """Emit Fortran parameter matrix (column-major). Chunks if >_MAX_CONTINUATION_VALS to avoid continuation limit."""
    flat = arr.flatten(order="F")  # Fortran column-major
    if len(flat) <= _MAX_CONTINUATION_VALS:
        lines = []
        lines.append("  real(WP), parameter, dimension(" + str(nrows) + "," + str(ncols) + ") :: " + name + " = reshape([ &")
        for i in range(0, len(flat), max_per_line):
            chunk = flat[i : i + max_per_line]
            vals = ", ".join(f90_real(v) for v in chunk)
            if i + max_per_line < len(flat):
                lines.append("       " + vals + ", &")
            else:
                lines.append("       " + vals + " ], shape=[" + str(nrows) + ", " + str(ncols) + "])")
        return lines
    # Chunk to stay under continuation limit; use [ ] (F2003) for parser compatibility
    lines = []
    chunk_size = _MAX_CONTINUATION_VALS
    n_chunks = (len(flat) + chunk_size - 1) // chunk_size
    for c in range(n_chunks):
        start = c * chunk_size
        end = min(start + chunk_size, len(flat))
        chunk_flat = flat[start:end]
        chunk_name = name + "_c" + str(c + 1)
        lines.append("  real(WP), parameter, dimension(" + str(len(chunk_flat)) + ") :: " + chunk_name + " = [ &")
        lines.extend(_emit_array_constructor(chunk_flat, max_per_line, bracket_style=True))
        lines.append("")
    # Concatenate chunks and reshape; use [ ] for concatenation
    concat = ", ".join(name + "_c" + str(c + 1) for c in range(n_chunks))
    lines.append("  real(WP), parameter, dimension(" + str(nrows) + "," + str(ncols) + ") :: " + name + " = reshape([ " + concat + " ], shape=[" + str(nrows) + ", " + str(ncols) + "])")
    return lines


def cgs_to_si_A(A_cgs, n_reactants, length_unit):
    """Convert pre-exponential A from YAML (cm, mol) to SI (m, mol) for fcmech."""
    if length_unit != "cm":
        return A_cgs
    # k in (conc)^(1-n)/s; conc: mol/cm³ -> mol/m³ gives factor 1e-6 per conc power
    return A_cgs * (CM_TO_M3_FACTOR ** (n_reactants - 1))


def _falloff_n_fwd(left):
    """
    Falloff forward molecularity (n0, ninf) for cgs_to_si.
    n_left = sum of reactant stoichiometric coefficients.
    Low-P: M participates -> order = n_left + 1.
    High-P: no M -> order = n_left.
    """
    n_left = sum(c for _, c in left)
    return (n_left + 1, n_left)


def _species_f90_name(sp_name, species_names):
    """Get Fortran sXXX parameter name for species."""
    nm_str = str(sp_name)
    return "s" + nm_str.replace("(", "").replace(")", "").replace("+", "").replace("-", "").replace(" ", "")[:24]


def _resolve_species_name(sp_name, species_names):
    """
    Resolve equation species name to mechanism species name.
    Equation names may have suffixes (e.g. 1C7H15O2-C7H15O2 -> 1C7H15O2).
    Handles case differences (e.g. p-C4H9 -> P-C4H9, n-C3H7 -> N-C3H7).
    Returns matched species name or None.
    """
    if sp_name in species_names:
        return sp_name
    # Try longest match: sp_name starts with species + "-"
    for s in sorted(species_names, key=len, reverse=True):
        if sp_name.startswith(s + "-"):
            return s
    # Case-insensitive fallback (e.g. p-C4H9 -> P-C4H9, n-C3H7 -> N-C3H7)
    sp_upper = sp_name.upper()
    for s in species_names:
        if sp_upper == s.upper():
            return s
    for s in sorted(species_names, key=len, reverse=True):
        if sp_upper.startswith((s + "-").upper()):
            return s
    return None


def _build_reactant_expr(reactants, species_names):
    """Build Fortran expression: product of c(sXXX)**coeff for each reactant."""
    terms = []
    for sp_name, coeff in reactants:
        resolved = _resolve_species_name(sp_name, species_names)
        # Leading digit may have been parsed as coeff (e.g. "1C7H15O2-C7H15O2" -> coeff=1, sp="C7H15O2-C7H15O2")
        if resolved is None and 1 <= coeff <= 9:
            resolved = _resolve_species_name(str(coeff) + sp_name, species_names)
            if resolved is not None:
                coeff = 1
        if resolved is None:
            continue
        sn = _species_f90_name(resolved, species_names)
        if coeff == 1:
            terms.append("c(" + sn + ")")
        else:
            terms.append("c(" + sn + ")**" + str(float(coeff)) + "_WP")
    return " * ".join(terms) if terms else "1.0_WP"


def _format_reaction_side(side):
    """Format reactants/products: no coeff when 1; space between coeff and name when coeff != 1."""
    parts = []
    for sp, c in side:
        if c == 1:
            parts.append(sp)
        else:
            parts.append(f"{c} {sp}")
    return " + ".join(parts) if parts else ""


def parse_equation(eq):
    """Parse reaction equation 'aA + bB <=> cC + dD' or 'aA + bB => cC'."""
    eq = eq.split("#")[0].strip()
    rev = "<=>" in eq
    if rev:
        left, right = eq.split("<=>")
    else:
        left, right = eq.split("=>") if "=>" in eq else (eq, "")
    left = left.strip()
    right = right.strip()

    def parse_side(s):
        if not s:
            return []
        # Remove (+ M) / (+M) modifiers before splitting (e.g. "CO2 (+ M)" -> "CO2")
        s = re.sub(r"\s*\(\+?\s*M\s*\)\s*", " ", s)
        terms = re.split(r"\s+\+\s+", s)
        species = []
        for t in terms:
            t = t.strip()
            if not t:
                continue
            m = re.match(r"^(\d*)\s*([A-Za-z0-9_()\-]+)\s*$", t)
            if m:
                coeff_str = m.group(1)
                sp = m.group(2)
                # Leading digit + hyphen: digit is part of species name (e.g. "3-CH2" = triplet methylene)
                if sp.startswith("-") and coeff_str:
                    sp = coeff_str + sp
                    coeff = 1
                else:
                    coeff = int(coeff_str) if coeff_str else 1
                if sp.upper() != "M" and sp != "(+M)" and sp != "(+ M)":
                    species.append((sp, coeff))
        return species

    return parse_side(left), parse_side(right), rev


def load_mechanism(path):
    """Load and normalize YAML mechanism."""
    with open(path) as f:
        data = yaml.safe_load(f)

    # Get units
    units = data.get("units", {})
    act_energy_unit = units.get("activation-energy", "cal/mol")
    if "cal" in act_energy_unit.lower():
        Ea_scale = CAL_TO_J
    else:
        Ea_scale = 1.0

    # Get species from top-level or first phase
    all_species = data.get("species", [])
    name_to_sp = {}
    for s in all_species:
        if isinstance(s, dict):
            n = s.get("name", "X")
            name_to_sp[str(n)] = s
            if n is True:
                name_to_sp["yes"] = s
            elif n is False:
                name_to_sp["no"] = s

    species_list = []
    species_names = []
    # YAML parses NO/yes as bool; heuristics: False->NO, True->yes (common in chem)
    def canon_name(n):
        return "NO" if n is False else ("yes" if n is True else str(n))

    if "phases" in data:
        for phase in data["phases"]:
            if "species" in phase:
                sp_names = phase["species"]
                if isinstance(sp_names, list):
                    for n in sp_names:
                        key = "no" if n is False else ("yes" if n is True else str(n))
                        if key in name_to_sp or (n is False and "no" in name_to_sp) or (n is True and "yes" in name_to_sp):
                            k = key if key in name_to_sp else ("no" if n is False else "yes")
                            species_list.append(name_to_sp[k])
                            species_names.append(canon_name(n))
                break
    if not species_list and "species" in data:
        species_list = data["species"]
        species_names = [canon_name(sp.get("name", "X")) for sp in species_list if isinstance(sp, dict)]

    if not species_list:
        raise ValueError("No species found in mechanism")

    name_to_idx = {nm: i + 1 for i, nm in enumerate(species_names)}

    # Atoms: from phases elements or from union of species compositions
    atom_masses_kg = {
        "H": 0.00100794, "C": 0.0120107, "O": 0.0159994, "N": 0.0140067,
        "AR": 0.039948, "Ar": 0.039948,
        "HE": 0.004002602, "He": 0.004002602,
    }
    atom_names = []
    if "phases" in data:
        for phase in data["phases"]:
            if "elements" in phase:
                atom_names = [str(e) if e not in (True, False) else ("NO" if e is False else "yes") for e in phase["elements"]]
                break
    if not atom_names:
        # Collect from species compositions
        all_atoms = set()
        for sp in species_list:
            if isinstance(sp, dict):
                comp = sp.get("composition", {})
                for elem in comp:
                    all_atoms.add(str(elem).upper() if str(elem) != "Ar" else "Ar")
        # Order: O, H, C, N, Ar, He (common convention)
        for a in ("O", "H", "C", "N", "Ar", "He"):
            if a in all_atoms or (a == "Ar" and "AR" in all_atoms) or (a == "He" and "HE" in all_atoms):
                atom_names.append(a if a != "Ar" else "Ar")
        for a in sorted(all_atoms):
            if a not in ("O", "H", "C", "N", "Ar", "AR", "He", "HE"):
                atom_names.append(a)
    nA = len(atom_names)
    nS = len(species_list)

    # Composition matrix: comp_matrix[a][s] = count of atom a in species s
    comp_matrix = np.zeros((nA, nS), dtype=int)
    for j, sp in enumerate(species_list):
        if isinstance(sp, dict):
            comp = sp.get("composition", {})
            for elem, count in comp.items():
                ekey = str(elem).upper()
                for i, an in enumerate(atom_names):
                    if an.upper() == ekey:
                        comp_matrix[i, j] = int(count)
                        break

    atom_masses = np.array([
        atom_masses_kg.get(a.upper(), atom_masses_kg.get(a, 0.03))
        for a in atom_names
    ], dtype=float)

    # Reactions
    reactions = data.get("reactions", [])

    return {
        "species": species_list,
        "species_names": species_names,
        "name_to_idx": name_to_idx,
        "reactions": reactions,
        "Ea_scale": Ea_scale,
        "units": units,
        "atom_names": atom_names,
        "atom_masses": atom_masses,
        "comp_matrix": comp_matrix,
        "nA": nA,
    }


def extract_yaml_labels(yaml_path):
    """Extract YAML reaction labels (e.g. # Reaction 1) from raw file. YAML strips # as comments."""
    labels = {}
    with open(yaml_path) as f:
        rid = 0
        for line in f:
            if "equation:" in line:
                rid += 1
                if "#" in line:
                    parts = line.split("#", 1)
                    labels[rid] = "# " + parts[1].strip()
    return labels


def get_thermo_coeffs(spec):
    """Extract NASA7 coefficients and T_mid from species."""
    thermo = spec.get("thermo", {})
    if isinstance(thermo, list):
        thermo = thermo[0]
    data = thermo.get("data", [])
    if not data:
        return None, None, 1000.0
    # First row = low T, second = high T (Cantera order)
    low = data[0]
    high = data[1] if len(data) > 1 else data[0]
    # T_mid from temperature-ranges: [T_low, T_mid, T_high]
    tr = thermo.get("temperature-ranges", [200.0, 1000.0, 5000.0])
    T_mid = float(tr[1]) if len(tr) > 1 else 1000.0
    return low, high, T_mid


def get_transport(spec):
    """Extract transport data."""
    trans = spec.get("transport", {})
    if isinstance(trans, list):
        trans = trans[0]
    return {
        "diameter": trans.get("diameter", 3.0),
        "well_depth": trans.get("well-depth", 100.0),
        "geometry": trans.get("geometry", "nonlinear"),
    }


def compute_molar_mass(spec):
    """Compute molar mass from composition (kg/mol)."""
    comp = spec.get("composition", {})
    masses = {"H": 0.00100794, "C": 0.0120107, "O": 0.0159994, "N": 0.0140067, "AR": 0.039948}
    m = 0.0
    for elem, count in comp.items():
        m += masses.get(elem.upper(), 0.03) * count  # fallback for unknown
    return m if m > 0 else 0.03


def collision_integral_omega_mu(T_star):
    """Approximate omega_mu for viscosity (Lennard-Jones)."""
    if T_star < 0.3:
        return 100.0
    if T_star > 100:
        return 0.5
    # Simplified fit
    a, b, c, d = 1.16145, -0.14874, 0.52487, -0.77320
    e, f, g, h = 2.16178, -2.43787, -6.435e-4, 0.14874
    ln_T = np.log(T_star)
    n = 2.0 - a * np.exp(b * ln_T) + c * np.exp(d * ln_T) + e * np.exp(f * ln_T) + g * np.exp(h * ln_T)
    return n


def collision_integral_omega_D(T_star):
    """Approximate omega_D for diffusion."""
    if T_star < 0.3:
        return 100.0
    if T_star > 400:
        return 0.5
    a, b, c, d = 1.06036, 0.15610, 0.19300, 0.47635
    e, f, g, h = 1.03587, 1.52996, 1.76474, 3.89411
    ln_T = np.log(T_star)
    n = 1.0 / (a * T_star**b) + c / np.exp(d * T_star) + e / np.exp(f * T_star) + g / np.exp(h * T_star)
    return n


def nasa7_hort_SR(T, a, T_switch):
    """
    NASA7 H and S at temperature T (evaluation at T, not at T_switch).
    T_switch: threshold to select low-T (T < T_switch) vs high-T (T >= T_switch) coefficient set.
    h/(RT) = a0 + a1*T/2 + a2*T^2/3 + a3*T^3/4 + a4*T^4/5 + a5/T
    S/R = a0*ln(T) + a1*T + a2*T^2/2 + a3*T^3/3 + a4*T^4/4 + a6
    """
    if T < T_switch:
        c = a[0:7]
    else:
        c = a[7:14]
    T2, T3, T4 = T * T, T * T * T, T * T * T * T
    hort = c[0] + c[1] * T / 2 + c[2] * T2 / 3 + c[3] * T3 / 4 + c[4] * T4 / 5 + c[5] / T
    s_R = c[0] * np.log(T) + c[1] * T + c[2] * T2 / 2 + c[3] * T3 / 3 + c[4] * T4 / 4 + c[6]
    return hort, s_R


def compute_keq_from_thermo(T_arr, nu, thermo_coeffs, T_switch, length_unit="m"):
    """
    Compute K_c(T) from NASA7 thermo for use in k_b = k_f/K_c.
    H and S are evaluated at each T in T_arr (varying temperature).
    K_p = exp(-Delta_G^0/(R*T)), K_c = K_p * (P^0/(R*T))^(sum nu_i).
    nu[i] = stoichiometric coefficient (positive products, negative reactants).
    thermo_coeffs[i] = [a0..a6 low, a0..a6 high], T_switch[i] = threshold (K) for low/high range.
    length_unit: "m" -> K_c in (mol/m^3)^dnu (SI); "cm" -> K_c in (mol/cm^3)^dnu (matches YAML)
    """
    P0 = 101325.0  # Pa (1 atm)
    dnu = np.sum(nu)
    # P0/(R*T) in mol/m^3; for cm we need mol/cm^3 = 1e6 * mol/m^3
    conc_factor = 1e6 if length_unit == "cm" else 1.0
    K_eq = np.ones_like(T_arr)
    for it, T in enumerate(T_arr):
        dG_RT = 0.0
        for i, nui in enumerate(nu):
            if nui == 0:
                continue
            hort, s_R = nasa7_hort_SR(T, thermo_coeffs[i], T_switch[i])
            dG_RT += nui * (hort - s_R)  # Delta_G^0/(R*T) = sum(nu_i * (H_i/(RT) - S_i/R))
        K_p = np.exp(-dG_RT)
        K_eq[it] = K_p * (conc_factor * P0 / (R_UNIV * T)) ** dnu
    return K_eq


def fit_backward_arrhenius(kf_A, kf_b, kf_Ea, K_eq_arr=None, delta_S0=None, delta_H0=None, T_range=(300, 3000), n_pts=50, kf_A_si=None):
    """
    Fit k_b = k_f / K_eq to Arrhenius form over T_range.
    K_eq can be provided as:
    - K_eq_arr: array of length n_pts (thermo-derived, temperature-dependent)
    - delta_S0, delta_H0: constant (K_eq = exp(delta_S0/R - delta_H0/(R*T)))
    - If none: K_eq=1 (k_b=k_f).
    kf_A_si: if provided, k_f is already in SI; K_eq_arr must be in SI (mol/m^3)^dnu.
      Returns A_b in SI (no cgs_to_si needed). Use for fcmech consistency.
    """
    T_arr = np.linspace(T_range[0], T_range[1], n_pts)
    kf_A_use = kf_A_si if kf_A_si is not None else kf_A
    kf = kf_A_use * (T_arr**kf_b) * np.exp(-kf_Ea / (R_UNIV * T_arr))
    if K_eq_arr is not None:
        K_eq = np.maximum(K_eq_arr, 1e-300)
    elif delta_S0 is not None and delta_H0 is not None:
        K_eq = np.exp(delta_S0 / R_UNIV - delta_H0 / (R_UNIV * T_arr))
    else:
        K_eq = np.ones_like(T_arr)
    kb = kf / K_eq
    # Avoid log(0): use minimal floor only for numerical stability (no artificial clipping)
    kb = np.maximum(kb, 1e-300)

    # Fit log(kb) = log(A) + b*log(T) - Ea/(R*T)
    # X[:,2] = 1/(R*T), so coeffs[2] = -Ea => Ea = -coeffs[2] (coeff multiplies 1/(R*T))
    log_kb = np.log(kb)
    log_T = np.log(T_arr)
    inv_RT = 1.0 / (R_UNIV * T_arr)
    X = np.column_stack([np.ones_like(T_arr), log_T, inv_RT])
    coeffs, _, _, _ = np.linalg.lstsq(X, log_kb, rcond=None)
    A_fit = np.exp(coeffs[0])
    b_fit = coeffs[1]
    Ea_fit = -coeffs[2]  # coeffs[2] * inv_RT = -Ea/(R*T), so coeffs[2] = -Ea
    return A_fit, b_fit, Ea_fit


def yaml2nga(mech, output_path, yaml_path=None, collect_viz=False):
    """Generate chem_data_fc.f90 from mechanism data. If collect_viz=True, returns viz_data for --viz."""
    nS = len(mech["species_names"])
    reactions = mech["reactions"]
    name_to_idx = mech["name_to_idx"]
    species_names = mech["species_names"]
    Ea_scale = mech["Ea_scale"]

    # Build thermo_coeffs and T_mid early (needed for backward rate fitting)
    T_mid_sp = np.zeros(nS)
    thermo_coeffs = np.zeros((nS, 14))
    default_coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for i, sp in enumerate(mech["species"]):
        if isinstance(sp, dict):
            low, high, T_mid = get_thermo_coeffs(sp)
            T_mid_sp[i] = T_mid
            if low and high:
                thermo_coeffs[i, 0:7] = low[0:7]
                thermo_coeffs[i, 7:14] = high[0:7]
            else:
                thermo_coeffs[i, 0:7] = default_coeffs
                thermo_coeffs[i, 7:14] = default_coeffs
        else:
            T_mid_sp[i] = 1000.0
            thermo_coeffs[i, 0:7] = default_coeffs
            thermo_coeffs[i, 7:14] = default_coeffs

    # Classify reactions - forward and backward (pre-fitted Arrhenius for reversible)
    arrhenius_fwd = []
    arrhenius_bwd = []  # (rid, A, b, Ea) - pre-fitted from K_eq
    viz_data = [] if collect_viz else None
    three_body = []
    falloff = []
    plog = []  # (rid, left, right, rev, [(P, A, b, Ea), ...]) - P in Pa

    # Length unit from mechanism: "cm" -> K_c in (mol/cm^3)^dnu to match YAML rates
    length_unit = mech.get("units", {}).get("length", "m")
    if length_unit not in ("cm", "m"):
        length_unit = "m"

    # Parse reactions
    for i, r in enumerate(reactions):
        if isinstance(r, str):
            continue
        eq = r.get("equation", "")
        if not eq:
            continue
        if r.get("negative-A", False):
            continue  # Skip: negative-A not supported (e.g. Nakamura Reaction 68)
        left, right, rev = parse_equation(eq)
        rtype = r.get("type", "elementary")
        rc = r.get("rate-constant", r.get("rate-constant", {}))
        if isinstance(rc, dict):
            A = float(rc.get("A", 1.0))
            b = float(rc.get("b", 0.0))
            Ea = float(rc.get("Ea", 0.0)) * Ea_scale
        else:
            A, b, Ea = 1.0, 0.0, 0.0

        if rtype == "three-body":
            A_b, b_b, Ea_b = None, None, None
            if rev:
                nu = np.zeros(nS)
                for sp_name, coeff in left:
                    resolved = _resolve_species_name(sp_name, species_names)
                    if resolved is not None:
                        nu[name_to_idx[resolved] - 1] -= coeff
                for sp_name, coeff in right:
                    resolved = _resolve_species_name(sp_name, species_names)
                    if resolved is not None:
                        nu[name_to_idx[resolved] - 1] += coeff
                T_range = (300.0, 3000.0)
                n_pts = 50
                T_arr = np.linspace(T_range[0], T_range[1], n_pts)
                # Use SI (m) for K_eq and k_f so fit returns A_b in SI (fcmech uses mol/m^3)
                K_eq_arr = compute_keq_from_thermo(T_arr, nu, thermo_coeffs, T_mid_sp, length_unit="m")
                # Three-body: n = n_left + 1 (dissociation HNCO+M: n=2; association O+O+M: n=3)
                n_fwd = sum(c for _, c in left) + 1
                kf_A_si = cgs_to_si_A(A, n_fwd, length_unit)
                A_b, b_b, Ea_b = fit_backward_arrhenius(A, b, Ea, K_eq_arr=K_eq_arr, T_range=T_range, n_pts=n_pts, kf_A_si=kf_A_si)
                if collect_viz:
                    kf_si = kf_A_si * (T_arr**b) * np.exp(-Ea / (R_UNIV * T_arr))
                    k_b_keq = np.maximum(kf_si / np.maximum(K_eq_arr, 1e-300), 1e-300)
                    k_b_fit = A_b * (T_arr**b_b) * np.exp(-Ea_b / (R_UNIV * T_arr))
                    viz_data.append((eq.split("#")[0].strip() + " (3-body)", i + 1, T_arr, k_b_keq, k_b_fit))
            default_eff = float(r.get("default-efficiency", 1.0))
            three_body.append((i + 1, left, right, rev, A, b, Ea, r.get("efficiencies", {}), default_eff, A_b, b_b, Ea_b))
        elif rtype in ("pressure-dependent-Arrhenius", "plog"):
            # Plog: rate-constants = [{P, A, b, Ea}, ...], P in Pa
            rcs = r.get("rate-constants", r.get("rate-constants", []))
            if not rcs:
                continue
            pts = []
            for rd in rcs:
                if isinstance(rd, dict):
                    P = parse_pressure(rd.get("P", rd.get("p", 101325.0)))
                    Ai = float(rd.get("A", 1.0))
                    bi = float(rd.get("b", 0.0))
                    Eai = float(rd.get("Ea", 0.0)) * Ea_scale
                    pts.append((P, Ai, bi, Eai))
            pts.sort(key=lambda x: x[0])  # sort by P ascending
            if pts:
                # For reversible Plog: fit backward Arrhenius at each pressure point (k_b = k_f/K_eq)
                pts_b = None
                if rev:
                    nu = np.zeros(nS)
                    for sp_name, coeff in left:
                        resolved = _resolve_species_name(sp_name, species_names)
                        if resolved is not None:
                            nu[name_to_idx[resolved] - 1] -= coeff
                    for sp_name, coeff in right:
                        resolved = _resolve_species_name(sp_name, species_names)
                        if resolved is not None:
                            nu[name_to_idx[resolved] - 1] += coeff
                    T_range = (300.0, 3000.0)
                    n_pts = 50
                    T_arr = np.linspace(T_range[0], T_range[1], n_pts)
                    K_eq_arr = compute_keq_from_thermo(T_arr, nu, thermo_coeffs, T_mid_sp, length_unit="m")
                    n_fwd = sum(c for _, c in left)
                    pts_b = []
                    for P, Ai, bi, Eai in pts:
                        kf_A_si = cgs_to_si_A(Ai, n_fwd, length_unit)
                        Ab, bb, Eab = fit_backward_arrhenius(Ai, bi, Eai, K_eq_arr=K_eq_arr, T_range=T_range, n_pts=n_pts, kf_A_si=kf_A_si)
                        pts_b.append((Ab, bb, Eab))
                    if collect_viz and pts_b:
                        P0, Ai, bi, Eai = pts[0]
                        Ab, bb, Eab = pts_b[0]
                        kf_A_si = cgs_to_si_A(Ai, n_fwd, length_unit)
                        T_arr = np.linspace(T_range[0], T_range[1], n_pts)
                        kf_si = kf_A_si * (T_arr**bi) * np.exp(-Eai / (R_UNIV * T_arr))
                        k_b_keq = np.maximum(kf_si / np.maximum(K_eq_arr, 1e-300), 1e-300)
                        k_b_fit = Ab * (T_arr**bb) * np.exp(-Eab / (R_UNIV * T_arr))
                        viz_data.append((eq.split("#")[0].strip() + f" (plog P={P0:.0f}Pa)", i + 1, T_arr, k_b_keq, k_b_fit))
                plog.append((i + 1, left, right, rev, pts, pts_b))
        elif rtype == "falloff":
            low = r.get("low-P-rate-constant", rc)
            high = r.get("high-P-rate-constant", rc)
            if isinstance(low, dict):
                A0 = float(low.get("A", 1.0))
                b0 = float(low.get("b", 0.0))
                Ea0 = float(low.get("Ea", 0.0)) * Ea_scale
            else:
                A0, b0, Ea0 = 1.0, 0.0, 0.0
            if isinstance(high, dict):
                Ainf = float(high.get("A", 1.0))
                binf = float(high.get("b", 0.0))
                Eainf = float(high.get("Ea", 0.0)) * Ea_scale
            else:
                Ainf, binf, Eainf = 1.0, 0.0, 0.0
            troe = r.get("Troe", r.get("troe", []))
            # For reversible falloff: fit backward Arrhenius for low-P and high-P limits.
            # Forward molecularity: n0=n_left+1, ninf=n_left (by reactant stoichiometry).
            # k_b = k_f/K_eq; K_eq in SI (mol/m^3)^dnu; fit returns A_b in SI for correct units.
            # Troe coefficients and efficiencies are shared (same fc_fo, M) for forward and reverse.
            A0_b, b0_b, Ea0_b = None, None, None
            Ainf_b, binf_b, Eainf_b = None, None, None
            if rev:
                nu = np.zeros(nS)
                for sp_name, coeff in left:
                    resolved = _resolve_species_name(sp_name, species_names)
                    if resolved is not None:
                        nu[name_to_idx[resolved] - 1] -= coeff
                for sp_name, coeff in right:
                    resolved = _resolve_species_name(sp_name, species_names)
                    if resolved is not None:
                        nu[name_to_idx[resolved] - 1] += coeff
                T_range = (300.0, 3000.0)
                n_pts = 50
                T_arr = np.linspace(T_range[0], T_range[1], n_pts)
                K_eq_arr = compute_keq_from_thermo(T_arr, nu, thermo_coeffs, T_mid_sp, length_unit="m")
                n_fwd_low, n_fwd_high = _falloff_n_fwd(left)
                kf0_si = cgs_to_si_A(A0, n_fwd_low, length_unit)
                kfinf_si = cgs_to_si_A(Ainf, n_fwd_high, length_unit)
                A0_b, b0_b, Ea0_b = fit_backward_arrhenius(A0, b0, Ea0, K_eq_arr=K_eq_arr, T_range=T_range, n_pts=n_pts, kf_A_si=kf0_si)
                Ainf_b, binf_b, Eainf_b = fit_backward_arrhenius(Ainf, binf, Eainf, K_eq_arr=K_eq_arr, T_range=T_range, n_pts=n_pts, kf_A_si=kfinf_si)
                if collect_viz:
                    kf0_si_arr = kf0_si * (T_arr**b0) * np.exp(-Ea0 / (R_UNIV * T_arr))
                    k_b0_keq = np.maximum(kf0_si_arr / np.maximum(K_eq_arr, 1e-300), 1e-300)
                    k_b0_fit = A0_b * (T_arr**b0_b) * np.exp(-Ea0_b / (R_UNIV * T_arr))
                    viz_data.append((eq.split("#")[0].strip() + " (falloff low-P)", i + 1, T_arr, k_b0_keq, k_b0_fit))
                    kfinf_si_arr = kfinf_si * (T_arr**binf) * np.exp(-Eainf / (R_UNIV * T_arr))
                    k_binf_keq = np.maximum(kfinf_si_arr / np.maximum(K_eq_arr, 1e-300), 1e-300)
                    k_binf_fit = Ainf_b * (T_arr**binf_b) * np.exp(-Eainf_b / (R_UNIV * T_arr))
                    viz_data.append((eq.split("#")[0].strip() + " (falloff high-P)", i + 1, T_arr, k_binf_keq, k_binf_fit))
            default_eff = float(r.get("default-efficiency", 1.0))
            falloff.append((i + 1, left, right, rev, A0, b0, Ea0, Ainf, binf, Eainf, troe, r.get("efficiencies", {}),
                           default_eff, A0_b, b0_b, Ea0_b, Ainf_b, binf_b, Eainf_b))
        else:
            arrhenius_fwd.append((i + 1, left, right, A, b, Ea))
            if rev:
                # Pre-fit backward Arrhenius: k_b = k_f/K_eq with thermo-derived K_eq (SI for fcmech)
                nu = np.zeros(nS)
                for sp_name, coeff in left:
                    resolved = _resolve_species_name(sp_name, species_names)
                    if resolved is not None:
                        nu[name_to_idx[resolved] - 1] -= coeff
                for sp_name, coeff in right:
                    resolved = _resolve_species_name(sp_name, species_names)
                    if resolved is not None:
                        nu[name_to_idx[resolved] - 1] += coeff
                T_range = (300.0, 3000.0)
                n_pts = 50
                T_arr = np.linspace(T_range[0], T_range[1], n_pts)
                K_eq_arr = compute_keq_from_thermo(T_arr, nu, thermo_coeffs, T_mid_sp, length_unit="m")
                n_fwd = sum(c for _, c in left)
                kf_A_si = cgs_to_si_A(A, n_fwd, length_unit)
                Ab, bb, Eab = fit_backward_arrhenius(A, b, Ea, K_eq_arr=K_eq_arr, T_range=T_range, n_pts=n_pts, kf_A_si=kf_A_si)
                n_bwd = sum(c for _, c in right)  # backward reactants = forward products
                arrhenius_bwd.append((i + 1, Ab, bb, Eab, n_bwd))
                if collect_viz:
                    kf_si = kf_A_si * (T_arr**b) * np.exp(-Ea / (R_UNIV * T_arr))
                    k_b_keq = np.maximum(kf_si / np.maximum(K_eq_arr, 1e-300), 1e-300)
                    k_b_fit = Ab * (T_arr**bb) * np.exp(-Eab / (R_UNIV * T_arr))
                    viz_data.append((eq.split("#")[0].strip(), i + 1, T_arr, k_b_keq, k_b_fit))

    n_plog = len(plog)
    n_tb_rev = sum(1 for r in three_body if r[3])  # reversible three-body count
    n_fo_rev = sum(1 for r in falloff if r[3])  # reversible falloff count
    n_plog_rev = sum(1 for r in plog if r[3])  # reversible plog count
    nR = len(arrhenius_fwd) + len(arrhenius_bwd) + len(three_body) + n_tb_rev + len(falloff) + n_fo_rev + n_plog + n_plog_rev
    n_arr = len(arrhenius_fwd) + len(arrhenius_bwd)
    n_tb = len(three_body)
    n_fo = len(falloff)

    # Build Arrhenius parameter arrays for Fortran (convert to SI when YAML uses cm)
    # Store ln(A) and E/R: ln(k) = ln(A) + b*ln(T) - E/(R*T), then k = exp(ln(k))
    A_fwd = [cgs_to_si_A(r[3], sum(c for _, c in r[1]), length_unit) for r in arrhenius_fwd]
    A_bwd = [r[1] for r in arrhenius_bwd]  # r[1] = A_b already in SI (fit used kf_A_si + K_eq in m)
    A_arr = np.array(A_fwd + A_bwd)
    b_arr = np.array([r[4] for r in arrhenius_fwd] + [r[2] for r in arrhenius_bwd])
    Ea_arr = np.array([r[5] for r in arrhenius_fwd] + [r[3] for r in arrhenius_bwd])
    ln_A_arr = np.log(np.maximum(A_arr, 1e-300))
    E_R_arr = Ea_arr / R_UNIV

    # Three-body: n = n_left + 1 (e.g. O+O+M: n_left=2 -> n=3)
    A_tb = np.array([cgs_to_si_A(r[4], sum(c for _, c in r[1]) + 1, length_unit) for r in three_body]) if n_tb else np.array([])
    b_tb = np.array([r[5] for r in three_body]) if n_tb else np.array([])
    Ea_tb = np.array([r[6] for r in three_body]) if n_tb else np.array([])
    ln_A_tb = np.log(np.maximum(A_tb, 1e-300)) if n_tb else np.array([])
    E_R_tb = Ea_tb / R_UNIV if n_tb else np.array([])
    # Three-body backward (reversible only): bimolecular (n=2) e.g. O2 + M => 2O + M
    # A_b from fit is already in SI (K_eq and k_f used SI)
    tb_rev_idx = []
    if n_tb_rev > 0:
        tb_rev_list = [(i, r) for i, r in enumerate(three_body) if r[3]]
        A_tb_b = np.array([r[9] for _, r in tb_rev_list])
        b_tb_b = np.array([r[10] for _, r in tb_rev_list])
        Ea_tb_b = np.array([r[11] for _, r in tb_rev_list])
        ln_A_tb_b = np.log(np.maximum(A_tb_b, 1e-300))
        E_R_tb_b = Ea_tb_b / R_UNIV
        tb_rev_idx = [-1] * n_tb
        for j, (i, _) in enumerate(tb_rev_list):
            tb_rev_idx[i] = j
    else:
        ln_A_tb_b = b_tb_b = E_R_tb_b = np.array([])
        tb_rev_idx = [-1] * n_tb if n_tb else []

    # Falloff: n0 = n_left+1, ninf = n_left (by reactant count)
    A0_fo = np.array([cgs_to_si_A(r[4], _falloff_n_fwd(r[1])[0], length_unit) for r in falloff]) if n_fo else np.array([])
    b0_fo = np.array([r[5] for r in falloff]) if n_fo else np.array([])
    Ea0_fo = np.array([r[6] for r in falloff]) if n_fo else np.array([])
    Ainf_fo = np.array([cgs_to_si_A(r[7], _falloff_n_fwd(r[1])[1], length_unit) for r in falloff]) if n_fo else np.array([])
    binf_fo = np.array([r[8] for r in falloff]) if n_fo else np.array([])
    Eainf_fo = np.array([r[9] for r in falloff]) if n_fo else np.array([])
    ln_A0_fo = np.log(np.maximum(A0_fo, 1e-300)) if n_fo else np.array([])
    E_R0_fo = Ea0_fo / R_UNIV if n_fo else np.array([])
    ln_Ainf_fo = np.log(np.maximum(Ainf_fo, 1e-300)) if n_fo else np.array([])
    E_Rinf_fo = Eainf_fo / R_UNIV if n_fo else np.array([])
    # Falloff backward: fit gives k_b=k_f/K_eq already in SI. A0_b,Ainf_b in correct units for getlindratecoeff.
    # Troe and efficiencies: same as forward (shared fc_fo, M(n_tb+i) in getlindratecoeff).
    fo_rev_idx = []  # fo_rev_idx[i] = backward index (0-based) for falloff i, -1 if irreversible
    if n_fo_rev > 0:
        fo_rev_list = [(i, r) for i, r in enumerate(falloff) if r[3]]
        A0_b_fo = np.array([r[13] for _, r in fo_rev_list])  # [CO2][M]
        b0_b_fo = np.array([r[14] for _, r in fo_rev_list])
        Ea0_b_fo = np.array([r[15] for _, r in fo_rev_list])
        Ainf_b_fo = np.array([r[16] for _, r in fo_rev_list])  # unimolecular 1/s
        binf_b_fo = np.array([r[17] for _, r in fo_rev_list])
        Eainf_b_fo = np.array([r[18] for _, r in fo_rev_list])
        ln_A0_b_fo = np.log(np.maximum(A0_b_fo, 1e-300))
        E_R0_b_fo = Ea0_b_fo / R_UNIV
        ln_Ainf_b_fo = np.log(np.maximum(Ainf_b_fo, 1e-300))
        E_Rinf_b_fo = Eainf_b_fo / R_UNIV
        fo_rev_idx = [-1] * n_fo
        for j, (i, _) in enumerate(fo_rev_list):
            fo_rev_idx[i] = j
    else:
        A0_b_fo = b0_b_fo = Ea0_b_fo = Ainf_b_fo = binf_b_fo = Eainf_b_fo = np.array([])
        ln_A0_b_fo = E_R0_b_fo = ln_Ainf_b_fo = E_Rinf_b_fo = np.array([])
        fo_rev_idx = [-1] * n_fo if n_fo else []
    # Troe: [a, T3, T1, T2?] or {A/a, T3, T1, T2}; F_cent = (1-a)*exp(-T/T3) + a*exp(-T/T1) + exp(-T2/T)
    # When T2 is 0 or omitted, Cantera omits the exp(-T2/T) term; use 1e30 so exp(-1e30/T)≈0
    def parse_troe(t):
        if isinstance(t, (list, tuple)) and len(t) >= 3:
            a, T3, T1 = float(t[0]), float(t[1]), float(t[2])
            T2 = float(t[3]) if len(t) > 3 and t[3] is not None else 1.0e30
            if T2 == 0.0:
                T2 = 1.0e30
            return a, T3, T1, T2
        if isinstance(t, dict):
            a = float(t.get("A", t.get("a", 1.0)))
            T3 = float(t.get("T3", 1.0e30))
            T1 = float(t.get("T1", 1.0e30))
            T2 = float(t.get("T2", 1.0e30))
            if T2 == 0.0:
                T2 = 1.0e30
            return a, T3, T1, T2
        return 1.0, 1.0e30, 1.0e30, 1.0e30  # Lindemann: FC=1
    troe_pars = [parse_troe(r[10]) for r in falloff] if n_fo else []
    Troe_a = np.array([p[0] for p in troe_pars]) if n_fo else np.array([])
    Troe_T3 = np.array([p[1] for p in troe_pars]) if n_fo else np.array([])
    Troe_T1 = np.array([p[2] for p in troe_pars]) if n_fo else np.array([])
    Troe_T2 = np.array([p[3] for p in troe_pars]) if n_fo else np.array([])

    # Plog: flattened arrays + start/length per reaction (Fortran 1-based indices)
    if n_plog > 0:
        P_plog = []
        ln_A_plog = []
        b_plog = []
        E_R_plog = []
        ln_A_b_plog = []
        b_b_plog = []
        E_R_b_plog = []
        plog_reac_start = []
        plog_reac_len = []
        plog_rev_idx = []
        idx = 1  # Fortran 1-based
        plog_rev_j = 0
        for j, reac in enumerate(plog):
            pts = reac[4]
            pts_b = reac[5]
            left_plog = reac[1]
            n_reactants = sum(c for _, c in left_plog)
            n_bwd = sum(c for _, c in reac[2])  # backward reactants = forward products
            plog_reac_start.append(idx)
            plog_reac_len.append(len(pts))
            plog_rev_idx.append(plog_rev_j + 1 if reac[3] else 0)
            if reac[3]:
                plog_rev_j += 1
            for ipt, (P, A, b, Ea) in enumerate(pts):
                P_plog.append(P)
                A_si = cgs_to_si_A(A, n_reactants, length_unit)
                ln_A_plog.append(np.log(np.maximum(A_si, 1e-300)))
                b_plog.append(b)
                E_R_plog.append(Ea / R_UNIV)
                if pts_b is not None:
                    Ab, bb, Eab = pts_b[ipt]
                    # Ab from fit is already in SI (K_eq and k_f used SI)
                    ln_A_b_plog.append(np.log(np.maximum(Ab, 1e-300)))
                    b_b_plog.append(bb)
                    E_R_b_plog.append(Eab / R_UNIV)
                else:
                    ln_A_b_plog.append(0.0)
                    b_b_plog.append(0.0)
                    E_R_b_plog.append(0.0)
            idx += len(pts)
        P_plog = np.array(P_plog)
        ln_A_plog = np.array(ln_A_plog)
        b_plog = np.array(b_plog)
        E_R_plog = np.array(E_R_plog)
        ln_A_b_plog = np.array(ln_A_b_plog)
        b_b_plog = np.array(b_b_plog)
        E_R_b_plog = np.array(E_R_b_plog)
        plog_reac_start = np.array(plog_reac_start)
        plog_reac_len = np.array(plog_reac_len)
        plog_rev_idx = np.array(plog_rev_idx)
        n_plog_pts = len(P_plog)
    else:
        n_plog_pts = 0
        ln_A_b_plog = b_b_plog = E_R_b_plog = np.array([])
        plog_rev_idx = np.array([])

    # Third-body efficiency: reduced matrix for species with non-default efficiencies only
    # M(r) = sum_j eff_matrix_reduced(r,j)*c(eff_species_idx(j)) + sum(c for species not in eff set)
    n_m = max(1, n_tb + n_fo)
    species_names = mech["species_names"]
    species_in_efficiencies = set()
    for r in range(n_tb):
        for sp_name in three_body[r][7]:  # eff_dict
            species_in_efficiencies.add(sp_name)
    for r in range(n_fo):
        for sp_name in falloff[r][11]:  # eff_dict
            species_in_efficiencies.add(sp_name)
    # Only species that exist in mechanism and appear in at least one efficiency dict
    eff_species_names = sorted([sp for sp in species_in_efficiencies if sp in species_names])
    n_eff = len(eff_species_names)
    eff_species_idx = np.array([species_names.index(sp) + 1 for sp in eff_species_names])  # 1-based
    eff_matrix_reduced = np.zeros((n_m, n_eff))
    for r in range(n_tb):
        eff_dict = three_body[r][7]
        default_eff = three_body[r][8]
        for j, sp in enumerate(eff_species_names):
            eff_matrix_reduced[r, j] = float(eff_dict.get(sp, default_eff))
    for r in range(n_fo):
        eff_dict = falloff[r][11]
        default_eff = falloff[r][12]
        for j, sp in enumerate(eff_species_names):
            eff_matrix_reduced[n_tb + r, j] = float(eff_dict.get(sp, default_eff))

    # Per-reaction default efficiency for species not in eff_species_idx (Cantera default-efficiency)
    default_eff_arr = np.array(
        [three_body[r][8] for r in range(n_tb)] + [falloff[r][12] for r in range(n_fo)]
    ) if n_m > 0 else np.array([1.0])

    # Molar masses (kg/mol) - needed for mass conservation check
    W_sp = []
    for sp in mech["species"]:
        if isinstance(sp, dict):
            W_sp.append(compute_molar_mass(sp))
        else:
            W_sp.append(0.03)
    W_sp = np.array(W_sp)

    # Index boundaries for k array (1-based)
    i_tb_end = n_arr + n_tb
    i_tb_bwd_end = i_tb_end + n_tb_rev
    i_fo_end = i_tb_bwd_end + n_fo
    i_fo_bwd_end = i_fo_end + n_fo_rev
    i_plog_start = i_fo_bwd_end + 1
    i_plog_end = i_plog_start + n_plog - 1 if n_plog > 0 else i_plog_start - 1
    i_plog_bwd_start = i_plog_end + 1

    # Build reaction rate expressions w(i) = k(i) * product(reactants) [* M(r) for 3-body/falloff]
    # k layout: Arrhenius(1..n_arr) | 3-body | falloff_fwd | falloff_bwd | plog_fwd | plog_bwd
    w_exprs = [""] * nR
    idx = 0
    # Arrhenius forward
    for r in arrhenius_fwd:
        left = r[1]
        expr = _build_reactant_expr(left, species_names)
        w_exprs[idx] = expr
        idx += 1
    # Arrhenius backward
    for r in arrhenius_bwd:
        rid = r[0]
        right = None
        for af in arrhenius_fwd:
            if af[0] == rid:
                right = af[2]  # right = products = backward reactants
                break
        if right is None:
            right = []
        expr = _build_reactant_expr(right, species_names)
        w_exprs[idx] = expr
        idx += 1
    # Three-body: left * M(r)
    for i, r in enumerate(three_body):
        left = r[1]
        expr = _build_reactant_expr(left, species_names)
        m_term = "m(" + str(i + 1) + ")"
        w_exprs[idx] = expr + " * " + m_term if expr != "1.0_WP" else m_term
        idx += 1
    # Three-body backward: right * M(r) for same three-body
    for i, r in enumerate(three_body):
        if not r[3]:
            continue
        right = r[2]
        expr = _build_reactant_expr(right, species_names)
        m_term = "m(" + str(i + 1) + ")"
        tb_rev_j = tb_rev_idx[i]
        w_exprs[i_tb_end + tb_rev_j] = expr + " * " + m_term if expr != "1.0_WP" else m_term
    # Falloff forward: k already includes [M] via Pr=k0*M/k_inf; rate = k * product(reactants excl. M)
    idx = i_tb_bwd_end  # falloff forward starts after three-body backward block
    for i, r in enumerate(falloff):
        left = r[1]
        expr = _build_reactant_expr(left, species_names)
        w_exprs[idx] = expr
        idx += 1
    # Falloff backward: same, k includes [M]; rate = k * product(reactants excl. M)
    for i, r in enumerate(falloff):
        if not r[3]:
            continue
        right = r[2]
        expr = _build_reactant_expr(right, species_names)
        fo_rev_j = fo_rev_idx[i]
        w_exprs[i_fo_end + fo_rev_j] = expr
    idx = i_fo_bwd_end
    # Plog forward
    for j, r in enumerate(plog):
        left = r[1]
        expr = _build_reactant_expr(left, species_names)
        w_exprs[idx] = expr
        idx += 1
    # Plog backward
    for j, r in enumerate(plog):
        if not r[3]:
            continue
        right = r[2]
        expr = _build_reactant_expr(right, species_names)
        plog_rev_j = plog_rev_idx[j]
        w_exprs[i_plog_bwd_start - 1 + plog_rev_j - 1] = expr
    # Build stoichiometric contributions for production rates: cdot(i) = sum_r nu(i,r)*w(r)
    # nu(i,r) = (coeff of species i in products) - (coeff of species i in reactants)
    # Build list of (w_idx, left, right) for all reactions in w order
    # Also build (w_idx, left, right, rev, is_backward, rid) for reaction index file
    reactions_for_nu = []  # (w_idx, left, right)
    reactions_for_index = []  # (w_idx_1based, left, right, rev, is_backward, rid)
    idx = 0
    for r in arrhenius_fwd:
        rid = r[0]
        reactions_for_nu.append((idx, r[1], r[2]))
        reactions_for_index.append((idx + 1, r[1], r[2], False, False, rid))
        idx += 1
    for r in arrhenius_bwd:
        rid = r[0]
        fwd_left, fwd_right = None, None
        for af in arrhenius_fwd:
            if af[0] == rid:
                fwd_left, fwd_right = af[1], af[2]
                break
        # Backward: reactants = forward products, products = forward reactants
        reactions_for_nu.append((idx, fwd_right if fwd_right else [], fwd_left if fwd_left else []))
        reactions_for_index.append((idx + 1, fwd_right if fwd_right else [], fwd_left if fwd_left else [], True, True, rid))
        idx += 1
    for r in three_body:
        reactions_for_nu.append((idx, r[1], r[2]))
        reactions_for_index.append((idx + 1, r[1], r[2], r[3], False, r[0]))
        idx += 1
    for i, r in enumerate(three_body):
        if r[3]:
            reactions_for_nu.append((i_tb_end + tb_rev_idx[i], r[2], r[1]))
            reactions_for_index.append((i_tb_end + tb_rev_idx[i] + 1, r[2], r[1], True, True, r[0]))
    idx = i_tb_bwd_end
    for r in falloff:
        reactions_for_nu.append((idx, r[1], r[2]))
        reactions_for_index.append((idx + 1, r[1], r[2], r[3], False, r[0]))
        idx += 1
    for i, r in enumerate(falloff):
        if r[3]:
            reactions_for_nu.append((i_fo_end + fo_rev_idx[i], r[2], r[1]))
            reactions_for_index.append((i_fo_end + fo_rev_idx[i] + 1, r[2], r[1], True, True, r[0]))
    idx = i_fo_bwd_end
    for r in plog:
        reactions_for_nu.append((idx, r[1], r[2]))
        reactions_for_index.append((idx + 1, r[1], r[2], r[3], False, r[0]))
        idx += 1
    for j, r in enumerate(plog):
        if r[3]:
            reactions_for_nu.append((i_plog_bwd_start - 1 + plog_rev_idx[j] - 1, r[2], r[1]))
            reactions_for_index.append((i_plog_bwd_start + plog_rev_idx[j] - 1, r[2], r[1], True, True, r[0]))

    # For each species, collect (w_idx, nu) with nu != 0
    def nu_coeff(left, right, sp_name):
        def resolve(sp, c):
            r = _resolve_species_name(sp, species_names)
            if r is None and 1 <= c <= 9:
                r = _resolve_species_name(str(c) + sp, species_names)
                if r is not None:
                    c = 1
            return r, c

        n = 0
        for sp, c in right:
            r, c = resolve(sp, c)
            if r == sp_name:
                n += c
        for sp, c in left:
            r, c = resolve(sp, c)
            if r == sp_name:
                n -= c
        return n

    cdot_terms = [[] for _ in range(nS)]  # cdot_terms[i] = [(w_idx, nu), ...]
    for w_idx, left, right in reactions_for_nu:
        for sp_name in species_names:
            nu = nu_coeff(left, right, sp_name)
            if nu != 0:
                i = species_names.index(sp_name)
                cdot_terms[i].append((w_idx, nu))

    # Mass conservation check: for each reaction, sum(nu_i * W_i) = 0
    # Mass-based production: sum_i (W_i * cdot_i) = sum_r w_r * sum_i (nu_i,r * W_i) = 0
    mass_cons_tol = 1e-10
    mass_cons_violations = []
    def resolve_for_mass(sp, c):
        r = _resolve_species_name(sp, species_names)
        if r is None and 1 <= c <= 9:
            r = _resolve_species_name(str(c) + sp, species_names)
            if r is not None:
                c = 1
        return r, c

    for w_idx, left, right in reactions_for_nu:
        dm = 0.0
        for sp_name, coeff in right:
            resolved, coeff = resolve_for_mass(sp_name, coeff)
            if resolved is not None:
                i = species_names.index(resolved)
                dm += coeff * W_sp[i]
        for sp_name, coeff in left:
            resolved, coeff = resolve_for_mass(sp_name, coeff)
            if resolved is not None:
                i = species_names.index(resolved)
                dm -= coeff * W_sp[i]
        if abs(dm) > mass_cons_tol:
            mass_cons_violations.append((w_idx + 1, dm, left, right))
    if mass_cons_violations:
        print("WARNING: Mass conservation violations (reaction, delta_mass kg/mol):", file=sys.stderr)
        for w_idx, dm, left, right in mass_cons_violations[:10]:
            print(f"  w({w_idx}): delta={dm:.2e}  {' + '.join(f'{c}{s}' for s,c in left)} => {' + '.join(f'{c}{s}' for s,c in right)}", file=sys.stderr)
        if len(mass_cons_violations) > 10:
            print(f"  ... and {len(mass_cons_violations)-10} more", file=sys.stderr)

    # Write reaction index file: fcmech index, equation, YAML label
    # Also build eq_strs for inline comments in get_reaction_rates
    out_stem = Path(output_path).stem
    if out_stem == "chem_data_fc":
        index_path = Path(output_path).parent / "chem_data_reactions.txt"
    else:
        index_path = Path(output_path).parent / (out_stem + "_reactions.txt")
    reactions_for_index.sort(key=lambda x: x[0])
    idx_width = len(str(nR))
    eq_width = 60
    eq_strs = [""] * nR  # equation string for each w(i), for inline comment
    with open(index_path, "w") as fidx:
        fidx.write("chem_data reaction index\n" if out_stem == "chem_data_fc" else "fcmech reaction index\n")
        fidx.write("=" * (idx_width + eq_width + 30) + "\n")
        fidx.write(f"{'w_idx':>{idx_width}}  {'Equation':<{eq_width}}  YAML label\n")
        fidx.write("-" * (idx_width + eq_width + 30) + "\n")
        yaml_labels = extract_yaml_labels(yaml_path) if yaml_path else {}
        for w_idx, left, right, rev, is_backward, rid in reactions_for_index:
            arrow = "=>" if is_backward else ("<=>" if rev else "=>")
            left_str = _format_reaction_side(left)
            right_str = _format_reaction_side(right)
            eq_str = f"{left_str} {arrow} {right_str}"
            if w_idx <= nR:
                eq_strs[w_idx - 1] = eq_str
            label = yaml_labels.get(rid, "")
            fidx.write(f"{w_idx:>{idx_width}}  {eq_str:<{eq_width}}  {label}\n")
    print(f"  Reaction index: {index_path}")

    lines = []
    yaml_name = Path(yaml_path).name if yaml_path else "(unknown)"
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append(f"!  FILE {Path(output_path).name}")
    lines.append("!  Module for chemical kinetics in NGA2 (finite chemistry)")
    lines.append("!  Generated by yaml2nga.py from YAML mechanism. Do not edit manually.")
    lines.append(f"!  Source: {yaml_name}")
    lines.append("!--------------------------------------------------------------------------------------------------")
    lines.append("")
    lines.append("module fcmech")
    lines.append("  use precision")
    lines.append("  use string")
    lines.append("  implicit none")
    lines.append("")
    lines.append("  !--------------------------------------------------------------------------------------------------")
    lines.append("  !  VARIABLE DEFINITIONS")
    lines.append("  !--------------------------------------------------------------------------------------------------")
    lines.append("  !  Rcst           : universal gas constant (J/(mol·K))")
    lines.append("  !  nS             : number of species")
    lines.append("  !  nR             : number of reactions")
    lines.append("  !  sXXX           : species index (1..nS) for species XXX")
    lines.append("  !  W_sp           : molar mass (kg/mol), dimension(nS)")
    lines.append("  !  koveps         : k/epsilon for Lennard-Jones (dimensionless), dimension(nS)")
    lines.append("  !  mucoeff        : viscosity coefficient for Chapman-Enskog, dimension(nS)")
    lines.append("  !  Dcoeffs        : binary diffusion coefficients, dimension(nS,nS)")
    lines.append("  !  Ocoeffs        : sqrt(k/eps_i * k/eps_j) for reduced T* in collision integral, dimension(nS,nS)")
    lines.append("  !  Cpsp, hsp      : species Cp (J/(mol·K)) and enthalpy (J/mol), dimension(nS), set by fcmech_thermodata")
    lines.append("  !  T_mid          : NASA7 mid temperature (K), dimension(nS)")
    lines.append("  !  thermo_coeffs  : NASA7 polynomial coefficients, dimension(nS,14), cols 1-7 low T, 8-14 high T")
    lines.append("  !  n_arr          : number of Arrhenius reactions (forward + backward)")
    lines.append("  !  n_tb           : number of three-body reactions")
    lines.append("  !  n_fo           : number of falloff reactions")
    lines.append("  !  n_plog         : number of pressure-dependent (plog) reactions")
    lines.append("  !  n_plog_pts     : total number of plog pressure points")
    lines.append("  !  i_arr_end      : last index of Arrhenius block (1..i_arr_end)")
    lines.append("  !  i_tb_end       : last index of three-body block (i_arr_end+1..i_tb_end)")
    lines.append("  !  i_tb_bwd_end   : last index of three-body backward block (i_tb_end+1..i_tb_bwd_end)")
    lines.append("  !  i_fo_end       : last index of falloff block (i_tb_bwd_end+1..i_fo_end)")
    lines.append("  !  i_plog_start   : first index of plog block")
    lines.append("  !  ln_A_arr       : ln(A) for Arrhenius, dimension(n_arr)")
    lines.append("  !  b_arr          : temperature exponent for Arrhenius, dimension(n_arr)")
    lines.append("  !  E_R_arr        : Ea/R (K) for Arrhenius, dimension(n_arr)")
    lines.append("  !  ln_A_tb,b_tb,E_R_tb : same for three-body, dimension(n_tb)")
    lines.append("  !  ln_A0_fo,b0_fo,E_R0_fo   : low-P falloff params, dimension(n_fo)")
    lines.append("  !  ln_Ainf_fo,binf_fo,E_Rinf_fo : high-P falloff params, dimension(n_fo)")
    lines.append("  !  Troe_a,Troe_T3,Troe_T1,Troe_T2 : Troe falloff centering params, dimension(n_fo)")
    lines.append("  !  n_eff          : number of species with non-default third-body efficiencies")
    lines.append("  !  eff_species_idx: species indices (1..nS) for efficiency matrix columns, dimension(n_eff)")
    lines.append("  !  eff_matrix_reduced : efficiency matrix (reaction,eff_species), dimension(n_m,n_eff)")
    lines.append("  !  default_eff    : per-reaction default efficiency for species not in eff_species (Cantera default-efficiency)")
    lines.append("  !  P_plog         : pressure points (Pa) for plog, dimension(n_plog_pts)")
    lines.append("  !  ln_A_plog,b_plog,E_R_plog : plog Arrhenius params, dimension(n_plog_pts)")
    lines.append("  !  plog_reac_start, plog_reac_len : plog indexing per reaction")
    lines.append("  !  nA             : number of atom types")
    lines.append("  !  atom_names     : atom names (character len=2), dimension(nA)")
    lines.append("  !  atom_masses    : atom molar masses (kg/mol), dimension(nA)")
    lines.append("  !  comp_matrix    : comp_matrix(a,s) = count of atom a in species s, dimension(nA,nS)")
    lines.append("  !--------------------------------------------------------------------------------------------------")
    lines.append("")
    nA = mech.get("nA", 0)
    lines.append(f"  real(WP), parameter :: Rcst = {R_UNIV}_WP")
    lines.append("  integer, parameter :: nS = " + str(nS))
    lines.append("  integer, parameter :: nR = " + str(nR))
    lines.append("  integer, parameter :: nA = " + str(nA))
    lines.append("")

    # Species indices (sXXX = 1..nS)
    lines.append("  ! --- Species indices ---")
    for i, nm in enumerate(mech["species_names"]):
        nm_str = str(nm)  # ensure str (YAML may parse NO/yes as bool)
        sn = "s" + nm_str.replace("(", "").replace(")", "").replace("+", "").replace("-", "").replace(" ", "")[:24]
        lines.append("  integer, parameter :: " + sn + " = " + str(i + 1))
    lines.append("")

    # Molar masses
    W_sp = []
    for sp in mech["species"]:
        if isinstance(sp, dict):
            W_sp.append(compute_molar_mass(sp))
        else:
            W_sp.append(0.03)
    lines.append("  ! --- Molar masses (kg/mol) ---")
    lines.append("  real(WP), parameter, dimension(nS) :: W_sp = (/ &")
    for i, w in enumerate(W_sp):
        if i < len(W_sp) - 1:
            lines.append("       " + f90_real(w) + ", & ! " + mech["species_names"][i] + " &")
        else:
            lines.append("       " + f90_real(w) + " /)  ! " + mech["species_names"][i])
    lines.append("")

    # Atom data: names, masses, composition matrix (nA x nS)
    if nA > 0:
        lines.append("  ! --- Atom names (character len=2) ---")
        lines.append("  character(len=2), parameter, dimension(nA) :: atom_names = (/ &")
        for i, an in enumerate(mech["atom_names"]):
            an2 = str(an)[:2].ljust(2)
            lines.append("       '" + an2 + "'" + (", &" if i < nA - 1 else " /)"))
        lines.append("")
        lines.append("  ! --- Atom molar masses (kg/mol) ---")
        lines.append("  real(WP), parameter, dimension(nA) :: atom_masses = (/ &")
        for i, am in enumerate(mech["atom_masses"]):
            lines.append("       " + f90_real(float(am)) + (", &" if i < nA - 1 else " /)"))
        lines.append("")
        lines.append("  ! --- Species composition: comp_matrix(a,s) = count of atom a in species s ---")
        lines.extend(f90_int_2d_array_parameter("comp_matrix", mech["comp_matrix"], nA, nS))
        lines.append("")

    # Transport: k/eps, mucoeff, diameter for diffusion
    koveps = []
    mucoeff = []
    diameters = []
    for sp in mech["species"]:
        if isinstance(sp, dict):
            trans = get_transport(sp)
            eps_k = trans["well_depth"]
            d = trans["diameter"]
            w = compute_molar_mass(sp)
            koveps.append(1.0 / eps_k if eps_k > 0 else 0.01)
            mucoeff.append(2.67e-6 * np.sqrt(w * 1000) / (d**2) if d > 0 else 1e-5)
            diameters.append(d)
        else:
            koveps.append(0.01)
            mucoeff.append(1e-5)
            diameters.append(3.0)

    # Precompute Dcoeffs, Ocoeffs for diffusion
    Dcoeffs = np.zeros((nS, nS))
    Ocoeffs_arr = np.zeros((nS, nS))
    for i in range(nS):
        for j in range(nS):
            if i != j:
                sigma_ij = 0.5 * (diameters[i] + diameters[j])
                M_ij = 2.0 * W_sp[i] * W_sp[j] * 1000 / (W_sp[i] * 1000 + W_sp[j] * 1000)
                Dcoeffs[i, j] = 1.858e-7 * np.sqrt(300.0) * np.sqrt(1.0 / M_ij) / (sigma_ij**2)
                # Ocoeffs = sqrt(k/eps_i * k/eps_j) for reduced T* in collision integral (ARCANE)
                Ocoeffs_arr[i, j] = np.sqrt(koveps[i] * koveps[j])
            else:
                Dcoeffs[i, j] = 1.0
                Ocoeffs_arr[i, j] = 1.0
    lines.append("  ! --- Transport: k/epsilon (Lennard-Jones) ---")
    lines.append("  real(WP), parameter, dimension(nS) :: koveps = (/ &")
    for i in range(nS):
        if i < nS - 1:
            lines.append("       " + f90_real(koveps[i]) + ", & ! " + mech["species_names"][i] + " &")
        else:
            lines.append("       " + f90_real(koveps[i]) + " /)  ! " + mech["species_names"][i])
    lines.append("")
    lines.append("  ! --- Transport: viscosity coefficients ---")
    lines.append("  real(WP), parameter, dimension(nS) :: mucoeff = (/ &")
    for i in range(nS):
        if i < nS - 1:
            lines.append("       " + f90_real(mucoeff[i]) + ", & ! " + mech["species_names"][i] + " &")
        else:
            lines.append("       " + f90_real(mucoeff[i]) + " /)  ! " + mech["species_names"][i])
    lines.append("")

    lines.append("  ! --- Binary diffusion coefficients ---")
    lines.extend(f90_2d_array_parameter("Dcoeffs", Dcoeffs, nS, nS))
    lines.append("  ! --- Collision integral coefficients (Ocoeffs) ---")
    lines.extend(f90_2d_array_parameter("Ocoeffs", Ocoeffs_arr, nS, nS))
    lines.append("")

    lines.append("  ! --- Thermodynamic: Cp and enthalpy (module variables, set by fcmech_thermodata) ---")
    lines.append("  real(WP), dimension(nS) :: Cpsp, hsp")
    lines.append("")

    # Thermo parameter arrays (already built above for backward rate fitting)
    lines.append("  ! --- NASA7 thermo: mid temperature (K) ---")
    lines.append("  real(WP), parameter, dimension(nS) :: T_mid = (/ &")
    for i in range(nS):
        if i < nS - 1:
            lines.append("       " + f90_real(T_mid_sp[i]) + ", & ! " + mech["species_names"][i] + " &")
        else:
            lines.append("       " + f90_real(T_mid_sp[i]) + " /)  ! " + mech["species_names"][i])
    lines.append("")
    lines.append("  ! --- NASA7 thermo: polynomial coefficients (low cols 1-7, high cols 8-14) ---")
    lines.extend(f90_2d_array_parameter("thermo_coeffs", thermo_coeffs, nS, 14))
    lines.append("")

    # Arrhenius coefficient parameters: ln(A), b, E/R
    lines.append("  ! --- Reaction counts and index boundaries ---")
    # ln(k) = ln(A) + b*ln(T) - E/(R*T), then k = exp(ln(k))
    lines.append("  integer, parameter :: n_arr = " + str(n_arr))
    lines.append("  integer, parameter :: n_tb = " + str(n_tb))
    lines.append("  integer, parameter :: n_tb_rev = " + str(n_tb_rev))
    lines.append("  integer, parameter :: n_fo = " + str(n_fo))
    lines.append("  integer, parameter :: n_fo_rev = " + str(n_fo_rev))
    lines.append("  integer, parameter :: n_plog = " + str(n_plog))
    lines.append("  integer, parameter :: n_plog_rev = " + str(n_plog_rev))
    lines.append("  integer, parameter :: n_plog_pts = " + str(n_plog_pts))
    lines.append("  integer, parameter :: n_m = " + str(n_m) + "  ! third-body/falloff M groups")
    i_tb_end = n_arr + n_tb
    i_tb_bwd_end = i_tb_end + n_tb_rev
    i_fo_end = i_tb_bwd_end + n_fo
    i_fo_bwd_end = i_fo_end + n_fo_rev
    i_plog_start = i_fo_bwd_end + 1
    i_plog_end = i_plog_start + n_plog - 1 if n_plog > 0 else i_plog_start - 1
    i_plog_bwd_start = i_plog_end + 1
    i_plog_bwd_end = i_plog_bwd_start + n_plog_rev - 1 if n_plog_rev > 0 else i_plog_end
    lines.append("  integer, parameter :: i_arr_end = " + str(n_arr) + ", i_tb_end = " + str(n_arr + n_tb) + ", i_tb_bwd_end = " + str(i_tb_bwd_end) + ", i_fo_end = " + str(i_fo_end) + ", i_fo_bwd_end = " + str(i_fo_bwd_end) + ", i_plog_start = " + str(i_plog_start) + ", i_plog_end = " + str(i_plog_end) + ", i_plog_bwd_start = " + str(i_plog_bwd_start) + ", i_plog_bwd_end = " + str(i_plog_bwd_end))
    lines.append("")
    if n_arr > 0:
        lines.append("  ! --- Arrhenius: ln(A), b, E/R ---")
        lines.extend(f90_array_parameter("ln_A_arr", ln_A_arr))
        lines.append("")
        lines.extend(f90_array_parameter("b_arr", b_arr))
        lines.append("")
        lines.extend(f90_array_parameter("E_R_arr", E_R_arr))
        lines.append("")
    if n_tb > 0:
        lines.append("  ! --- Three-body: ln(A), b, E/R ---")
        lines.extend(f90_array_parameter("ln_A_tb", ln_A_tb))
        lines.append("")
        lines.extend(f90_array_parameter("b_tb", b_tb))
        lines.append("")
        lines.extend(f90_array_parameter("E_R_tb", E_R_tb))
        lines.append("")
    if n_tb_rev > 0:
        lines.append("  ! --- Three-body backward (reversible): ln(A), b, E/R ---")
        lines.extend(f90_array_parameter("ln_A_tb_b", ln_A_tb_b))
        lines.append("")
        lines.extend(f90_array_parameter("b_tb_b", b_tb_b))
        lines.append("")
        lines.extend(f90_array_parameter("E_R_tb_b", E_R_tb_b))
        lines.append("")
        lines.append("  ! --- Three-body reversible index: tb_rev_idx(i)=backward index for 3-body i (1-based), 0 if irreversible ---")
        tb_rev_idx_f90 = [j + 1 if j >= 0 else 0 for j in tb_rev_idx]
        lines.extend(f90_int_array_parameter("tb_rev_idx", tb_rev_idx_f90))
        lines.append("")
    if n_fo > 0:
        lines.append("  ! --- Falloff: low-P (ln_A0,b0,E_R0) and high-P (ln_Ainf,binf,E_Rinf) ---")
        lines.extend(f90_array_parameter("ln_A0_fo", ln_A0_fo))
        lines.append("")
        lines.extend(f90_array_parameter("b0_fo", b0_fo))
        lines.append("")
        lines.extend(f90_array_parameter("E_R0_fo", E_R0_fo))
        lines.append("")
        lines.extend(f90_array_parameter("ln_Ainf_fo", ln_Ainf_fo))
        lines.append("")
        lines.extend(f90_array_parameter("binf_fo", binf_fo))
        lines.append("")
        lines.extend(f90_array_parameter("E_Rinf_fo", E_Rinf_fo))
        lines.append("")
        lines.extend(f90_array_parameter("Troe_a", Troe_a))
        lines.append("")
        lines.extend(f90_array_parameter("Troe_T3", Troe_T3))
        lines.append("")
        lines.extend(f90_array_parameter("Troe_T1", Troe_T1))
        lines.append("")
        lines.extend(f90_array_parameter("Troe_T2", Troe_T2))
        lines.append("")
    if n_fo_rev > 0:
        lines.append("  ! --- Falloff backward (reversible): low-P and high-P ---")
        lines.extend(f90_array_parameter("ln_A0_b_fo", ln_A0_b_fo))
        lines.append("")
        lines.extend(f90_array_parameter("b0_b_fo", b0_b_fo))
        lines.append("")
        lines.extend(f90_array_parameter("E_R0_b_fo", E_R0_b_fo))
        lines.append("")
        lines.extend(f90_array_parameter("ln_Ainf_b_fo", ln_Ainf_b_fo))
        lines.append("")
        lines.extend(f90_array_parameter("binf_b_fo", binf_b_fo))
        lines.append("")
        lines.extend(f90_array_parameter("E_Rinf_b_fo", E_Rinf_b_fo))
        lines.append("")
        lines.append("  ! --- Falloff reversible index: fo_rev_idx(i)=backward index for falloff i (1-based), 0 if irreversible ---")
        fo_rev_idx_f90 = [j + 1 if j >= 0 else 0 for j in fo_rev_idx]
        lines.extend(f90_int_array_parameter("fo_rev_idx", fo_rev_idx_f90))
        lines.append("")
    # Third-body efficiency: reduced matrix (n_m x n_eff) + sum(c) for species with default 1.0
    # M(r) = sum_j eff_matrix_reduced(r,j)*c(eff_species_idx(j)) + (sum(c) - sum over eff species of c)
    if n_m > 0:
        lines.append("  ! --- Third-body efficiencies: species with non-default eff ---")
        lines.append("  integer, parameter :: n_eff = " + str(n_eff))
        lines.append("  ! Species with non-default efficiencies: " + ", ".join(eff_species_names))
        if n_eff > 0:
            lines.extend(f90_int_array_parameter("eff_species_idx", eff_species_idx))
            lines.append("")
            lines.extend(f90_2d_array_parameter("eff_matrix_reduced", eff_matrix_reduced, n_m, n_eff))
            lines.append("")
        lines.extend(f90_array_parameter("default_eff", default_eff_arr))
        lines.append("")
        lines.append("")
    if n_plog > 0:
        lines.append("  ! --- Plog: pressure points (Pa) and Arrhenius params ---")
        lines.extend(f90_array_parameter("P_plog", P_plog))
        lines.append("")
        lines.extend(f90_array_parameter("ln_A_plog", ln_A_plog))
        lines.append("")
        lines.extend(f90_array_parameter("b_plog", b_plog))
        lines.append("")
        lines.extend(f90_array_parameter("E_R_plog", E_R_plog))
        lines.append("")
        lines.extend(f90_int_array_parameter("plog_reac_start", plog_reac_start))
        lines.append("")
        lines.extend(f90_int_array_parameter("plog_reac_len", plog_reac_len))
        lines.append("")
        lines.append("  ! --- Plog reversible index: plog_rev_idx(j)=backward index for plog j (1-based), 0 if irreversible ---")
        lines.extend(f90_int_array_parameter("plog_rev_idx", plog_rev_idx))
        lines.append("")
        if n_plog_rev > 0:
            lines.append("  ! --- Plog backward (reversible): ln_A, b, E/R at each pressure point ---")
            lines.extend(f90_array_parameter("ln_A_b_plog", ln_A_b_plog))
            lines.append("")
            lines.extend(f90_array_parameter("b_b_plog", b_b_plog))
            lines.append("")
            lines.extend(f90_array_parameter("E_R_b_plog", E_R_b_plog))
            lines.append("")
    lines.append("  contains")
    lines.append("")

    # fcmech_omegamu, fcmech_omegaD (Neufeld collision integrals, from ARCANE)
    lines.append("  ! Collision integral Omega^(2,2)* for viscosity (Neufeld et al.)")
    lines.append("  real(WP) function fcmech_omegamu(T_)")
    lines.append("    implicit none")
    lines.append("    ! T_ : reduced temperature T* = T*k_B/epsilon for collision integral Omega^(2,2)*")
    lines.append("    real(WP), intent(in) :: T_")
    lines.append("    real(WP), parameter, dimension(9) :: mArray = (/ &")
    lines.append("         3.3530622607_WP, 2.53272006_WP, 2.9024238575_WP, &")
    lines.append("         0.11186138893_WP, 0.8662326188_WP, 1.3913958626_WP, &")
    lines.append("         3.158490576_WP, 0.18973411754_WP, 0.00018682962894_WP/)")
    lines.append("    integer :: arrIndex")
    lines.append("    real(WP) :: omegamu_Nr, omegamu_Dr")
    lines.append("    omegamu_Dr = mArray(9)")
    lines.append("    do arrIndex = 1, 4")
    lines.append("      omegamu_Dr = mArray(9-arrIndex) + T_*omegamu_Dr")
    lines.append("    end do")
    lines.append("    omegamu_Nr = mArray(4)")
    lines.append("    do arrIndex = 1, 3")
    lines.append("      omegamu_Nr = mArray(4-arrIndex) + T_*omegamu_Nr")
    lines.append("    end do")
    lines.append("    fcmech_omegamu = omegamu_Nr/omegamu_Dr")
    lines.append("  end function fcmech_omegamu")
    lines.append("")
    lines.append("  ! Collision integral Omega^(1,1)* for diffusion (Neufeld et al., from ARCANE)")
    lines.append("  real(WP) function fcmech_omegaD(T_)")
    lines.append("    implicit none")
    lines.append("    ! T_ : reduced temperature T* for collision integral Omega^(1,1)* (diffusion)")
    lines.append("    real(WP), intent(in) :: T_")
    lines.append("    real(WP), parameter, dimension(9) :: mArray = (/ &")
    lines.append("         6.8728271691_WP, 9.4122316321_WP, 7.7442359037_WP, &")
    lines.append("         0.23424661229_WP, 1.45337701568_WP, 5.2269794238_WP, &")
    lines.append("         9.7108519575_WP, 0.46539437353_WP, 0.00041908394781_WP/)")
    lines.append("    integer :: arrIndex")
    lines.append("    real(WP) :: omegaD_Nr, omegaD_Dr")
    lines.append("    omegaD_Dr = mArray(9)")
    lines.append("    do arrIndex = 1, 4")
    lines.append("      omegaD_Dr = mArray(9-arrIndex) + T_*omegaD_Dr")
    lines.append("    end do")
    lines.append("    omegaD_Nr = mArray(4)")
    lines.append("    do arrIndex = 1, 3")
    lines.append("      omegaD_Nr = mArray(4-arrIndex) + T_*omegaD_Nr")
    lines.append("    end do")
    lines.append("    fcmech_omegaD = omegaD_Nr/omegaD_Dr")
    lines.append("  end function fcmech_omegaD")
    lines.append("")

    # getlindratecoeff for falloff (Troe/Lindemann, from ARCANE custom_printing)
    lines.append("  !--------------------------------------------------------------------------------------------------")
    lines.append("  !  function getlindratecoeff")
    lines.append("  !  Pressure-dependent rate coeff (Troe/Lindemann)")
    lines.append("  !--------------------------------------------------------------------------------------------------")
    lines.append("  real(WP) function getlindratecoeff(Tloc,k0,kinf,fc,concin,Ploc)")
    lines.append("    implicit none")
    lines.append("    ! Tloc,Ploc : temperature (K), pressure (Pa)")
    lines.append("    ! k0, kinf : low-P and high-P rate coefficients")
    lines.append("    ! fc      : Troe centering factor F_cent")
    lines.append("    ! concin  : third-body concentration (M); if<=0 uses P/(R*T)")
    lines.append("    real(WP), intent(in) :: Tloc,k0,kinf,fc,concin,Ploc")
    lines.append("    real(WP) :: conc,redP,ntmp,ccoeff,dcoeff,lgknull,f")
    lines.append("    if (concin > 0.0_WP) then")
    lines.append("      conc = concin")
    lines.append("    else")
    lines.append("      conc = Ploc / (Rcst * Tloc)")
    lines.append("    end if")
    lines.append("    redP = ABS(k0) * conc / MAX(ABS(kinf), TINY(1.0_WP)) + TINY(1.0_WP)")
    lines.append("    ntmp = 0.75_WP - 1.27_WP * LOG10(MAX(fc, TINY(1.0_WP)))")
    lines.append("    ccoeff = -0.4_WP - 0.67_WP * LOG10(MAX(fc, TINY(1.0_WP)))")
    lines.append("    dcoeff = 0.14_WP")
    lines.append("    lgknull = LOG10(redP)")
    lines.append("    f = (lgknull+ccoeff)/(ntmp-dcoeff*(lgknull+ccoeff))")
    lines.append("    f = fc**(1.0_WP / (f*f + 1.0_WP))")
    lines.append("    getlindratecoeff = kinf * f * redP / (1.0_WP + redP)")
    lines.append("  end function getlindratecoeff")
    lines.append("")

    # get_thirdbodies: M(r) = default_eff(r)*sum_c_rest + sum_j eff_matrix_reduced(r,j)*c_eff(j)
    # sum_c_rest = sum(c) for species not in eff_species_idx; default_eff from Cantera default-efficiency
    lines.append("  subroutine get_thirdbodies(M,c)")
    lines.append("    implicit none")
    lines.append("    ! c(nS)      : species concentrations (mol/m^3)")
    lines.append("    ! M(n_m)     : effective third-body concentration per reaction (mol/m^3)")
    lines.append("    real(WP), dimension(nS), intent(in) :: c")
    lines.append("    real(WP), dimension(" + str(n_m) + "), intent(out) :: M")
    if n_m > 0:
        if n_eff > 0:
            lines.append("    ! c_eff     : c(eff_species_idx), concentration of species with explicit efficiencies")
            lines.append("    ! sum_c_rest: sum(c) for species not in eff_species_idx")
            lines.append("    real(WP), dimension(n_eff) :: c_eff")
            lines.append("    real(WP) :: sum_c_rest")
            lines.append("    integer :: j")
            lines.append("    do j = 1, n_eff")
            lines.append("      c_eff(j) = c(eff_species_idx(j))")
            lines.append("    end do")
            lines.append("    sum_c_rest = sum(c) - sum(c_eff)")
            lines.append("    M(1:" + str(n_m) + ") = default_eff(1:" + str(n_m) + ") * sum_c_rest + matmul(eff_matrix_reduced, c_eff)")
        else:
            lines.append("    M(1:" + str(n_m) + ") = default_eff(1:" + str(n_m) + ") * sum(c)")
    else:
        lines.append("    M(1:1) = sum(c)")
    lines.append("  end subroutine get_thirdbodies")
    lines.append("")

    # get_rate_coefficients - VECTORIZED using parameter arrays
    lines.append("  !--------------------------------------------------------------------------------------------------")
    lines.append("  !  subroutine get_rate_coefficients")
    lines.append("  !  Evaluate rate coefficients (Arrhenius, 3-body, falloff)")
    lines.append("  !--------------------------------------------------------------------------------------------------")
    lines.append("  subroutine get_rate_coefficients(k,M,Tloc,Ploc)")
    lines.append("    implicit none")
    lines.append("    ! k(nR)     : rate coefficients (m^3/(mol·s) or 1/s)")
    lines.append("    ! M(n_m)    : effective third-body concentration from get_thirdbodies")
    lines.append("    ! Tloc, Ploc: temperature (K), pressure (Pa)")
    lines.append("    real(WP), dimension(nR), intent(out) :: k")
    lines.append("    real(WP), dimension(:), intent(in) :: M")
    lines.append("    real(WP), intent(in) :: Tloc, Ploc")
    lines.append("    real(WP) :: T_log, T_inv, fc_fo, ln_k0, ln_kinf")
    lines.append("    integer :: i")
    lines.append("    T_log = log(Tloc)")
    lines.append("    T_inv = 1.0_WP/Tloc")
    lines.append("    ! Arrhenius: ln(k) = ln(A) + b*ln(T) - E/(R*T), then k = exp(ln(k))")
    if n_arr > 0:
        lines.append("    k(1:n_arr) = exp(ln_A_arr(:) + b_arr(:)*T_log - E_R_arr(:)*T_inv)")
    lines.append("    ! Three-body: k = A*T^b*exp(-E/RT); M is applied in get_reaction_rates, not here")
    if n_tb > 0:
        lines.append("    k(i_arr_end+1:i_tb_end) = exp(ln_A_tb(:) + b_tb(:)*T_log - E_R_tb(:)*T_inv)")
    if n_tb_rev > 0:
        lines.append("    do i = 1, n_tb")
        lines.append("      if (tb_rev_idx(i) > 0) then")
        lines.append("        k(i_tb_end+tb_rev_idx(i)) = exp(ln_A_tb_b(tb_rev_idx(i)) + b_tb_b(tb_rev_idx(i))*T_log - E_R_tb_b(tb_rev_idx(i))*T_inv)")
        lines.append("      end if")
        lines.append("    end do")
    lines.append("    ! Falloff (Troe/Lindemann)")
    if n_fo > 0:
        lines.append("    do i = 1, n_fo")
        lines.append("      fc_fo = (1.0_WP - Troe_a(i))*EXP(-Tloc/Troe_T3(i)) + Troe_a(i)*EXP(-Tloc/Troe_T1(i))")
        lines.append("      if (Troe_T2(i) < 1.0e29_WP) fc_fo = fc_fo + EXP(-Troe_T2(i)/Tloc)")
        lines.append("      ln_k0 = ln_A0_fo(i) + b0_fo(i)*T_log - E_R0_fo(i)*T_inv")
        lines.append("      ln_kinf = ln_Ainf_fo(i) + binf_fo(i)*T_log - E_Rinf_fo(i)*T_inv")
        lines.append("      k(i_tb_bwd_end+i) = getlindratecoeff(Tloc, &")
        lines.append("        exp(ln_k0), exp(ln_kinf), fc_fo, M(n_tb+i), Ploc)")
        if n_fo_rev > 0:
            lines.append("      if (fo_rev_idx(i) > 0) then")
            lines.append("        ln_k0 = ln_A0_b_fo(fo_rev_idx(i)) + b0_b_fo(fo_rev_idx(i))*T_log - E_R0_b_fo(fo_rev_idx(i))*T_inv")
            lines.append("        ln_kinf = ln_Ainf_b_fo(fo_rev_idx(i)) + binf_b_fo(fo_rev_idx(i))*T_log - E_Rinf_b_fo(fo_rev_idx(i))*T_inv")
            lines.append("        k(i_fo_end+fo_rev_idx(i)) = getlindratecoeff(Tloc, &")
            lines.append("          exp(ln_k0), exp(ln_kinf), fc_fo, M(n_tb+i), Ploc)")
            lines.append("      end if")
        lines.append("    end do")
    if n_plog > 0:
        lines.append("    call get_pdep_rate_coefficients(k, Tloc, Ploc)")
    lines.append("  end subroutine get_rate_coefficients")
    lines.append("")

    # get_pdep_rate_coefficients - Plog (logarithmic interpolation)
    if n_plog > 0:
        lines.append("  !--------------------------------------------------------------------------------------------------")
        lines.append("  !  subroutine get_pdep_rate_coefficients")
        lines.append("  !  Pressure-dependent Arrhenius (Plog) with log-log interpolation")
        lines.append("  !--------------------------------------------------------------------------------------------------")
        lines.append("  subroutine get_pdep_rate_coefficients(k, Tloc, Ploc)")
        lines.append("    implicit none")
        lines.append("    ! k(nR)     : rate coefficients (inout; plog block overwritten)")
        lines.append("    real(WP), dimension(nR), intent(inout) :: k")
        lines.append("    real(WP), intent(in) :: Tloc, Ploc")
        lines.append("    real(WP), dimension(n_plog_pts) :: ln_k_pdep")
        lines.append("    real(WP) :: T_log, T_inv, P_log, P_inf, P_sup, a, b, log_P_sup, log_P_inf")
        lines.append("    integer :: j, is, ilen, idx, index_inf, index_sup, ipt")
        lines.append("    T_log = log(Tloc)")
        lines.append("    T_inv = 1.0_WP/Tloc")
        lines.append("    P_log = log(Ploc)")
        lines.append("    do j = 1, n_plog")
        lines.append("      is = plog_reac_start(j)")
        lines.append("      ilen = plog_reac_len(j)")
        lines.append("      do ipt = 1, ilen")
        lines.append("        idx = is + ipt - 1")
        lines.append("        ln_k_pdep(idx) = ln_A_plog(idx) + b_plog(idx)*T_log - E_R_plog(idx)*T_inv")
        lines.append("      end do")
        lines.append("      P_inf = 0.0_WP")
        lines.append("      P_sup = 1.0e10_WP")
        lines.append("      index_inf = is")
        lines.append("      index_sup = is + ilen - 1")
        lines.append("      do ipt = 1, ilen")
        lines.append("        idx = is + ipt - 1")
        lines.append("        if (Ploc > P_plog(idx)) then")
        lines.append("          if (P_plog(idx) > P_inf) then")
        lines.append("            P_inf = P_plog(idx)")
        lines.append("            index_inf = idx")
        lines.append("          end if")
        lines.append("        end if")
        lines.append("        if (Ploc < P_plog(idx)) then")
        lines.append("          if (P_plog(idx) < P_sup) then")
        lines.append("            P_sup = P_plog(idx)")
        lines.append("            index_sup = idx")
        lines.append("          end if")
        lines.append("        end if")
        lines.append("      end do")
        lines.append("      log_P_sup = log(P_plog(index_sup))")
        lines.append("      log_P_inf = log(P_plog(index_inf))")
        lines.append("      if (abs(log_P_sup - log_P_inf) < 1e-30_WP) then")
        lines.append("        k(i_plog_start + j - 1) = exp(ln_k_pdep(index_inf))")
        lines.append("      else")
        lines.append("        a = (ln_k_pdep(index_sup) - ln_k_pdep(index_inf)) / (log_P_sup - log_P_inf)")
        lines.append("        b = ln_k_pdep(index_inf) - a * log_P_inf")
        lines.append("        k(i_plog_start + j - 1) = exp(a * P_log + b)")
        lines.append("      end if")
        if n_plog_rev > 0:
            lines.append("      if (plog_rev_idx(j) > 0) then")
            lines.append("        do ipt = 1, ilen")
            lines.append("          idx = is + ipt - 1")
            lines.append("          ln_k_pdep(idx) = ln_A_b_plog(idx) + b_b_plog(idx)*T_log - E_R_b_plog(idx)*T_inv")
            lines.append("        end do")
            lines.append("        if (abs(log_P_sup - log_P_inf) < 1e-30_WP) then")
            lines.append("          k(i_plog_bwd_start + plog_rev_idx(j) - 1) = exp(ln_k_pdep(index_inf))")
            lines.append("        else")
            lines.append("          a = (ln_k_pdep(index_sup) - ln_k_pdep(index_inf)) / (log_P_sup - log_P_inf)")
            lines.append("          b = ln_k_pdep(index_inf) - a * log_P_inf")
            lines.append("          k(i_plog_bwd_start + plog_rev_idx(j) - 1) = exp(a * P_log + b)")
            lines.append("        end if")
            lines.append("      end if")
        lines.append("    end do")
        lines.append("  end subroutine get_pdep_rate_coefficients")
        lines.append("")

    # get_reaction_rates - explicit w(i) = k(i) * product(reactants) [* M(r)]
    lines.append("  subroutine get_reaction_rates(w,k,m,c)")
    lines.append("    implicit none")
    lines.append("    ! c(nS) : species concentrations (mol/m^3)")
    lines.append("    ! k(nR) : rate coefficients")
    lines.append("    ! m(n_m): effective third-body concentrations from get_thirdbodies")
    lines.append("    ! w(nR): reaction rates (mol/(m^3·s))")
    lines.append("    real(WP), dimension(nS), intent(in) :: c")
    lines.append("    real(WP), dimension(nR), intent(in) :: k")
    lines.append("    real(WP), dimension(:), intent(in) :: m")
    lines.append("    real(WP), dimension(nR), intent(out) :: w")
    for i, expr in enumerate(w_exprs):
        if expr:
            comment = f"  ! {eq_strs[i]}" if i < len(eq_strs) and eq_strs[i] else ""
            lines.append("    w(" + str(i + 1) + ") = k(" + str(i + 1) + ") * " + expr + comment)
    lines.append("  end subroutine get_reaction_rates")
    lines.append("")

    # get_production_rates: cdot(i) = sum_r nu(i,r)*w(r)
    lines.append("  subroutine get_production_rates(cdot,w)")
    lines.append("    implicit none")
    lines.append("    ! w(nR)    : reaction rates (mol/(m^3·s))")
    lines.append("    ! cdot(nS) : species production rates (mol/(m^3·s)) = sum_r nu(i,r)*w(r)")
    lines.append("    real(WP), dimension(nR), intent(in) :: w")
    lines.append("    real(WP), dimension(nS), intent(out) :: cdot")
    max_terms_per_line = 5
    for i in range(nS):
        terms = cdot_terms[i]
        sn = _species_f90_name(species_names[i], species_names)
        if not terms:
            lines.append("    cdot(" + sn + ") = 0.0_WP")
            continue
        parts = []
        for w_idx, nu in terms:
            if nu == 1:
                parts.append("+ w(" + str(w_idx + 1) + ")")
            elif nu == -1:
                parts.append("- w(" + str(w_idx + 1) + ")")
            elif nu > 0:
                parts.append("+ " + str(float(nu)) + "_WP * w(" + str(w_idx + 1) + ")")
            else:
                parts.append("- " + str(float(-nu)) + "_WP * w(" + str(w_idx + 1) + ")")
        # Format: cdot(sXXX) = 0.0_WP & + term1 & + term2 & ... (continuation with &)
        line = "    cdot(" + sn + ") = 0.0_WP"
        for j, p in enumerate(parts):
            line += " &"
            lines.append(line)
            line = "         " + p
        lines.append(line)
    lines.append("  end subroutine get_production_rates")
    lines.append("")

    # y2c
    lines.append("  subroutine y2c(y, W_sp, P, T, c)")
    lines.append("    implicit none")
    lines.append("    ! y(nS)     : mass fractions")
    lines.append("    ! W_sp(nS)  : molar masses (kg/mol)")
    lines.append("    ! P, T      : pressure (Pa), temperature (K)")
    lines.append("    ! c(nS)     : molar concentrations (mol/m^3) = rho*y/W_sp, rho=P*W_mix/(R*T)")
    lines.append("    real(WP), dimension(nS), intent(in) :: y")
    lines.append("    real(WP), dimension(nS), intent(in) :: W_sp")
    lines.append("    real(WP), intent(in) :: P, T")
    lines.append("    real(WP), dimension(nS), intent(out) :: c")
    lines.append("    real(WP) :: W_mix, rho")
    lines.append("    W_mix = 1.0_WP / sum(y(1:nS) / W_sp(1:nS))")
    lines.append("    rho = P * W_mix / (Rcst * T)")
    lines.append("    c(1:nS) = rho * y(1:nS) / W_sp(1:nS)")
    lines.append("  end subroutine y2c")
    lines.append("")

    # fcmech_thermodata - NASA7 from parameter arrays (chem_state_class style)
    # thermo_coeffs(i,1:7)=low a0..a6, thermo_coeffs(i,8:14)=high a0..a6
    # h/(RT) = dot_product(a0..a5, th) with th=[1, T/2, T^2/3, T^3/4, T^4/5, 1/T]
    # Cp/R = dot_product(a0..a4, tc) with tc=[1, T, T^2, T^3, T^4]
    lines.append("  subroutine fcmech_thermodata(T)")
    lines.append("    implicit none")
    lines.append("    ! T(K)     : temperature")
    lines.append("    ! th(6)    : coefficient multipliers for h/(RT): [1, T/2, T^2/3, T^3/4, T^4/5, 1/T]")
    lines.append("    ! tc(5)    : coefficient multipliers for Cp/R: [1, T, T^2, T^3, T^4]")
    lines.append("    ! hort     : h/(RT) from dot_product")
    lines.append("    ! cpor     : Cp/R from dot_product")
    lines.append("    real(WP), intent(in) :: T")
    lines.append("    real(WP) :: th(6), tc(5), Tpnm1, hort, cpor")
    lines.append("    integer :: i, n")
    lines.append("    th(1) = 1.0_WP")
    lines.append("    th(6) = 1.0_WP / T")
    lines.append("    Tpnm1 = 1.0_WP")
    lines.append("    do n = 2, 5")
    lines.append("      Tpnm1 = Tpnm1 * T")
    lines.append("      th(n) = Tpnm1 / real(n, WP)")
    lines.append("    end do")
    lines.append("    tc(1) = 1.0_WP")
    lines.append("    do n = 2, 5")
    lines.append("      tc(n) = tc(n-1) * T")
    lines.append("    end do")
    lines.append("    do i = 1, nS")
    lines.append("      if (T < T_mid(i)) then")
    lines.append("        hort = thermo_coeffs(i,1)*th(1) + thermo_coeffs(i,2)*th(2) + thermo_coeffs(i,3)*th(3) + thermo_coeffs(i,4)*th(4) + thermo_coeffs(i,5)*th(5) + thermo_coeffs(i,6)*th(6)")
    lines.append("        cpor = thermo_coeffs(i,1)*tc(1) + thermo_coeffs(i,2)*tc(2) + thermo_coeffs(i,3)*tc(3) + thermo_coeffs(i,4)*tc(4) + thermo_coeffs(i,5)*tc(5)")
    lines.append("      else")
    lines.append("        hort = thermo_coeffs(i,8)*th(1) + thermo_coeffs(i,9)*th(2) + thermo_coeffs(i,10)*th(3) + thermo_coeffs(i,11)*th(4) + thermo_coeffs(i,12)*th(5) + thermo_coeffs(i,13)*th(6)")
    lines.append("        cpor = thermo_coeffs(i,8)*tc(1) + thermo_coeffs(i,9)*tc(2) + thermo_coeffs(i,10)*tc(3) + thermo_coeffs(i,11)*tc(4) + thermo_coeffs(i,12)*tc(5)")
    lines.append("      end if")
    lines.append("      hsp(i) = Rcst * T * hort")
    lines.append("      Cpsp(i) = Rcst * cpor")
    lines.append("    end do")
    lines.append("  end subroutine fcmech_thermodata")
    lines.append("")

    # fcmech_get_thermodata
    lines.append("  subroutine fcmech_get_thermodata(h,cp,T)")
    lines.append("    implicit none")
    lines.append("    ! h(nS), cp(nS) : enthalpy (J/mol), heat capacity (J/(mol·K))")
    lines.append("    ! T(K)          : temperature")
    lines.append("    real(WP), dimension(nS), intent(out) :: h, cp")
    lines.append("    real(WP), intent(in) :: T")
    lines.append("    call fcmech_thermodata(T)")
    lines.append("    h = hsp")
    lines.append("    cp = Cpsp")
    lines.append("  end subroutine fcmech_get_thermodata")
    lines.append("")

    # fcmech_get_speciesnames
    lines.append("  subroutine fcmech_get_speciesnames(names)")
    lines.append("    implicit none")
    lines.append("    ! names(nS) : species names (character)")
    lines.append("    character(len=*), dimension(nS), intent(out) :: names")
    for i, nm in enumerate(mech["species_names"]):
        lines.append("    names(" + str(i + 1) + ") = '" + nm[:30] + "'")
    lines.append("  end subroutine fcmech_get_speciesnames")
    lines.append("")
    if nA > 0:
        lines.append("  subroutine fcmech_get_atomnames(names)")
        lines.append("    implicit none")
        lines.append("    ! names(nA) : atom names (character)")
        lines.append("    character(len=*), dimension(nA), intent(out) :: names")
        for i, an in enumerate(mech["atom_names"]):
            lines.append("    names(" + str(i + 1) + ") = '" + str(an)[:2].ljust(2) + "'")
        lines.append("  end subroutine fcmech_get_atomnames")
        lines.append("")
        lines.append("  subroutine fcmech_get_atommasses(masses)")
        lines.append("    implicit none")
        lines.append("    ! masses(nA) : atom molar masses (kg/mol)")
        lines.append("    real(WP), dimension(nA), intent(out) :: masses")
        lines.append("    masses(1:nA) = atom_masses(1:nA)")
        lines.append("  end subroutine fcmech_get_atommasses")
        lines.append("")
        lines.append("  subroutine fcmech_get_composition(comp)")
        lines.append("    implicit none")
        lines.append("    ! comp(nA,nS) : comp(a,s) = count of atom a in species s")
        lines.append("    integer, dimension(nA,nS), intent(out) :: comp")
        lines.append("    comp(1:nA,1:nS) = comp_matrix(1:nA,1:nS)")
        lines.append("  end subroutine fcmech_get_composition")
        lines.append("")

    # fcmech_get_viscosity
    lines.append("  subroutine fcmech_get_viscosity(mu, T)")
    lines.append("    implicit none")
    lines.append("    ! mu(nS) : pure species viscosity (Pa·s)")
    lines.append("    ! T(K)   : temperature")
    lines.append("    real(WP), dimension(nS), intent(out) :: mu")
    lines.append("    real(WP), intent(in) :: T")
    lines.append("    integer :: i")
    lines.append("    do i = 1, nS")
    lines.append("      mu(i) = mucoeff(i)*sqrt(T)/fcmech_omegamu(T*koveps(i))")
    lines.append("    end do")
    lines.append("  end subroutine fcmech_get_viscosity")
    lines.append("")

    # fcmech_get_conductivity
    lines.append("  subroutine fcmech_get_conductivity(lambda, T, mu)")
    lines.append("    implicit none")
    lines.append("    ! lambda(nS) : thermal conductivity (W/(m·K))")
    lines.append("    ! T(K)      : temperature")
    lines.append("    ! mu(nS)    : viscosity from fcmech_get_viscosity")
    lines.append("    real(WP), dimension(nS), intent(out) :: lambda")
    lines.append("    real(WP), intent(in) :: T")
    lines.append("    real(WP), dimension(nS), intent(in) :: mu")
    lines.append("    lambda = mu*(Cpsp + 1.2_WP*Rcst/W_sp)")
    lines.append("  end subroutine fcmech_get_conductivity")
    lines.append("")

    # fcmech_get_invDij
    lines.append("  subroutine fcmech_get_invDij(invDij, T, P)")
    lines.append("    implicit none")
    lines.append("    ! invDij(nS,nS): inverse binary diffusion coefficients (s/m^2)")
    lines.append("    ! T(K), P(Pa) : temperature, pressure")
    lines.append("    real(WP), dimension(nS,nS), intent(out) :: invDij")
    lines.append("    real(WP), intent(in) :: T, P")
    lines.append("    real(WP) :: TPterm")
    lines.append("    integer :: i, j")
    lines.append("    TPterm = P/(T*sqrt(T))")
    lines.append("    do i = 1, nS")
    lines.append("      do j = 1, i-1")
    lines.append("        invDij(i,j) = TPterm*fcmech_omegaD(T*Ocoeffs(i,j))/Dcoeffs(i,j)")
    lines.append("        invDij(j,i) = invDij(i,j)")
    lines.append("      end do")
    lines.append("      invDij(i,i) = 0.0_WP")
    lines.append("    end do")
    lines.append("  end subroutine fcmech_get_invDij")
    lines.append("")

    # fcmech_get_ydot - main entry: mass-based production rates dY_i/dt
    lines.append("  subroutine fcmech_get_ydot(P, T, Y, ydot)")
    lines.append("    implicit none")
    lines.append("    ! P, T(Pa,K)    : pressure, temperature")
    lines.append("    ! Y(nS)         : mass fractions")
    lines.append("    ! ydot(nS)      : mass-based production rates dY_i/dt (1/s)")
    lines.append("    !   ydot_i = W_i * cdot_i / rho, with cdot_i = molar production rate (mol/(m^3·s))")
    lines.append("    ! c(nS)         : molar concentrations (mol/m^3)")
    lines.append("    ! k(nR), w(nR) : rate coefficients, reaction rates")
    lines.append("    ! M(n_m)        : effective third-body concentrations")
    lines.append("    real(WP), intent(in) :: P, T")
    lines.append("    real(WP), dimension(nS), intent(in) :: Y")
    lines.append("    real(WP), dimension(nS), intent(out) :: ydot")
    lines.append("    real(WP), dimension(nS) :: c, wdot")
    lines.append("    real(WP), dimension(nR) :: k, w")
    lines.append("    real(WP), dimension(" + str(n_m) + ") :: M")
    lines.append("    real(WP) :: W_mix, rho")
    lines.append("    call y2c(Y, W_sp, P, T, c)")
    lines.append("    call get_thirdbodies(M, c)")
    lines.append("    call get_rate_coefficients(k, M, T, P)")
    lines.append("    call get_reaction_rates(w, k, M, c)")
    lines.append("    call get_production_rates(wdot, w)")
    lines.append("    W_mix = 1.0_WP / sum(Y(1:nS) / W_sp(1:nS))")
    lines.append("    rho = P * W_mix / (Rcst * T)")
    lines.append("    ydot(1:nS) = W_sp(1:nS) * wdot(1:nS) / rho")
    lines.append("  end subroutine fcmech_get_ydot")
    lines.append("")

    lines.append("end module fcmech")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))

    if collect_viz:
        return viz_data
    return None


def run_viz(viz_data, output_prefix="kc_backward"):
    """Plot k_b from K_eq vs Arrhenius fit in 5x5 subplot grids."""
    if not HAS_MATPLOTLIB:
        print("Error: matplotlib required for --viz. Install with: pip install matplotlib", file=sys.stderr)
        sys.exit(1)
    n_reac = len(viz_data)
    if n_reac == 0:
        print("No reversible reactions to visualize.")
        return
    n_per_fig = 25  # 5x5
    n_figs = (n_reac + n_per_fig - 1) // n_per_fig
    for fig_idx in range(n_figs):
        start = fig_idx * n_per_fig
        end = min(start + n_per_fig, n_reac)
        n_sub = end - start
        fig, axes = plt.subplots(5, 5, figsize=(14, 14))
        axes = axes.flatten()
        for k in range(n_per_fig):
            ax = axes[k]
            if k < n_sub:
                idx = start + k
                eq, rid, T_arr, k_b_keq, k_b_fit = viz_data[idx]
                ax.semilogy(T_arr, k_b_keq, linewidth=2)
                ax.semilogy(T_arr, k_b_fit, linewidth=2)
                ax.set_xlabel("T (K)")
                ax.set_ylabel("k_b (mol/(m³·s) or 1/s)")
                ax.set_title(f"R{rid}: {eq[:40]}{'...' if len(eq) > 40 else ''}", fontsize=8)
                ax.grid(True, alpha=0.3)
            else:
                ax.axis("off")
        plt.tight_layout()
        outpath = f"{output_prefix}_fig{fig_idx + 1}.png"
        plt.savefig(outpath, dpi=150)
        plt.close()
        print(f"  Saved {outpath}")


def main():
    parser = argparse.ArgumentParser(description="Generate chem_data_fc.f90 from YAML mechanism (yaml2nga)")
    parser.add_argument("input", help="Input YAML mechanism file")
    parser.add_argument("output", nargs="?", default="chem_data_fc.f90", help="Output Fortran file (default: chem_data_fc.f90)")
    parser.add_argument("--viz", action="store_true",
                        help="Compare k_b from K_eq vs Arrhenius fit for each reversible reaction")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading mechanism from {input_path}")
    mech = load_mechanism(input_path)
    print(f"  Species: {len(mech['species_names'])}")
    print(f"  Reactions: {len(mech['reactions'])}")

    print(f"Writing chem_data_fc.f90 to {output_path}")
    viz_data = yaml2nga(mech, output_path, yaml_path=input_path, collect_viz=args.viz)
    if args.viz and viz_data:
        print(f"Generating backward rate comparison plots ({len(viz_data)} reversible reactions)...")
        run_viz(viz_data, output_prefix="kc_backward")
    print("Done.")


if __name__ == "__main__":
    main()
