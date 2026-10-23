#!/usr/bin/env python3
"""
Generate chem_data_fc.f90 (NGA2 finite-rate chemistry module ``fcmech``) from a Cantera YAML mechanism. Cantera is not
needed: the YAML is read following Cantera 3.0's rules and the transport fits are computed here, so the module holds
the same rates, thermo and transport as Cantera would use.

Usage: python yaml2nga_lib_cf.py <mechanism.yaml> [output.f90] [--phase NAME] [--no-transport] [--custom-f90 FILE]
                                 [--backend yaml|cantera] [--check]
       Default output: chem_data_fc.f90 (+ chem_data_reactions.txt, or <stem>_reactions.txt)

--backend yaml (default): read the YAML here (units: block or per-value units, Cantera's SI defaults when absent;
       reaction types, colliders and efficiencies, NASA7 thermo and atomic weights as Cantera 3.0 interprets them) and
       compute the mixture-averaged transport fits here (port of Cantera's GasTransport / MMCollisionInt).
       Needs numpy and PyYAML (or ruamel.yaml).
--backend cantera: take the mechanism and the transport fits from Cantera (>= 3.0) instead.
--check: read the mechanism with both backends (needs Cantera) and report the largest differences. Rates, thermo and
       species data agree exactly; transport properties to ~1e-12 (least-squares round-off).
--custom-f90 FILE: ARCANE reduced (QSS) mechanism. Species, thermo and transport come from the YAML phase (the
       transported species); fcmech_get_ydot calls ARCANE's `customkinetics(P, T, Y, wdot)` from FILE (wdot in
       mol/(m^3 s)), which must be compiled together with the generated module. The species order and molar masses
       of FILE are checked against the YAML. The rate routines of fcmech (get_rate_coefficients, ...) are empty.

Same module interface as yaml2nga.py (fcmech_get_ydot, fcmech_get_thermodata, fcmech_get_speciesnames,
fcmech_get_atomnames/atommasses/composition, fcmech_get_viscosity, fcmech_get_conductivity, fcmech_get_invDij,
y2c, get_thirdbodies, get_rate_coefficients, get_reaction_rates, get_production_rates; c in mol/m^3, h in J/mol,
cp in J/(mol K), rates in mol/(m^3 s)). Differences to yaml2nga.py:

- Units, reaction types, colliders and efficiencies exactly as Cantera interprets the YAML.
- Reverse rates are k_b = k_f/K_c with K_c from the NASA7 thermo at run time (as Cantera does); no Arrhenius fit.
- PLOG as in Cantera: entries at the same pressure are summed, a pressure on a node uses that node, log-log
  interpolation between nodes, nearest node outside the table.
- Supported: NASA7 thermo (one or two ranges); Arrhenius (negative A allowed), three-body (M or explicit collider),
  falloff Lindemann/Troe/SRI (also chemically activated, explicit collider), PLOG, Chebyshev, explicit reaction
  orders; ideal-gas phases with a species list, 'all' or sections of the same or another file. Anything else stops
  with an error; no reaction is skipped silently.
- Fortran species identifiers: 's' + name with ( ) + - and blanks removed (as in yaml2nga.py), '*' -> 'x', any
  other character that is not alphanumeric -> '_'; a clash (Fortran is case-insensitive) gets a suffix _2, _3, ...
- Transport: Cantera's mixture-averaged species fits (4th-order polynomials in ln T):
  mu_k = sqrt(T)*poly^2, lambda_k = sqrt(T)*poly, D_ij = T^1.5*poly/P. Pure-species values and binary
  diffusivities are those of Cantera's "mixture-averaged" model (mixing rules are left to the caller).
- Real constants are written with 17 significant digits.
- Reaction numbering: w(1:nRf) are the forward rates in YAML order, the reverse rates follow; the file
  chem_data_reactions.txt lists the mapping. fcmech_get_conductivity keeps its mu argument but does not use it.

Requires: numpy, PyYAML or ruamel.yaml (Cantera >= 3.0 only for --backend cantera and --check).
"""

import argparse
import math
import re
import sys
import warnings
from pathlib import Path

try:
    import numpy as np
except ImportError:
    print("Error: NumPy required. Install with: pip install numpy", file=sys.stderr)
    sys.exit(1)

# Physical constants as in Cantera 3.0 (ct_defs.h)
AVOGADRO = 6.02214076e26             # 1/kmol
BOLTZMANN = 1.380649e-23             # J/K
GAS_CONSTANT = AVOGADRO * BOLTZMANN  # J/(kmol K)
ONE_ATM = 1.01325e5                  # Pa
ELECTRON_CHARGE = 1.602176634e-19    # C
SMALL_NUMBER = 1.e-300
KMOL = 1.0e3                         # mol per kmol (Cantera works in kmol)
R_UNIV = GAS_CONSTANT / KMOL         # J/(mol K)

PER_LINE = 3                       # real values per line in data statements
PER_STMT = 60                      # values per data statement (keeps continuation lines short)
LINES_PER_STMT = 20                # continuation lines per production-rate statement


class UnsupportedMechanism(Exception):
    pass


# --------------------------------------------------------------------------------------------------
#  Fortran formatting helpers
# --------------------------------------------------------------------------------------------------
def fr(x):
    """Fortran real(WP) literal with full double precision."""
    x = float(x)
    if not np.isfinite(x):
        raise ValueError(f"non-finite constant {x}")
    s = repr(x)
    if "e" not in s and "." not in s:
        s += ".0"
    return s + "_WP"


def fi(x):
    return str(int(x))


def declare(lines, kind, name, shape):
    dims = ",".join(str(s) for s in shape)
    lines.append(f"  {kind}, dimension({dims}) :: {name}")


def data_1d(lines, name, values, fmt, per_line=PER_LINE, per_stmt=PER_STMT):
    """data name(i1:i2) / ... / statements in chunks (no continuation-limit issues for large arrays)."""
    values = list(values)
    for s in range(0, len(values), per_stmt):
        chunk = values[s:s + per_stmt]
        lines.append(f"  data {name}({s + 1}:{s + len(chunk)}) / &")
        for i in range(0, len(chunk), per_line):
            sub = ", ".join(fmt(v) for v in chunk[i:i + per_line])
            lines.append("       " + sub + (", &" if i + per_line < len(chunk) else " /"))


def data_cols(lines, name, arr, fmt, per_line=PER_LINE, per_stmt=PER_STMT):
    """data name(:,j1:j2) / ... / for a 2D array arr(nrow, ncol), column-major."""
    arr = np.asarray(arr)
    nrow, ncol = arr.shape
    cols_per_stmt = max(1, per_stmt // nrow)
    for j0 in range(0, ncol, cols_per_stmt):
        j1 = min(ncol, j0 + cols_per_stmt)
        flat = arr[:, j0:j1].flatten(order="F")
        lines.append(f"  data {name}(:,{j0 + 1}:{j1}) / &")
        for i in range(0, len(flat), per_line):
            sub = ", ".join(fmt(v) for v in flat[i:i + per_line])
            lines.append("       " + sub + (", &" if i + per_line < len(flat) else " /"))


def wrap_sum(lhs, terms, first=None, width=110):
    """lhs = first + t1 + t2 ... as continuation lines of <= width characters, LINES_PER_STMT lines per statement."""
    rows, row = [], ""
    for t in terms:
        if row and len(row) + 1 + len(t) > width:
            rows.append(row)
            row = t
        else:
            row = (row + " " + t) if row else t
    if row:
        rows.append(row)
    out = []
    if not rows:
        return [f"    {lhs} = {first if first is not None else '0.0_WP'}"]
    for s in range(0, len(rows), LINES_PER_STMT):
        chunk = rows[s:s + LINES_PER_STMT]
        head = (first if first is not None else "0.0_WP") if s == 0 else lhs
        out.append(f"    {lhs} = {head} &")
        for i, r in enumerate(chunk):
            out.append("      " + r + (" &" if i < len(chunk) - 1 else ""))
    return out


def comment_safe(s):
    return s.replace("!", "")


# --------------------------------------------------------------------------------------------------
#  Species identifiers
# --------------------------------------------------------------------------------------------------
_RESERVED = {"s", "sum", "sqrt", "sign", "size", "sin", "sinh", "shape", "spread", "scale", "spacing"}


def f90_species_ids(names):
    ids, seen = [], set()
    for nm in names:
        base = "s" + re.sub(r"[()+\- ]", "", nm).replace("*", "x")
        base = re.sub(r"[^A-Za-z0-9_]", "_", base)[:60]
        cand, n = base, 2
        while cand.lower() in seen or cand.lower() in _RESERVED:
            cand, n = f"{base}_{n}", n + 1
        seen.add(cand.lower())
        ids.append(cand)
    return ids



def lnA_sgn(A):
    if A == 0.0:
        return -700.0, 0.0
    return float(np.log(abs(A))), (1.0 if A > 0 else -1.0)


def arr_expr(A, b, E_R):
    """Fortran expression for A*T^b*exp(-E_R/T) using T_log, T_inv."""
    if A == 0.0:
        return "0.0_WP"
    lnA, sgn = lnA_sgn(A)
    e = f"exp({fr(lnA)} {'-' if b < 0 else '+'} {fr(abs(b))}*T_log {'+' if E_R < 0 else '-'} {fr(abs(E_R))}*T_inv)"
    return e if sgn > 0 else "(-" + e + ")"


# --------------------------------------------------------------------------------------------------
#  Units (port of Cantera 3.0 Units / UnitSystem: same factors, same arithmetic)
# --------------------------------------------------------------------------------------------------
class Units:
    """factor and dimensions (mass, length, time, temperature, current, quantity); kmol-based like Cantera."""

    def __init__(self, factor=1.0, mass=0.0, length=0.0, time=0.0, temperature=0.0, current=0.0, quantity=0.0):
        self.factor = factor
        self.mass, self.length, self.time = mass, length, time
        self.temperature, self.current, self.quantity = temperature, current, quantity
        self.pressure = self.energy = 0.0
        if (mass != 0 and length == -mass and time == -2 * mass and temperature == 0 and current == 0
                and quantity == 0):
            self.pressure = mass
        elif (mass != 0 and length == 2 * mass and time == -2 * mass and temperature == 0 and current == 0
              and quantity == 0):
            self.energy = mass

    def dims(self):
        return (self.mass, self.length, self.time, self.temperature, self.current, self.quantity)

    def convertible(self, other):
        return self.dims() == other.dims()

    def imul(self, o):
        self.factor *= o.factor
        self.mass += o.mass; self.length += o.length; self.time += o.time
        self.temperature += o.temperature; self.current += o.current; self.quantity += o.quantity
        self.pressure += o.pressure; self.energy += o.energy
        return self

    def pow(self, e):
        return Units(math.pow(self.factor, e), self.mass * e, self.length * e, self.time * e, self.temperature * e,
                     self.current * e, self.quantity * e)


_KNOWN_UNITS = {
    "": Units(1.0), "1": Units(1.0),
    "kg": Units(1.0, 1, 0, 0), "g": Units(1e-3, 1, 0, 0),
    "m": Units(1.0, 0, 1, 0), "micron": Units(1e-6, 0, 1, 0), "angstrom": Units(1e-10, 0, 1, 0),
    "Å": Units(1e-10, 0, 1, 0),
    "s": Units(1.0, 0, 0, 1), "min": Units(60, 0, 0, 1), "hr": Units(3600, 0, 0, 1),
    "K": Units(1.0, 0, 0, 0, 1), "C": Units(1.0, 0, 0, 0, 1),
    "A": Units(1.0, 0, 0, 0, 0, 1),
    "mol": Units(1e-3, 0, 0, 0, 0, 0, 1), "gmol": Units(1e-3, 0, 0, 0, 0, 0, 1),
    "mole": Units(1e-3, 0, 0, 0, 0, 0, 1), "kmol": Units(1.0, 0, 0, 0, 0, 0, 1),
    "kgmol": Units(1.0, 0, 0, 0, 0, 0, 1), "molec": Units(1.0 / AVOGADRO, 0, 0, 0, 0, 0, 1),
    "J": Units(1.0, 1, 2, -2), "cal": Units(4.184, 1, 2, -2), "erg": Units(1e-7, 1, 2, -2),
    "eV": Units(ELECTRON_CHARGE, 1, 2, -2),
    "N": Units(1.0, 1, 1, -2), "dyn": Units(1e-5, 1, 1, -2),
    "Pa": Units(1.0, 1, -1, -2), "atm": Units(ONE_ATM, 1, -1, -2), "bar": Units(1.0e5, 1, -1, -2),
    "dyn/cm^2": Units(0.1, 1, -1, -2),
    "m^3": Units(1.0, 0, 3, 0), "liter": Units(0.001, 0, 3, 0), "L": Units(0.001, 0, 3, 0),
    "l": Units(0.001, 0, 3, 0), "cc": Units(1.0e-6, 0, 3, 0),
    "ohm": Units(1.0, 1, 2, -3, 0, -2), "V": Units(1.0, 1, 2, -3, 0, -1), "coulomb": Units(1.0, 0, 0, 1, 0, 1),
    "J/kmol": Units(1.0, 1, 2, -2, 0, 0, -1),
}
_PREFIXES = {"Y": 1e24, "Z": 1e21, "E": 1e18, "P": 1e15, "T": 1e12, "G": 1e9, "M": 1e6, "k": 1e3, "h": 1e2,
             "d": 1e-1, "c": 1e-2, "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15, "a": 1e-18,
             "z": 1e-21, "y": 1e-24}
_FACTOR_RE = re.compile(r"[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)")


def parse_units(name):
    """Units(string) as Cantera: optional leading factor, then 'unit^exp' groups separated by * and /."""
    u, start = Units(1.0), 0
    m = _FACTOR_RE.search(name)
    if m and name.find(m.group(0)) == 0:
        u.factor = float(m.group(0))
        start = len(m.group(0))
    while True:
        stops = [i for i in (name.find("*", start), name.find("/", start)) if i >= 0]
        stop = min(stops) if stops else -1
        caret = name.find("^", start)
        if caret >= 0 and stop >= 0 and caret > stop:
            caret = -1
        end = len(name) if stop < 0 else stop
        unit = name[start:(caret if caret >= 0 else end)].strip()
        exponent = float(name[caret + 1:end]) if caret >= 0 else 1.0
        if start != 0 and name[start - 1] == "/":
            exponent = -exponent
        if unit in _KNOWN_UNITS:
            u.imul(_KNOWN_UNITS[unit].pow(exponent))
        elif unit[:1] in _PREFIXES and unit[1:] in _KNOWN_UNITS:
            v = _KNOWN_UNITS[unit[1:]]
            s = Units(v.factor, *v.dims())
            s.factor *= _PREFIXES[unit[:1]]
            u.imul(s.pow(exponent))
        else:
            raise UnsupportedMechanism(f"unknown unit '{unit}' in '{name}'")
        if stop < 0:
            break
        start = stop + 1
    return u


class UnitSystem:
    """Default units of a YAML file (Cantera UnitSystem): SI with kmol, unless a 'units' map says otherwise."""

    def __init__(self, parent=None):
        if parent is None:
            self.mass_f = self.length_f = self.time_f = self.quantity_f = 1.0
            self.pressure_f = self.energy_f = self.act_f = 1.0
            self.explicit_act = False
        else:
            self.__dict__.update(parent.__dict__)

    def set_defaults(self, units):
        for dim in sorted(units):                        # std::map order, as Cantera
            name = str(units[dim])
            u = parse_units(name)
            if dim == "mass" and u.convertible(_KNOWN_UNITS["kg"]):
                self.mass_f = u.factor
            elif dim == "length" and u.convertible(_KNOWN_UNITS["m"]):
                self.length_f = u.factor
            elif dim == "time" and u.convertible(_KNOWN_UNITS["s"]):
                self.time_f = u.factor
            elif dim == "temperature" and u.convertible(_KNOWN_UNITS["K"]):
                if u.factor != 1.0:
                    raise UnsupportedMechanism("temperature units other than K are not supported")
            elif dim == "current" and u.convertible(_KNOWN_UNITS["A"]):
                pass
            elif dim == "quantity" and u.convertible(_KNOWN_UNITS["kmol"]):
                self.quantity_f = u.factor
            elif dim == "pressure" and u.convertible(_KNOWN_UNITS["Pa"]):
                self.pressure_f = u.factor
            elif dim == "energy" and u.convertible(_KNOWN_UNITS["J"]):
                self.energy_f = u.factor
            elif dim == "activation-energy":
                pass
            else:
                raise UnsupportedMechanism(f"cannot set default unit for '{dim}' to '{name}'")
        if "activation-energy" in units:
            u = parse_units(str(units["activation-energy"]))
            if u.convertible(parse_units("J/kmol")):
                self.act_f = u.factor
            elif u.convertible(_KNOWN_UNITS["K"]):
                self.act_f = GAS_CONSTANT
            elif u.convertible(_KNOWN_UNITS["eV"]):
                self.act_f = u.factor * AVOGADRO
            else:
                raise UnsupportedMechanism(f"'{units['activation-energy']}' is not an activation-energy unit")
            self.explicit_act = True
        elif not self.explicit_act:
            self.act_f = self.energy_f / self.quantity_f

    def convert_to(self, value, dest):
        return (value / dest.factor
                * math.pow(self.mass_f, dest.mass - dest.pressure - dest.energy)
                * math.pow(self.length_f, dest.length + dest.pressure - 2 * dest.energy)
                * math.pow(self.time_f, dest.time + 2 * dest.pressure + 2 * dest.energy)
                * math.pow(self.quantity_f, dest.quantity)
                * math.pow(self.pressure_f, dest.pressure)
                * math.pow(self.energy_f, dest.energy))

    def convert(self, v, dest):
        """A YAML value (number, or 'number unit' string) in the units 'dest' (a Units object or a unit string)."""
        dest = parse_units(dest) if isinstance(dest, str) else dest
        if isinstance(v, str):
            val, _, src = v.strip().partition(" ")
            if not src:
                raise UnsupportedMechanism(f"cannot parse '{v}' as a value and a unit")
            src = parse_units(src)
            if not src.convertible(dest):
                raise UnsupportedMechanism(f"incompatible units in '{v}'")
            return float(val) * src.factor / dest.factor
        return self.convert_to(float(v), dest)

    def activation_energy_K(self, v):
        """Activation energy as E/R in K (Cantera convertActivationEnergy(v, 'K'))."""
        if isinstance(v, str):
            val, _, src = v.strip().partition(" ")
            value, u = float(val), parse_units(src)
            if u.convertible(parse_units("J/kmol")):
                value *= u.factor
            elif u.convertible(parse_units("K")):
                value *= GAS_CONSTANT * u.factor
            elif u.convertible(parse_units("eV")):
                value *= AVOGADRO * u.factor
            else:
                raise UnsupportedMechanism(f"'{v}' is not an activation energy")
            return value / (GAS_CONSTANT * 1.0)
        return float(v) * self.act_f / GAS_CONSTANT


def rate_units(e):
    """Units of a rate coefficient: (kmol/m^3)^e / s with e = 1 - reaction order (factor 1, as Cantera)."""
    return Units(1.0, 0, -3.0 * e, -1, 0, 0, e)


# --------------------------------------------------------------------------------------------------
#  Elements: Cantera's atomic-weight table (kg/kmol; CIAAW 2017, IUPAC 2018)
# --------------------------------------------------------------------------------------------------
_ATOMIC_WEIGHTS = {
    "H": (1.008, "hydrogen"), "He": (4.002602, "helium"), "Li": (6.94, "lithium"), "Be": (9.0121831, "beryllium"),
    "B": (10.81, "boron"), "C": (12.011, "carbon"), "N": (14.007, "nitrogen"), "O": (15.999, "oxygen"),
    "F": (18.998403163, "fluorine"), "Ne": (20.1797, "neon"), "Na": (22.98976928, "sodium"),
    "Mg": (24.305, "magnesium"), "Al": (26.9815384, "aluminum"), "Si": (28.085, "silicon"),
    "P": (30.973761998, "phosphorus"), "S": (32.06, "sulfur"), "Cl": (35.45, "chlorine"), "Ar": (39.95, "argon"),
    "K": (39.0983, "potassium"), "Ca": (40.078, "calcium"), "Sc": (44.955908, "scandium"),
    "Ti": (47.867, "titanium"), "V": (50.9415, "vanadium"), "Cr": (51.9961, "chromium"),
    "Mn": (54.938043, "manganese"), "Fe": (55.845, "iron"), "Co": (58.933194, "cobalt"), "Ni": (58.6934, "nickel"),
    "Cu": (63.546, "copper"), "Zn": (65.38, "zinc"), "Ga": (69.723, "gallium"), "Ge": (72.630, "germanium"),
    "As": (74.921595, "arsenic"), "Se": (78.971, "selenium"), "Br": (79.904, "bromine"), "Kr": (83.798, "krypton"),
    "Rb": (85.4678, "rubidium"), "Sr": (87.62, "strontium"), "Y": (88.90584, "yttrium"), "Zr": (91.224, "zirconium"),
    "Nb": (92.90637, "nobelium"), "Mo": (95.95, "molybdenum"), "Tc": (-1.0, "technetium"),
    "Ru": (101.07, "ruthenium"), "Rh": (102.90549, "rhodium"), "Pd": (106.42, "palladium"),
    "Ag": (107.8682, "silver"), "Cd": (112.414, "cadmium"), "In": (114.818, "indium"), "Sn": (118.710, "tin"),
    "Sb": (121.760, "antimony"), "I": (126.90447, "iodine"), "Xe": (131.293, "xenon"),
    "Cs": (132.90545196, "cesium"), "Ba": (137.327, "barium"), "La": (138.90547, "lanthanum"),
    "Ce": (140.116, "cerium"), "Pr": (140.90766, "praseodymium"), "Nd": (144.242, "neodymium"),
    "Pm": (-1.0, "promethium"), "Sm": (150.36, "samarium"), "Eu": (151.964, "europium"),
    "Gd": (157.25, "gadolinium"), "Tb": (158.925354, "terbium"), "Dy": (162.500, "dysprosium"),
    "Ho": (164.930328, "holmium"), "Er": (167.259, "erbium"), "Tm": (168.934218, "thulium"),
    "Yb": (173.045, "ytterbium"), "Lu": (174.9668, "lutetium"), "Hf": (178.49, "hafnium"),
    "Ta": (180.94788, "tantalum"), "W": (183.84, "tungsten"), "Re": (186.207, "rhenium"), "Ir": (192.217, "iridium"),
    "Pt": (195.084, "platinum"), "Au": (196.966570, "gold"), "Hg": (200.592, "mercury"), "Tl": (204.38, "thallium"),
    "Bi": (208.98040, "bismuth"), "Po": (-1.0, "polonium"), "At": (-1.0, "astatine"), "Rn": (-1.0, "radon"),
    "Fr": (-1.0, "francium"), "Ra": (-1.0, "radium"), "Ac": (-1.0, "actinium"), "Th": (232.0377, "thorium"),
    "Pa": (231.03588, "protactinium"), "U": (238.02891, "uranium"), "Np": (-1.0, "neptunium"),
    "Pu": (-1.0, "plutonium"), "Am": (-1.0, "americium"), "Cm": (-1.0, "curium"), "Bk": (-1.0, "berkelium"),
    "Cf": (-1.0, "californium"), "Es": (-1.0, "einsteinium"), "Fm": (-1.0, "fermium"), "Md": (-1.0, "mendelevium"),
    "No": (-1.0, "nobelium"), "Lr": (-1.0, "lawrencium"), "Rf": (-1.0, "rutherfordium"), "Db": (-1.0, "dubnium"),
    "Sg": (-1.0, "seaborgium"), "Bh": (-1.0, "bohrium"), "Hs": (-1.0, "hassium"), "Mt": (-1.0, "meitnerium"),
    "Ds": (-1.0, "darmstadtium"), "Rg": (-1.0, "roentgenium"), "Cn": (-1.0, "copernicium"), "Nh": (-1.0, "nihonium"),
    "Gl": (-1.0, "flerovium"), "Mc": (-1.0, "moscovium"), "Lv": (-1.0, "livermorium"), "Ts": (-1.0, "tennessine"),
    "Og": (-1.0, "oganesson"), "D": (2.0141017781, "deuterium"), "Tr": (3.0160492820, "tritium"),
}


def element_weight(sym, custom):
    if sym in custom:
        return custom[sym]
    if sym.strip() in _ATOMIC_WEIGHTS:
        w = _ATOMIC_WEIGHTS[sym.strip()][0]
    else:
        w = next((v[0] for v in _ATOMIC_WEIGHTS.values() if v[1] == sym.strip().lower()), 0.0)
    if w > 0.0:
        return w
    if w < 0.0:
        raise UnsupportedMechanism(f"element '{sym}' has no stable isotopes")
    raise UnsupportedMechanism(f"element not found: {sym}")


# --------------------------------------------------------------------------------------------------
#  YAML input (YAML 1.2 scalars, as Cantera's yaml-cpp: NO/Yes/On stay strings, 1e5 is a float)
# --------------------------------------------------------------------------------------------------
_YAML_CACHE = {}


def load_yaml(path):
    path = Path(path).resolve()
    if path in _YAML_CACHE:
        return _YAML_CACHE[path]
    text = path.read_text()
    try:
        import yaml as _yaml

        class _Loader(_yaml.SafeLoader):
            pass

        _Loader.yaml_implicit_resolvers = {
            k: [(tag, rx) for tag, rx in v if tag not in ("tag:yaml.org,2002:bool", "tag:yaml.org,2002:float")]
            for k, v in _yaml.SafeLoader.yaml_implicit_resolvers.items()}
        _Loader.add_implicit_resolver("tag:yaml.org,2002:bool", re.compile(r"^(?:true|True|TRUE|false|False|FALSE)$"),
                                      list("tTfF"))
        _Loader.add_implicit_resolver(
            "tag:yaml.org,2002:float",
            re.compile(r"^(?:[-+]?(?:[0-9][0-9_]*)\.[0-9_]*(?:[eE][-+]?[0-9]+)?|[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)"
                       r"|[-+]?\.[0-9_]+(?:[eE][-+]?[0-9]+)?|[-+]?\.(?:inf|Inf|INF)|\.(?:nan|NaN|NAN))$"),
            list("-+0123456789."))
        data = _yaml.load(text, Loader=_Loader)
    except ImportError:
        from ruamel.yaml import YAML
        data = YAML(typ="safe").load(text)
    _YAML_CACHE[path] = data
    return data


def _section(ref, data, base):
    """'name' (section of this file) or 'file.yaml/name' -> (list of entries, file data, file dir)."""
    if "/" in ref:
        fname, sec = ref.rsplit("/", 1)
        f = (base / fname)
        d = load_yaml(f)
        return d.get(sec, []), d, f.parent
    return data.get(ref, []), data, base


# --------------------------------------------------------------------------------------------------
#  Reactions (port of Cantera 3.0 parseReactionEquation / Reaction::setEquation / ThirdBody)
# --------------------------------------------------------------------------------------------------
def parse_equation(eq):
    tokens = eq.split() + ["+"]
    reac, prod, reversible = {}, {}, None
    last, reactants = None, True
    for i in range(1, len(tokens)):
        t = tokens[i]
        if t == "+" or t.startswith("(+") or t in ("<=>", "=", "=>"):
            sp, stoich, mass_action = tokens[i - 1], 1.0, True
            if last is not None and tokens[last] == "(+" and sp.endswith(")"):
                mass_action = False
                sp = "(+" + sp
            elif last == i - 1 and sp.startswith("(+") and sp.endswith(")"):
                mass_action = False
            elif last == (i - 2 if i >= 2 else None) or (last is None and i == 1):
                pass
            elif last == i - 3 or (last is None and i == 2):
                stoich = float(tokens[i - 2])
            else:
                raise UnsupportedMechanism(f"cannot parse reaction equation '{eq}'")
            side = reac if reactants else prod
            side[sp] = side.get(sp, 0.0) + stoich
            last = i
        if t in ("<=>", "="):
            reversible, reactants = True, False
        elif t == "=>":
            reversible, reactants = False, False
    if reversible is None:
        raise UnsupportedMechanism(f"no '<=>', '=' or '=>' in reaction equation '{eq}'")
    return reac, prod, reversible


def _third_body(node, reac, prod, rate_type):
    """Cantera Reaction::setEquation: detect the collider and remove it from reactants and products.
    Returns None or (name, mass_action, default_efficiency, efficiencies)."""
    if rate_type == "elementary":
        return None
    tb, count, countM = None, 0, 0
    for name, st in sorted(reac.items()):              # std::map order, as Cantera
        if name in prod:
            tb = name
            generic = tb in ("M", "(+M)", "(+ M)")
            count += 1; countM += generic
            if st > 1 and prod[tb] > 1:
                count += 1; countM += generic
    explicit_effs = None
    if count == 0:
        if rate_type.startswith("three-body"):
            raise UnsupportedMechanism(f"reaction '{node['equation']}': no valid third-body collider")
        return None
    if countM > 1:
        raise UnsupportedMechanism(f"reaction '{node['equation']}': several generic colliders M")
    if count > 1:
        if countM:
            tb = "M"
        elif "efficiencies" in node:
            effs = {str(k): v for k, v in node["efficiencies"].items()}
            if len(effs) != 1 or next(iter(effs)) not in reac:
                raise UnsupportedMechanism(f"reaction '{node['equation']}': ambiguous third-body colliders")
            tb = str(next(iter(effs)))
            explicit_effs = True
        elif "default-efficiency" in node or rate_type.startswith("three-body"):
            raise UnsupportedMechanism(f"reaction '{node['equation']}': ambiguous third-body colliders")
        else:
            return None
    elif tb != "M" and not rate_type.startswith("three-body") and not tb.startswith("(+"):
        if any(math.trunc(v) != v for v in reac.values()) or any(math.trunc(v) != v for v in prod.values()):
            return None
        if int(sum(reac.values())) != 3 and int(sum(prod.values())) != 3:
            return None
    for side in (reac, prod):
        if math.trunc(side[tb]) != 1:
            side[tb] -= 1.0
        else:
            del side[tb]
    name, mass_action = tb, True
    if tb.startswith("(+ "):
        name, mass_action = tb[3:-1], False
    elif tb.startswith("(+"):
        name, mass_action = tb[2:-1], False
    default, effs = (1.0, {}) if name == "M" else (0.0, {name: 1.0})
    if "default-efficiency" in node:
        value = float(node["default-efficiency"])
        if name != "M" and value != 0.0:
            raise UnsupportedMechanism(f"reaction '{node['equation']}': default efficiency of explicit collider")
        default = value
    if "efficiencies" in node:
        effs = {str(k): float(v) for k, v in node["efficiencies"].items()}
    if name != "M" and (len(effs) != 1 or next(iter(effs)) != name):
        raise UnsupportedMechanism(f"reaction '{node['equation']}': incompatible third-body definitions")
    return name, mass_action, default, effs


def _arrhenius(rc, units, e):
    """Cantera ArrheniusBase::setRateParameters -> (A in kmol units, b, Ea in J/kmol as Cantera reports it)."""
    if isinstance(rc, dict):
        A, b, Ea = rc["A"], rc["b"], rc.get("Ea", 0.0)
    else:
        A, b, Ea = rc[0], rc[1], (rc[2] if len(rc) > 2 else 0.0)
    A = units.convert(A, rate_units(e))
    Ea_R = units.activation_energy_K(Ea) if Ea != 0.0 else 0.0
    return dict(A=A, b=float(b), Ea=Ea_R * GAS_CONSTANT)


def _reaction(node, units, sidx, skip_undeclared_tb):
    """One YAML reaction -> dict (Cantera conventions: kmol-based A, Ea in J/kmol, compositions in name order).
    Returns None when a species is not declared (handled by the caller)."""
    eq = node["equation"]
    reac, prod, reversible = parse_equation(eq)
    rate_type = str(node.get("type", "Arrhenius"))
    tbody = _third_body(node, reac, prod, rate_type)
    for s in list(reac) + list(prod):
        if s not in sidx:
            return None
    orders = {str(k): float(v) for k, v in node.get("orders", {}).items()}
    reac = dict(sorted(reac.items())); prod = dict(sorted(prod.items())); orders = dict(sorted(orders.items()))
    effs_default, effs = None, None
    if tbody:
        name, mass_action, effs_default, effs = tbody
        undeclared = [s for s in effs if s not in sidx]
        if undeclared:
            if not skip_undeclared_tb:
                raise UnsupportedMechanism(f"reaction '{eq}': third-body efficiencies for undeclared species "
                                           f"{undeclared}")
            if name != "M":
                return "skip"
            effs = {s: v for s, v in effs.items() if s in sidx}
    # exponent of the concentration units in the rate coefficient (Cantera calculateRateCoeffUnits)
    e = 1.0
    for s, o in orders.items():
        e -= o
    for s, v in reac.items():
        if s not in orders:
            e -= v
    if tbody and tbody[1]:
        e -= 1.0
    r = dict(equation=eq, reactants=reac, products=prod, orders=orders, reversible=reversible,
             third_body=None if not tbody else dict(name=tbody[0], mass_action=tbody[1], default=effs_default,
                                                     efficiencies=dict(sorted(effs.items()))))
    if rate_type in ("Arrhenius", "elementary", "three-body"):
        if rate_type == "three-body" and not (tbody and tbody[1]):
            raise UnsupportedMechanism(f"reaction '{eq}': three-body reaction without collider")
        r["rate"] = dict(type="Arrhenius", **_arrhenius(node["rate-constant"], units, e))
    elif rate_type in ("falloff", "chemically-activated"):
        if not tbody or tbody[1]:
            raise UnsupportedMechanism(f"reaction '{eq}': falloff without (+M) collider")
        chem = rate_type == "chemically-activated"
        low = _arrhenius(node["low-P-rate-constant"], units, e if chem else e - 1.0)
        high = _arrhenius(node["high-P-rate-constant"], units, e + 1.0 if chem else e)
        if "Troe" in node:
            f = node["Troe"]
            c = [float(f["A"]), float(f["T3"]), float(f["T1"])] + ([float(f["T2"])] if "T2" in f else [])
            rt3 = math.inf if abs(c[1]) < SMALL_NUMBER else 1.0 / c[1]
            rt1 = math.inf if abs(c[2]) < SMALL_NUMBER else 1.0 / c[2]
            t2 = c[3] if len(c) > 3 else 0.0
            coeffs = [c[0], 1.0 / rt3, 1.0 / rt1] + ([t2] if abs(t2) >= SMALL_NUMBER else [])
            ftype = "troe"
        elif "SRI" in node:
            f = node["SRI"]
            coeffs = [float(f["A"]), float(f["B"]), float(f["C"])]
            coeffs += [float(f["D"]), float(f["E"])] if "D" in f or "E" in f else [1.0, 0.0]
            if "D" in f and "E" not in f or "E" in f and "D" not in f:
                raise UnsupportedMechanism(f"reaction '{eq}': SRI needs both D and E")
            ftype = "sri"
        elif "Tsang" in node:
            raise UnsupportedMechanism(f"reaction '{eq}': Tsang falloff not supported")
        else:
            coeffs, ftype = [], "lindemann"
        r["rate"] = dict(type="falloff", ftype=ftype, chem=chem, low=low, high=high, coeffs=coeffs)
    elif rate_type == "pressure-dependent-Arrhenius":
        if tbody:
            raise UnsupportedMechanism(f"reaction '{eq}': PLOG with third body")
        rates = []
        for rc in node["rate-constants"]:
            rates.append((math.exp(math.log(units.convert(rc["P"], "Pa"))), _arrhenius(rc, units, e)))  # as Cantera
        r["rate"] = dict(type="plog", rates=sorted(rates, key=lambda x: x[0]))   # stable: multimap order
    elif rate_type == "Chebyshev":
        if tbody:
            raise UnsupportedMechanism(f"reaction '{eq}': Chebyshev with third body")
        data = np.array(node["data"], dtype=float)
        data[0, 0] += math.log10(units.convert(1.0, rate_units(e)))
        Tr, Pr = node["temperature-range"], node["pressure-range"]
        r["rate"] = dict(type="chebyshev", Tmin=units.convert(Tr[0], "K"), Tmax=units.convert(Tr[1], "K"),
                         Pmin=units.convert(Pr[0], "Pa"), Pmax=units.convert(Pr[1], "Pa"), data=data)
    else:
        raise UnsupportedMechanism(f"reaction '{eq}': rate type {rate_type} not supported")
    if r["rate"]["type"] == "Arrhenius" and not node.get("negative-A", False) and r["rate"]["A"] < 0:
        raise UnsupportedMechanism(f"reaction '{eq}': negative A needs 'negative-A: true'")
    return r


# --------------------------------------------------------------------------------------------------
#  Mechanism readers -> mech dict (species, elements, thermo, reactions in Cantera conventions, transport fits)
# --------------------------------------------------------------------------------------------------
def read_yaml(yaml_path, phase=None, transport=True, custom=False):
    """Cantera-free reader (Cantera 3.0 rules)."""
    yaml_path = Path(yaml_path).resolve()
    data, base = load_yaml(yaml_path), yaml_path.parent
    top_units = UnitSystem()
    if "units" in data:
        top_units.set_defaults(data["units"])
    phases = data.get("phases", [])
    if not phases:
        raise UnsupportedMechanism("no phases in the YAML file")
    ph = phases[0] if phase is None else next((p for p in phases if p.get("name") == phase), None)
    if ph is None:
        raise UnsupportedMechanism(f"phase '{phase}' not found")
    if ph.get("thermo") not in ("ideal-gas", None):
        raise UnsupportedMechanism(f"phase thermo '{ph.get('thermo')}' not supported (ideal-gas only)")

    # species (in phase order)
    spec, species = ph.get("species", "all"), []
    if spec == "all":
        species = list(data.get("species", []))
    elif isinstance(spec, list) and all(isinstance(s, str) for s in spec):
        allsp = {str(s["name"]): s for s in data.get("species", [])}
        missing = [s for s in spec if s not in allsp]
        if missing:
            raise UnsupportedMechanism(f"species not defined: {missing}")
        species = [allsp[s] for s in spec]
    elif isinstance(spec, list):
        for item in spec:
            for ref, names in item.items():
                sec, _, _ = _section(ref, data, base)
                if names == "all":
                    species += list(sec)
                else:
                    allsp = {str(s["name"]): s for s in sec}
                    species += [allsp[n] for n in names]
    else:
        raise UnsupportedMechanism(f"phase species '{spec}' not supported")
    names = [str(s["name"]) for s in species]
    if len(set(names)) != len(names):
        raise UnsupportedMechanism("duplicate species names")

    # elements
    custom_w = {str(e["symbol"]): float(e["atomic-weight"]) for e in data.get("elements", [])
                if isinstance(e, dict) and "symbol" in e}
    els = ph.get("elements")
    if els is None:                                       # added as met (composition keys in name order)
        elements = []
        for s in species:
            for el in sorted(str(k) for k in s["composition"]):
                if el not in elements:
                    elements.append(el)
    elif all(isinstance(e, str) for e in els):
        elements = list(els)
    else:
        elements = []
        for item in els:
            for ref, syms in item.items():
                sec, _, _ = _section(ref, data, base)
                custom_w.update({str(e["symbol"]): float(e["atomic-weight"]) for e in sec})
                elements += list(syms)
    aw = [element_weight(e, custom_w) for e in elements]
    comp = np.zeros((len(elements), len(species)))
    W, W_kmol = np.zeros(len(species)), np.zeros(len(species))
    for k, s in enumerate(species):
        if s.get("charge", 0) or "E" in s["composition"]:
            raise UnsupportedMechanism(f"species {names[k]}: charged species not supported")
        c = [0.0] * len(elements)
        for el, n in s["composition"].items():
            if str(el) not in elements:
                raise UnsupportedMechanism(f"species {names[k]} contains undeclared element {el}")
            c[elements.index(str(el))] = float(n)
        wt = 0.0
        for m in range(len(elements)):
            wt += c[m] * aw[m]
        W_kmol[k] = max(wt, 1e-20)
        W[k] = W_kmol[k] / KMOL
        comp[:, k] = c

    # thermo
    T_mid, th, Pref = np.zeros(len(species)), np.zeros((14, len(species))), None
    Tlim = np.zeros((2, len(species)))
    for k, s in enumerate(species):
        t = s.get("thermo", {})
        tr = [top_units.convert(v, "K") for v in t.get("temperature-ranges", [])]
        if t.get("model") != "NASA7" or len(tr) not in (2, 3) or len(t["data"]) != len(tr) - 1:
            raise UnsupportedMechanism(f"species {names[k]}: thermo {t.get('model')} (only NASA7 supported)")
        T_mid[k], Tlim[:, k] = tr[1], (tr[0], tr[-1])        # one range: repeated on both sides (as Cantera)
        th[0:7, k] = [float(v) for v in t["data"][0]]
        th[7:14, k] = [float(v) for v in t["data"][-1]]
        p = top_units.convert(t["reference-pressure"], "Pa") if "reference-pressure" in t else ONE_ATM
        if Pref is not None and p != Pref:
            raise UnsupportedMechanism("species with different thermo reference pressures")
        Pref = p

    mech = dict(source=yaml_path, species=names, W=W, W_kmol=W_kmol, elements=elements,
                atomic_weights=np.array(aw) / KMOL, comp=np.rint(comp).astype(int), T_mid=T_mid, thermo=th,
                Pref=Pref if Pref is not None else ONE_ATM, Tlim=Tlim, reactions=[], transport=None)

    # reactions
    sidx = {s: i for i, s in enumerate(names)}
    if not custom:
        kin = ph.get("kinetics")
        rspec = ph.get("reactions", "all" if kin else "none")
        if kin not in (None, "gas", "bulk"):
            raise UnsupportedMechanism(f"phase kinetics '{kin}' not supported (gas only; use --custom-f90 for "
                                       "ARCANE custom kinetics)")
        if isinstance(rspec, str):
            rspec = [{"reactions": rspec}] if rspec != "none" else []
        elif all(isinstance(x, str) for x in rspec):
            rspec = [{x: "all"} for x in rspec]
        skip_tb = bool(ph.get("skip-undeclared-third-bodies", False))
        for item in rspec:
            for ref, rule in item.items():
                if rule == "none":
                    continue
                sec, fdata, _ = _section(ref, data, base)
                u = UnitSystem()
                if "units" in fdata:
                    u.set_defaults(fdata["units"])
                for node in sec:
                    ru = u
                    if "units" in node:
                        ru = UnitSystem(u); ru.set_defaults(node["units"])
                    rx = _reaction(node, ru, sidx, skip_tb)
                    if rx is None:
                        if rule == "declared-species":
                            continue
                        raise UnsupportedMechanism(f"reaction '{node['equation']}' has undeclared species")
                    if rx == "skip":
                        continue
                    mech["reactions"].append(rx)
    if transport:
        try:
            mech["transport"] = transport_fits(species, mech)
        except UnsupportedMechanism as e:
            print(f"WARNING: no mixture-averaged transport ({e}); transport routines not generated", file=sys.stderr)
    return mech


def read_cantera(yaml_path, phase=None, transport=True, custom=False):
    """Reference reader through Cantera (--backend cantera, --check)."""
    import cantera as ct
    warnings.filterwarnings("ignore", category=UserWarning)
    ct.suppress_thermo_warnings()
    yaml_path = Path(yaml_path).resolve()
    if custom:
        data = load_yaml(yaml_path)
        ph = data["phases"][0] if phase is None else next(p for p in data["phases"] if p["name"] == phase)
        names = ph["species"]
        if not (isinstance(names, list) and all(isinstance(n, str) for n in names)):
            raise UnsupportedMechanism("custom mode needs an explicit species list in the YAML phase")
        allsp = {s.name: s for s in ct.Species.list_from_file(str(yaml_path))}
        sps = [allsp[n] for n in names]
        gas = ct.Solution(thermo="ideal-gas", species=sps, kinetics="gas", reactions=[])
        mk_tr = lambda: ct.Solution(thermo="ideal-gas", species=sps, kinetics="gas", reactions=[],
                                    transport_model="mixture-averaged")
    else:
        gas = ct.Solution(str(yaml_path), phase or "")
        mk_tr = lambda: ct.Solution(str(yaml_path), phase or "", transport_model="mixture-averaged")
    nS = gas.n_species
    T_mid, th = np.zeros(nS), np.zeros((14, nS))
    for k in range(nS):
        t = gas.species(k).thermo
        if type(t).__name__ != "NasaPoly2":
            raise UnsupportedMechanism(f"species {gas.species_names[k]}: thermo {type(t).__name__} "
                                       "(only NASA7 supported)")
        c = t.coeffs
        T_mid[k] = c[0]; th[0:7, k] = c[8:15]; th[7:14, k] = c[1:8]
    nA = gas.n_elements
    mech = dict(source=yaml_path, species=list(gas.species_names), W=gas.molecular_weights / KMOL,
                elements=[gas.element_name(m) for m in range(nA)], atomic_weights=gas.atomic_weights / KMOL,
                comp=np.array([[int(round(gas.n_atoms(k, m))) for k in range(nS)] for m in range(nA)], dtype=int),
                T_mid=T_mid, thermo=th, Pref=gas.reference_pressure, reactions=[], transport=None)
    sidx = {s: i for i, s in enumerate(gas.species_names)}
    arr = lambda a: dict(A=a.pre_exponential_factor, b=a.temperature_exponent, Ea=a.activation_energy)
    for ir, r in enumerate(gas.reactions() if not custom else []):
        rate = r.rate
        tb = None
        if r.third_body is not None:
            tb = dict(name=r.third_body.name, mass_action=r.third_body.mass_action,
                      default=float(r.third_body.default_efficiency),
                      efficiencies=dict(sorted((s, float(e)) for s, e in r.third_body.efficiencies.items()
                                               if s in sidx)))
        x = dict(equation=r.equation, reactants=dict(r.reactants), products=dict(r.products), orders=dict(r.orders),
                 reversible=r.reversible, third_body=tb)
        if isinstance(rate, ct.ArrheniusRate):
            x["rate"] = dict(type="Arrhenius", **arr(rate))
        elif isinstance(rate, (ct.LindemannRate, ct.TroeRate, ct.SriRate)):
            ftype = {"LindemannRate": "lindemann", "TroeRate": "troe", "SriRate": "sri"}[type(rate).__name__]
            x["rate"] = dict(type="falloff", ftype=ftype, chem=bool(rate.chemically_activated), low=arr(rate.low_rate),
                             high=arr(rate.high_rate), coeffs=[float(v) for v in rate.falloff_coeffs])
        elif isinstance(rate, ct.PlogRate):
            x["rate"] = dict(type="plog", rates=[(float(P), arr(a)) for P, a in rate.rates])
        elif isinstance(rate, ct.ChebyshevRate):
            (Tmin, Tmax), (Pmin, Pmax) = rate.temperature_range, rate.pressure_range
            x["rate"] = dict(type="chebyshev", Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax,
                             data=np.array(rate.data, dtype=float).copy())
        else:
            raise UnsupportedMechanism(f"reaction {ir + 1} ({r.equation}): rate type {type(rate).__name__} "
                                       "not supported")
        mech["reactions"].append(x)
    if transport:
        try:
            gtr = mk_tr()
        except ct.CanteraError as e:
            print(f"WARNING: no mixture-averaged transport ({str(e).strip().splitlines()[-1]}); "
                  "transport routines not generated", file=sys.stderr)
            gtr = None
        if gtr is not None:
            pairs = [(i, j) for j in range(1, nS) for i in range(j)]
            mech["transport"] = dict(
                visc=np.array([gtr.get_viscosity_polynomial(k) for k in range(nS)]).T,
                cond=np.array([gtr.get_thermal_conductivity_polynomial(k) for k in range(nS)]).T,
                diff=(np.array([gtr.get_binary_diff_coeffs_polynomial(i, j) for i, j in pairs]).T if pairs
                      else np.zeros((5, 0))))
    return mech


# --------------------------------------------------------------------------------------------------
#  Mechanism -> internal tables (all in mol, m, J, K)
# --------------------------------------------------------------------------------------------------
def arr_si(a, order):
    """Arrhenius (Cantera conventions, kmol units) -> (A in mol units, b, E/R in K)."""
    return a["A"] * KMOL ** (1.0 - order), a["b"], a["Ea"] / GAS_CONSTANT


def build_tables(mech):
    sp = mech["species"]
    sidx = {s: i for i, s in enumerate(sp)}
    nRf = len(mech["reactions"])

    arr, tb, fo, plog, cheb = [], [], [], [], []
    mgroups, mkey_to_idx = [], {}
    rate_conc, net_nu, rev_list, kinds = [], [], [], []

    for ir, r in enumerate(mech["reactions"]):
        rate, eq = r["rate"], r["equation"]
        conc = dict(r["reactants"])
        conc.update(r["orders"])                      # explicit orders replace stoichiometric ones
        order = float(sum(conc.values()))
        rate_conc.append(conc)
        nu = {}
        for s, v in r["products"].items():
            nu[s] = nu.get(s, 0.0) + v
        for s, v in r["reactants"].items():
            nu[s] = nu.get(s, 0.0) - v
        net_nu.append({s: v for s, v in nu.items() if v != 0.0})

        im, t3 = 0, r["third_body"]
        if t3 is not None:
            effs = {s: float(e) for s, e in t3["efficiencies"].items() if s in sidx}
            key = (float(t3["default"]), tuple(sorted(effs.items())))
            if key not in mkey_to_idx:
                mgroups.append(key)
                mkey_to_idx[key] = len(mgroups)
            im = mkey_to_idx[key]

        if rate["type"] == "Arrhenius":
            if t3 is None:
                arr.append((ir + 1,) + arr_si(rate, order))
                kinds.append("Arrhenius")
            else:
                tb.append((ir + 1,) + arr_si(rate, order + 1.0) + (im,))
                kinds.append("three-body")
        elif rate["type"] == "falloff":
            if t3 is None:
                raise UnsupportedMechanism(f"reaction {ir + 1} ({eq}): falloff without third body")
            chem = bool(rate["chem"])                    # Pr = k0*[M]/kinf dimensionless:
            low = arr_si(rate["low"], order if chem else order + 1.0)      # falloff k0: n+1, kinf: n
            high = arr_si(rate["high"], order - 1.0 if chem else order)    # chem.-act. k0: n, kinf: n-1
            if rate["ftype"] == "lindemann":
                ftype, pars = "lindemann", []
            elif rate["ftype"] == "troe":
                c = list(rate["coeffs"])
                a, T3, T1 = c[0], c[1], c[2]
                T2 = c[3] if len(c) > 3 else 0.0
                rT3 = 1.0 / T3 if abs(T3) > 1e-300 else 1e300
                rT1 = 1.0 / T1 if abs(T1) > 1e-300 else 1e300
                ftype, pars = "troe", [a, rT3, rT1, T2]
            else:
                c = list(rate["coeffs"]) + [1.0, 0.0][len(rate["coeffs"]) - 3:]
                ftype, pars = "sri", c[:5]
            fo.append((ir + 1, ftype, chem, low, high, pars, im))
            kinds.append(("chem.-activated " if chem else "falloff ") + ftype)
        elif rate["type"] == "plog":
            if t3 is not None:
                raise UnsupportedMechanism(f"reaction {ir + 1} ({eq}): PLOG with third body")
            nodes = {}
            for P, a in rate["rates"]:
                nodes.setdefault(float(P), []).append(arr_si(a, order))
            plog.append((ir + 1, sorted(nodes.items())))
            kinds.append("plog")
        elif rate["type"] == "chebyshev":
            if t3 is not None:
                raise UnsupportedMechanism(f"reaction {ir + 1} ({eq}): Chebyshev with third body")
            data = np.array(rate["data"], dtype=float).copy()
            data[0, 0] += 3.0 * (1.0 - order)          # log10(kmol -> mol) unit shift
            cheb.append((ir + 1, rate["Tmin"], rate["Tmax"], rate["Pmin"], rate["Pmax"], data))
            kinds.append("chebyshev")
        else:
            raise UnsupportedMechanism(f"reaction {ir + 1} ({eq}): rate type {rate['type']} not supported")

        if r["reversible"]:
            if r["orders"]:
                raise UnsupportedMechanism(f"reaction {ir + 1} ({eq}): reversible with explicit orders")
            rev_list.append(ir + 1)

    return dict(nRf=nRf, arr=arr, tb=tb, fo=fo, plog=plog, cheb=cheb, mgroups=mgroups,
                rate_conc=rate_conc, net_nu=net_nu, rev_list=rev_list, kinds=kinds)


# --------------------------------------------------------------------------------------------------
#  Transport fits: port of Cantera 3.0 GasTransport (mixture-averaged species fits) and MMCollisionInt
#  (Monchick and Mason collision integrals). Same formulas and order of operations; least-squares fits
#  by Householder QR (Cantera uses Eigen's column-pivoting QR: coefficients agree to round-off).
# --------------------------------------------------------------------------------------------------
PI = 3.14159265358979323846
LIGHT_SPEED = 299792458.0
PLANCK = 6.62607015e-34
FINE_STRUCTURE = 7.2973525693e-3
PERMEABILITY_0 = 2 * FINE_STRUCTURE * PLANCK / (ELECTRON_CHARGE * ELECTRON_CHARGE * LIGHT_SPEED)
EPSILON_0 = 1.0 / (LIGHT_SPEED * LIGHT_SPEED * PERMEABILITY_0)

_MM_DELTA = np.array([
    0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5,
])
_MM_TSTAR22 = np.array([
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0,
    8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0,
])
_MM_OMEGA22_TABLE = np.array([
    4.1005, 4.266, 4.833, 5.742, 6.729, 8.624, 10.34, 11.89, 3.2626, 3.305, 3.516, 3.914, 4.433, 5.57, 6.637,
    7.618, 2.8399, 2.836, 2.936, 3.168, 3.511, 4.329, 5.126, 5.874, 2.531, 2.522, 2.586, 2.749, 3.004, 3.64, 4.282,
    4.895, 2.2837, 2.277, 2.329, 2.46, 2.665, 3.187, 3.727, 4.249, 2.0838, 2.081, 2.13, 2.243, 2.417, 2.862, 3.329,
    3.786, 1.922, 1.924, 1.97, 2.072, 2.225, 2.614, 3.028, 3.435, 1.7902, 1.795, 1.84, 1.934, 2.07, 2.417, 2.788,
    3.156, 1.6823, 1.689, 1.733, 1.82, 1.944, 2.258, 2.596, 2.933, 1.5929, 1.601, 1.644, 1.725, 1.838, 2.124,
    2.435, 2.746, 1.4551, 1.465, 1.504, 1.574, 1.67, 1.913, 2.181, 2.451, 1.3551, 1.365, 1.4, 1.461, 1.544, 1.754,
    1.989, 2.228, 1.28, 1.289, 1.321, 1.374, 1.447, 1.63, 1.838, 2.053, 1.2219, 1.231, 1.259, 1.306, 1.37, 1.532,
    1.718, 1.912, 1.1757, 1.184, 1.209, 1.251, 1.307, 1.451, 1.618, 1.795, 1.0933, 1.1, 1.119, 1.15, 1.193, 1.304,
    1.435, 1.578, 1.0388, 1.044, 1.059, 1.083, 1.117, 1.204, 1.31, 1.428, 0.99963, 1.004, 1.016, 1.035, 1.062,
    1.133, 1.22, 1.319, 0.96988, 0.9732, 0.983, 0.9991, 1.021, 1.079, 1.153, 1.236, 0.92676, 0.9291, 0.936, 0.9473,
    0.9628, 1.005, 1.058, 1.121, 0.89616, 0.8979, 0.903, 0.9114, 0.923, 0.9545, 0.9955, 1.044, 0.87272, 0.8741,
    0.878, 0.8845, 0.8935, 0.9181, 0.9505, 0.9893, 0.85379, 0.8549, 0.858, 0.8632, 0.8703, 0.8901, 0.9164, 0.9482,
    0.83795, 0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916, 0.82435, 0.8251, 0.8273, 0.8308, 0.8356,
    0.8493, 0.8676, 0.8901, 0.80184, 0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504, 0.78363, 0.784,
    0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212, 0.76834, 0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878,
    0.7983, 0.75518, 0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797, 0.74364, 0.7438, 0.7445, 0.7455,
    0.747, 0.7512, 0.7569, 0.7642, 0.71982, 0.72, 0.7204, 0.7211, 0.7221, 0.725, 0.7289, 0.7339, 0.70097, 0.7011,
    0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112, 0.68545, 0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905,
    0.6932, 0.67232, 0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784, 0.65099, 0.651, 0.6512, 0.6513,
    0.6516, 0.6524, 0.6534, 0.6546, 0.61397, 0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147, 0.5887,
    0.5889, 0.5894, 0.59, 0.5903, 0.5901, 0.5895, 0.5885,
])
_MM_TSTAR = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0,
    7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0, 500.0,
])
_MM_ASTAR_TABLE = np.array([
    1.0065, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.084, 1.0231, 1.066, 1.038, 1.04, 1.043, 1.05, 1.052, 1.051,
    1.0424, 1.045, 1.048, 1.052, 1.056, 1.065, 1.066, 1.064, 1.0719, 1.067, 1.06, 1.055, 1.058, 1.068, 1.071,
    1.071, 1.0936, 1.087, 1.077, 1.069, 1.068, 1.075, 1.078, 1.078, 1.1053, 1.098, 1.088, 1.08, 1.078, 1.082,
    1.084, 1.084, 1.1104, 1.104, 1.096, 1.089, 1.086, 1.089, 1.09, 1.09, 1.1114, 1.107, 1.1, 1.095, 1.093, 1.095,
    1.096, 1.095, 1.1104, 1.107, 1.102, 1.099, 1.098, 1.1, 1.1, 1.099, 1.1086, 1.106, 1.102, 1.101, 1.101, 1.105,
    1.105, 1.104, 1.1063, 1.104, 1.103, 1.103, 1.104, 1.108, 1.109, 1.108, 1.102, 1.102, 1.103, 1.105, 1.107,
    1.112, 1.115, 1.115, 1.0985, 1.099, 1.101, 1.104, 1.108, 1.115, 1.119, 1.12, 1.096, 1.096, 1.099, 1.103, 1.108,
    1.116, 1.121, 1.124, 1.0943, 1.095, 1.099, 1.102, 1.108, 1.117, 1.123, 1.126, 1.0934, 1.094, 1.097, 1.102,
    1.107, 1.116, 1.123, 1.128, 1.0926, 1.094, 1.097, 1.099, 1.105, 1.115, 1.123, 1.13, 1.0934, 1.095, 1.097,
    1.099, 1.104, 1.113, 1.122, 1.129, 1.0948, 1.096, 1.098, 1.1, 1.103, 1.112, 1.119, 1.127, 1.0965, 1.097, 1.099,
    1.101, 1.104, 1.11, 1.118, 1.126, 1.0997, 1.1, 1.101, 1.102, 1.105, 1.11, 1.116, 1.123, 1.1025, 1.103, 1.104,
    1.105, 1.106, 1.11, 1.115, 1.121, 1.105, 1.105, 1.106, 1.107, 1.108, 1.111, 1.115, 1.12, 1.1072, 1.107, 1.108,
    1.108, 1.109, 1.112, 1.115, 1.119, 1.1091, 1.109, 1.109, 1.11, 1.111, 1.113, 1.115, 1.119, 1.1107, 1.111,
    1.111, 1.111, 1.112, 1.114, 1.116, 1.119, 1.1133, 1.114, 1.113, 1.114, 1.114, 1.115, 1.117, 1.119, 1.1154,
    1.115, 1.116, 1.116, 1.116, 1.117, 1.118, 1.12, 1.1172, 1.117, 1.117, 1.118, 1.118, 1.118, 1.119, 1.12, 1.1186,
    1.119, 1.119, 1.119, 1.119, 1.119, 1.12, 1.121, 1.1199, 1.12, 1.12, 1.12, 1.12, 1.121, 1.121, 1.122, 1.1223,
    1.122, 1.122, 1.122, 1.122, 1.123, 1.123, 1.124, 1.1243, 1.124, 1.124, 1.124, 1.124, 1.124, 1.125, 1.125,
    1.1259, 1.126, 1.126, 1.126, 1.126, 1.126, 1.126, 1.126, 1.1273, 1.127, 1.127, 1.127, 1.127, 1.127, 1.127,
    1.128, 1.1297, 1.13, 1.13, 1.13, 1.13, 1.13, 1.13, 1.129, 1.1339, 1.134, 1.134, 1.135, 1.135, 1.134, 1.134,
    1.132, 1.1364, 1.137, 1.137, 1.138, 1.139, 1.138, 1.137, 1.135, 1.14187, 1.14187, 1.14187, 1.14187, 1.14187,
    1.14187, 1.14187, 1.14187,
])
_MM_BSTAR_TABLE = np.array([
    1.1852, 1.2963, 1.2963, 1.2963, 1.2963, 1.2963, 1.2963, 1.2963, 1.196, 1.216, 1.237, 1.269, 1.285, 1.29, 1.297,
    1.294, 1.2451, 1.257, 1.34, 1.389, 1.366, 1.327, 1.314, 1.278, 1.29, 1.294, 1.272, 1.258, 1.262, 1.282, 1.29,
    1.299, 1.2986, 1.291, 1.284, 1.278, 1.277, 1.288, 1.294, 1.297, 1.2865, 1.281, 1.276, 1.272, 1.277, 1.286,
    1.292, 1.298, 1.2665, 1.264, 1.261, 1.263, 1.269, 1.284, 1.292, 1.298, 1.2455, 1.244, 1.248, 1.255, 1.262,
    1.278, 1.289, 1.296, 1.2253, 1.225, 1.234, 1.24, 1.252, 1.271, 1.284, 1.295, 1.2078, 1.21, 1.216, 1.227, 1.242,
    1.264, 1.281, 1.292, 1.1919, 1.192, 1.205, 1.216, 1.23, 1.256, 1.273, 1.287, 1.1678, 1.172, 1.181, 1.195,
    1.209, 1.237, 1.261, 1.277, 1.1496, 1.155, 1.161, 1.174, 1.189, 1.221, 1.246, 1.266, 1.1366, 1.141, 1.147,
    1.159, 1.174, 1.202, 1.231, 1.256, 1.127, 1.13, 1.138, 1.148, 1.162, 1.191, 1.218, 1.242, 1.1197, 1.122, 1.129,
    1.14, 1.149, 1.178, 1.205, 1.231, 1.108, 1.11, 1.116, 1.122, 1.132, 1.154, 1.18, 1.205, 1.1016, 1.103, 1.107,
    1.112, 1.12, 1.138, 1.16, 1.183, 1.098, 1.099, 1.102, 1.106, 1.112, 1.127, 1.145, 1.165, 1.0958, 1.097, 1.099,
    1.102, 1.107, 1.119, 1.135, 1.153, 1.0935, 1.094, 1.095, 1.097, 1.1, 1.109, 1.12, 1.134, 1.0925, 1.092, 1.094,
    1.095, 1.098, 1.104, 1.112, 1.122, 1.0922, 1.092, 1.093, 1.094, 1.096, 1.1, 1.106, 1.115, 1.0922, 1.092, 1.093,
    1.093, 1.095, 1.098, 1.103, 1.11, 1.0923, 1.092, 1.093, 1.093, 1.094, 1.097, 1.101, 1.106, 1.0923, 1.092,
    1.092, 1.093, 1.094, 1.096, 1.099, 1.103, 1.0927, 1.093, 1.093, 1.093, 1.094, 1.095, 1.098, 1.101, 1.093,
    1.093, 1.093, 1.093, 1.094, 1.094, 1.096, 1.099, 1.0933, 1.094, 1.093, 1.094, 1.094, 1.095, 1.096, 1.098,
    1.0937, 1.093, 1.094, 1.094, 1.094, 1.094, 1.096, 1.097, 1.0939, 1.094, 1.094, 1.094, 1.094, 1.095, 1.095,
    1.097, 1.0943, 1.094, 1.094, 1.094, 1.095, 1.095, 1.096, 1.096, 1.0944, 1.095, 1.094, 1.094, 1.094, 1.095,
    1.095, 1.096, 1.0944, 1.094, 1.095, 1.094, 1.094, 1.095, 1.096, 1.096, 1.0943, 1.095, 1.094, 1.094, 1.095,
    1.095, 1.095, 1.095, 1.0941, 1.094, 1.094, 1.094, 1.094, 1.094, 1.094, 1.096, 1.0947, 1.095, 1.094, 1.094,
    1.093, 1.093, 1.094, 1.095, 1.0957, 1.095, 1.094, 1.093, 1.092, 1.093, 1.093, 1.094, 1.10185, 1.10185, 1.10185,
    1.10185, 1.10185, 1.10185, 1.10185, 1.10185,
])
_MM_CSTAR_TABLE = np.array([
    0.8889, 0.77778, 0.77778, 0.77778, 0.77778, 0.77778, 0.77778, 0.77778, 0.88575, 0.8988, 0.8378, 0.8029, 0.7876,
    0.7805, 0.7799, 0.7801, 0.87268, 0.8692, 0.8647, 0.8479, 0.8237, 0.7975, 0.7881, 0.7784, 0.85182, 0.8525,
    0.8366, 0.8198, 0.8054, 0.7903, 0.7839, 0.782, 0.83542, 0.8362, 0.8306, 0.8196, 0.8076, 0.7918, 0.7842, 0.7806,
    0.82629, 0.8278, 0.8252, 0.8169, 0.8074, 0.7916, 0.7838, 0.7802, 0.82299, 0.8249, 0.823, 0.8165, 0.8072,
    0.7922, 0.7839, 0.7798, 0.82357, 0.8257, 0.8241, 0.8178, 0.8084, 0.7927, 0.7839, 0.7794, 0.82657, 0.828,
    0.8264, 0.8199, 0.8107, 0.7939, 0.7842, 0.7796, 0.8311, 0.8234, 0.8295, 0.8228, 0.8136, 0.796, 0.7854, 0.7798,
    0.8363, 0.8366, 0.8342, 0.8267, 0.8168, 0.7986, 0.7864, 0.7805, 0.84762, 0.8474, 0.8438, 0.8358, 0.825, 0.8041,
    0.7904, 0.7822, 0.85846, 0.8583, 0.853, 0.8444, 0.8336, 0.8118, 0.7957, 0.7854, 0.8684, 0.8674, 0.8619, 0.8531,
    0.8423, 0.8186, 0.8011, 0.7898, 0.87713, 0.8755, 0.8709, 0.8616, 0.8504, 0.8265, 0.8072, 0.7939, 0.88479,
    0.8831, 0.8779, 0.8695, 0.8578, 0.8338, 0.8133, 0.799, 0.89972, 0.8986, 0.8936, 0.8846, 0.8742, 0.8504, 0.8294,
    0.8125, 0.91028, 0.9089, 0.9043, 0.8967, 0.8869, 0.8649, 0.8438, 0.8253, 0.91793, 0.9166, 0.9125, 0.9058,
    0.897, 0.8768, 0.8557, 0.8372, 0.92371, 0.9226, 0.9189, 0.9128, 0.905, 0.8861, 0.8664, 0.8484, 0.93135, 0.9304,
    0.9274, 0.9226, 0.9164, 0.9006, 0.8833, 0.8662, 0.93607, 0.9353, 0.9329, 0.9291, 0.924, 0.9109, 0.8958, 0.8802,
    0.93927, 0.9387, 0.9366, 0.9334, 0.9292, 0.9162, 0.905, 0.8911, 0.94149, 0.9409, 0.9393, 0.9366, 0.9331,
    0.9236, 0.9122, 0.8997, 0.94306, 0.9426, 0.9412, 0.9388, 0.9357, 0.9276, 0.9175, 0.9065, 0.94419, 0.9437,
    0.9425, 0.9406, 0.938, 0.9308, 0.9219, 0.9119, 0.94571, 0.9455, 0.9445, 0.943, 0.9409, 0.9353, 0.9283, 0.9201,
    0.94662, 0.9464, 0.9456, 0.9444, 0.9428, 0.9382, 0.9325, 0.9258, 0.94723, 0.9471, 0.9464, 0.9455, 0.9442,
    0.9405, 0.9355, 0.9298, 0.94764, 0.9474, 0.9469, 0.9462, 0.945, 0.9418, 0.9378, 0.9328, 0.9479, 0.9478, 0.9474,
    0.9465, 0.9457, 0.943, 0.9394, 0.9352, 0.94827, 0.9481, 0.948, 0.9472, 0.9467, 0.9447, 0.9422, 0.9391, 0.94842,
    0.9484, 0.9481, 0.9478, 0.9472, 0.9458, 0.9437, 0.9415, 0.94852, 0.9484, 0.9483, 0.948, 0.9475, 0.9465, 0.9449,
    0.943, 0.94861, 0.9487, 0.9484, 0.9481, 0.9479, 0.9468, 0.9455, 0.943, 0.94872, 0.9486, 0.9486, 0.9483, 0.9482,
    0.9475, 0.9464, 0.9452, 0.94881, 0.9488, 0.9489, 0.949, 0.9487, 0.9482, 0.9476, 0.9468, 0.94863, 0.9487,
    0.9489, 0.9491, 0.9493, 0.9491, 0.9483, 0.9476, 0.94444, 0.94444, 0.94444, 0.94444, 0.94444, 0.94444, 0.94444,
    0.94444,
])


def polyfit(x, y, w, deg):
    """Cantera polyfit: least squares for p_0 + p_1 x + ... + p_deg x^deg; w (squared weights) used if w[0] > 0."""
    x = np.asarray(x, dtype=float)
    A = np.empty((len(x), deg + 1))
    A[:, 0] = 1.0
    if deg > 0:
        A[:, 1] = x
    for i in range(1, deg):
        A[:, i + 1] = A[:, i] * x
    y = np.array(y, dtype=float)
    if w is not None and w[0] > 0:
        sw = np.sqrt(np.asarray(w, dtype=float))
        A = sw[:, None] * A
        y = y * sw
    Q, R = np.linalg.qr(A)
    return np.linalg.solve(R, Q.T @ y)


def _poly(x, c):
    """Horner, highest power first (Cantera poly4 / poly6)."""
    r = c[-1]
    for ck in c[-2::-1]:
        r = r * x + ck
    return r


def _quad_interp(x0, x, y):
    dx21 = x[1] - x[0]
    dx32 = x[2] - x[1]
    dx31 = dx21 + dx32
    dy32 = y[2] - y[1]
    dy21 = y[1] - y[0]
    a = (dx21 * dy32 - dy21 * dx32) / (dx21 * dx31 * dx32)
    return a * (x0 - x[0]) * (x0 - x[1]) + (dy21 / dx21) * (x0 - x[1]) + y[1]


class MMCollisionInt:
    """Monchick-Mason reduced collision integrals Omega22*, A* (Omega11* = Omega22*/A*) vs T* and delta*."""

    def __init__(self):
        self.logT = [math.log(_MM_TSTAR[i + 1]) for i in range(37)]
        self.o22poly = [polyfit(_MM_DELTA, _MM_OMEGA22_TABLE[8 * i:8 * i + 8], [-1.0], 6) for i in range(37)]
        self.apoly = [polyfit(_MM_DELTA, _MM_ASTAR_TABLE[8 * (i + 1):8 * (i + 1) + 8], [-1.0], 6) for i in range(37)]

    def _interp(self, ts, ds, table, polys, shift):
        i = next((i for i in range(37) if ts < _MM_TSTAR22[i]), 37)
        i1 = max(i - 1, 0)
        i2 = i1 + 3
        if i2 > 36:
            i2 = 36
            i1 = i2 - 3
        vals = [table[8 * (i + shift)] if ds == 0.0 else _poly(ds, polys[i]) for i in range(i1, i2)]
        return _quad_interp(math.log(ts), self.logT[i1:i1 + 3], vals)

    def omega22(self, ts, ds):
        return self._interp(ts, ds, _MM_OMEGA22_TABLE, self.o22poly, 0)

    def astar(self, ts, ds):
        return self._interp(ts, ds, _MM_ASTAR_TABLE, self.apoly, 1)

    def omega11(self, ts, ds):
        return self.omega22(ts, ds) / self.astar(ts, ds)


def _cp_R(th, T_mid, k, t):
    """NASA7 cp/R as Cantera NasaPoly1::updateProperties (low range for T <= T_mid)."""
    c = th[0:7, k] if t <= T_mid[k] else th[7:14, k]
    tt0, tt1 = t, t * t
    tt2 = tt1 * t
    tt3 = tt2 * t
    return c[0] + c[1] * tt0 + c[2] * tt1 + c[3] * tt2 + c[4] * tt3


def transport_fits(species, mech):
    """Mixture-averaged species fits (Cantera GasTransport::fitProperties / fitDiffCoeffs, degree 4 in ln T):
    visc (sqrt(mu/sqrt(T))), cond (lambda/sqrt(T)), diff (D_ij*P/T^1.5, pairs i < j)."""
    nsp = len(species)
    geo = {"atom": 0.0, "linear": 1.0, "nonlinear": 1.5}
    crot, sigma, eps, dip, polar, alpha, zrot = [], [], [], [], [], [], []
    for s in species:
        tr = s.get("transport")
        if not tr or tr.get("model") != "gas":
            raise UnsupportedMechanism(f"missing gas-phase transport data for species '{s['name']}'")
        crot.append(geo[tr["geometry"]])
        sigma.append(1e-10 * float(tr["diameter"]))
        eps.append(BOLTZMANN * float(tr["well-depth"]))
        dip.append(1e-21 / LIGHT_SPEED * float(tr.get("dipole", 0.0)))
        polar.append(dip[-1] > 0)
        alpha.append(1e-30 * float(tr.get("polarizability", 0.0)))
        zrot.append(float(tr.get("rotational-relaxation", 0.0)))
    mw = mech["W_kmol"]
    red = np.zeros((nsp, nsp)); diam = np.zeros((nsp, nsp)); epsil = np.zeros((nsp, nsp)); delta = np.zeros((nsp, nsp))
    for i in range(nsp):
        for j in range(i, nsp):
            red[i, j] = mw[i] * mw[j] / (AVOGADRO * (mw[i] + mw[j]))
            diam[i, j] = 0.5 * (sigma[i] + sigma[j])
            epsil[i, j] = math.sqrt(eps[i] * eps[j])
            dij = math.sqrt(dip[i] * dip[j])
            d = diam[i, j]
            delta[i, j] = 0.5 * dij * dij / (4 * PI * EPSILON_0 * epsil[i, j] * d * d * d)
            if polar[i] == polar[j]:
                f_eps = f_sigma = 1.0
            else:
                kp = i if polar[i] else j
                knp = j if kp == i else i
                d3np, d3p = math.pow(sigma[knp], 3), math.pow(sigma[kp], 3)
                alpha_star = alpha[knp] / d3np
                mu_p_star = dip[kp] / math.sqrt(4 * PI * EPSILON_0 * d3p * eps[kp])
                xi = 1.0 + 0.25 * alpha_star * mu_p_star * mu_p_star * math.sqrt(eps[kp] / eps[knp])
                f_sigma = math.pow(xi, -1.0 / 6.0)
                f_eps = xi * xi
            diam[i, j] *= f_sigma
            epsil[i, j] *= f_eps
            for a in (red, diam, epsil, delta):
                a[j, i] = a[i, j]
    mm = MMCollisionInt()
    npt, deg = 50, 4
    tmin, tmax = float(np.max(mech["Tlim"][0])), float(np.min(mech["Tlim"][1]))
    dt = (tmax - tmin) / (npt - 1)
    ts = [tmin + dt * n for n in range(npt)]
    tlog = [math.log(t) for t in ts]
    th, T_mid = mech["thermo"], mech["T_mid"]
    visc_c, cond_c = [], []
    for k in range(nsp):
        tstar = BOLTZMANN * 298.0 / eps[k]
        fz_298 = 1.0 + math.pow(PI, 1.5) / math.sqrt(tstar) * (0.5 + 1.0 / tstar) + (0.25 * PI * PI + 2) / tstar
        spvisc, spcond = [], []
        for t in ts:
            cp_R = _cp_R(th, T_mid, k, t)
            tstar = BOLTZMANN * t / eps[k]
            sqrt_T = math.sqrt(t)
            om22 = mm.omega22(tstar, delta[k, k])
            om11 = mm.omega11(tstar, delta[k, k])
            diffcoeff = (3.0 / 16.0 * math.sqrt(2.0 * PI / red[k, k]) * math.pow(BOLTZMANN * t, 1.5)
                         / (PI * sigma[k] * sigma[k] * om11))
            visc = 5.0 / 16.0 * math.sqrt(PI * mw[k] * BOLTZMANN * t / AVOGADRO) / (om22 * PI * sigma[k] * sigma[k])
            f_int = mw[k] / (GAS_CONSTANT * t) * diffcoeff / visc
            cv_rot = crot[k]
            A_factor = 2.5 - f_int
            fz_tstar = 1.0 + math.pow(PI, 1.5) / math.sqrt(tstar) * (0.5 + 1.0 / tstar) + (0.25 * PI * PI + 2) / tstar
            B_factor = zrot[k] * fz_298 / fz_tstar + 2.0 / PI * (5.0 / 3.0 * cv_rot + f_int)
            c1 = 2.0 / PI * A_factor / B_factor
            cv_int = cp_R - 2.5 - cv_rot
            f_rot = f_int * (1.0 + c1)
            f_trans = 2.5 * (1.0 - c1 * cv_rot / 1.5)
            cond = (visc / mw[k]) * GAS_CONSTANT * (f_trans * 1.5 + f_rot * cv_rot + f_int * cv_int)
            spvisc.append(math.sqrt(visc / sqrt_T))
            spcond.append(cond / sqrt_T)
        visc_c.append(polyfit(tlog, spvisc, [1.0 / (v * v) for v in spvisc], deg))
        cond_c.append(polyfit(tlog, spcond, [1.0 / (v * v) for v in spcond], deg))
    diff_c = []
    for j in range(1, nsp):                           # pairs (i < j) in the order of the generated module
        for i in range(j):
            vals = []
            for t in ts:
                tstar = BOLTZMANN * t / epsil[j, i]
                om11 = mm.omega11(tstar, delta[j, i])
                dc = (3.0 / 16.0 * math.sqrt(2.0 * PI / red[i, j]) * math.pow(BOLTZMANN * t, 1.5)
                      / (PI * diam[j, i] * diam[j, i] * om11))
                vals.append(dc / math.pow(t, 1.5))
            diff_c.append(polyfit(tlog, vals, [1.0 / (v * v) for v in vals], deg))
    return dict(visc=np.array(visc_c).T, cond=np.array(cond_c).T,
                diff=np.array(diff_c).T if diff_c else np.zeros((deg + 1, 0)))


# --------------------------------------------------------------------------------------------------
#  Code generation
# --------------------------------------------------------------------------------------------------
def conc_product(conc, ids, sidx):
    terms = []
    for s, o in conc.items():
        c = f"c({ids[sidx[s]]})"
        if o == 1.0:
            terms.append(c)
        elif float(o).is_integer() and o > 0:
            terms.append(f"{c}**{int(o)}")
        else:
            terms.append(f"{c}**{fr(o)}")
    return " * ".join(terms)



def write_fcmech(mech, output_path, custom_f90=None):
    yaml_path, output_path = Path(mech["source"]), Path(output_path)
    sp = list(mech["species"])
    nS = len(sp)
    sidx = {s: i for i, s in enumerate(sp)}
    ids = f90_species_ids(sp)
    T = build_tables(mech)
    nRf, rev = T["nRf"], T["rev_list"]
    nRr = len(rev)
    nR = nRf + nRr
    rev_w = {ir: nRf + j + 1 for j, ir in enumerate(rev)}
    mgroups = T["mgroups"]
    n_m = max(1, len(mgroups))
    rxn = mech["reactions"]

    # thermo (NASA7 only)
    T_mid = np.asarray(mech["T_mid"], dtype=float)
    th = np.asarray(mech["thermo"], dtype=float)  # rows 1-7 low, 8-14 high (Fortran thermo_coeffs(14,nS))
    W = np.asarray(mech["W"], dtype=float)
    nA = len(mech["elements"])
    comp = np.asarray(mech["comp"], dtype=int).reshape(nA, nS)

    L = []
    L += ["!--------------------------------------------------------------------------------------------------",
          f"!  FILE {output_path.name}",
          "!  Module for chemical kinetics in NGA2 (finite chemistry)",
          "!  Generated by yaml2nga_lib_cf.py from a Cantera YAML mechanism. Do not edit manually.",
          *([f"!  Source terms: ARCANE customkinetics in {Path(custom_f90).name} (compile it with this module)"]
            if custom_f90 else []),
          f"!  Source: {yaml_path.name}   ({nS} species, {nRf} reactions, {nRr} reversible)",
          "!  Units: c mol/m^3, k (m^3/mol)^(n-1)/s, w and cdot mol/(m^3 s), h J/mol, cp J/(mol K), P Pa, T K",
          "!  Reverse rates: k_b = k_f/K_c from NASA7 thermo at run time (Cantera convention)",
          "!--------------------------------------------------------------------------------------------------",
          "", "module fcmech", "  use precision", "  use string", "  implicit none", ""]
    L += ["  !  nS, nR, nA      : species, reaction directions (nRf forward + nRr reverse), elements",
          "  !  sXXX            : species index of species XXX",
          "  !  w(1:nRf)        : forward rate of Cantera reaction i (same order as the YAML)",
          "  !  w(nRf+j)        : reverse rate of the j-th reversible reaction (rev_idx(j) = its forward index)",
          "  !  M(1:n_m)        : effective third-body concentrations (one per distinct efficiency set)",
          ""]
    L.append(f"  real(WP), parameter :: Rcst = {fr(R_UNIV)}")
    L.append(f"  real(WP), parameter :: Pref = {fr(mech["Pref"])}  ! NASA7 reference pressure (Pa)")
    L.append(f"  integer, parameter :: nS = {nS}")
    L.append(f"  integer, parameter :: nRf = {nRf}")
    L.append(f"  integer, parameter :: nRr = {nRr}")
    L.append(f"  integer, parameter :: nR = {nR}")
    L.append(f"  integer, parameter :: nA = {nA}")
    L.append(f"  integer, parameter :: n_m = {n_m}")
    L.append("")
    L.append("  ! --- Species indices ---")
    for i, (s, sid) in enumerate(zip(sp, ids)):
        L.append(f"  integer, parameter :: {sid} = {i + 1}  ! {comment_safe(s)}")
    L.append("")

    decl, data = [], []
    declare(decl, "real(WP)", "W_sp", [nS]); data_1d(data, "W_sp", W, fr)
    if nA > 0:
        decl.append(f"  character(len=2), dimension(nA), parameter :: atom_names = (/ &")
        for m in range(nA):
            decl.append(f"       '{mech["elements"][m][:2].ljust(2)}'" + (", &" if m < nA - 1 else " /)"))
        declare(decl, "real(WP)", "atom_masses", [nA]); data_1d(data, "atom_masses", mech["atomic_weights"], fr)
        declare(decl, "integer", "comp_matrix", [nA, nS]); data_cols(data, "comp_matrix", comp, fi, per_line=12, per_stmt=240)
    declare(decl, "real(WP)", "T_mid", [nS]); data_1d(data, "T_mid", T_mid, fr)
    declare(decl, "real(WP)", "thermo_coeffs", [14, nS]); data_cols(data, "thermo_coeffs", th, fr)

    # Arrhenius / three-body arrays
    arr, tb = T["arr"], T["tb"]
    neg_arr = any(a[1] < 0 for a in arr)
    neg_tb = any(a[1] < 0 for a in tb)
    for tag, rows, neg in (("arr", arr, neg_arr), ("tb", tb, neg_tb)):
        n = len(rows)
        L.append(f"  integer, parameter :: n_{tag} = {n}")
        if n == 0:
            continue
        declare(decl, "integer", f"{tag}_iw", [n]); data_1d(data, f"{tag}_iw", [r[0] for r in rows], fi, per_line=12, per_stmt=240)
        declare(decl, "real(WP)", f"{tag}_lnA", [n]); data_1d(data, f"{tag}_lnA", [lnA_sgn(r[1])[0] for r in rows], fr)
        declare(decl, "real(WP)", f"{tag}_b", [n]); data_1d(data, f"{tag}_b", [r[2] for r in rows], fr)
        declare(decl, "real(WP)", f"{tag}_ER", [n]); data_1d(data, f"{tag}_ER", [r[3] for r in rows], fr)
        if neg:
            declare(decl, "real(WP)", f"{tag}_sgn", [n]); data_1d(data, f"{tag}_sgn", [lnA_sgn(r[1])[1] for r in rows], fr)

    # PLOG tables
    plog = T["plog"]
    n_plog = len(plog)
    L.append(f"  integer, parameter :: n_plog = {n_plog}")
    if n_plog:
        p_iw, p_n0, p_nn, n_lnP, n_e0, n_ne, e_lnA, e_b, e_ER, e_sgn = [], [], [], [], [], [], [], [], [], []
        for iw, nodes in plog:
            p_iw.append(iw); p_n0.append(len(n_lnP) + 1); p_nn.append(len(nodes))
            for P, ents in nodes:
                n_lnP.append(np.log(P)); n_e0.append(len(e_lnA) + 1); n_ne.append(len(ents))
                for A, b, E_R in ents:
                    lnA, sg = lnA_sgn(A)
                    e_lnA.append(lnA); e_b.append(b); e_ER.append(E_R); e_sgn.append(sg)
        L.append(f"  integer, parameter :: n_plog_nodes = {len(n_lnP)}, n_plog_ent = {len(e_lnA)}")
        for nm, v, f in (("plog_iw", p_iw, fi), ("plog_n0", p_n0, fi), ("plog_nn", p_nn, fi)):
            declare(decl, "integer", nm, [n_plog]); data_1d(data, nm, v, f, per_line=12, per_stmt=240)
        declare(decl, "real(WP)", "plog_node_lnP", ["n_plog_nodes"]); data_1d(data, "plog_node_lnP", n_lnP, fr)
        for nm, v in (("plog_node_e0", n_e0), ("plog_node_ne", n_ne)):
            declare(decl, "integer", nm, ["n_plog_nodes"]); data_1d(data, nm, v, fi, per_line=12, per_stmt=240)
        for nm, v in (("plog_lnA", e_lnA), ("plog_b", e_b), ("plog_ER", e_ER), ("plog_sgn", e_sgn)):
            declare(decl, "real(WP)", nm, ["n_plog_ent"]); data_1d(data, nm, v, fr)

    # Chebyshev tables
    cheb = T["cheb"]
    n_cheb = len(cheb)
    L.append(f"  integer, parameter :: n_cheb = {n_cheb}")
    if n_cheb:
        c_iw, c_nT, c_nP, c_c0, c_rTmin, c_rTmax, c_lPmin, c_lPmax, c_a = [], [], [], [], [], [], [], [], []
        for iw, Tmin, Tmax, Pmin, Pmax, a in cheb:
            c_iw.append(iw); c_nT.append(a.shape[0]); c_nP.append(a.shape[1]); c_c0.append(len(c_a) + 1)
            c_rTmin.append(1.0 / Tmin); c_rTmax.append(1.0 / Tmax)
            c_lPmin.append(np.log10(Pmin)); c_lPmax.append(np.log10(Pmax))
            c_a.extend(a.flatten(order="F"))
        L.append(f"  integer, parameter :: n_cheb_coef = {len(c_a)}, cheb_max = {max(max(c_nT), max(c_nP))}")
        for nm, v in (("cheb_iw", c_iw), ("cheb_nT", c_nT), ("cheb_nP", c_nP), ("cheb_c0", c_c0)):
            declare(decl, "integer", nm, [n_cheb]); data_1d(data, nm, v, fi, per_line=12, per_stmt=240)
        for nm, v in (("cheb_rTmin", c_rTmin), ("cheb_rTmax", c_rTmax), ("cheb_lPmin", c_lPmin), ("cheb_lPmax", c_lPmax)):
            declare(decl, "real(WP)", nm, [n_cheb]); data_1d(data, nm, v, fr)
        declare(decl, "real(WP)", "cheb_a", ["n_cheb_coef"]); data_1d(data, "cheb_a", c_a, fr)

    if nRr:
        declare(decl, "integer", "rev_idx", [nRr]); data_1d(data, "rev_idx", rev, fi, per_line=12, per_stmt=240)

    # Transport fits
    npair = nS * (nS - 1) // 2
    tr = mech.get("transport")
    if tr is not None:
        vf, cf, df = tr["visc"], tr["cond"], tr["diff"]
        nfit = vf.shape[0]
        L.append(f"  integer, parameter :: n_fit = {nfit}, n_pair = {npair}")
        declare(decl, "real(WP)", "visc_fit", ["n_fit", nS]); data_cols(data, "visc_fit", vf, fr)
        declare(decl, "real(WP)", "cond_fit", ["n_fit", nS]); data_cols(data, "cond_fit", cf, fr)
        if npair:
            declare(decl, "real(WP)", "diff_fit", ["n_fit", "n_pair"]); data_cols(data, "diff_fit", df, fr)

    L.append("")
    L.append("  ! --- Thermodynamic: Cp and enthalpy (module variables, set by fcmech_thermodata) ---")
    L.append("  real(WP), dimension(nS) :: Cpsp, hsp")
    L.append("")
    L.append("  ! --- Constant tables (data statements below) ---")
    L += decl
    L.append("")
    L += data
    L += ["", "contains", ""]

    # ---------------------------------------------------------------- thermo
    L += ["  ! NASA7: h/RT, cp/R, s/R; low range for T <= T_mid (Cantera convention)",
          "  subroutine fcmech_nasa7(T, hort, cpor, sor)",
          "    implicit none",
          "    real(WP), intent(in) :: T",
          "    real(WP), dimension(nS), intent(out) :: hort, cpor, sor",
          "    real(WP) :: T2, T3, T4, T_inv, T_log",
          "    integer :: i, o",
          "    T2 = T*T; T3 = T2*T; T4 = T3*T; T_inv = 1.0_WP/T; T_log = log(T)",
          "    do i = 1, nS",
          "      o = 7",
          "      if (T <= T_mid(i)) o = 0",
          "      hort(i) = thermo_coeffs(o+1,i) + thermo_coeffs(o+2,i)*T/2.0_WP + thermo_coeffs(o+3,i)*T2/3.0_WP &",
          "              + thermo_coeffs(o+4,i)*T3/4.0_WP + thermo_coeffs(o+5,i)*T4/5.0_WP + thermo_coeffs(o+6,i)*T_inv",
          "      cpor(i) = thermo_coeffs(o+1,i) + thermo_coeffs(o+2,i)*T + thermo_coeffs(o+3,i)*T2 &",
          "              + thermo_coeffs(o+4,i)*T3 + thermo_coeffs(o+5,i)*T4",
          "      sor(i)  = thermo_coeffs(o+1,i)*T_log + thermo_coeffs(o+2,i)*T + thermo_coeffs(o+3,i)*T2/2.0_WP &",
          "              + thermo_coeffs(o+4,i)*T3/3.0_WP + thermo_coeffs(o+5,i)*T4/4.0_WP + thermo_coeffs(o+7,i)",
          "    end do",
          "  end subroutine fcmech_nasa7", "",
          "  subroutine fcmech_thermodata(T)",
          "    implicit none",
          "    real(WP), intent(in) :: T",
          "    real(WP), dimension(nS) :: hort, cpor, sor",
          "    call fcmech_nasa7(T, hort, cpor, sor)",
          "    hsp = Rcst*T*hort",
          "    Cpsp = Rcst*cpor",
          "  end subroutine fcmech_thermodata", "",
          "  subroutine fcmech_get_thermodata(h, cp, T)",
          "    implicit none",
          "    ! h(nS), cp(nS) : enthalpy (J/mol), heat capacity (J/(mol K))",
          "    real(WP), dimension(nS), intent(out) :: h, cp",
          "    real(WP), intent(in) :: T",
          "    call fcmech_thermodata(T)",
          "    h = hsp",
          "    cp = Cpsp",
          "  end subroutine fcmech_get_thermodata", "",
          "  ! Standard-state g/RT of every species",
          "  subroutine fcmech_get_gort(gort, T)",
          "    implicit none",
          "    real(WP), dimension(nS), intent(out) :: gort",
          "    real(WP), intent(in) :: T",
          "    real(WP), dimension(nS) :: hort, cpor, sor",
          "    call fcmech_nasa7(T, hort, cpor, sor)",
          "    gort = hort - sor",
          "  end subroutine fcmech_get_gort", ""]

    # ---------------------------------------------------------------- third bodies
    L += ["  subroutine get_thirdbodies(M, c)",
          "    implicit none",
          "    ! c(nS): concentrations (mol/m^3); M(n_m): effective third-body concentrations (mol/m^3)",
          "    real(WP), dimension(nS), intent(in) :: c",
          "    real(WP), dimension(n_m), intent(out) :: M",
          "    real(WP) :: ctot",
          "    ctot = sum(c)"]
    if not mgroups:
        L.append("    M = ctot")
    for im, (dflt, effs) in enumerate(mgroups, start=1):
        terms = []
        for s, e in effs:
            d = e - dflt
            if d != 0.0:
                terms.append(("+ " if d > 0 else "- ") + f"{fr(abs(d))}*c({ids[sidx[s]]})")
        first = "ctot" if dflt == 1.0 else ("0.0_WP" if dflt == 0.0 else f"{fr(dflt)}*ctot")
        L += wrap_sum(f"M({im})", terms, first=first)
    L += ["  end subroutine get_thirdbodies", ""]

    # ---------------------------------------------------------------- falloff helpers
    L += ["  ! Falloff blending: F = 1 (Lindemann), Troe, SRI; chem = .true. for chemically activated (as Cantera)",
          "  real(WP) function falloff_k(k0, kinf, Mc, F, chem)",
          "    implicit none",
          "    real(WP), intent(in) :: k0, kinf, Mc, F",
          "    logical, intent(in) :: chem",
          "    real(WP) :: Pr",
          "    Pr = Mc*k0/(kinf + tiny(1.0_WP))",
          "    if (chem) then",
          "      falloff_k = k0*F/(1.0_WP + Pr)",
          "    else",
          "      falloff_k = kinf*Pr*F/(1.0_WP + Pr)",
          "    end if",
          "  end function falloff_k", "",
          "  real(WP) function troe_F(Pr, T, a, rT3, rT1, T2)",
          "    implicit none",
          "    real(WP), intent(in) :: Pr, T, a, rT3, rT1, T2",
          "    real(WP) :: Fc, lgFc, lpr, cc, nn, f1",
          "    Fc = (1.0_WP - a)*exp(-min(T*rT3, 700.0_WP)) + a*exp(-min(T*rT1, 700.0_WP))",
          "    if (T2 /= 0.0_WP) Fc = Fc + exp(-T2/T)",
          "    lgFc = log10(max(Fc, 1.0e-300_WP))",
          "    lpr = log10(max(Pr, 1.0e-300_WP))",
          "    cc = -0.4_WP - 0.67_WP*lgFc",
          "    nn = 0.75_WP - 1.27_WP*lgFc",
          "    f1 = (lpr + cc)/(nn - 0.14_WP*(lpr + cc))",
          "    troe_F = 10.0_WP**(lgFc/(1.0_WP + f1*f1))",
          "  end function troe_F", "",
          "  real(WP) function sri_F(Pr, T, a, b, cc, d, e)",
          "    implicit none",
          "    real(WP), intent(in) :: Pr, T, a, b, cc, d, e",
          "    real(WP) :: lpr, x",
          "    lpr = log10(max(Pr, 1.0e-300_WP))",
          "    x = 1.0_WP/(1.0_WP + lpr*lpr)",
          "    sri_F = d*(a*exp(-b/T) + exp(-T/cc))**x * T**e",
          "  end function sri_F", ""]

    # ---------------------------------------------------------------- rate coefficients
    L += ["  !--------------------------------------------------------------------------------------------------",
          "  !  subroutine get_rate_coefficients: k(1:nRf) forward, k(nRf+1:nR) reverse = k_f/K_c",
          "  !--------------------------------------------------------------------------------------------------",
          "  subroutine get_rate_coefficients(k, M, Tloc, Ploc)",
          "    implicit none",
          "    real(WP), dimension(nR), intent(out) :: k",
          "    real(WP), dimension(:), intent(in) :: M",
          "    real(WP), intent(in) :: Tloc, Ploc",
          "    real(WP) :: T_log, T_inv, P_log, ln_c0, k0, kinf, Pr",
          "    real(WP), dimension(nS) :: gort"]
    if n_plog:
        L.append("    integer :: j, n1, n2, i1, i2")
        L.append("    real(WP) :: k1, k2")
    if n_cheb:
        L.append("    integer :: jc, it, ip, c0, nT, nP")
        L.append("    real(WP) :: Tr, Prr, lk, phT(cheb_max), phP(cheb_max)")
    L += ["    T_log = log(Tloc)", "    T_inv = 1.0_WP/Tloc", "    P_log = log(Ploc)"]
    if arr:
        L.append("    k(arr_iw) = " + ("arr_sgn*" if neg_arr else "") + "exp(arr_lnA + arr_b*T_log - arr_ER*T_inv)")
    if tb:
        L.append("    k(tb_iw) = " + ("tb_sgn*" if neg_tb else "") + "exp(tb_lnA + tb_b*T_log - tb_ER*T_inv)")
    for iw, ftype, chem, low, high, pars, im in T["fo"]:
        L.append(f"    ! {comment_safe(rxn[iw - 1]["equation"])}  ({ftype})")
        L.append(f"    k0 = {arr_expr(*low)}")
        L.append(f"    kinf = {arr_expr(*high)}")
        L.append(f"    Pr = M({im})*k0/(kinf + tiny(1.0_WP))")
        if ftype == "lindemann":
            F = "1.0_WP"
        elif ftype == "troe":
            F = "troe_F(Pr, Tloc, " + ", ".join(fr(p) for p in pars) + ")"
        else:
            F = "sri_F(Pr, Tloc, " + ", ".join(fr(p) for p in pars) + ")"
        L.append(f"    k({iw}) = falloff_k(k0, kinf, M({im}), &")
        L.append(f"      {F}, {'.true.' if chem else '.false.'})")
    if n_plog:
        L += ["    ! PLOG: same-pressure entries summed, exact at nodes, log-log interpolation, clamped outside",
              "    do j = 1, n_plog",
              "      n1 = plog_n0(j)",
              "      n2 = n1 + plog_nn(j) - 1",
              "      if (P_log <= plog_node_lnP(n1)) then",
              "        i1 = n1; i2 = n1",
              "      else if (P_log >= plog_node_lnP(n2)) then",
              "        i1 = n2; i2 = n2",
              "      else",
              "        i1 = n1",
              "        do while (plog_node_lnP(i1+1) <= P_log)",
              "          i1 = i1 + 1",
              "        end do",
              "        i2 = i1 + 1",
              "      end if",
              "      k1 = max(plog_node_k(i1, T_log, T_inv), tiny(1.0_WP))",
              "      if (i2 == i1) then",
              "        k(plog_iw(j)) = k1",
              "      else",
              "        k2 = max(plog_node_k(i2, T_log, T_inv), tiny(1.0_WP))",
              "        k(plog_iw(j)) = exp(log(k1) + (log(k2) - log(k1))*(P_log - plog_node_lnP(i1)) &",
              "                                     /(plog_node_lnP(i2) - plog_node_lnP(i1)))",
              "      end if",
              "    end do"]
    if n_cheb:
        L += ["    ! Chebyshev: log10 k = sum_t sum_p a(t,p) T_t(Tr) T_p(Pr)",
              "    do jc = 1, n_cheb",
              "      nT = cheb_nT(jc); nP = cheb_nP(jc); c0 = cheb_c0(jc)",
              "      Tr = (2.0_WP*T_inv - cheb_rTmin(jc) - cheb_rTmax(jc))/(cheb_rTmax(jc) - cheb_rTmin(jc))",
              "      Prr = (2.0_WP*P_log/log(10.0_WP) - cheb_lPmin(jc) - cheb_lPmax(jc))/(cheb_lPmax(jc) - cheb_lPmin(jc))",
              "      phT(1) = 1.0_WP; phP(1) = 1.0_WP",
              "      if (nT > 1) phT(2) = Tr",
              "      if (nP > 1) phP(2) = Prr",
              "      do it = 3, nT",
              "        phT(it) = 2.0_WP*Tr*phT(it-1) - phT(it-2)",
              "      end do",
              "      do ip = 3, nP",
              "        phP(ip) = 2.0_WP*Prr*phP(ip-1) - phP(ip-2)",
              "      end do",
              "      lk = 0.0_WP",
              "      do ip = 1, nP",
              "        do it = 1, nT",
              "          lk = lk + cheb_a(c0 + (ip-1)*nT + it - 1)*phT(it)*phP(ip)",
              "        end do",
              "      end do",
              "      k(cheb_iw(jc)) = 10.0_WP**lk",
              "    end do"]
    if nRr:
        L += ["    ! Reverse: k_b = k_f/K_c, K_c = exp(-sum(nu g/RT))*(Pref/(R T))^dnu",
              "    call fcmech_get_gort(gort, Tloc)",
              "    ln_c0 = log(Pref/(Rcst*Tloc))"]
        for j, ir in enumerate(rev):
            nu = T["net_nu"][ir - 1]
            dnu = sum(nu.values())
            terms = []
            for s, v in nu.items():
                g = f"gort({ids[sidx[s]]})"
                if v == 1.0:
                    terms.append(f"+ {g}")
                elif v == -1.0:
                    terms.append(f"- {g}")
                else:
                    terms.append(("+ " if v > 0 else "- ") + f"{fr(abs(v))}*{g}")
            if dnu != 0.0:
                terms.append(("- " if dnu > 0 else "+ ") + f"{fr(abs(dnu))}*ln_c0")
            stmt = f"    k({nRf + j + 1}) = k({ir})*exp({' '.join(terms)})"
            if len(stmt) <= 130:
                L.append(stmt)
            else:
                L.append(f"    k({nRf + j + 1}) = k({ir})*exp( &")
                rows = [" ".join(terms[i:i + 4]) for i in range(0, len(terms), 4)]
                for i, row in enumerate(rows):
                    L.append("      " + row + (" &" if i < len(rows) - 1 else ")"))
    L += ["  end subroutine get_rate_coefficients", ""]

    if n_plog:
        L += ["  real(WP) function plog_node_k(n, T_log, T_inv)",
              "    implicit none",
              "    integer, intent(in) :: n",
              "    real(WP), intent(in) :: T_log, T_inv",
              "    integer :: e",
              "    plog_node_k = 0.0_WP",
              "    do e = plog_node_e0(n), plog_node_e0(n) + plog_node_ne(n) - 1",
              "      plog_node_k = plog_node_k + plog_sgn(e)*exp(plog_lnA(e) + plog_b(e)*T_log - plog_ER(e)*T_inv)",
              "    end do",
              "  end function plog_node_k", ""]

    # ---------------------------------------------------------------- reaction rates
    tb_m = {r[0]: r[4] for r in tb}
    L += ["  subroutine get_reaction_rates(w, k, m, c)",
          "    implicit none",
          "    ! w(nR): rates of the nR reaction directions (mol/(m^3 s))",
          "    real(WP), dimension(nS), intent(in) :: c",
          "    real(WP), dimension(nR), intent(in) :: k",
          "    real(WP), dimension(:), intent(in) :: m",
          "    real(WP), dimension(nR), intent(out) :: w"]
    eqs = [""] * nR
    for ir in range(1, nRf + 1):
        r = rxn[ir - 1]
        expr = conc_product(T["rate_conc"][ir - 1], ids, sidx)
        if ir in tb_m:
            expr = (expr + " * " if expr else "") + f"m({tb_m[ir]})"
        eqs[ir - 1] = r["equation"]
        L.append(f"    w({ir}) = k({ir})" + (f" * {expr}" if expr else "") + f"  ! {comment_safe(r["equation"])}")
    for j, ir in enumerate(rev):
        r = rxn[ir - 1]
        expr = conc_product(r["products"], ids, sidx)
        if ir in tb_m:
            expr = (expr + " * " if expr else "") + f"m({tb_m[ir]})"
        iw = nRf + j + 1
        eqs[iw - 1] = f"reverse of {ir}: " + r["equation"]
        L.append(f"    w({iw}) = k({iw})" + (f" * {expr}" if expr else "") + f"  ! reverse {ir}")
    L += ["  end subroutine get_reaction_rates", ""]

    # ---------------------------------------------------------------- production rates
    contrib = [[] for _ in range(nS)]
    for ir in range(1, nRf + 1):
        for s, v in T["net_nu"][ir - 1].items():
            contrib[sidx[s]].append((ir, v))
    for j, ir in enumerate(rev):
        for s, v in T["net_nu"][ir - 1].items():
            contrib[sidx[s]].append((nRf + j + 1, -v))
    L += ["  subroutine get_production_rates(cdot, w)",
          "    implicit none",
          "    ! cdot(nS): molar production rates (mol/(m^3 s)) = sum_r nu(i,r)*w(r)",
          "    real(WP), dimension(nR), intent(in) :: w",
          "    real(WP), dimension(nS), intent(out) :: cdot"]
    for k in range(nS):
        terms = []
        for iw, v in sorted(contrib[k]):
            if v == 1.0:
                terms.append(f"+ w({iw})")
            elif v == -1.0:
                terms.append(f"- w({iw})")
            else:
                terms.append(("+ " if v > 0 else "- ") + f"{fr(abs(v))}*w({iw})")
        L += wrap_sum(f"cdot({ids[k]})", terms)
    L += ["  end subroutine get_production_rates", ""]

    # ---------------------------------------------------------------- misc
    L += ["  subroutine y2c(y, W_sp, P, T, c)",
          "    implicit none",
          "    ! c = rho*y/W_sp (mol/m^3), rho = P*W_mix/(R*T)",
          "    real(WP), dimension(nS), intent(in) :: y",
          "    real(WP), dimension(nS), intent(in) :: W_sp",
          "    real(WP), intent(in) :: P, T",
          "    real(WP), dimension(nS), intent(out) :: c",
          "    real(WP) :: W_mix, rho",
          "    W_mix = 1.0_WP/sum(y(1:nS)/W_sp(1:nS))",
          "    rho = P*W_mix/(Rcst*T)",
          "    c(1:nS) = rho*y(1:nS)/W_sp(1:nS)",
          "  end subroutine y2c", "",
          "  subroutine fcmech_get_speciesnames(names)",
          "    implicit none",
          "    character(len=*), dimension(nS), intent(out) :: names"]
    for i, s in enumerate(sp):
        L.append(f"    names({i + 1}) = '" + s.replace("'", "''") + "'")
    L += ["  end subroutine fcmech_get_speciesnames", ""]
    if nA > 0:
        L += ["  subroutine fcmech_get_atomnames(names)",
              "    implicit none",
              "    character(len=*), dimension(nA), intent(out) :: names",
              "    names(1:nA) = atom_names(1:nA)",
              "  end subroutine fcmech_get_atomnames", "",
              "  subroutine fcmech_get_atommasses(masses)",
              "    implicit none",
              "    real(WP), dimension(nA), intent(out) :: masses",
              "    masses(1:nA) = atom_masses(1:nA)",
              "  end subroutine fcmech_get_atommasses", "",
              "  subroutine fcmech_get_composition(comp)",
              "    implicit none",
              "    integer, dimension(nA,nS), intent(out) :: comp",
              "    comp(1:nA,1:nS) = comp_matrix(1:nA,1:nS)",
              "  end subroutine fcmech_get_composition", ""]

    if tr is not None:
        horner = lambda arr_, col: "".join(["("] * (nfit - 1)) + f"{arr_}({nfit},{col})" + "".join(
            f"*lnT + {arr_}({p},{col}))" for p in range(nfit - 1, 0, -1))
        L += ["  ! Transport: Cantera mixture-averaged species fits, polynomials in ln T",
              "  subroutine fcmech_get_viscosity(mu, T)",
              "    implicit none",
              "    ! mu(nS): pure-species viscosity (Pa s) = sqrt(T)*poly(lnT)**2",
              "    real(WP), dimension(nS), intent(out) :: mu",
              "    real(WP), intent(in) :: T",
              "    real(WP) :: lnT, p",
              "    integer :: i",
              "    lnT = log(T)",
              "    do i = 1, nS",
              "      p = " + horner("visc_fit", "i"),
              "      mu(i) = sqrt(T)*p*p",
              "    end do",
              "  end subroutine fcmech_get_viscosity", "",
              "  subroutine fcmech_get_conductivity(lambda, T, mu)",
              "    implicit none",
              "    ! lambda(nS): pure-species thermal conductivity (W/(m K)) = sqrt(T)*poly(lnT)",
              "    ! mu(nS)    : not used (kept for interface compatibility with yaml2nga.py)",
              "    real(WP), dimension(nS), intent(out) :: lambda",
              "    real(WP), intent(in) :: T",
              "    real(WP), dimension(nS), intent(in) :: mu",
              "    real(WP) :: lnT",
              "    integer :: i",
              "    lnT = log(T)",
              "    do i = 1, nS",
              "      lambda(i) = sqrt(T)*" + horner("cond_fit", "i"),
              "    end do",
              "  end subroutine fcmech_get_conductivity", "",
              "  subroutine fcmech_get_invDij(invDij, T, P)",
              "    implicit none",
              "    ! invDij(nS,nS): inverse binary diffusion coefficients (s/m^2), D_ij = T**1.5*poly(lnT)/P",
              "    real(WP), dimension(nS,nS), intent(out) :: invDij",
              "    real(WP), intent(in) :: T, P",
              "    real(WP) :: lnT, TPterm, poly",
              "    integer :: i, j, ip",
              "    lnT = log(T)",
              "    TPterm = P/(T*sqrt(T))",
              "    ip = 0",
              "    do j = 1, nS",
              "      invDij(j,j) = 0.0_WP",
              "      do i = 1, j-1",
              "        ip = ip + 1",
              "        poly = " + horner("diff_fit", "ip"),
              "        invDij(i,j) = TPterm/poly",
              "        invDij(j,i) = invDij(i,j)",
              "      end do",
              "    end do",
              "  end subroutine fcmech_get_invDij", ""]

    if custom_f90:
        L += ["  subroutine fcmech_get_ydot(P, T, Y, ydot)",
              "    implicit none",
              "    ! ydot(nS): dY_i/dt = W_i*cdot_i/rho (1/s); cdot from ARCANE customkinetics (mol/(m^3 s))",
              "    real(WP), intent(in) :: P, T",
              "    real(WP), dimension(nS), intent(in) :: Y",
              "    real(WP), dimension(nS), intent(out) :: ydot",
              "    real(WP), dimension(nS) :: wdot, Yloc",
              "    real(WP) :: Ploc, Tloc, W_mix, rho",
              "    external :: customkinetics",
              "    Ploc = P; Tloc = T; Yloc = Y",
              "    call customkinetics(Ploc, Tloc, Yloc, wdot)",
              "    W_mix = 1.0_WP/sum(Y(1:nS)/W_sp(1:nS))",
              "    rho = P*W_mix/(Rcst*T)",
              "    ydot(1:nS) = W_sp(1:nS)*wdot(1:nS)/rho",
              "  end subroutine fcmech_get_ydot", "",
              "end module fcmech", ""]
    else:
        L += ["  subroutine fcmech_get_ydot(P, T, Y, ydot)",
              "    implicit none",
              "    ! ydot(nS): dY_i/dt = W_i*cdot_i/rho (1/s)",
              "    real(WP), intent(in) :: P, T",
              "    real(WP), dimension(nS), intent(in) :: Y",
              "    real(WP), dimension(nS), intent(out) :: ydot",
              "    real(WP), dimension(nS) :: c, wdot",
              "    real(WP), dimension(nR) :: k, w",
              "    real(WP), dimension(n_m) :: M",
              "    real(WP) :: W_mix, rho",
              "    call y2c(Y, W_sp, P, T, c)",
              "    call get_thirdbodies(M, c)",
              "    call get_rate_coefficients(k, M, T, P)",
              "    call get_reaction_rates(w, k, M, c)",
              "    call get_production_rates(wdot, w)",
              "    W_mix = 1.0_WP/sum(Y(1:nS)/W_sp(1:nS))",
              "    rho = P*W_mix/(Rcst*T)",
              "    ydot(1:nS) = W_sp(1:nS)*wdot(1:nS)/rho",
              "  end subroutine fcmech_get_ydot", "",
              "end module fcmech", ""]

    long = [l for l in L if len(l) > 132]
    if long:
        raise RuntimeError(f"internal: {len(long)} generated lines exceed 132 characters, e.g.\n{long[0]}")
    output_path.write_text("\n".join(L))

    # reaction index file
    stem = output_path.stem
    idx_path = output_path.parent / ("chem_data_reactions.txt" if stem == "chem_data_fc" else stem + "_reactions.txt")
    if custom_f90:
        idx_path = None
    else:
        with open(idx_path, "w") as f:
            f.write(f"fcmech reaction index (yaml2nga_lib_cf.py, source {yaml_path.name})\n")
            f.write("w(1:nRf) = forward rate of YAML/Cantera reaction i; w(nRf+j) = reverse of a reversible reaction\n")
            f.write(f"{'w_idx':>6}  {'yaml':>5}  {'dir':<3}  {'type':<26}  equation\n")
            for ir in range(1, nRf + 1):
                r = rxn[ir - 1]
                f.write(f"{ir:>6}  {ir:>5}  {'fwd':<3}  {T['kinds'][ir - 1]:<26}  {r["equation"]}\n")
            for j, ir in enumerate(rev):
                r = rxn[ir - 1]
                f.write(f"{nRf + j + 1:>6}  {ir:>5}  {'rev':<3}  {T['kinds'][ir - 1]:<26}  {r["equation"]}\n")

    counts = {}
    for kd in T["kinds"]:
        counts[kd] = counts.get(kd, 0) + 1
    return dict(nS=nS, nRf=nRf, nRr=nRr, n_m=n_m, counts=counts, index=idx_path,
                renamed=[(s, i) for s, i in zip(sp, ids) if "s" + s != i], transport=tr is not None)


def check_custom_f90(path, mech):
    """ARCANE custom kinetics file: species order (W_sp comments) and molar masses must match the YAML phase."""
    src = Path(path).read_text()
    nS = len(mech["species"])
    m = re.search(r"INTEGER,\s*PARAMETER\s*::\s*nspec\s*=\s*(\d+)", src, re.I)
    if not m or "subroutine customkinetics" not in src.lower():
        raise UnsupportedMechanism(f"{path}: not an ARCANE custom kinetics file (nspec / customkinetics not found)")
    if int(m.group(1)) != nS:
        raise UnsupportedMechanism(f"{path}: nspec = {m.group(1)} but the YAML phase has {nS} species")
    blk = re.search(r"DIMENSION\(1:nspec\)\s*::\s*W_sp\s*=\s*\(/(.*?)/\)", src, re.S | re.I).group(1)
    rows = [m.groups() for ln in blk.splitlines() if not ln.strip().startswith("!")   # '!' rows: QSS species
            for m in [re.search(r"([-+0-9.eEdD]+)_pr\s*,?\s*&?\s*!\s*(\S+)", ln)] if m]
    names = [n for _, n in rows]
    # ARCANE writes mangled names in these comments ('-' -> X, '(' ')' '*' -> G): compare alphanumerics only
    same = lambda a, b: len(a) == len(b) and all(x == y or not y.isalnum() for x, y in zip(a, b))
    if len(names) != nS or not all(same(a, b) for a, b in zip(names, mech["species"])):
        bad = next((i for i, (a, b) in enumerate(zip(names, mech["species"])) if not same(a, b)), len(names))
        raise UnsupportedMechanism(f"{path}: species order differs from the YAML at index {bad + 1}")
    W = np.array([float(v.replace("d", "e").replace("D", "e")) for v, _ in rows])
    err = np.abs(W - mech["W"]) / mech["W"]
    if err.max() > 1e-3:
        raise UnsupportedMechanism(f"{path}: molar mass of {names[int(err.argmax())]} differs from the YAML")


def all_rates_zero(mech):
    """True for placeholder mechanisms (no reactions, or every pre-exponential factor zero): their rates come from a
    separate routine."""
    def zero(rate):
        if rate["type"] == "Arrhenius":
            return rate["A"] == 0.0
        if rate["type"] == "falloff":
            return rate["low"]["A"] == 0.0 and rate["high"]["A"] == 0.0
        if rate["type"] == "plog":
            return all(a["A"] == 0.0 for _, a in rate["rates"])
        return False
    return all(zero(r["rate"]) for r in mech["reactions"])


def read_mechanism(yaml_path, phase=None, transport=True, custom=False, backend="yaml"):
    if backend == "cantera":
        return read_cantera(yaml_path, phase, transport, custom)
    return read_yaml(yaml_path, phase, transport, custom)


def yaml2nga_lib(yaml_path, output_path, phase=None, transport=True, custom_f90=None, backend="yaml"):
    mech = read_mechanism(yaml_path, phase, transport, bool(custom_f90), backend)
    if custom_f90:
        check_custom_f90(custom_f90, mech)
    return write_fcmech(mech, output_path, custom_f90)


def compare_backends(yaml_path, phase=None, transport=True, custom=False):
    """Largest differences between the Cantera-free reader and Cantera (relative, per table)."""
    a = read_mechanism(yaml_path, phase, transport, custom, "yaml")
    b = read_mechanism(yaml_path, phase, transport, custom, "cantera")
    rel = lambda x, y: float(np.max(np.abs(np.asarray(x, float) - np.asarray(y, float))
                                    / np.maximum(np.abs(np.asarray(y, float)), 1e-300))) if np.size(y) else 0.0
    out = {}
    for k in ("species", "elements"):
        if list(a[k]) != list(b[k]):
            raise RuntimeError(f"{k} differ: {a[k][:5]} ... vs {b[k][:5]} ...")
    out["W"] = rel(a["W"], b["W"]); out["atomic_weights"] = rel(a["atomic_weights"], b["atomic_weights"])
    out["comp"] = float(np.max(np.abs(a["comp"] - b["comp"]))) if a["comp"].size else 0.0
    out["thermo"] = max(rel(a["thermo"], b["thermo"]), rel(a["T_mid"], b["T_mid"]), rel(a["Pref"], b["Pref"]))
    ta, tb = build_tables(a), build_tables(b)
    for key in ("rev_list", "kinds", "mgroups", "net_nu", "rate_conc"):
        if ta[key] != tb[key]:
            i = next(i for i, (x, y) in enumerate(zip(ta[key], tb[key])) if x != y) if len(ta[key]) == len(tb[key]) \
                else -1
            raise RuntimeError(f"reaction tables '{key}' differ (first at {i}: {ta[key][i]} vs {tb[key][i]})")
    num = lambda t: [v for row in t for v in _flatten(row)]
    out["rates"] = max([rel(num(ta[k]), num(tb[k])) for k in ("arr", "tb", "fo", "plog", "cheb")] + [0.0])
    if a["transport"] is not None and b["transport"] is not None:
        out["transport"] = max(rel(a["transport"][k], b["transport"][k]) for k in ("visc", "cond", "diff"))
    return out


def _flatten(x):
    if isinstance(x, (list, tuple)):
        for y in x:
            yield from _flatten(y)
    elif isinstance(x, np.ndarray):
        yield from x.ravel().tolist()
    elif isinstance(x, (int, float, np.floating, np.integer)) and not isinstance(x, bool):
        yield float(x)


def main():
    ap = argparse.ArgumentParser(description="Generate NGA2 chem_data_fc.f90 (module fcmech) from a Cantera YAML "
                                 "mechanism, without Cantera")
    ap.add_argument("input", help="Input YAML mechanism")
    ap.add_argument("output", nargs="?", default="chem_data_fc.f90", help="Output Fortran file (default: chem_data_fc.f90)")
    ap.add_argument("--phase", default=None, help="Phase name in the YAML (default: first phase)")
    ap.add_argument("--no-transport", action="store_true", help="Do not generate the transport routines")
    ap.add_argument("--custom-f90", default=None, help="ARCANE QSS custom kinetics file providing the source terms")
    ap.add_argument("--backend", choices=["yaml", "cantera"], default="yaml",
                    help="read the mechanism here (yaml, default) or through Cantera")
    ap.add_argument("--check", action="store_true", help="compare both backends (needs Cantera) and exit")
    a = ap.parse_args()
    if not Path(a.input).exists():
        print(f"Error: Input file not found: {a.input}", file=sys.stderr)
        sys.exit(1)
    try:
        if a.check:
            d = compare_backends(a.input, a.phase, not a.no_transport, bool(a.custom_f90))
            print("Largest relative difference, Cantera-free reader vs Cantera:")
            for k, v in d.items():
                print(f"  {k:15s} {v:.3e}")
            return
        print(f"Loading mechanism from {a.input} ({'Cantera' if a.backend == 'cantera' else 'Cantera-free reader'})")
        info = yaml2nga_lib(a.input, a.output, phase=a.phase, transport=not a.no_transport, custom_f90=a.custom_f90,
                            backend=a.backend)
    except UnsupportedMechanism as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"  Species: {info['nS']}   Reactions: {info['nRf']} ({info['nRr']} reversible) -> nR = {info['nRf'] + info['nRr']}")
    print("  Types: " + ", ".join(f"{k} {v}" for k, v in sorted(info["counts"].items())))
    print(f"  Third-body groups: {info['n_m']}   Transport: {'yes' if info['transport'] else 'no'}")
    if info["renamed"]:
        print("  Fortran identifiers changed by sanitizing: " + ", ".join(f"{s}->{i}" for s, i in info["renamed"][:12])
              + (" ..." if len(info["renamed"]) > 12 else ""))
    print(f"  Wrote {a.output}" + (f" and {info['index']}" if info["index"] else ""))


if __name__ == "__main__":
    main()
