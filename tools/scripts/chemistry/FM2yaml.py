#!/usr/bin/env python3
"""
FM2yaml - FlameMaster to YAML Converter

Converts chemical mechanism files from FlameMaster format (.mech, .thermo, .trans)
to YAML format for NGA2.

Reaction labels:
- Number only (e.g. 25, 26): Irreversible
- Number + f (e.g. 1f, i66f): Forward; reversible via K_eq if no matching b
- Number + b (e.g. i66b): Backward
- If both f and b exist for same base id: use explicit forward and backward rates
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

KJ_TO_CAL = 239.006  # kJ/mol to cal/mol


def strip_tabs_and_trim(s: str) -> str:
    """Replace tabs with spaces and trim. YAML does not allow tabs."""
    return s.replace("\t", " ").strip()


def uppercase_species_in_equation(eq_str: str) -> str:
    """
    Convert all species names in a reaction equation to uppercase.
    Preserves stoichiometric coefficients and structure (e.g. 2 H2 + O2 => 2 H2O).
    """
    if "<=>" in eq_str:
        left, right = eq_str.split("<=>", 1)
        arrow = "<=>"
    elif "=>" in eq_str:
        left, right = eq_str.split("=>", 1)
        arrow = "=>"
    else:
        return eq_str
    left = left.strip()
    right = right.strip()

    def process_side(s: str) -> str:
        terms = [t.strip() for t in s.split("+") if t.strip()]
        result = []
        for t in terms:
            parts = t.split()
            if len(parts) >= 2 and parts[0].isdigit():
                result.append(f"{parts[0]} {parts[1].upper()}")
            else:
                result.append(t.upper())
        return " + ".join(result)

    return f"{process_side(left)} {arrow} {process_side(right)}"


def format_eng_number(val: float) -> str:
    """Format number in engineering notation, strip trailing zeros from mantissa."""
    if val == 0:
        return "0.0"
    s = f"{val:.15e}"
    if "e" in s.lower():
        idx = s.lower().index("e")
        mantissa = s[:idx].strip()
        exp_part = s[idx:]
        try:
            exp_val = int(s[idx + 1 :].replace("+", ""))
            if abs(exp_val) < 100 and "E" not in exp_part:
                s = f"{val:.15e}"
                if "e" in s:
                    idx = s.lower().index("e")
                    mantissa = s[:idx].strip()
                    exp_part = s[idx:]
        except ValueError:
            pass
        dot_pos = mantissa.find(".")
        if dot_pos >= 0 and dot_pos < len(mantissa) - 1:
            last_nonzero = dot_pos
            for i in range(len(mantissa) - 1, dot_pos, -1):
                if mantissa[i] != "0":
                    last_nonzero = i
                    break
            if last_nonzero < len(mantissa) - 1:
                if last_nonzero == dot_pos and len(mantissa) > dot_pos:
                    mantissa = mantissa[: dot_pos + 1] + "0"
                else:
                    mantissa = mantissa[: last_nonzero + 1]
        return mantissa + exp_part
    return s


def parse_rate_params(params_str: str) -> tuple[dict, bool]:
    """
    Parse rate parameters from reaction string.
    Returns (dict with a, n, E, ai, ni, Ei, is_pressure_dep), parse_success.
    For falloff: ai, ni, Ei = low-P; a, n, E = high-P.
    """
    result = {
        "a": 0.0,
        "n": 0.0,
        "E": 0.0,
        "ai": 0.0,
        "ni": 0.0,
        "Ei": 0.0,
        "is_pressure_dep": False,
    }
    found = {k: False for k in ["a", "n", "E", "ai", "ni", "Ei"]}

    remaining = params_str.strip()
    if "ai" in remaining or "Ai" in remaining:
        result["is_pressure_dep"] = True

    # Parse key=value pairs
    while remaining:
        eq_pos = remaining.find("=")
        if eq_pos < 0:
            break

        key_part = remaining[:eq_pos].strip()
        # Extract last word (letters only) as key
        key_match = re.search(r"([a-zA-Z]+)\s*$", key_part)
        key = key_match.group(1) if key_match else key_part.strip()

        remaining = remaining[eq_pos + 1 :].lstrip()

        # Find next key or end
        next_key = re.search(
            r"\b([a-zA-Z]+)\s*=", remaining
        )
        tab_pos = remaining.find("\t")

        if next_key and next_key.start() > 1:
            value = remaining[: next_key.start()].strip()
            remaining = remaining[next_key.start() :]
        elif tab_pos > 1:
            value = remaining[:tab_pos].strip()
            remaining = remaining[tab_pos + 1 :]
        else:
            value = remaining.strip()
            remaining = ""

        value = value.rstrip("\t ")

        try:
            if key in ("a", "A"):
                result["a"] = float(value)
                found["a"] = True
            elif key in ("ai", "Ai"):
                result["ai"] = float(value)
                found["ai"] = True
            elif key in ("n", "N"):
                result["n"] = float(value)
                found["n"] = True
            elif key in ("ni", "Ni"):
                result["ni"] = float(value)
                found["ni"] = True
            elif key in ("E",):
                result["E"] = float(value)
                found["E"] = True
            elif key in ("Ei",):
                result["Ei"] = float(value)
                found["Ei"] = True
        except ValueError:
            pass

    if result["is_pressure_dep"]:
        if not found["a"] and found["ai"]:
            result["a"] = result["ai"]
            found["a"] = True
        if not found["n"] and found["ni"]:
            result["n"] = result["ni"]
            found["n"] = True
        if not found["E"] and found["Ei"]:
            result["E"] = result["Ei"]
            found["E"] = True
        if not found["ai"]:
            result["ai"] = result["a"]
        if not found["ni"]:
            result["ni"] = result["n"]
        if not found["Ei"]:
            result["Ei"] = result["E"]

    parse_success = found["a"] and found["n"] and found["E"]
    return result, parse_success


def parse_reaction_to_store(line: str) -> dict:
    """Parse one reaction line. Extract base_id and suffix from reaction_id."""
    colon_pos = line.find(":")
    reaction_id = line[:colon_pos].strip()

    if len(reaction_id) >= 1:
        if reaction_id[-1] == "f":
            suffix = "f"
            base_id = reaction_id[:-1]
        elif reaction_id[-1] == "b":
            suffix = "b"
            base_id = reaction_id[:-1]
        else:
            suffix = " "
            base_id = reaction_id
    else:
        suffix = " "
        base_id = ""

    arrow_pos = line.find("->")
    if arrow_pos < 0:
        arrow_pos = line.find("<=>")

    reactants_str = line[colon_pos + 1 : arrow_pos].strip()
    brace_pos = line.find("{")
    if brace_pos >= 0:
        products_str = line[arrow_pos + 2 : brace_pos].strip()
        rate_params = line[brace_pos + 1 :].strip()
        eq_pos = rate_params.find("}")
        if eq_pos >= 0:
            rate_params = rate_params[:eq_pos]
    else:
        products_str = line[arrow_pos + 2 :].strip()
        rate_params = ""

    reactants_str = strip_tabs_and_trim(reactants_str)
    products_str = strip_tabs_and_trim(products_str)

    rate_dict, parse_success = parse_rate_params(rate_params)
    if not parse_success and rate_params:
        print(f"  Warning: Could not parse rate parameters for reaction")

    return {
        "id": reaction_id,
        "base_id": base_id,
        "suffix": suffix,
        "reactants": reactants_str,
        "products": products_str,
        **rate_dict,
    }


def write_reaction_yaml(
    f,
    rxn_num: int,
    reactants: str,
    products: str,
    is_reversible: bool,
    a: float,
    n: float,
    E: float,
    ai: float,
    ni: float,
    Ei: float,
    is_pressure_dep: bool,
) -> None:
    """Write one reaction to YAML output. Species names are converted to uppercase."""
    Ea_cal = E * KJ_TO_CAL
    arrow = "<=>" if is_reversible else "=>"
    eq_str = uppercase_species_in_equation(f"{reactants} {arrow} {products}")
    f.write(f"- equation: {eq_str}  # Reaction {rxn_num}\n")

    if "M'" in reactants or "M'" in products:
        f.write("  type: three-body\n")
    elif is_pressure_dep:
        f.write("  type: falloff\n")

    if is_pressure_dep:
        Ea_low_cal = Ei * KJ_TO_CAL
        f.write(
            f"  high-P-rate-constant: {{A: {format_eng_number(a)}, "
            f"b: {format_eng_number(n)}, Ea: {format_eng_number(Ea_cal)}}}\n"
        )
        f.write(
            f"  low-P-rate-constant: {{A: {format_eng_number(ai)}, "
            f"b: {format_eng_number(ni)}, Ea: {format_eng_number(Ea_low_cal)}}}\n"
        )
    else:
        f.write(
            f"  rate-constant: {{A: {format_eng_number(a)}, "
            f"b: {format_eng_number(n)}, Ea: {format_eng_number(Ea_cal)}}}\n"
        )


def parse_mech_file(mech_file: Path, f) -> None:
    """Parse FlameMaster mechanism file and write reactions to YAML."""
    max_reactions = 1000
    reactions = []

    with open(mech_file) as fp:
        lines = fp.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        trimmed = line.strip()
        i += 1

        if not trimmed or trimmed.startswith("#"):
            continue

        if ":" in trimmed and "->" in trimmed:
            full_line = trimmed
            while "{" in full_line and "}" not in full_line:
                if i < len(lines):
                    full_line += " " + lines[i].strip()
                    i += 1
                else:
                    break

            if len(reactions) >= max_reactions:
                print(f"ERROR: Too many reactions (max {max_reactions})")
                sys.exit(1)

            reactions.append(parse_reaction_to_store(full_line))

    f.write("reactions:\n")

    nreactions = 0
    for i, rxn in enumerate(reactions):
        base_id = rxn["base_id"]
        has_f = any(r["base_id"] == base_id and r["suffix"] == "f" for r in reactions)
        has_b = any(r["base_id"] == base_id and r["suffix"] == "b" for r in reactions)

        if rxn["suffix"] == " ":
            nreactions += 1
            write_reaction_yaml(
                f,
                nreactions,
                rxn["reactants"],
                rxn["products"],
                False,
                rxn["a"],
                rxn["n"],
                rxn["E"],
                rxn["ai"],
                rxn["ni"],
                rxn["Ei"],
                rxn["is_pressure_dep"],
            )
        elif rxn["suffix"] == "f":
            nreactions += 1
            write_reaction_yaml(
                f,
                nreactions,
                rxn["reactants"],
                rxn["products"],
                not has_b,
                rxn["a"],
                rxn["n"],
                rxn["E"],
                rxn["ai"],
                rxn["ni"],
                rxn["Ei"],
                rxn["is_pressure_dep"],
            )
        elif rxn["suffix"] == "b":
            nreactions += 1
            write_reaction_yaml(
                f,
                nreactions,
                rxn["reactants"],
                rxn["products"],
                False,
                rxn["a"],
                rxn["n"],
                rxn["E"],
                rxn["ai"],
                rxn["ni"],
                rxn["Ei"],
                rxn["is_pressure_dep"],
            )

    print(f"  Parsed {len(reactions)} reactions -> {nreactions} YAML reactions")


def parse_transport_file(trans_file: Path) -> tuple[list[str], list[float], list[float]]:
    """Parse transport file. Returns (species_names, eps_over_k, sigma)."""
    species_names = []
    eps_over_k = []
    sigma = []

    if not trans_file.exists():
        return species_names, eps_over_k, sigma

    with open(trans_file) as fp:
        lines = fp.readlines()

    skip_header = False
    for line in lines:
        if not line.strip() or line.strip().startswith("#"):
            continue
        if "SpecName" in line or "specname" in line:
            skip_header = True
            continue

        if len(line.strip()) < 30:
            continue

        parts = line.split()
        if not parts:
            continue
        species_name = parts[0]
        if not re.search(r"[a-zA-Z]", species_name):
            continue

        try:
            if len(parts) >= 4:
                eps_val = float(parts[2])
                sigma_val = float(parts[3])
            elif len(parts) >= 3:
                eps_val = float(parts[1])
                sigma_val = float(parts[2])
            else:
                continue
        except (ValueError, IndexError):
            try:
                eps_val = float(parts[1])
                sigma_val = float(parts[2])
            except (ValueError, IndexError):
                continue

        species_names.append(species_name)
        eps_over_k.append(eps_val)
        sigma.append(sigma_val)

    print(f"  Parsed {len(species_names)} transport entries")
    return species_names, eps_over_k, sigma


def parse_thermo_composition(name_line: str) -> tuple[dict, bool]:
    """
    Parse composition from thermo file line (fixed width format).
    Format: columns 29+ have N_count O, O_count H, H_count C, C_count G
    Returns (dict with N, O, H, C, Ar), success.
    """
    comp = {"N": 0, "O": 0, "H": 0, "C": 0, "Ar": 0}
    if len(name_line) < 28:
        return comp, False

    pos = 28
    name_line = name_line.ljust(50)

    def parse_count():
        nonlocal pos
        count = 0
        if pos < len(name_line) and name_line[pos].isdigit():
            count = int(name_line[pos])
            pos += 1
            if pos < len(name_line) and name_line[pos].isdigit():
                count = count * 10 + int(name_line[pos])
                pos += 1
        pos += 1
        return count

    try:
        comp["N"] = parse_count()
        while pos < min(45, len(name_line)) and name_line[pos] == " ":
            pos += 1
        comp["O"] = parse_count()
        while pos < min(45, len(name_line)) and name_line[pos] == " ":
            pos += 1
        comp["H"] = parse_count()
        while pos < min(45, len(name_line)) and name_line[pos] == " ":
            pos += 1
        comp["C"] = parse_count()
        return comp, True
    except (ValueError, IndexError):
        return comp, False


def parse_composition_from_name(species_name: str) -> dict:
    """Parse composition from species name (e.g. N2, H2O, N-C7H16)."""
    comp = {"N": 0, "O": 0, "H": 0, "C": 0, "Ar": 0}
    i = 0
    name_len = len(species_name)

    while i < name_len:
        while i < name_len and not species_name[i].isalpha():
            i += 1
        if i >= name_len:
            break

        if i + 1 < name_len and species_name[i : i + 2] in ("Ar", "AR"):
            i += 2
            num = 0
            while i < name_len and species_name[i].isdigit():
                num = num * 10 + int(species_name[i])
                i += 1
            comp["Ar"] += num or 1
        else:
            elem = species_name[i]
            i += 1
            num = 0
            while i < name_len and species_name[i].isdigit():
                num = num * 10 + int(species_name[i])
                i += 1
            num = num or 1
            if elem in ("N", "n"):
                comp["N"] += num
            elif elem in ("O", "o"):
                comp["O"] += num
            elif elem in ("H", "h"):
                comp["H"] += num
            elif elem in ("C", "c"):
                comp["C"] += num

    return comp


def parse_thermo_species(
    fp,
    name_line: str,
    f,
    sp_num: int,
    trans_names: list[str],
    trans_eps: list[float],
    trans_sigma: list[float],
    has_trans_data: bool,
) -> None:
    """Parse a single species from thermo file. Species names are written in uppercase."""
    species_name = name_line[:18].strip().upper()

    comp, has_composition = parse_thermo_composition(name_line)
    if not has_composition:
        comp = {"N": 0, "O": 0, "H": 0, "C": 0, "Ar": 0}
        comp = parse_composition_from_name(species_name)
        has_composition = any(comp.values())

    if len(name_line) >= 77:
        try:
            T_low = float(name_line[44:55])
            T_high = float(name_line[55:66])
            T_mid = float(name_line[66:77])
        except ValueError:
            T_low, T_mid, T_high = 300.0, 1000.0, 5000.0
    else:
        T_low, T_mid, T_high = 300.0, 1000.0, 5000.0

    low_coeffs = [0.0] * 7
    high_coeffs = [0.0] * 7

    for i in range(3):
        coeff_line = fp.readline()
        if not coeff_line:
            break
        coeff_line = coeff_line.rstrip("\n").ljust(105)
        def safe_float(s: str) -> float:
            s = s.strip()
            return float(s) if s else 0.0

        try:
            if i == 0:
                for j in range(5):
                    low_coeffs[j] = safe_float(coeff_line[j * 15 : (j + 1) * 15])
            elif i == 1:
                # Line 2: low[5], low[6], high[0], high[1], high[2] (5 values; last 2 cols are line number)
                low_coeffs[5] = safe_float(coeff_line[0:15])
                low_coeffs[6] = safe_float(coeff_line[15:30])
                high_coeffs[0] = safe_float(coeff_line[30:45])
                high_coeffs[1] = safe_float(coeff_line[45:60])
                high_coeffs[2] = safe_float(coeff_line[60:75])
            elif i == 2:
                # Line 3: high[3], high[4], high[5], high[6] (4 values)
                high_coeffs[3] = safe_float(coeff_line[0:15])
                high_coeffs[4] = safe_float(coeff_line[15:30])
                high_coeffs[5] = safe_float(coeff_line[30:45])
                high_coeffs[6] = safe_float(coeff_line[45:60])
        except (ValueError, IndexError) as e:
            if sp_num <= 1:
                print(f"  Warning: Thermo parse error for {species_name}: {e}")
            low_coeffs = [0.0] * 7
            high_coeffs = [0.0] * 7

    f.write(f"- name: {species_name}\n")

    if has_composition:
        comp_parts = [f"{k}: {v}" for k, v in comp.items() if v > 0]
        f.write(f"  composition: {{{', '.join(comp_parts)}}}\n")
    else:
        f.write("  composition: {}\n")

    f.write("  thermo:\n")
    f.write("    model: NASA7\n")
    f.write(f"    temperature-ranges: [{T_low}, {T_mid}, {T_high}]\n")
    f.write("    data:\n")
    low_str = ", ".join(format_eng_number(c) for c in low_coeffs)
    high_str = ", ".join(format_eng_number(c) for c in high_coeffs)
    f.write(f"    - [{low_str}]\n")
    f.write(f"    - [{high_str}]\n")

    if has_trans_data:
        found_trans = False
        for i, tname in enumerate(trans_names):
            if tname.strip().upper() == species_name.strip().upper():
                found_trans = True
                total_atoms = sum(comp.values())
                if total_atoms == 1:
                    geometry = "atom"
                elif total_atoms == 2:
                    geometry = "linear"
                elif total_atoms == 3:
                    if (
                        (comp["C"] == 1 and comp["O"] == 2)
                        or (comp["N"] == 2 and comp["O"] == 1)
                        or (comp["C"] == 2 and comp["H"] == 1)
                        or (comp["C"] == 1 and comp["H"] == 2)
                        or (comp["H"] == 1 and comp["C"] == 1 and comp["N"] == 1)
                        or (comp["N"] == 1 and comp["C"] == 1 and comp["O"] == 1)
                    ):
                        geometry = "linear"
                    else:
                        geometry = "nonlinear"
                else:
                    geometry = "nonlinear"

                f.write("  transport:\n")
                f.write("    model: gas\n")
                f.write(f"    geometry: {geometry}\n")
                f.write(f"    diameter: {format_eng_number(trans_sigma[i])}\n")
                f.write(f"    well-depth: {format_eng_number(trans_eps[i])}\n")
                break

        if not found_trans:
            print(f'  Warning: No transport data found for species "{species_name}"')


def parse_thermo_file(
    thermo_file: Path,
    f,
    trans_names: list[str],
    trans_eps: list[float],
    trans_sigma: list[float],
    has_trans_data: bool,
) -> None:
    """Parse FlameMaster thermo file."""
    if not thermo_file.exists():
        print(f"  Warning: Cannot open thermo file: {thermo_file}")
        return

    with open(thermo_file) as fp:
        fp.readline()
        fp.readline()

        nspecies = 0
        while True:
            line = fp.readline()
            if not line:
                break
            if "G" in line or "L" in line or "S" in line:
                nspecies += 1
                parse_thermo_species(
                    fp,
                    line,
                    f,
                    nspecies,
                    trans_names,
                    trans_eps,
                    trans_sigma,
                    has_trans_data,
                )

    print(f"  Parsed {nspecies} species")


def get_date_string() -> str:
    """Return current date string for YAML metadata."""
    now = datetime.now()
    return now.strftime("%a %d %b %Y %H:%M:%S")


def convert_mechanism(
    mech_file: Path,
    thermo_file: Path,
    trans_file: Path,
    yaml_file: Path,
    has_thermo: bool,
    has_trans: bool,
) -> None:
    """Convert mechanism, thermo, and transport files to YAML."""
    trans_names, trans_eps, trans_sigma = [], [], []
    if has_trans and trans_file.exists():
        trans_names, trans_eps, trans_sigma = parse_transport_file(trans_file)
    has_trans_data = len(trans_names) > 0

    with open(yaml_file, "w") as f:
        f.write("description: |-\n")
        f.write("  Converted from FlameMaster format\n")
        f.write(f"  Mechanism file: {mech_file}\n")
        if has_thermo:
            f.write(f"  Thermo file: {thermo_file}\n")
        if has_trans:
            f.write(f"  Transport file: {trans_file}\n")
        f.write("\n")
        f.write("generator: FM2yaml\n")
        f.write(f"date: {get_date_string()}\n")
        f.write(f"input-files: [{mech_file}\n")
        if has_thermo:
            f.write(f"  , {thermo_file}\n")
        if has_trans:
            f.write(f"  , {trans_file}\n")
        f.write("]\n")
        f.write("\n")
        f.write("units: {length: cm, quantity: mol, activation-energy: cal/mol}\n")
        f.write("\n")
        f.write("phases:\n")
        f.write("- name: gas\n")
        f.write("  thermo: ideal-gas\n")
        f.write("  elements: [O, H, C, N, Ar]  # TODO: parse from mechanism\n")
        f.write("  species: all\n")
        f.write("  kinetics: gas\n")
        f.write("  reactions: all\n")
        if has_trans:
            f.write("  transport: mixture-averaged\n")
        f.write("  state:\n")
        f.write("    T: 300.0\n")
        f.write("    P: 1.01325e+05\n")
        f.write("\n")

        f.write("species:\n")
        if has_thermo:
            parse_thermo_file(
                thermo_file, f, trans_names, trans_eps, trans_sigma, has_trans_data
            )

        f.write("\n")
        parse_mech_file(mech_file, f)


def main():
    parser = argparse.ArgumentParser(
        description="FM2yaml - FlameMaster to YAML Converter"
    )
    parser.add_argument(
        "mech_file",
        type=Path,
        help="Mechanism file (.mech)",
    )
    parser.add_argument(
        "thermo_file",
        type=Path,
        nargs="?",
        default=None,
        help="Thermo file (.thermo). Default: <mech_base>.thermo",
    )
    parser.add_argument(
        "trans_file",
        type=Path,
        nargs="?",
        default=None,
        help="Transport file (.trans). Default: <mech_base>.trans",
    )
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        default=None,
        help="Output YAML file. Default: <mech_base>.yaml",
    )
    args = parser.parse_args()

    mech_file = args.mech_file
    if not mech_file.exists():
        print(f"ERROR: Mechanism file not found: {mech_file}")
        sys.exit(1)

    base = mech_file.parent / mech_file.stem

    thermo_file = args.thermo_file or base.with_suffix(".thermo")
    trans_file = args.trans_file or base.with_suffix(".trans")
    yaml_file = args.output or base.with_suffix(".yaml")

    has_thermo = thermo_file.exists()
    has_trans = trans_file.exists()

    print("FM2yaml - FlameMaster to YAML Converter")
    print("======================================")
    print(f"Mechanism file: {mech_file}")
    print(f"Thermo file:    {thermo_file}" if has_thermo else "Thermo file:    (not found, skipping)")
    print(f"Transport file: {trans_file}" if has_trans else "Transport file: (not found, skipping)")
    print(f"Output file:    {yaml_file}")
    print()

    convert_mechanism(mech_file, thermo_file, trans_file, yaml_file, has_thermo, has_trans)

    print("Conversion complete!")


if __name__ == "__main__":
    main()
