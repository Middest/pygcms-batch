"""
Py-GC-MS Raw-to-Analysis Pipeline
==================================
Complete workflow from raw TXT/QGD files to analysis-ready data.

Steps:
  1. PARSE    - Read NIST export TXT files
  2. FILTER   - Remove contaminants, SI threshold
  3. ALIGN    - Cross-treatment RT alignment
  4. RESOLVE  - ID conflicts via EI spectrum or corrections file
  5. CLASSIFY - Shahriar 2026 12-class scheme
  6. VALIDATE - Mass balance, SI stats, batch check
  7. EXPORT   - Clean matrix CSV + verification report

Usage:
  python pipeline.py --input <TXT_dir> --output <output_dir>
       [--sample_map <mapping.json>]
       [--corrections <corrections.json>]
       [--si_threshold 70] [--keep_contaminants]
       [--keep_tmah] [--no_renormalize]

Cleaning notes:
  - TMAH thermochemolysis artifacts (trimethylamine, trimethyltriazine,
    tetramethylmethanediamine, methenamine, etc.) are removed by default.
    These are reagent-derived N compounds, NOT soil N — they inflate Other_N
    and the R_MP (microbial/plant) ratio. Use --keep_tmah to keep them as a
    control. See TMAH_ARTIFACTS below.
  - Class composition is renormalized to 100% of kept-peak total per treatment
    by default. Use --no_renormalize for raw Conc% sums.
"""
import os, sys, re, json, csv, argparse
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime

# ============================================================
# CONFIGURATION
# ============================================================

# Known contaminants to flag (not necessarily remove)
CONTAMINANTS = [
    "silane", "siloxane", "phthalate", "phthalic",
    "trichlorodocosyl", "pentafluoropropionic",
    "bis(2-ethylhexyl)", "1,4-benzenedicarboxylic",
    "column bleed", "septum",
]

# TMAH (tetramethylammonium hydroxide) thermochemolysis artifacts.
# At 500 °C TMAH undergoes Hofmann elimination to trimethylamine, which reacts
# with formaldehyde (from carbohydrates) to form N-methyl condensation products.
# These are REAGENT chemistry, not soil-derived N — they must be removed before
# interpreting Other_N / N-MAH as microbial N signal.
# Reference: TMAH-Py-GC-MS thermochemolysis artifacts (e.g., Fabbri et al.;
# del Río & Hatcher; Estournel-Pelardy et al.). Case-insensitive substrings.
TMAH_ARTIFACTS = [
    "methylamine, n,n-dimethyl",          # trimethylamine (TMAH Hofmann elimination)
    "methanediamine",                     # N,N,N',N'-tetramethylmethanediamine
    "triazine, hexahydro",                # trimethyltriazine / triazine dione (TMA+formaldehyde cyclization)
    "methenamine",                        # hexamethylenetetramine
    "acetonitrile, (dimethylamino)",      # dimethylaminoacetonitrile
    # Spectrum-confirmed reagent peaks (2026-08-01, paddy soil). NIST may name
    # the trimethylamine peak / TMAH cation anything; these names were observed.
    # Robust detection is spectral (base m/z 58) — see SKILL.md "Spectral caveat".
    "tetramethylammonium",                # TMAH reagent cation (e.g. "...acetate")
    "butanoic acid, 4-(dimethylamino)",   # trimethylamine misID (POC BC15)
    "1-methyldodecylamine",               # trimethylamine misID (MAOC BC7.5; RT2.47 impossible for C13 amine)
]

# GC derivatization reagent artifacts (acylation/silylation): fluoro-acyl esters,
# trifluoroacetates, TMS derivatives. These come from TFAA/PFPA/HFBA/BSTFA-type
# reagents reacting with alcohols/amines, NOT from soil OM. Same status as TMAH
# artifacts — reagent chemistry to remove before class interpretation.
DERIVATIZATION_ARTIFACTS = [
    "pentafluoroprop",                    # pentafluoropropionic acid esters
    "trifluoroacetate", "trifluoroacet",  # trifluoroacetate esters
    "heptafluorobuty",                    # heptafluorobutyric acid esters
    "pentafluorooctano", "perfluorooctanoic",
    "trimethylsilyl", "tert-butyldimethylsilyl", "tms derivative",
    "bis(trimethylsilyl)", "silyl",
]


def is_tmah_artifact(name):
    """Return the matching reagent-artifact keyword, or None.

    Checks both TMAH thermochemolysis artifacts and GC derivatization
    reagent artifacts (fluoro-acyl esters / TMS). Both are reagent
    chemistry, not soil-derived compounds.
    """
    if not name:
        return None
    n = name.lower()
    for kw in TMAH_ARTIFACTS:
        if kw in n:
            return f"tmah:{kw}"
    for kw in DERIVATIZATION_ARTIFACTS:
        if kw in n:
            return f"deriv:{kw}"
    return None

# Shahriar 2026 class map (full name -> short code)
CLASS_MAP = {
    "Monocyclic aromatic hydrocarbons (MAH)": "MAH",
    "N-containing monocyclic aromatic hydrocarbons (N-MAH)": "N-MAH",
    "Polycyclic aromatic hydrocarbons (PAH)": "PAH",
    "Lignin": "Lignin", "Phenols": "Phenols",
    "Degraded saccharides": "Sugars",
    "Fatty acids, alcohols and esters   ": "Fatty_acids_lipids",
    "Fatty acids, alcohols and esters": "Fatty_acids_lipids",
    "Alkenes                            ": "Alkenes", "Alkenes": "Alkenes",
    "Short alkanes                      ": "Short_alkanes", "Short alkanes": "Short_alkanes",
    "Long alkanes                       ": "Long_alkanes", "Long alkanes": "Long_alkanes",
    "Other hydrocarbons                 ": "Other_hydrocarbons",
    "Other hydrocarbons": "Other_hydrocarbons",
    "Other N-containing compounds (N-containing)": "Other_N",
    "Other N-containing compounds": "Other_N",
}

PLANT = {"Lignin", "Long_alkanes", "Alkenes"}
MICROBIAL = {"N-MAH", "Other_N"}


# ============================================================
# STEP 1: PARSE
# ============================================================

def parse_txt(filepath):
    """Parse Shimadzu GCMSsolution NIST export TXT file.

    Returns list of dicts with keys: rt, area, conc, height, name, si, mark
    """
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = [l.strip("﻿").rstrip("\n") for l in f.readlines()]

    peaks, sim = [], {}
    section, peak_start = None, False

    for line in lines:
        parts = line.split("\t")
        if "[MC Peak Table]" in line:
            section, peak_start = "mc", False; continue
        elif "[MS Similarity Search Results" in line:
            section, peak_start = "sim", False; continue
        elif line.startswith("[") and line.endswith("]"):
            section, peak_start = None, False; continue

        if section == "mc":
            if "Peak#" in line and "Ret.Time" in line:
                peak_start = True; continue
            if peak_start and len(parts) >= 11:
                try:
                    peaks.append({
                        "rt": float(parts[1]),
                        "area": int(parts[5]) if parts[5].strip() else 0,
                        "height": int(parts[6]) if parts[6].strip() else 0,
                        "conc": float(parts[8]) if parts[8].strip() else 0.0,
                        "mark": parts[9].strip(),
                        "name": parts[10].strip() if len(parts) > 10 else "",
                    })
                except (ValueError, IndexError):
                    continue

        if section == "sim":
            if "Spectrum#" in line and "SI" in line: continue
            if len(parts) >= 6:
                try:
                    if int(parts[1]) == 1:  # only best hit
                        sim[int(parts[0])] = {
                            "si": int(parts[2]), "cas": parts[3].strip(),
                            "name": parts[4].strip(),
                        }
                except (ValueError, IndexError):
                    continue

    # Merge SI into peaks
    for i, p in enumerate(peaks):
        match = sim.get(i + 1, {})
        p["si"] = match.get("si", 0)
        if not p["name"] or len(p["name"]) < 3:
            p["name"] = match.get("name", "")

    return peaks


# ============================================================
# STEP 2: FILTER
# ============================================================

def filter_peaks(peaks, si_threshold=70, remove_contaminants=True, remove_tmah=True):
    """Apply quality filters to peak list.

    Returns (kept_peaks, removed_peaks) with filtering reasons.
    """
    kept, removed = [], []

    for p in peaks:
        reasons = []

        # SI check
        if p["si"] > 0 and p["si"] < si_threshold:
            reasons.append(f"SI={p['si']}<{si_threshold}")

        # Contaminant check
        is_contam = False
        if remove_contaminants:
            name_lower = p["name"].lower()
            for kw in CONTAMINANTS:
                if kw in name_lower:
                    reasons.append(f"contaminant:{kw}")
                    is_contam = True
                    break

        # TMAH thermochemolysis artifact check (reagent-derived, not soil N)
        if remove_tmah:
            kw = is_tmah_artifact(p["name"])
            if kw:
                reasons.append(f"reagent_artifact:{kw}")
                is_contam = True

        # Zero/negative area
        if p["area"] <= 0:
            reasons.append("zero_area")

        if reasons:
            p["filter_reason"] = "; ".join(reasons)
            p["filtered"] = True
            removed.append(p)
        else:
            p["filtered"] = False
            kept.append(p)

    return kept, removed


# ============================================================
# STEP 3: ALIGN
# ============================================================

def align_peaks(all_peaks, reference_idx=0, rt_tolerance=0.08):
    """Align peaks across treatments by retention time.

    Args:
        all_peaks: list of (treatment_name, peaks_list)
        reference_idx: index of reference treatment (default 0 = CK)
        rt_tolerance: max RT difference for matching

    Returns aligned list of dicts, one per aligned feature.
    """
    ref_name, ref_peaks = all_peaks[reference_idx]
    other_list = [(n, p) for i, (n, p) in enumerate(all_peaks) if i != reference_idx]

    matrix = []
    for ref_p in ref_peaks:
        row = {
            "rt_ref": ref_p["rt"],
            f"name_{ref_name}": ref_p["name"],
            f"conc_{ref_name}": ref_p["conc"],
            f"area_{ref_name}": ref_p["area"],
            f"si_{ref_name}": ref_p["si"],
            f"filtered_{ref_name}": ref_p.get("filtered", False),
        }
        conflicts = []

        for other_name, other_peaks in other_list:
            best, best_d = None, rt_tolerance + 0.01
            for op in other_peaks:
                d = abs(op["rt"] - ref_p["rt"])
                if d < best_d:
                    best, best_d = op, d

            if best:
                row[f"rt_{other_name}"] = best["rt"]
                row[f"delta_rt_{other_name}"] = round(best_d, 4)
                row[f"name_{other_name}"] = best["name"]
                row[f"conc_{other_name}"] = best["conc"]
                row[f"area_{other_name}"] = best["area"]
                row[f"si_{other_name}"] = best["si"]
                row[f"filtered_{other_name}"] = best.get("filtered", False)

                if best["name"] != ref_p["name"] and best["name"] and ref_p["name"]:
                    conflicts.append((other_name, ref_p["name"], best["name"]))
            else:
                for k in ["rt", "delta_rt", "name", "conc", "area", "si", "filtered"]:
                    row[f"{k}_{other_name}"] = None

        row["id_conflicts"] = conflicts
        row["has_conflict"] = len(conflicts) > 0
        matrix.append(row)

    return matrix


# ============================================================
# STEP 4: RESOLVE
# ============================================================

def load_shahriar_library():
    """Load Shahriar 2026 compound classification library."""
    lib_path = os.path.join(os.path.dirname(__file__), "..", "data", "shahriar_library.json")
    with open(lib_path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return {k.lower().strip(): CLASS_MAP.get(v, v) for k, v in raw.items()}


def classify_compound(name, library):
    """Classify using Shahriar library with fuzzy fallback."""
    if not name:
        return "Unknown"

    n = name.lower().strip()

    # 1. Exact library match
    if n in library:
        return library[n]

    # 2. Remove CAS/suffix
    base = n.split("$$")[0].strip()
    if base in library:
        return library[base]

    # 3. First part before comma
    first = re.split(r"[,;]", n)[0].strip()
    if first in library:
        return library[first]

    # 4. Remove stereochemistry
    clean = re.sub(r"\(r\)|\(s\)|\(e\)|\(z\)|\[r-\]|\[s-\]", "", first, flags=re.I)
    clean = re.sub(r"\s+", " ", clean).strip()
    if clean in library:
        return library[clean]

    # 5. Keyword fallback
    return _keyword_classify(n)


def _keyword_classify(n):
    """Keyword-based classification. Order matters: N-containing BEFORE hydrocarbons."""
    # N-containing (must precede alkanes/alkenes)
    if any(k in n for k in ["nitrile", "cyanide", "butanenitrile", "pentanenitrile",
            "heptadecanenitrile", "pentadecanenitrile", "octadecenenitrile",
            "dodecanenitrile", "tetradecanenitrile"]):
        return "Other_N"
    if any(k in n for k in ["pyrrole", "pyridine", "indole", "trimethylindole",
            "tetramethylindole", "pyrrolidine", "piperidinamine", "methylindole"]):
        return "N-MAH"
    if any(k in n for k in ["amine", "amide", "methanediamine", "triazine",
            "carbamate", "tetramethylammonium", "tetramethyl-"]):
        return "Other_N"

    # PAH (before MAH)
    if any(k in n for k in ["naphthalene", "fluorene", "anthracene", "pyrene",
            "fluoranthene", "dibenzofuran"]):
        return "PAH"

    # Specific markers
    if any(k in n for k in ["furfural", "levoglucosan", "furan", "methylfurfural"]):
        return "Sugars"
    if any(k in n for k in ["guaiacol", "syringol", "vanillin", "eugenol", "cinnamyl"]):
        return "Lignin"
    if any(k in n for k in ["phenol", "cresol", "di-tert-butylphenol"]):
        return "Phenols"

    # MAH
    if any(k in n for k in ["benzene", "toluene", "xylene", "ethylbenzene",
            "mesitylene", "cymene", "indene", "indane", "styrene", "phenyl"]):
        return "MAH"

    # Lipids
    if any(k in n for k in ["dodecanol", "undecanol", "tetradecanol", "hexadecanol",
            "alcohol", "methyl ester", "ethyl ester", "acetic acid",
            "hexadecanoic", "octadecanoic", "dodecanoic", "tetradecanoic",
            "tridecanoic", "octadecenoic", "succinic"]):
        return "Fatty_acids_lipids"

    # Alkenes
    if any(k in n for k in ["alkene", "octene", "nonene", "decene", "undecene",
            "dodecene", "tridecene", "tetradecene", "pentadecene", "hexadecene",
            "heptadecene", "octadecene", "nonadecene", "eicosene",
            "cyclooctatetraene", "diene", "triene", "cyclopentadiene"]):
        return "Alkenes"

    # Long alkanes (C19+)
    if any(k in n for k in ["nonadecane", "eicosane", "heneicosane", "docosane",
            "tricosane", "tetracosane", "pentacosane", "hexacosane", "heptacosane",
            "octacosane", "nonacosane", "triacontane", "tetratriacontane",
            "pentatriacontane", "hexatriacontane", "tetratetracontane",
            "tritetracontane", "squalane"]):
        return "Long_alkanes"

    # Short alkanes
    if any(k in n for k in ["nonane", "decane", "undecane", "dodecane", "tridecane",
            "tetradecane", "pentadecane", "hexadecane", "heptadecane",
            "octadecane", "octane", "heptane", "hexane",
            "cyclopropane", "cyclohexane", "cyclododecane", "cycloundecane"]):
        return "Short_alkanes"

    # Other hydrocarbons
    if any(k in n for k in ["yne", "acetylene", "dioxane", "bicyclo", "spiro", "metheno"]):
        return "Other_hydrocarbons"

    return "Unknown"


def apply_corrections(peaks, corrections):
    """Apply manual compound name corrections.

    Args:
        peaks: list of peak dicts
        corrections: dict like {"3.215": "Toluene", "18.737": "Pentadecanenitrile"}
    """
    for rt_str, correct_name in corrections.items():
        rt_targ = float(rt_str)
        for p in peaks:
            if abs(p["rt"] - rt_targ) < 0.04:
                p["name_original"] = p.get("name_original", p["name"])
                p["name"] = correct_name
                p["corrected"] = True


def resolve_conflicts(matrix, library):
    """Flag and categorize identification conflicts."""
    for row in matrix:
        if not row["has_conflict"]:
            row["conflict_severity"] = "none"
            continue

        # Determine severity
        conflicts = row["id_conflicts"]
        ref_name = row.get("name_CK", row.get(f"name_{list(row.keys())[0].split('_')[1]}", ""))

        ref_cat = classify_compound(ref_name, library)
        cross_class = False
        for other_treat, _, other_name in conflicts:
            other_cat = classify_compound(other_name, library)
            if other_cat != ref_cat and other_cat != "Unknown" and ref_cat != "Unknown":
                cross_class = True
                break

        if cross_class:
            row["conflict_severity"] = "HIGH"
        else:
            row["conflict_severity"] = "LOW"


# ============================================================
# STEP 6: VALIDATE
# ============================================================

def validate(peaks_by_treatment, matrix, library, renormalize=True):
    """Run validation checks and return results dict.

    Args:
        renormalize: If True, class conc is renormalized to relative % of
            kept-peak total (sums to 100 per treatment). If False, raw
            Conc% sums are reported (do not sum to 100 when peaks filtered).
    """
    results = {}

    # Mass balance
    for tname, peaks in peaks_by_treatment:
        total = sum(p["conc"] for p in peaks if not p.get("filtered"))
        results[f"mass_balance_{tname}"] = round(total, 1)

    # SI statistics
    for tname, peaks in peaks_by_treatment:
        si_vals = [p["si"] for p in peaks if p["si"] > 0]
        si_ok = sum(1 for s in si_vals if s >= 80)
        si_hi = sum(1 for s in si_vals if s >= 90)
        si_lo = sum(1 for s in si_vals if s < 70)
        results[f"si_{tname}"] = {
            "total": len(si_vals), "ok": si_ok, "high": si_hi, "low": si_lo
        }

    # Conflict summary
    high = sum(1 for r in matrix if r.get("conflict_severity") == "HIGH")
    low = sum(1 for r in matrix if r.get("conflict_severity") == "LOW")
    results["conflicts"] = {"high": high, "low": low, "total": high + low}

    # Class composition per treatment
    for tname, peaks in peaks_by_treatment:
        cats = defaultdict(lambda: {"conc": 0.0, "count": 0})
        total_conc = 0.0
        for p in peaks:
            if p.get("filtered"): continue
            cat = classify_compound(p["name"], library)
            cats[cat]["conc"] += p["conc"]
            cats[cat]["count"] += 1
            total_conc += p["conc"]

        if renormalize and total_conc > 0:
            scale = 100.0 / total_conc
            for c in cats:
                cats[c]["conc"] *= scale

        plant = sum(cats[c]["conc"] for c in PLANT)
        microb = sum(cats[c]["conc"] for c in MICROBIAL)
        unk = cats["Unknown"]["conc"]
        r_mp = microb / plant if plant > 0 else 0
        results[f"source_{tname}"] = {
            "plant": round(plant, 1), "microb": round(microb, 1),
            "unknown": round(unk, 1), "R_MP": round(r_mp, 2),
            "renormalized": bool(renormalize),
            "classes": {c: {"conc": round(d["conc"], 1), "count": d["count"]}
                       for c, d in sorted(cats.items(), key=lambda x: x[1]["conc"], reverse=True)},
        }

    return results


# ============================================================
# STEP 7: EXPORT
# ============================================================

def export_results(matrix, validation, treatment_names, output_dir, library, peaks_by_treatment):
    """Export analysis-ready CSV and verification report."""
    os.makedirs(output_dir, exist_ok=True)

    # ---- CSV: analysis-ready matrix ----
    csv_path = os.path.join(output_dir, "analysis_ready_matrix.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        header = ["Feature_ID", "RT_ref", "Name_ref", "Class_ref", "Conc%_ref", "SI_ref", "MSI_Level"]
        for t in treatment_names[1:]:
            header += [f"RT_{t}", f"Name_{t}", f"Class_{t}", f"Conc%_{t}",
                       f"SI_{t}", f"DeltaRT_{t}"]
        header += ["Conflict_Severity", "QC_Flag"]
        writer.writerow(header)

        ref_name = treatment_names[0]
        for row in matrix:
            ref_cat = classify_compound(row.get(f"name_{ref_name}", ""), library)
            line = [
                row.get("feature_id", ""), row["rt_ref"],
                row.get(f"name_{ref_name}", ""), ref_cat,
                row.get(f"conc_{ref_name}", ""), row.get(f"si_{ref_name}", ""),
                row.get("msi_level", ""),
            ]
            for t in treatment_names[1:]:
                t_cat = classify_compound(row.get(f"name_{t}", ""), library)
                line += [
                    row.get(f"rt_{t}", ""), row.get(f"name_{t}", ""), t_cat,
                    row.get(f"conc_{t}", ""), row.get(f"si_{t}", ""),
                    row.get(f"delta_rt_{t}", ""),
                ]

            # QC flags
            flags = []
            if row.get("conflict_severity") == "HIGH":
                flags.append("CROSS_CLASS_CONFLICT")
            if row.get("conflict_severity") == "LOW":
                flags.append("naming_conflict")
            if any(row.get(f"filtered_{t}") for t in treatment_names):
                flags.append("filtered_in_some")

            line.append(row.get("conflict_severity", ""))
            line.append("; ".join(flags) if flags else "PASS")
            writer.writerow(line)

    print(f"  Matrix: {csv_path}")

    # ---- CSV: class composition ----
    comp_path = os.path.join(output_dir, "class_composition.csv")
    with open(comp_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        all_cats = set()
        for t in treatment_names:
            all_cats.update(validation[f"source_{t}"]["classes"].keys())
        writer.writerow(["Class", "Source"] + treatment_names)
        for cat in sorted(all_cats, key=lambda c:
                validation[f"source_{treatment_names[0]}"]["classes"].get(c, {}).get("conc", 0),
                reverse=True):
            src = "Plant" if cat in PLANT else ("Microbial" if cat in MICROBIAL else "Mixed")
            row = [cat, src]
            for t in treatment_names:
                row.append(validation[f"source_{t}"]["classes"].get(cat, {}).get("conc", 0))
            writer.writerow(row)

    print(f"  Composition: {comp_path}")

    # ---- MD: verification report ----
    report_path = os.path.join(output_dir, "verification_report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(f"# Py-GC-MS Data Quality Report\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")

        f.write("## 1. Filtering Summary\n\n")
        f.write("| Treatment | Raw Peaks | Kept | Filtered | SI>=80 | SI>=90 | SI<70 | TMAH_artifacts |\n")
        f.write("|-----------|-----------|------|----------|--------|--------|-------|----------------|\n")
        for tname, peaks in peaks_by_treatment:
            raw = len(peaks)
            kept = sum(1 for p in peaks if not p.get("filtered"))
            filt = raw - kept
            n_tmah = sum(1 for p in peaks if p.get("filtered") and "artifact" in p.get("filter_reason", ""))
            si = validation[f"si_{tname}"]
            f.write(f"| {tname} | {raw} | {kept} | {filt} | {si['ok']} | {si['high']} | {si['low']} | {n_tmah} |\n")

        f.write("\n_TMAH artifacts: reagent-derived thermochemolysis products (trimethylamine,\n")
        f.write("trimethyltriazine, tetramethylmethanediamine, methenamine, dimethylaminoacetonitrile).\n")
        f.write("Plus GC derivatization reagent artifacts (fluoro-acyl esters, TMS derivatives).\n")
        f.write("Removed before class composition — they are reagent chemistry, not soil OM._\n\n")

        renorm = all(validation.get(f"source_{t}", {}).get("renormalized") for t in treatment_names)
        f.write("_Class composition values are " +
                ("renormalized to 100% of kept-peak total." if renorm else "raw Conc% sums (not renormalized).") +
                "_\n\n")

        f.write("\n## 2. Compound Class Composition\n\n")
        f.write("| Class | Source | " + " | ".join(treatment_names) + " |\n")
        f.write("|-------|--------|" + "|".join(["-------"] * len(treatment_names)) + "|\n")
        all_cats = set()
        for t in treatment_names:
            all_cats.update(validation[f"source_{t}"]["classes"].keys())
        for cat in sorted(all_cats, key=lambda c:
                validation[f"source_{treatment_names[0]}"]["classes"].get(c, {}).get("conc", 0),
                reverse=True):
            src = "Plant" if cat in PLANT else ("Microb" if cat in MICROBIAL else "Mixed")
            vals = " | ".join(
                str(validation[f"source_{t}"]["classes"].get(cat, {}).get("conc", 0))
                for t in treatment_names)
            f.write(f"| {cat} | {src} | {vals} |\n")

        f.write("\n## 3. Source Attribution\n\n")
        f.write("| Treatment | Plant% | Microbial% | R_MP |\n")
        f.write("|-----------|--------|------------|------|\n")
        for t in treatment_names:
            s = validation[f"source_{t}"]
            f.write(f"| {t} | {s['plant']} | {s['microb']} | {s['R_MP']} |\n")
        f.write(f"\nR_MP = (N-MAH + Other_N) / (Lignin + Long_alkanes + Alkenes)\n")

        f.write(f"\n## 4. ID Conflicts\n\n")
        f.write(f"- HIGH (cross-class): {validation['conflicts']['high']} peaks\n")
        f.write(f"- LOW (naming only): {validation['conflicts']['low']} peaks\n")

        high_conflicts = [r for r in matrix if r.get("conflict_severity") == "HIGH"]
        if high_conflicts:
            f.write("\n### High-severity conflicts (cross-class)\n\n")
            f.write("| RT | Ref Name | Ref Class | Conflicting Treatments |\n")
            f.write("|----|----------|-----------|----------------------|\n")
            for row in sorted(high_conflicts, key=lambda r:
                    max(r.get(f"conc_{t}", 0) or 0 for t in treatment_names), reverse=True)[:20]:
                ref_name = row.get(f"name_{treatment_names[0]}", "")[:40]
                ref_cat = classify_compound(ref_name, library)
                items = []
                for t in treatment_names[1:]:
                    other_name = row.get(f"name_{t}", "")
                    if other_name and other_name != ref_name:
                        other_cat = classify_compound(other_name, library)
                        if other_cat != ref_cat:
                            items.append(f"{t}: {other_name[:25]}[{other_cat}]")
                f.write(f"| {row['rt_ref']:.3f} | {ref_name} | {ref_cat} | {'; '.join(items[:3])} |\n")

        f.write(f"\n## 5. Decision\n\n")
        issues = []
        if any(validation[f"mass_balance_{t}"] < 95 for t in treatment_names):
            issues.append("Mass balance < 95% in some samples")
        if validation["conflicts"]["high"] > 10:
            issues.append(f"{validation['conflicts']['high']} cross-class conflicts need manual review")
        for t in treatment_names:
            if validation[f"si_{t}"]["low"] > validation[f"si_{t}"]["total"] * 0.3:
                issues.append(f"{t}: >30% peaks with SI<70")

        if issues:
            f.write("**Issues requiring attention:**\n\n")
            for i in issues:
                f.write(f"- {i}\n")
        else:
            f.write("**Data quality: ACCEPTABLE for analysis.**\n\n")
            f.write("The filtered, aligned, and classified dataset is ready for use.\n")

    print(f"  Report: {report_path}")

    return csv_path, report_path


# ============================================================
# MAIN PIPELINE
# ============================================================

def run_pipeline(input_dir, output_dir, sample_map=None, corrections=None,
                 si_threshold=70, remove_contaminants=True, remove_tmah=True,
                 renormalize=True, reference_treatment="CK"):
    """Execute the complete raw-to-analysis pipeline.

    Args:
        input_dir: Directory with NIST export .TXT files
        output_dir: Output directory for results
        sample_map: Dict mapping filename prefix -> treatment name
        corrections: Dict mapping treatment -> {RT: correct_name}
        si_threshold: Minimum SI for retention (0 = keep all)
        remove_contaminants: Whether to flag/remove known contaminants
        remove_tmah: Whether to remove TMAH thermochemolysis artifacts
        renormalize: Whether to renormalize class composition to 100%
        reference_treatment: Treatment to use as RT alignment reference

    Returns:
        (matrix_path, report_path)
    """
    # Default sample map
    if sample_map is None:
        sample_map = {"5": "CK", "6": "BC7.5", "7": "BC15", "8": "BC30"}
    if corrections is None:
        corrections = {}

    library = load_shahriar_library()
    print(f"Loaded Shahriar 2026 library: {len(library)} compounds\n")

    # ---- Step 1: Parse ----
    print("[1/7] Parsing TXT files...")
    peaks_by_treatment = []
    for prefix in sorted(sample_map.keys()):
        path = os.path.join(input_dir, f"{prefix}.txt")
        if not os.path.exists(path):
            print(f"  SKIP: {path} not found")
            continue
        peaks = parse_txt(path)
        tname = sample_map[prefix]
        print(f"  {tname}: {len(peaks)} peaks")

        # Apply corrections
        if tname in corrections:
            apply_corrections(peaks, corrections[tname])
            n_corr = sum(1 for p in peaks if p.get("corrected"))
            if n_corr:
                print(f"    {n_corr} peaks corrected")

        peaks_by_treatment.append((tname, peaks))

    treatment_names = [t for t, _ in peaks_by_treatment]
    print(f"  Treatments: {', '.join(treatment_names)}")

    # ---- Step 2: Filter ----
    print(f"\n[2/7] Filtering (SI>={si_threshold}, contaminants={'removed' if remove_contaminants else 'kept'}, "
          f"TMAH artifacts={'removed' if remove_tmah else 'kept'})...")
    for i, (tname, peaks) in enumerate(peaks_by_treatment):
        kept, removed = filter_peaks(peaks, si_threshold, remove_contaminants, remove_tmah)
        peaks_by_treatment[i] = (tname, kept + removed)  # keep all, just flag
        n_filt = len(removed)
        print(f"  {tname}: {n_filt} peaks flagged ({len(kept)} pass)")

    # ---- Step 3: Align ----
    print(f"\n[3/7] RT alignment (reference={reference_treatment}, tolerance=0.08 min)...")
    ref_idx = next((i for i, (n, _) in enumerate(peaks_by_treatment) if n == reference_treatment), 0)
    matrix = align_peaks(peaks_by_treatment, ref_idx)
    print(f"  {len(matrix)} aligned features")

    # ---- Step 4: Resolve ----
    print("\n[4/7] Resolving identification conflicts...")
    resolve_conflicts(matrix, library)
    high = sum(1 for r in matrix if r.get("conflict_severity") == "HIGH")
    low = sum(1 for r in matrix if r.get("conflict_severity") == "LOW")
    print(f"  HIGH (cross-class): {high}, LOW (naming): {low}")

    if high > 0:
        print("  Top cross-class conflicts:")
        for row in sorted(matrix, key=lambda r:
                max(r.get(f"conc_{t}", 0) or 0 for t in treatment_names), reverse=True):
            if row.get("conflict_severity") == "HIGH":
                ref_name = row.get(f"name_{treatment_names[0]}", "")[:35]
                ref_cat = classify_compound(ref_name, library)
                print(f"    RT {row['rt_ref']:.3f}: {ref_name} [{ref_cat}]")
                for t in treatment_names[1:]:
                    oname = row.get(f"name_{t}", "")
                    if oname and oname != ref_name:
                        ocat = classify_compound(oname, library)
                        if ocat != ref_cat:
                            print(f"      vs {t}: {oname[:35]} [{ocat}]")
                high -= 1
                if high == 0: break

    # ---- Step 5: Feature IDs + MSI levels ----
    print("\n[5/7] Assigning Feature IDs and MSI confidence levels...")
    matrix = assign_feature_ids(matrix)

    # Add classification and MSI level to each row
    for row in matrix:
        for t in treatment_names:
            name = row.get(f"name_{t}", "")
            row[f"class_{t}"] = classify_compound(name, library)
        row["msi_level"] = assign_msi_level(row, library)

    levels = Counter(r["msi_level"] for r in matrix)
    print(f"  Level_2: {levels.get('Level_2',0)}, Level_3: {levels.get('Level_3',0)}, Level_4: {levels.get('Level_4',0)}")

    # Closed-sum detection
    closed_sum = detect_closed_sum_effect(matrix, treatment_names)
    if closed_sum:
        print(f"  Closed-sum warnings: {len(closed_sum)} peaks >15% in at least one treatment")
        for w in closed_sum[:5]:
            print(f"    {w['treatment']} {w['feature_id']} ({w['name'][:30]}): {w['conc']:.1f}%")

    # ---- Step 6: Classify (already done inline, summarize) ----
    print("\n[6/7] Classifying compounds (Shahriar 2026)...")
    for tname, peaks in peaks_by_treatment:
        cats = Counter()
        for p in peaks:
            if not p.get("filtered"):
                cats[classify_compound(p["name"], library)] += 1
        total = sum(cats.values())
        top3 = cats.most_common(3)
        print(f"  {tname}: {total} classified, top: {', '.join(f'{c}({n})' for c,n in top3)}")

    # ---- Step 7: Validate ----
    print("\n[7/8] Validating...")
    validation = validate(peaks_by_treatment, matrix, library, renormalize=renormalize)

    for t in treatment_names:
        s = validation[f"source_{t}"]
        print(f"  {t}: Plant={s['plant']:.1f}%, Microb={s['microb']:.1f}%, R_MP={s['R_MP']:.2f}")

    # ---- Step 8: Export ----
    print("\n[8/8] Exporting results...")
    # Generate filter log
    generate_filter_log(peaks_by_treatment, output_dir)
    csv_path, report_path = export_results(
        matrix, validation, treatment_names, output_dir, library, peaks_by_treatment
    )

    print(f"\n{'='*60}")
    print(f"PIPELINE COMPLETE")
    print(f"  Analysis-ready data: {csv_path}")
    print(f"  Quality report: {report_path}")
    print(f"{'='*60}")

    return csv_path, report_path


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Py-GC-MS Raw-to-Analysis Pipeline")
    parser.add_argument("--input", required=True,
                        help="Directory with NIST export .TXT files")
    parser.add_argument("--output", required=True,
                        help="Output directory for results")
    parser.add_argument("--sample_map",
                        help="JSON: {'5':'CK','6':'BC7.5',...}")
    parser.add_argument("--corrections",
                        help="JSON: {'BC15':{'3.215':'Toluene'}}")
    parser.add_argument("--si_threshold", type=int, default=70,
                        help="Minimum SI (default: 70)")
    parser.add_argument("--keep_contaminants", action="store_true",
                        help="Keep known contaminants in output")
    parser.add_argument("--keep_tmah", action="store_true",
                        help="Keep TMAH thermochemolysis artifacts in output")
    parser.add_argument("--no_renormalize", action="store_true",
                        help="Report raw Conc%% sums instead of renormalizing to 100%%")
    parser.add_argument("--reference", default="CK",
                        help="Reference treatment for RT alignment")
    args = parser.parse_args()

    # Load sample map
    sample_map = {"5": "CK", "6": "BC7.5", "7": "BC15", "8": "BC30"}
    if args.sample_map and os.path.exists(args.sample_map):
        with open(args.sample_map) as f:
            loaded = json.load(f)
            loaded.pop("_notes", None)
            sample_map = loaded

    # Load corrections
    corrections = {}
    if args.corrections and os.path.exists(args.corrections):
        with open(args.corrections) as f:
            corrections = json.load(f)

    run_pipeline(
        args.input, args.output,
        sample_map=sample_map,
        corrections=corrections,
        si_threshold=args.si_threshold,
        remove_contaminants=not args.keep_contaminants,
        remove_tmah=not args.keep_tmah,
        renormalize=not args.no_renormalize,
        reference_treatment=args.reference,
    )


# ============================================================
# ENHANCEMENT MODULE: MSI Confidence, Feature IDs, Closed-Sum
# ============================================================

def assign_msi_level(row, library):
    """Assign Metabolomics Standards Initiative (MSI) confidence level.

    Level 1: Standard confirmed (RT+RI+MS with authentic standard)
    Level 2: Putative annotation (SI>=90 + consistent RT across samples)
    Level 3: Compound class only (SI>=80 or library match)
    Level 4: Unknown feature (SI<80 or no reliable classification)
    """
    # Check SI values across all treatments
    si_vals = []
    for key in row:
        if key.startswith("si_") and row[key] and row[key] > 0:
            si_vals.append(row[key])

    avg_si = sum(si_vals) / len(si_vals) if si_vals else 0

    # Check if all treatments name this consistently
    ref_name = None
    names_consistent = True
    for key in row:
        if key.startswith("name_") and row[key]:
            if ref_name is None:
                ref_name = row[key]
            elif row[key] != ref_name:
                names_consistent = False

    # Level 1: would need standard - not available from TXT alone
    # Level 2 proxy: SI>=90 AND consistent naming across all
    if avg_si >= 90 and names_consistent and not row.get("has_conflict"):
        return "Level_2"

    # Level 3: SI>=80 OR classified into a defined category
    ref_cat = None
    for key in row:
        if key.startswith("class_") and row[key] and row[key] != "Unknown":
            ref_cat = row[key]; break
    cat = classify_compound(row.get("name_CK", ref_name or ""), library) if not ref_cat else ref_cat

    if avg_si >= 80 or (cat and cat != "Unknown"):
        return "Level_3"

    # Level 4: everything else
    return "Level_4"


def assign_feature_ids(matrix, prefix="F"):
    """Assign sequential Feature IDs (F001, F002...) based on RT order."""
    sorted_matrix = sorted(matrix, key=lambda r: r.get("rt_ref", 999))
    for i, row in enumerate(sorted_matrix):
        row["feature_id"] = f"{prefix}{i+1:03d}"
    return sorted_matrix


def detect_closed_sum_effect(matrix, treatment_names, threshold=15.0):
    """Detect closed-sum effect: when a single peak >threshold% inflates
    the total area denominator, artificially suppressing other peaks' relative%.

    Returns list of (treatment, feature_id, peak_conc%) warnings.
    """
    warnings = []
    for row in matrix:
        for t in treatment_names:
            conc_key = f"conc_{t}"
            if conc_key in row and row[conc_key] and row[conc_key] > threshold:
                warnings.append({
                    "treatment": t,
                    "feature_id": row.get("feature_id", f"RT{row['rt_ref']:.2f}"),
                    "rt": row["rt_ref"],
                    "conc": row[conc_key],
                    "name": row.get(f"name_{t}", ""),
                })
    return warnings


def generate_filter_log(peaks_by_treatment, output_dir):
    """Generate comprehensive filter log documenting why each peak was kept/removed."""
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "filter_log.csv")
    with open(log_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Treatment", "RT", "Name", "Area", "Conc%", "SI",
            "Status", "Reason", "Action_Recommended"
        ])
        for tname, peaks in peaks_by_treatment:
            for p in peaks:
                status = "FILTERED" if p.get("filtered") else "PASS"
                reason = p.get("filter_reason", "")
                action = ""
                if "SI" in reason:
                    action = "Verify spectrum manually before use"
                elif "contaminant" in reason:
                    action = "Exclude from ecological interpretation"
                elif "artifact" in reason:
                    action = "Reagent artifact (TMAH/derivatization) - exclude from interpretation"
                elif "zero_area" in reason:
                    action = "Check integration"
                if p.get("corrected"):
                    reason += "; name_corrected"
                    action = "Corrected per manual review"
                writer.writerow([
                    tname, p["rt"], p["name"], p["area"], p["conc"], p["si"],
                    status, reason, action
                ])
    print(f"  Filter log: {log_path}")
    return log_path


# ============================================================
# MAIN (updated)
# ============================================================

if __name__ == "__main__":
    main()
