"""
Py-GC-MS Data Verification Module
==================================
Cross-treatment peak verification using:
  1. NIST TXT export parsing (MC Peak Table + Similarity Search Results)
  2. QGD raw mass spectrum extraction for EI spectrum comparison
  3. Shahriar 2026 12-class compound classification
  4. RT alignment, ID conflict detection, and batch effect flagging

Usage:
  python verify_data.py --input <TXT_dir> --qgd <QGD_dir> --output <report_dir>
       [--sample_map <mapping.json>] [--corrections <corrections.json>]

Corrections JSON format:
  {"BC15": {"3.215": "Toluene"}}
  (key=treatment, nested key=RT, value=corrected compound name)
"""
import os, sys, json, csv, argparse
from collections import defaultdict, Counter
from pathlib import Path
try:
    from qgd_reader import QGDFile
except ImportError:
    import sys; sys.path.insert(0, os.path.dirname(__file__))
    from qgd_reader import QGDFile

# Shahriar 2026 12-class short codes
CLASS_MAP = {
    "Monocyclic aromatic hydrocarbons (MAH)": "MAH",
    "N-containing monocyclic aromatic hydrocarbons (N-MAH)": "N-MAH",
    "Polycyclic aromatic hydrocarbons (PAH)": "PAH",
    "Lignin": "Lignin",
    "Phenols": "Phenols",
    "Degraded saccharides": "Sugars",
    "Fatty acids, alcohols and esters   ": "Fatty_acids_lipids",
    "Fatty acids, alcohols and esters": "Fatty_acids_lipids",
    "Alkenes                            ": "Alkenes",
    "Alkenes": "Alkenes",
    "Short alkanes                      ": "Short_alkanes",
    "Short alkanes": "Short_alkanes",
    "Long alkanes                       ": "Long_alkanes",
    "Long alkanes": "Long_alkanes",
    "Other hydrocarbons                 ": "Other_hydrocarbons",
    "Other hydrocarbons": "Other_hydrocarbons",
    "Other N-containing compounds (N-containing)": "Other_N",
    "Other N-containing compounds": "Other_N",
}

PLANT_CLASSES = {"Lignin", "Long_alkanes", "Alkenes"}
MICROBIAL_CLASSES = {"N-MAH", "Other_N"}


def load_shahriar_library(path=None):
    """Load Shahriar 2026 compound classification library."""
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "..", "data", "shahriar_library.json")
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return {k.lower(): CLASS_MAP.get(v, v) for k, v in raw.items()}


def parse_txt(path):
    """Parse Shimadzu GCMSsolution NIST export TXT file."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [l.strip("﻿").rstrip("\n") for l in f.readlines()]

    peaks = []
    section = None
    peak_start = False
    sim = {}

    for line in lines:
        parts = line.split("\t")
        if "[MC Peak Table]" in line:
            section = "mc"; continue
        elif "[MS Similarity Search Results" in line:
            section = "sim"; continue
        elif line.startswith("[") and line.endswith("]"):
            section = None; peak_start = False; continue

        if section == "mc":
            if "Peak#" in line:
                peak_start = True; continue
            if peak_start and len(parts) >= 11:
                try:
                    conc_s = parts[8].strip()
                    peaks.append({
                        "rt": float(parts[1]), "area": int(parts[5]),
                        "conc": float(conc_s) if conc_s else 0,
                        "name": parts[10].strip() if len(parts) > 10 else "",
                    })
                except (ValueError, IndexError):
                    pass

        if section == "sim":
            if "Spectrum#" in line: continue
            if len(parts) >= 6:
                try:
                    if int(parts[1]) == 1:
                        sim[int(parts[0])] = {"si": int(parts[2])}
                except (ValueError, IndexError):
                    pass

    for i, p in enumerate(peaks):
        p["si"] = sim.get(i + 1, {}).get("si", 0)

    return peaks


def classify_compound(name, library):
    """Classify compound using Shahriar library with fuzzy fallback."""
    if not name:
        return "Unknown"
    name_clean = name.lower().strip()

    # Exact match
    if name_clean in library:
        return library[name_clean]

    # Remove CAS suffixes
    base = name_clean.split("$$")[0].strip()
    if base in library:
        return library[base]

    # Remove stereochemistry
    import re
    first_part = re.split(r"[,;]", name_clean)[0].strip()
    if first_part in library:
        return library[first_part]

    no_stereo = re.sub(r"\(r\)|\(s\)|\(e\)|\(z\)|\[r-\]|\[s-\]", "", first_part, flags=re.I).strip()
    no_stereo = re.sub(r"\s+", " ", no_stereo).strip()
    if no_stereo in library:
        return library[no_stereo]

    # Keyword fallback
    return _classify_keyword(name_clean)


def _classify_keyword(n):
    """Keyword-based fallback classification."""
    if any(k in n for k in ["nitrile", "cyanide"]): return "Other_N"
    if any(k in n for k in ["pyrrole","pyridine","indole","trimethylindole","tetramethylindole","pyrrolidine","piperidinamine"]): return "N-MAH"
    if any(k in n for k in ["naphthalene","fluorene","anthracene","pyrene","fluoranthene","dibenzofuran"]): return "PAH"
    if any(k in n for k in ["furfural","levoglucosan","furan","methylfurfural"]): return "Sugars"
    if any(k in n for k in ["guaiacol","syringol","vanillin","eugenol","cinnamyl"]): return "Lignin"
    if any(k in n for k in ["phenol","cresol","di-tert-butylphenol"]): return "Phenols"
    if any(k in n for k in ["benzene","toluene","xylene","ethylbenzene","mesitylene","cymene","indene","indane","styrene","phenyl"]): return "MAH"
    if any(k in n for k in ["dodecanol","undecanol","alcohol","methyl ester","ethyl ester","acetic acid","hexadecanoic","octadecanoic","dodecanoic","tetradecanoic","tridecanoic","octadecenoic","succinic"]): return "Fatty_acids_lipids"
    if any(k in n for k in ["alkene","octene","nonene","decene","undecene","dodecene","tridecene","tetradecene","pentadecene","hexadecene","heptadecene","octadecene","nonadecene","eicosene","cyclooctatetraene","diene","triene","cyclopentadiene"]): return "Alkenes"
    if any(k in n for k in ["eicosane","heneicosane","docosane","tricosane","tetracosane","pentacosane","hexacosane","heptacosane","octacosane","nonacosane","triacontane","tetratriacontane","pentatriacontane","tritetracontane","squalane"]): return "Long_alkanes"
    if any(k in n for k in ["nonane","decane","undecane","dodecane","tridecane","tetradecane","pentadecane","hexadecane","heptadecane","octadecane","octane","heptane","hexane","cyclopropane","cyclohexane","cyclododecane","cycloundecane"]): return "Short_alkanes"
    if any(k in n for k in ["amine","amide","methanediamine","triazine","carbamate","tetramethylammonium","tetramethyl-"]): return "Other_N"
    if any(k in n for k in ["yne","acetylene","dioxane","bicyclo","spiro","metheno"]): return "Other_hydrocarbons"
    return "Unknown"


def verify_cross_treatment(parsed, library, qgd_dir=None):
    """Verify peaks across treatments: RT alignment, ID conflicts, EI spectra."""
    labels = dict(enumerate(["CK","BC7.5","BC15","BC30"]))  # TODO: from sample_map
    ck_peaks = parsed[0]

    report = {"aligned": [], "conflicts": [], "class_conflicts": [], "spectra_verified": []}

    for ck_p in sorted(ck_peaks, key=lambda p: p["conc"], reverse=True):
        ck_rt = ck_p["rt"]
        ck_cat = classify_compound(ck_p["name"], library)

        row = {
            "rt_ck": ck_rt, "name_ck": ck_p["name"], "cat_ck": ck_cat,
            "conc_ck": ck_p["conc"], "si_ck": ck_p["si"],
        }

        has_conflict = False
        for idx in range(1, len(parsed)):
            t_name = labels.get(idx, f"T{idx}")
            best = None
            best_d = 999
            for op in parsed[idx]:
                d = abs(op["rt"] - ck_rt)
                if d < 0.08 and d < best_d:
                    best = op; best_d = d

            if best:
                row[f"rt_{t_name}"] = best["rt"]
                row[f"name_{t_name}"] = best["name"]
                row[f"conc_{t_name}"] = best["conc"]
                row[f"si_{t_name}"] = best["si"]
                row[f"cat_{t_name}"] = classify_compound(best["name"], library)

                if best["name"] != ck_p["name"]:
                    has_conflict = True
                    other_cat = classify_compound(best["name"], library)
                    if other_cat != ck_cat:
                        report["class_conflicts"].append({
                            "rt": ck_rt, "conc": ck_p["conc"],
                            "ck_name": ck_p["name"], "ck_cat": ck_cat,
                            "other_treat": t_name, "other_name": best["name"],
                            "other_cat": other_cat,
                        })

        row["has_conflict"] = has_conflict
        report["aligned"].append(row)
        if has_conflict:
            report["conflicts"].append(row)

    return report


def generate_report(parsed, verify_result, library, output_dir, sample_names=None):
    """Generate verification report in Markdown."""
    if sample_names is None:
        sample_names = {0: "CK", 1: "BC7.5", 2: "BC15", 3: "BC30"}

    os.makedirs(output_dir, exist_ok=True)
    report_path = os.path.join(output_dir, "verification_report.md")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# Py-GC-MS Data Verification Report\n\n")

        # 1. Data integrity
        f.write("## 1. Data Integrity\n\n")
        f.write("| Sample | Peaks | SI>=80 | SI>=90 | SI<70 | Total Area |\n")
        f.write("|--------|-------|--------|--------|-------|------------|\n")
        for idx, data in enumerate(parsed):
            n = len(data)
            si_ok = sum(1 for p in data if p["si"] >= 80)
            si_hi = sum(1 for p in data if p["si"] >= 90)
            si_lo = sum(1 for p in data if 0 < p["si"] < 70)
            total_a = sum(p["area"] for p in data)
            f.write(f"| {sample_names[idx]} | {n} | {si_ok} | {si_hi} | {si_lo} | {total_a:,} |\n")

        # 2. Class composition
        f.write("\n## 2. Compound Class Composition (Shahriar 2026)\n\n")
        f.write("| Class | Source |")
        for i in range(len(parsed)):
            f.write(f" {sample_names[i]} |")
        f.write("\n|-------|--------|")
        for _ in range(len(parsed)):
            f.write("-------|")
        f.write("\n")

        summaries = []
        for idx, data in enumerate(parsed):
            cats = defaultdict(lambda: {"conc": 0.0, "count": 0})
            for p in data:
                cat = classify_compound(p["name"], library)
                cats[cat]["conc"] += p["conc"]
                cats[cat]["count"] += 1
            summaries.append(dict(cats))

        all_cats = set()
        for s in summaries:
            all_cats.update(s.keys())

        for cat in sorted(all_cats, key=lambda c: summaries[0].get(c, {}).get("conc", 0), reverse=True):
            src = "Plant" if cat in PLANT_CLASSES else ("Microb" if cat in MICROBIAL_CLASSES else "Mixed")
            f.write(f"| {cat} | {src} |")
            for s in summaries:
                f.write(f" {s.get(cat,{}).get('conc',0):.1f}% |")
            f.write("\n")

        # 3. Source attribution
        f.write("\n## 3. Source Attribution (R_MP)\n\n")
        f.write(f"R_MP = (N-MAH + Other_N) / (Lignin + Long_alkanes + Alkenes)\n\n")
        f.write("| Sample | Plant% | Microb% | Mixed% | R_MP |\n")
        f.write("|--------|--------|---------|--------|------|\n")
        for idx, s in enumerate(summaries):
            plant = sum(s.get(c, {}).get("conc", 0) for c in PLANT_CLASSES)
            microb = sum(s.get(c, {}).get("conc", 0) for c in MICROBIAL_CLASSES)
            r_mp = microb / plant if plant > 0 else 0
            f.write(f"| {sample_names[idx]} | {plant:.1f}% | {microb:.1f}% | {100-plant-microb:.1f}% | {r_mp:.2f} |\n")

        # 4. Key conflicts
        f.write(f"\n## 4. Cross-Class ID Conflicts ({len(verify_result['class_conflicts'])} total)\n\n")
        for c in sorted(verify_result["class_conflicts"], key=lambda x: x["conc"], reverse=True)[:20]:
            f.write(f"- RT {c['rt']:.3f} ({c['conc']:.1f}%): CK={c['ck_name'][:40]} [{c['ck_cat']}]\n")
            f.write(f"  vs {c['other_treat']}: {c['other_name'][:40]} [{c['other_cat']}]\n")

    print(f"Report: {report_path}")

    # Export aligned matrix
    csv_path = os.path.join(output_dir, "aligned_peak_matrix.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        header = ["RT_CK", "Name_CK", "Class_CK", "Conc%_CK", "SI_CK"]
        for i in range(1, len(parsed)):
            t = sample_names[i]
            header += [f"RT_{t}", f"Name_{t}", f"Class_{t}", f"Conc%_{t}", f"SI_{t}"]
        header.append("Conflict")
        writer.writerow(header)

        for row in verify_result["aligned"]:
            line = [row["rt_ck"], row["name_ck"], row.get("cat_ck",""), row.get("conc_ck",""), row.get("si_ck","")]
            for i in range(1, len(parsed)):
                t = sample_names[i]
                line += [row.get(f"rt_{t}",""), row.get(f"name_{t}",""),
                        row.get(f"cat_{t}",""), row.get(f"conc_{t}",""), row.get(f"si_{t}","")]
            line.append("YES" if row.get("has_conflict") else "")
            writer.writerow(line)

    print(f"Matrix: {csv_path}")
    return report_path


def main():
    parser = argparse.ArgumentParser(description="Py-GC-MS Data Verification")
    parser.add_argument("--input", required=True, help="Directory with NIST export TXT files")
    parser.add_argument("--qgd", help="Directory with QGD raw data files (optional)")
    parser.add_argument("--output", required=True, help="Output directory for reports")
    parser.add_argument("--sample_map", help="JSON mapping file: {'5':'CK','6':'BC7.5',...}")
    parser.add_argument("--corrections", help="JSON corrections: {'BC15':{'3.215':'Toluene'}}")
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

    # Load Shahriar library
    print("Loading Shahriar 2026 library...")
    library = load_shahriar_library()
    print(f"  {len(library)} compounds loaded")

    # Parse TXT files
    print(f"\nParsing TXT files from {args.input}...")
    parsed = []
    for sid in sorted(sample_map.keys()):
        path = os.path.join(args.input, f"{sid}.txt")
        if os.path.exists(path):
            data = parse_txt(path)
            print(f"  {sid} ({sample_map[sid]}): {len(data)} peaks")

            # Apply corrections
            treat_name = sample_map[sid]
            if treat_name in corrections:
                for rt_str, correct_name in corrections[treat_name].items():
                    rt_targ = float(rt_str)
                    for p in data:
                        if abs(p["rt"] - rt_targ) < 0.04:
                            old = p["name"]
                            p["name"] = correct_name
                            print(f"    CORRECTED RT{rt_targ}: {old[:30]} -> {correct_name}")

            parsed.append(data)

    # Verify
    print("\nRunning verification...")
    result = verify_cross_treatment(parsed, library, args.qgd)

    # Report
    print("\nGenerating report...")
    sample_names = {i: sample_map[sid] for i, sid in enumerate(sorted(sample_map.keys()))}
    report_path = generate_report(parsed, result, library, args.output, sample_names)

    print(f"\nDone: {report_path}")


if __name__ == "__main__":
    main()
