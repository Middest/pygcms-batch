"""
Py-GC-MS EI-Spectrum Conflict Resolution
=========================================
Resolve cross-treatment HIGH ID conflicts using raw QGD mass spectra.

For each HIGH-conflict aligned feature, extract the raw EI spectrum at the
matching RT from every treatment's QGD file and compute spectral similarity
(cosine over top-N ions). This provides an independent ground-truth:
  - If spectra across treatments are highly similar (cosine >= threshold),
    the same compound is present and the differing NIST names are
    identification noise → assign a UNIFIED class (majority vote among
    non-artifact names, tie-broken by highest SI).
  - If spectra genuinely differ, the treatments carry different compounds
    (co-elution / genuine difference) → keep per-treatment classes.

Usage:
  python resolve_conflicts_ei.py --matrix <aligned_matrix.csv> --qgd <QGD_dir>
       [--sample_map <mapping.json>] [--cosine 0.85]
       [--top_ions 12] [--output <out_dir>]

Outputs:
  conflict_resolution.csv  - per-feature verdict + unified class
  unified_class_composition.csv - class composition using unified classes
"""
import os, sys, json, csv, argparse
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
try:
    from qgd_reader import QGDFile
except ImportError:
    from scripts.qgd_reader import QGDFile


# Common background ions in Py-GC-MS (CO2=44, air/N2=28, and their fragments).
# These dominate weak-peak spectra and drag cosine down; excluded from matching.
BACKGROUND_MZ = {28.0, 28.1, 43.9, 44.0, 44.1, 45.0, 45.1}


def norm_intensity(peaks, exclude_bg=True):
    """Normalize peak intensities to unit vector (for cosine).

    Background ions (CO2 m/z 44 etc.) are excluded by default so weak peaks
    are compared on their actual fragment ions, not the common background.
    """
    m = max((i for _, i in peaks), default=0)
    if m <= 0:
        return {}
    out = {}
    for mz, i in peaks:
        if i <= 0:
            continue
        if exclude_bg and round(mz, 1) in BACKGROUND_MZ:
            continue
        out[mz] = i / m
    return out


def cosine_similarity(spec_a, spec_b, top_n=12):
    """Cosine similarity over top-N ions of the query spectrum.

    Uses top-N by intensity from spec_a (background ions excluded), matched
    against spec_b with +/-0.5 Da window (nominal mass tolerance).
    """
    if not spec_a or not spec_b:
        return 0.0
    a = [p for p in spec_a if round(p[0], 1) not in BACKGROUND_MZ]
    if len(a) < 3:
        a = spec_a  # nothing left after bg removal — use full spec
    a = sorted(a, key=lambda x: x[1], reverse=True)[:top_n]
    b = norm_intensity(spec_b)
    num = denom_a = denom_b = 0.0
    for mz_a, i_a in a:
        hit = None
        best_d = 0.7
        for mz_b, i_b in b.items():
            d = abs(mz_a - mz_b)
            if d < best_d:
                hit = i_b; best_d = d
        num += i_a * hit if hit else 0.0
        denom_a += i_a * i_a
        denom_b += (hit * hit if hit else 0.0)
    if denom_a <= 0 or denom_b <= 0:
        return 0.0
    return num / (denom_a ** 0.5 * denom_b ** 0.5)


def load_artifact_kws():
    """Return list of reagent-artifact name substrings."""
    return [
        # TMAH thermochemolysis
        "methylamine, n,n-dimethyl", "methanediamine", "triazine, hexahydro",
        "methenamine", "acetonitrile, (dimethylamino)",
        # derivatization reagents
        "pentafluoroprop", "trifluoroacetate", "trifluoroacet",
        "heptafluorobuty", "pentafluorooctano", "perfluorooctanoic",
        "trimethylsilyl", "tert-butyldimethylsilyl", "bis(trimethylsilyl)", "silyl",
    ]


def is_artifact(name):
    if not name:
        return False
    n = name.lower()
    return any(k in n for k in load_artifact_kws())


CLASS_MAP = {
    "Monocyclic aromatic hydrocarbons (MAH)": "MAH",
    "N-containing monocyclic aromatic hydrocarbons (N-MAH)": "N-MAH",
    "Polycyclic aromatic hydrocarbons (PAH)": "PAH",
    "Lignin": "Lignin", "Phenols": "Phenols",
    "Degraded saccharides": "Sugars",
    "Fatty acids, alcohols and esters": "Fatty_acids_lipids",
    "Alkenes": "Alkenes",
    "Short alkanes": "Short_alkanes",
    "Long alkanes": "Long_alkanes",
    "Other hydrocarbons": "Other_hydrocarbons",
    "Other N-containing compounds (N-containing)": "Other_N",
    "Other N-containing compounds": "Other_N",
}


def load_library(path=None):
    """Load Shahriar library and return lowercased name -> short class."""
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "..", "data", "shahriar_library.json")
    with open(path, encoding="utf-8") as f:
        raw = json.load(f)
    return {k.lower().strip(): CLASS_MAP.get(v, v) for k, v in raw.items()}


def classify_compound(name, library, prev_class=None):
    """Classifier matching pipeline.py behavior (library match + keyword fallback).

    Imported from pipeline.py when available so results are consistent with the
    main pipeline output; falls back to a local keyword classifier otherwise.
    """
    if not name:
        return "Unknown"
    try:
        from pipeline import classify_compound as _pipeline_cls
        return _pipeline_cls(name, library)
    except (ImportError, Exception):
        pass
    n = name.lower().strip()
    if n in library:
        return library[n]
    base = n.split("$$")[0].strip()
    if base in library:
        return library[base]
    import re
    first = re.split(r"[,;]", n)[0].strip()
    if first in library:
        return library[first]
    # keyword fallback (mirrors pipeline._keyword_classify)
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
    if any(k in n for k in ["naphthalene", "fluorene", "anthracene", "pyrene",
            "fluoranthene", "dibenzofuran", "azulene"]):
        return "PAH"
    if any(k in n for k in ["furfural", "levoglucosan", "furan", "methylfurfural"]):
        return "Sugars"
    if any(k in n for k in ["guaiacol", "syringol", "vanillin", "eugenol", "cinnamyl"]):
        return "Lignin"
    if any(k in n for k in ["phenol", "cresol", "di-tert-butylphenol"]):
        return "Phenols"
    if any(k in n for k in ["benzene", "toluene", "xylene", "ethylbenzene",
            "mesitylene", "cymene", "indene", "indane", "styrene", "phenyl",
            "anethole", "biphenyl"]):
        return "MAH"
    if any(k in n for k in ["dodecanol", "undecanol", "alcohol", "methyl ester",
            "ethyl ester", "acetic acid", "hexadecanoic", "octadecanoic",
            "dodecanoic", "tetradecanoic", "tridecanoic", "octadecenoic", "succinic"]):
        return "Fatty_acids_lipids"
    if any(k in n for k in ["alkene", "octene", "nonene", "decene", "undecene",
            "dodecene", "tridecene", "tetradecene", "pentadecene", "hexadecene",
            "heptadecene", "octadecene", "nonadecene", "eicosene", "cyclooctatetraene",
            "diene", "triene", "cyclopentadiene"]):
        return "Alkenes"
    if any(k in n for k in ["eicosane", "heneicosane", "docosane", "tricosane",
            "tetracosane", "pentacosane", "hexacosane", "heptacosane", "octacosane",
            "nonacosane", "triacontane", "tetratriacontane", "pentatriacontane",
            "tritetracontane", "squalane"]):
        return "Long_alkanes"
    if any(k in n for k in ["nonane", "decane", "undecane", "dodecane", "tridecane",
            "tetradecane", "pentadecane", "hexadecane", "heptadecane", "octadecane",
            "octane", "heptane", "hexane", "cyclopropane", "cyclohexane",
            "cyclododecane", "cycloundecane"]):
        return "Short_alkanes"
    if any(k in n for k in ["yne", "acetylene", "dioxane", "bicyclo", "spiro", "metheno"]):
        return "Other_hydrocarbons"
    return "Unknown"


def main():
    ap = argparse.ArgumentParser(description="EI-spectrum conflict resolution")
    ap.add_argument("--matrix", required=True, help="analysis_ready_matrix.csv")
    ap.add_argument("--qgd", required=True, help="Directory with QGD raw files")
    ap.add_argument("--sample_map", help="JSON: {'5':'CK','6':'BC7.5',...}")
    ap.add_argument("--cosine", type=float, default=0.85, help="Similarity threshold (default 0.85)")
    ap.add_argument("--top_ions", type=int, default=12, help="Top-N ions for cosine (default 12)")
    ap.add_argument("--output", default=".", help="Output directory")
    args = ap.parse_args()

    smap = {"5": "CK", "6": "BC7.5", "7": "BC15", "8": "BC30"}
    if args.sample_map and os.path.exists(args.sample_map):
        with open(args.sample_map) as f:
            smap = {k: v for k, v in json.load(f).items() if k != "_notes"}
    file_ids = {t: sid for sid, t in smap.items()}

    # read matrix
    with open(args.matrix, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    trts = list(smap.values())
    cmap = {t: (f"Conc%_ref" if t == "CK" else f"Conc%_{t}") for t in trts}
    nmap = {t: (f"Name_ref" if t == "CK" else f"Name_{t}") for t in trts}
    rmap = {t: (f"RT_ref" if t == "CK" else f"RT_{t}") for t in trts}

    lib = load_library()
    # pre-classify each treatment name
    prev_class = {}
    for r in rows:
        for t in trts:
            nm = r.get(nmap[t])
            if nm and t not in prev_class:
                prev_class[t] = classify_compound(nm, lib)

    # open QGD files lazily
    qgd_cache = {}
    def get_qgd(t):
        if t not in qgd_cache:
            sid = file_ids.get(t)
            if not sid:
                return None
            # try .qgd or .QGD
            for ext in (".qgd", ".QGD"):
                p = os.path.join(args.qgd, f"{sid}{ext}")
                if os.path.exists(p):
                    qgd_cache[t] = QGDFile(p)
                    return qgd_cache[t]
            qgd_cache[t] = None
        return qgd_cache[t]

    # collect HIGH conflicts
    conflicts = [r for r in rows if r.get("Conflict_Severity") == "HIGH"]
    out_rows = []
    unified = defaultdict(lambda: {t: 0.0 for t in trts})  # class -> trt -> conc

    print(f"Evaluating {len(conflicts)} HIGH-conflict features (cosine>={args.cosine})...")
    n_unified = n_kept_sep = n_no_spec = 0

    for r in conflicts:
        ref_name = r.get("Name_ref", "")
        ref_rt = float(r["RT_ref"])
        # collect per-treatment names + concs for this feature
        names = {}
        concs = {}
        rts = {}
        for t in trts:
            nm = r.get(nmap[t])
            if nm:
                names[t] = nm
                try:
                    concs[t] = float(r.get(cmap[t]) or 0)
                except:
                    concs[t] = 0.0
                try:
                    rts[t] = float(r.get(rmap[t]))
                except:
                    rts[t] = None

        # extract spectra from QGD at each treatment's RT
        specs = {}
        for t in names:
            qgd = get_qgd(t)
            rt = rts[t] or ref_rt
            if qgd is not None:
                s = qgd.get_spectrum_at_rt(rt, tolerance=0.08)
                if s:
                    specs[t] = s["peaks"]

        # reference spectrum: CK if available, else any
        ref_spec = specs.get("CK")
        if ref_spec is None and specs:
            ref_t = next(iter(specs))
            ref_spec = specs[ref_t]

        # verdict
        if ref_spec is None:
            verdict = "NO_SPECTRA"
            n_no_spec += 1
        else:
            sims = {}
            for t in names:
                if specs.get(t) and t != (next(iter(specs)) if ref_spec is not None else "CK"):
                    sims[t] = cosine_similarity(ref_spec, specs[t], args.top_ions)
            # decide: all available spectra similar?
            sim_vals = list(sims.values()) if sims else [1.0]
            all_similar = bool(sims) and min(sim_vals) >= args.cosine
            # (if only one spectrum, can't decide — treat as similar)
            if not sims and len(names) == 1:
                all_similar = True

            if all_similar:
                # UNIFIED class: majority vote among non-artifact names, tie -> CK
                verdict = "UNIFIED"
                n_unified += 1
                valid = {t: nm for t, nm in names.items() if not is_artifact(nm)}
                if not valid:
                    valid = names  # all artifact names (rare) — keep but flag
                classes = defaultdict(int)
                for t, nm in valid.items():
                    classes[classify_compound(nm, lib, prev_class.get(t))] += 1
                max_c = max(classes.values())
                unified_class = sorted(c for c, n in classes.items() if n == max_c)[0]
                # apply unified class to all treatments with conc
                for t in concs:
                    unified[unified_class][t] += concs[t]
            else:
                verdict = "GENUINE_DIFF"
                n_kept_sep += 1
                # keep per-treatment classes
                for t in concs:
                    cls = classify_compound(names.get(t, ""), lib, prev_class.get(t))
                    unified[cls][t] += concs[t]

        out_rows.append({
            "Feature": r.get("Feature_ID", ""),
            "RT_ref": ref_rt,
            "Name_ref": ref_name,
            "CK_class": r.get("Class_ref", ""),
            "Verdict": verdict,
            "N_treatments": len(names),
            "Cosine_min": round(min(sims.values()), 3) if sims else "",
            "Unified_class": unified_class if verdict == "UNIFIED" else r.get("Class_ref", ""),
        })
    print(f"  UNIFIED (same compound, ID noise): {n_unified}")
    print(f"  GENUINE_DIFF (different spectra): {n_kept_sep}")
    print(f"  NO_SPECTRA (unresolved): {n_no_spec}")

    os.makedirs(args.output, exist_ok=True)
    # export resolution table
    with open(os.path.join(args.output, "conflict_resolution.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()) if out_rows else [])
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    # export unified class composition for the CONFLICT contributions
    with open(os.path.join(args.output, "unified_conflict_class_composition.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["Class"] + trts)
        for cls in sorted(unified, key=lambda c: -max(unified[c].values())):
            w.writerow([cls] + [round(unified[cls][t], 2) for t in trts])

    # NEW: canonical ei_decisions.csv (Stage 3 output for apply_final.py)
    # decision: UNIFIED / GENUINE_DIFF / NO_SPECTRA (kept enum for compatibility)
    ei_path = os.path.join(args.output, "ei_decisions.csv")
    with open(ei_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=[
            "feature_id", "rt_ref", "name_ref", "verdict", "n_treatments",
            "cosine_min", "unified_class", "reason",
        ])
        w.writeheader()
        for r in out_rows:
            reason = ("EI_COSINE_GE_THRESHOLD" if r["Verdict"] == "UNIFIED"
                      else "EI_COSINE_LT_THRESHOLD" if r["Verdict"] == "GENUINE_DIFF"
                      else "NO_QGD_SPECTRUM")
            w.writerow({
                "feature_id": r["Feature"],
                "rt_ref": r["RT_ref"],
                "name_ref": r["Name_ref"],
                "verdict": r["Verdict"],
                "n_treatments": r["N_treatments"],
                "cosine_min": r["Cosine_min"],
                "unified_class": r["Unified_class"],
                "reason": reason,
            })

    print(f"\nWrote: {args.output}/conflict_resolution.csv")
    print(f"      {args.output}/unified_conflict_class_composition.csv")
    print(f"      {args.output}/ei_decisions.csv")
    # close QGD files
    for q in qgd_cache.values():
        if q:
            q.close()


if __name__ == "__main__":
    main()
