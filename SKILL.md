---
name: pygcms-batch
description: Batch process Py-GC-MS NIST export TXT files. Filter by SI>=80, remove TMAH thermochemolysis artifacts (reagent-derived N compounds), renormalize class composition to 100%, classify compounds by Chen 2023 (Carbon Research) and Kallenbach 2016 (Nature Communications) schemes, compute elemental atoms (C/H/O/N/P/S), NOSC, DG_COX, and source attribution (plant/microbial/mixed). Output multi-sheet Excel. Data verification with Shahriar 2026 library, cross-treatment RT alignment, and ID conflict detection.
---

# Py-GC-MS Batch Analysis

Batch process Py-GC-MS NIST search export TXT files for soil organic matter molecular characterization.

## Quick Start: Raw-to-Analysis Pipeline

One command to go from raw TXT files to analysis-ready data:

```
python scripts/pipeline.py --input <TXT_directory> --output <output_directory>
```

This runs: PARSE -> FILTER -> ALIGN -> RESOLVE -> CLASSIFY -> VALIDATE -> EXPORT.

1. Auto-detect all `.TXT` files in input directory (recursively)
2. Parse `[MC Peak Table]` for TIC peak areas and retention times
3. Parse `[MS Similarity Search Results for Spectrum Process Table]` for NIST library matches
4. Filter compounds: SI >= 80 (NIST Match Factor)
5. **Remove TMAH thermochemolysis artifacts** (reagent-derived N compounds that
   otherwise inflate the microbial signal — see "TMAH Thermochemolysis Artifact Removal" below)
6. Merge MC area data with NIST match data (SI, CAS, molecular formula, molecular weight)
7. Parse molecular formula → atom counts (C, H, O, N, P, S), detect heteroatoms (Cl, Br, F, I, Si, B)
8. Compute NOSC and ΔG_COX for each valid formula
9. **(Optional)** Query NIST Chemistry WebBook by CAS/name to fill missing molecular formulas, InChI, InChIKey
10. Classify compounds by TWO schemes:
   - **Chen 2023** (Carbon Research): 10 chemical groups + 3 source categories
   - **Kallenbach 2016** (Nature Communications): 8 origin-based classes
11. Classification uses 286-entry exact lookup table (built from manually curated data) with keyword fallback
12. Generate source attribution statistics (plant/microbial/mixed %) by relative peak area
13. Output multi-sheet Excel workbook

## Usage

### Basic
```
python scripts/pygcms_batch.py --input <TXT_directory> --output <output.xlsx> [--sample_map <mapping.json>]
```

### With NIST WebBook enrichment
```
python scripts/pygcms_batch.py --input <TXT_directory> --output <output.xlsx> --enrich [--enrich_delay 0.5]
```

### Sample mapping file (optional JSON)
```json
{
  "5": "CK",
  "6": "BC7.5",
  "7": "BC15",
  "8": "BC30"
}
```

## Output Excel sheets

| Sheet | Content |
|-------|---------|
| `筛选统计` | Peak filtering summary per sample: raw peaks, SI≥80 peaks, total TIC area, unique compounds, exotic element counts |
| `Chen_分类详细` | Compound-level Chen classification with **SI, CAS, formula, MW, InChIKey, NOSC, ΔG_COX** |
| `Chen_来源统计` | Source attribution statistics: plant/microbial/mixed % by relative area |
| `Chen_来源饼图` | Pie chart data: plant/microbial/mixed % per treatment |
| `Kallenbach_分类详细` | Compound-level Kallenbach classification with **SI, CAS, formula, MW, InChIKey** |
| `Kallenbach_来源统计` | Kallenbach source attribution statistics |
| `Kallenbach_来源饼图` | Kallenbach pie chart data |
| `说明_QC` | Metadata, classification references, QC notes |

## Per-Compound Information

Each compound in the output includes:

| Field | Source | Description |
|-------|--------|-------------|
| Peak# | MC Peak Table | Peak number in chromatogram |
| RT (min) | MC Peak Table | Retention time |
| Area | MC Peak Table | TIC peak area |
| Rel.Area% | Calculated | Relative abundance (% of total sample area) |
| **SI** | NIST MS Search | NIST Match Factor (0-100), only SI≥80 retained |
| CAS | NIST MS Search | CAS Registry Number |
| Name | NIST MS Search | Top Hit#1 compound name |
| **Formula** | NIST MS Search / WebBook | Molecular formula |
| **Mol.Weight** | NIST MS Search / WebBook | Molecular weight (g/mol) |
| **C, H, O, N, P, S** | Parsed from formula | Atom counts per molecule |
| **Other elements** | Parsed from formula | Cl, Br, F, I, Si, B etc. |
| InChIKey | NIST WebBook (--enrich) | Standard InChIKey for cross-referencing |
| NOSC | Calculated | Nominal Oxidation State of Carbon |
| ΔG_COX (kJ/mol C) | Calculated | Gibbs free energy of oxidation |
| Chen Class | Lookup | Chen 2023 chemical group |
| Chen Source | Lookup | plant / microbial / mixed |
| Kallenbach Class | Lookup | Kallenbach 2016 chemical group |
| Kallenbach Source | Lookup | plant / microbial / mixed |

## Data Verification

When a user needs to verify Py-GC-MS data across multiple treatments (e.g., biochar dose-response), use the verification module:

### Quick verification
```
python scripts/verify_data.py --input <TXT_dir> --output <report_dir> [--sample_map <mapping.json>] [--corrections <corrections.json>] [--qgd <QGD_dir>]
```

### What it checks

| Step | What | Output |
|------|------|--------|
| 1. Data integrity | Peak counts, SI distributions (stratified: SI≥90, SI≥80, 70≤SI<80, SI<70), A/H ratio, total area per sample | Report table |
| 2. Anomalous peaks | Wide peaks, over-integrated peaks, single-peak dominance >30%, RT-RI mismatch | Flagged peak list with NA assignments |
| 3. Class composition | Shahriar 2026 14-class scheme (with corrections for haloalkanes, Cymene, triazines, benzonitrile) | Report + CSV |
| 4. Source attribution | R_MP ratio (microbial/plant signal balance) with per-class breakdown | Report |
| 5. RT alignment | Cross-treatment peak matching (ΔRT ≤ 0.10 min) with EI spectral constraint | Aligned matrix CSV |
| 6. ID conflicts | Flagged cross-class naming discrepancies across treatments | Report + conflict list |
| 7. Classification corrections | Auto-detect known misclassification patterns | Corrections log |
| 8. Sensitivity analysis | Class composition under S0–S4 cleaning scenarios | Multi-scenario table |
| 9. Reconciliation | Peak-count and network-metric closure checks | Numeric verification |
| 10. EI spectrum verification | Raw QGD or TXT spectrum comparison at key RTs | Console output |

### Shahriar 2026 classification

Uses the official 1,678-compound library from Shahriar et al. (2026) ES&T
Supporting Information Table S2 (`data/shahriar_library.json`).

See `references/shahriar_2026_classification.md` for the full scheme.

### QGD raw data extraction

When QGD files are available, use `scripts/qgd_reader.py` to extract
raw mass spectra for EI spectrum comparison across treatments:

```python
from qgd_reader import QGDFile
with QGDFile("sample.qgd") as qgd:
    spec = qgd.get_spectrum_at_rt(18.73)  # Get spectrum at RT
    print(spec["peaks"][:10])             # Top m/z peaks
```

### Corrections file format

```json
{
  "BC15": {"3.215": "Toluene"},
  "BC30": {"18.729": "Pentadecanenitrile"}
}
```
Key=treatment name, nested key=RT, value=corrected compound name.

### When to verify

- **Always**: Before using Py-GC-MS data in a manuscript
- **Especially**: Multi-treatment comparisons where NIST auto-matching may differ
- **Required**: When same RT peak gets different names across treatments
- **Optional**: When all QGD raw files are available for spectrum-level validation

### EI-Spectrum Conflict Resolution (`resolve_conflicts_ei.py`)

When the aligned matrix has HIGH cross-class ID conflicts (same RT → different
NIST names in different treatments), resolve them against the QGD raw spectra
instead of trusting the names:

```
python scripts/resolve_conflicts_ei.py \
    --matrix <dir>/analysis_ready_matrix.csv \
    --qgd <QGD_dir> \
    --sample_map <mapping.json> \
    --cosine 0.85 --top_ions 12 \
    --output <dir>/conflict_resolution
```

For each HIGH-conflict feature, the script extracts the EI spectrum at the
matching RT from every treatment's QGD file and computes cosine similarity
(top-N ions, background m/z 44/28 excluded):
- **cosine ≥ threshold** → same compound present, differing NIST names are
  ID noise → **UNIFIED** class (majority vote among non-artifact names, tie →
  reference treatment).
- **cosine < threshold** → treatments genuinely carry different compounds
  (co-elution / real difference) → **GENUINE_DIFF**, keep per-treatment classes.

Outputs: `conflict_resolution.csv` (per-feature verdict + unified class) and
`unified_conflict_class_composition.csv`. Paddy-soil results (2026): Bulk 62
conflicts → 57 GENUINE_DIFF (batch-level spectral differences — see batch
caveat below); POC 34 UNIFIED / 28 DIFF; MAOC 42 UNIFIED / 22 DIFF.

**Batch caveat (Bulk fraction)**: if treatments were measured in different
analytical batches, spectra will differ across batches at the *same* RT and
the conflict verdicts will be dominated by GENUINE_DIFF. That does NOT mean
the compounds genuinely differ per treatment — it means the batches are not
directly comparable. Check file timestamps / total-ion-area before interpreting
Bulk cross-treatment composition (paddy-soil 2026: Bulk CK/BC30 re-measured in
a clean batch; BC7.5/BC15 from an older weak batch → Short_alkanes 44/36% is a
batch artifact, not a biochar effect).

### Reproducible FINAL composition

The final class-composition table (EI-resolved + artifact-removed + renormalized)
must be produced by a script, not interactively. Reference implementation:
`correct_final_tma.py` (paddy-soil session) — applies targeted removal of
spectrum-confirmed TMAH reagent peaks on top of the resolved composition.
See the "TMAH Thermochemolysis Artifact Removal" section for the spectral check.

## Advanced Data Cleaning & Verification

Beyond the basic SI≥80 pipeline, Py-GC-MS data for multi-treatment comparisons requires
systematic cleaning and reconciliation. See `references/data-cleaning-molecular-networking.md`
for the complete workflow.

### Anomalous Peak Detection

Run after initial parsing to flag peaks that distort quantitative results:

```
python scripts/verify_data.py --input <TXT_dir> --output <report_dir> --detect_anomalies
```

| Check | Criterion | Action |
|-------|-----------|--------|
| Wide peak (co-elution) | A/H ratio >20 or >5× sample median | Flag, set area to NA |
| Over-integration | Integration window >1.0 min or contains multiple apexes | Flag, set area to NA |
| Single-peak dominance | One peak >30% of total sample area | Flag, verify blank/derivatization |
| RT-RI mismatch | NIST RI deviates from estimated RI by >200 | Flag as possible misidentification |

**Critical rule for NA vs 0**: Anomalous peak areas must be set to **NA** (not 0).
NA = "peak exists but area is unreliable"; 0 = "compound not detected".

### Classification Correction

Common systematic misclassifications are detected automatically:

| Pattern | Incorrect | Correct | Detection |
|---------|-----------|---------|-----------|
| p-/o-Cymene | Alkenes | MAH | CAS-based pattern match |
| Haloalkanes (Br/I-substituted) | Alkanes | Other_hydrocarbons | Formula parsing (halogen detection) |
| Hexahydro-triazines | N-MAH | Other_N | Saturated N-heterocycle check |
| Benzonitrile and related | Other_N | N-MAH | Cyano+aromatic pattern |
| Tetratetracontane (C44) | Other_hydrocarbons | Long_alkanes | Long n-alkane pattern |
| Invalid CAS (0-00-0) | Specific name retained | Downgrade to spectral type | CAS validation |

After applying corrections, **all derived tables must be recomputed**:
class composition, source attribution, R_MP ratio, dose-response trends.

### TMAH Thermochemolysis Artifact Removal (Required for N-source Interpretation)

When samples are derivatized with TMAH (tetramethylammonium hydroxide) before
500°C pyrolysis, the reagent itself decomposes and its fragments condense with
formaldehyde (from carbohydrates) to form **N-methylated reagent artifacts**.
These are NOT soil-derived N compounds, but `classify_compound()`'s keyword
fallback routes them into `Other_N` → **Microbial** source, inflating the
microbial signal and R_MP ratio.

**Documented TMAH artifacts removed by the pipeline** (`TMAH_ARTIFACTS` in
`pipeline.py`):

| Compound | Formation mechanism |
|----------|---------------------|
| Methylamine, N,N-dimethyl- (trimethylamine) | TMAH Hofmann elimination |
| Methanediamine, N,N,N',N'-tetramethyl- | TMA + formaldehyde condensation |
| 1,3,5-Triazine, hexahydro-1,3,5-trimethyl- | TMA + formaldehyde cyclization |
| 1,2,4-triazine, hexahydro-3,5-dione | TMA + formaldehyde cyclization |
| Methenamine (hexamethylenetetramine) | TMA + formaldehyde condensation |
| Acetonitrile, (dimethylamino)- | Formaldehyde + cyanide |

**Default behavior**: TMAH artifacts are removed by default. Toggle with CLI flags:
```
python scripts/pipeline.py --input <dir> --output <dir> --si_threshold 80          # default: remove TMAH + renormalize
python scripts/pipeline.py --input <dir> --output <dir> --keep_tmah                # keep artifacts (control)
python scripts/pipeline.py --input <dir> --output <dir> --no_renormalize           # raw Conc% sums, no renorm
```

**Renormalization**: After removing TMAH artifacts, class composition is
renormalized to 100% of kept-peak total per treatment (relative abundance).
This makes treatments comparable even when raw Conc% denominators differ.
The verification report notes whether values are renormalized.

**Impact pattern observed (paddy soil biochar, 2026)**: MAOC was nearly
unaffected by the *name-based* removal (ΔOther_N < 9 pp) — but see the
spectral caveat below. Bulk and POC had large inflation in BC treatments
(Other_N dropped 35-56 pp after name-based removal), so R_MP patterns in
Bulk/POC were partly reagent artifacts. **Always report both control
(keep-TMAH) and cleaned (remove-TMAH) R_MP values when TMAH was used**, and
verify the N-source narrative survives removal.

#### ⚠️ Spectral caveat: name-based removal can MISS TMAH reagent peaks

The `TMAH_ARTIFACTS` list is **name-based**, so a TMAH reagent peak that NIST
misidentifies under a different name escapes removal and stays in `Other_N`
(or `Unknown`). Trimethylamine (base peak **m/z 58**, with m/z 59 = M+1 and
m/z 42) is the diagnostic reagent signature — *regardless of the NIST name*.

Confirmed cases (2026-08-01, QGD EI spectra, all base m/z 58):
- **POC BC15** RT 2.664 "Butanoic acid, 4-(dimethylamino)-3-hydroxy-" — 33%
  of the sample, actually trimethylamine; was inflating `Unknown` to 40%.
- **MAOC BC7.5** RT 2.470 "1-Methyldodecylamine" (6.8%) and RT 2.625
  "Tetramethylammonium acetate" (7.3%, the TMAH cation itself) — together
  14% of Other_N.
- POC BC30 RT 2.573 "Methylamine, N,N-dimethyl-" (30.5%) *was* caught by the
  name list — the same compound with a different NIST name in BC15 was not.

**Workflow**: when QGD raw files exist, scan each treatment's EI spectra for
peaks whose base peak is m/z 58 (trimethylamine/tetramethylammonium) and
remove them as TMAH reagent peaks, independent of the NIST-assigned name.
Use the generic CLI (`diag_trimethylamine.py` is now part of this skill):

```
python scripts/diag_trimethylamine.py --qgd <QGD_dir> --txt <TXT_dir> \
    --sample_map sample_map.json [--rt_tolerance 0.35] [--tic_min 50000] \
    [--output tmah_peaks.json]
```

The spectral test is: base peak m/z 58±0.5 AND m/z 59/base >0.10 AND
m/z 42/base >0.15 AND frac58 >0.12 → TMAH reagent peak, regardless of NIST
name. Feed the output JSON into `corrections.json` or use it for targeted
FINAL removal.
The distinction from *genuine* N-methylated soil compounds (e.g.
N,N,N'-trimethylethylenediamine, kept by user decision) is the spectrum: a
pure TMA peak shows only the m/z 58/59/42 cluster (plus CO₂ m/z 44
background), while a real derivatized amine also has higher-mass fragments.

### Sensitivity Scenarios

Report class composition under multiple cleaning levels:

| Scenario | Filter | Purpose |
|----------|--------|---------|
| S0 | Raw (all peaks) | Reference — shows problems |
| S1 | SI ≥ 70 | Remove low-confidence names |
| S2 | S1 + anomalous peak area → NA | Basic cleaning |
| **S3 (recommended)** | **S2 + re-normalize per sample** | **Main result** |
| S4 | S3 + all integration-risk peaks removed | Robustness check |

### Data Reconciliation

After cleaning, verify numerical closure:

1. **Peak count**: Raw spectra = 4×N₄ + 3×N₃ + 2×N₂ + N₁ (where Nₓ = nodes in x treatments)
2. **SI filter**: Raw peaks − SI<70 deleted = retained peaks
3. **Network metrics**: Mean degree k̄ = 2E/N (must match exactly)
4. **Annotation propagation**: New_classified = After_total − Before_total
5. **Cross-fraction pairs**: Σ(26 + 6 + 8 + 10 + 1) = 51 (all spectral pairs accounted for)

**Critical**: The number of peaks used for class composition (SI≥70, typically ~900)
is NOT the same as the number of EI spectra used for network nodes (all spectra, ~1050).
State both numbers explicitly and explain the difference in Methods.

## EI Molecular Networking (Shahriar 2026 Method)

Molecular networking based on deconvoluted EI fragmentation spectra, following
Shahriar et al. (2026) ES&T methodology adapted for Shimadzu GC-MS data.

### Workflow Summary

```
Raw TXT (peak table + EI spectra)
  → Parse m/z:intensity from [MS Spectrum] blocks
  → m/z binning (nominal mass)
  → Intensity filtering (≥1% relative)
  → Cross-treatment alignment (|ΔRT| ≤0.10 min, EI cosine ≥0.85, ≥6 matched ions)
  → Consensus EI spectra (MGF export)
  → Network construction (cosine ≥0.70, mutual TopK=10, ≥6 matched ions)
  → Class propagation to unannotated nodes
  → GraphML export for Cytoscape visualization
```

### Key Parameters

| Parameter | Main Value | Sensitivity Range |
|-----------|-----------|-------------------|
| Pair cosine threshold | 0.70 | 0.65, 0.70, 0.75, 0.80 |
| TopK | 10 | 5, 10, 15 |
| Minimum matched ions | 6 | 4, 6, 8 |
| Ion intensity filter | 1% relative | 1%, 2%, 5% |
| RT alignment tolerance | 0.10 min | 0.05, 0.10, 0.20 |
| Alignment cosine | 0.85 | 0.80, 0.85, 0.90 |

### Annotation Propagation

Unannotated nodes (SI<70 or no NIST match) receive chemical class from their
nearest annotated neighbor(s) in the molecular network:

- Only propagate **chemical class** (e.g., "lignin-related"), never specific compound names
- Require highest neighbor cosine ≥ 0.70; weighted majority class ≥ 60% of annotated neighbors
- Mark conflicts as "class ambiguous"
- Report as "putative [class]-related spectral type" in papers
- Always preserve original NIST annotation vs. propagated class in separate columns

### Cross-Fraction Spectral Matching (MAOC ↔ POC)

When comparing two carbon pools (e.g., MAOC vs POC):

```
Mutual best match, EI cosine ≥ 0.90, matched ions ≥ 6
```

Classify matched pairs into 5 categories:
- **Same candidate compound** (same CAS, name, class) — can discuss carbon pool partitioning
- **Isomer/naming variant** (same formula, different CAS) — merge as compound family
- **Homolog/broad class** (same class, different chain length) — same structural class
- **Cross-class high similarity** — requires manual verification (possible misannotation)
- **Same class, ambiguous identity** — specific structure pending

### Output Files

For each carbon pool and for cross-pool comparison:
- Consensus EI spectra (`.mgf`)
- Feature quantification matrix (`.csv`)
- Network node table (class, SI level, QC flags, per-treatment abundance)
- Network edge table (cosine, matched ions)
- GraphML for Cytoscape (`_network.graphml`)
- Annotation propagation log (before/after counts)

### When to Network

- **Always**: When comparing molecular composition across treatments or carbon pools
- **Recommended**: ≥3 treatments with ≥100 peaks each
- **Requires**: Per-peak EI m/z:intensity data (from TXT [MS Spectrum] blocks or QGD extraction)
- **Not possible with**: Peak table only (compound names + areas, no spectra)

## Cross-Interpretation with TG-DSC

When TG-DSC thermal data is available for the same samples, joint interpretation
reveals **the molecular basis of thermal stability.**

### Quick cross-reference

| Py-GC-MS signal | Check TG-DSC for | Interpretation |
|---|---|---|
| ↑ Condensed aromatics | ↑ TG-T50 (Air), ↑ ΔT50 | Aromatic compounds provide antioxidant shielding |
| ↑ Lignin-derived phenols | ↑ TG-T50 (moderate), ↑ stable pool % | Plant structural C, moderate thermal stability |
| ↓ Polysaccharides | ↓ Labile pool % | Loss of easily-pyrolyzed compounds confirmed |
| ↓ Lipids | ↓ Energy Density (ED) | Loss of energy-rich reduced compounds |
| ↑ N-heterocycles | ↑ MAOC-associated TG-T50 | Microbial-processed C stabilized on minerals |

### When to cross-interpret

- **Always**: When both Py-GC-MS and TG-DSC data exist for the same samples/treatments
- **Caution**: When TG-DSC is from single atmosphere only (no ΔT50)
- **Skip**: When Py-GC-MS compounds have low SI (<80) or bulk soil TG deconvolution R² < 0

### Full cross-interpretation guide
See `biochar-soc-knowledge/references/tgdsc-pygcms-crosswalk.md` for:
- Complete cross-interpretation matrix (both directions)
- Joint diagnostic patterns (4 common scenarios)
- Recommended analysis workflow
- Reporting templates for convergent findings

> **Cross-skill**: When the user has or mentions TG-DSC data alongside Py-GC-MS,
> invoke both `pygcms-batch` and `tg-dsc-analysis`, then use the crosswalk reference
> for joint interpretation.

## Required Python packages
openpyxl, olefile (stdlib: os, re, json, argparse, time, urllib, ssl, collections, struct)

## Reference papers
- Chen et al. 2023, Carbon Research 2:1. DOI: 10.1007/s44246-022-00034-0
- Kallenbach et al. 2016, Nature Communications 7:13630. DOI: 10.1038/ncomms13630
- Shahriar et al. 2026, Environ. Sci. Technol. 60, 9237-9249. DOI: 10.1021/acs.est.5c16415
- LaRowe & Van Cappellen 2011, Geochimica et Cosmochimica Acta (NOSC/DG_COX framework)

## Internal References

- `references/shahriar_2026_classification.md` — Shahriar 12-class scheme details
- `references/chen_2023_classification.md` — Chen 2023 10-class scheme details
- `references/kallenbach_2016_classification.md` — Kallenbach 2016 8-class scheme details
- **`references/data-cleaning-molecular-networking.md`** — Complete post-export workflow: anomalous peak detection, classification correction, **TMAH thermochemolysis artifact removal**, sensitivity scenarios, data reconciliation, EI molecular networking, cross-fraction spectral matching, annotation propagation rules, writing guidelines
