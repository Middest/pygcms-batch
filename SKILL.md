---
name: pygcms-batch
description: Batch process Py-GC-MS NIST export TXT files. Filter by SI>=80, classify compounds by Chen 2023 (Carbon Research) and Kallenbach 2016 (Nature Communications) schemes, compute elemental atoms (C/H/O/N/P/S), NOSC, ΔG_COX, and source attribution (plant/microbial/mixed). Output multi-sheet Excel.
---

# Py-GC-MS Batch Analysis

Batch process Py-GC-MS NIST search export TXT files for soil organic matter molecular characterization.

## Workflow

1. Auto-detect all `.TXT` files in input directory (recursively)
2. Parse `[MC Peak Table]` for TIC peak areas and retention times
3. Parse `[MS Similarity Search Results for Spectrum Process Table]` for NIST library matches
4. Filter compounds: SI >= 80 (NIST Match Factor)
5. Merge MC area data with NIST match data (SI, CAS, molecular formula, molecular weight)
6. Parse molecular formula → atom counts (C, H, O, N, P, S), detect heteroatoms (Cl, Br, F, I, Si, B)
7. Compute NOSC and ΔG_COX for each valid formula
8. **(Optional)** Query NIST Chemistry WebBook by CAS/name to fill missing molecular formulas, InChI, InChIKey
9. Classify compounds by TWO schemes:
   - **Chen 2023** (Carbon Research): 10 chemical groups + 3 source categories
   - **Kallenbach 2016** (Nature Communications): 8 origin-based classes
10. Classification uses 286-entry exact lookup table (built from manually curated data) with keyword fallback
11. Generate source attribution statistics (plant/microbial/mixed %) by relative peak area
12. Output multi-sheet Excel workbook

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

## Required Python packages
openpyxl (stdlib: os, re, json, argparse, time, urllib, ssl, collections)

## Reference papers
- Chen et al. 2023, Carbon Research 2:1. DOI: 10.1007/s44246-022-00034-0
- Kallenbach et al. 2016, Nature Communications 7:13630. DOI: 10.1038/ncomms13630
- LaRowe & Van Cappellen 2011, Geochimica et Cosmochimica Acta (NOSC/ΔG_COX framework)
