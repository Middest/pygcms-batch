---
name: pygcms-batch
description: Batch process Py-GC-MS NIST export TXT files. Filter by SI>=80, classify compounds by Chen 2023 (Carbon Research) and Kallenbach 2016 (Nature Communications) schemes, compute elemental atoms (C/H/O/N/P/S), NOSC, ΔG_COX, and source attribution (plant/microbial/mixed). Output multi-sheet Excel.
---

# Py-GC-MS Batch Analysis

Batch process Py-GC-MS NIST search export TXT files for soil organic matter molecular characterization.

## Workflow

1. Auto-detect all `.TXT` files in input directory (recursively), extract MC Peak Table
2. Filter compounds: SI >= 80
3. Parse molecular formula → atom counts (C, H, O, N, P, S)
4. Compute NOSC and ΔG_COX for each valid formula
5. Classify compounds by TWO schemes:
   - **Chen 2023** (Carbon Research): 10 chemical groups + 3 source categories
   - **Kallenbach 2016** (Nature Communications): 8 origin-based classes
6. Generate source attribution statistics (plant/microbial/mixed %)
7. Output multi-sheet Excel workbook

## Usage

```
python scripts/pygcms_batch.py --input <TXT_directory> --output <output.xlsx> [--sample_map <mapping.json>]
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
- `筛选统计` — peak filtering summary per sample
- `元素汇总` — weighted elemental composition
- `Chen_分类详细` — compound-level Chen classification
- `Chen_来源统计` — source attribution statistics (Chen)
- `Kallenbach_分类详细` — compound-level Kallenbach classification  
- `Kallenbach_来源统计` — source attribution statistics (Kallenbach)
- `分子公式汇总` — all molecular formulas with atom counts, NOSC, ΔG_COX

## Required Python packages
openpyxl, pandas (optional for large datasets)

## Reference papers
- Chen et al. 2023, Carbon Research 2:1. DOI: 10.1007/s44246-022-00034-0
- Kallenbach et al. 2016, Nature Communications 7:13630. DOI: 10.1038/ncomms13630
