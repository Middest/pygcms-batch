# Py-GC-MS Data Cleaning, Verification & Molecular Networking

## Overview

This document captures the complete post-export workflow developed through the paddy soil biochar Py-GC-MS project. It covers the full chain from raw TXT parsing through molecular network construction and cross-fraction spectral matching.

---

## 1. Data Cleaning Pipeline

### 1.1 SI Threshold Strategy (Dual-Level)

The Shahriar 2026 paper uses SI≥70 for broad coverage (~50.6% annotation rate), while GNPS
defaults recommend SI≥85. Our pipeline uses a **dual-threshold system**:

| Level | SI Range | Action |
|-------|----------|--------|
| High confidence | SI ≥ 85 | Retain compound name + class for reporting |
| Moderate confidence | 70 ≤ SI < 85 | Retain name as "candidate", flag for manual review |
| Low confidence (spectrum only) | SI < 70 | **Delete compound name**, retain EI spectrum as "Unannotated" node for molecular networking |

**Critical rule**: SI<70 peaks are NOT deleted from EI molecular network input — only their
compound names are removed. The experimental EI spectrum is real data and can still form
edges and receive class propagation from annotated neighbors.

### 1.2 Anomalous Peak Detection

Automated checks for peaks that may distort quantitative results:

| Check | Criterion | Action |
|-------|-----------|--------|
| Area/Height ratio | A/H > 20 (or >5× sample median) | Flag as "wide peak — possible co-elution" |
| Integration window | Proc.To − Proc.From > 1.0 min (or >10× sample median) | Flag as "over-integration — contains multiple peak apexes" |
| Single-peak dominance | Single peak >30% of total sample area | Flag as "dominant peak — verify blank, check for derivatization artifact" |
| RT-RI mismatch | NIST RI differs from experimental RT-based estimated RI by >200 units | Flag as "RT-RI mismatch — possible misidentification" |
| Cross-class at same RT | Same RT assigned to different chemical classes across treatments | Flag as "cross-class conflict — verify via EI spectrum" |

**Default action for flagged anomalous peaks**:
- Retain EI spectrum for molecular networking (as Unannotated or with "Candidate" tag)
- Set peak area to **NA** (NOT zero) for class composition calculations
- NA = "peak exists but area is unreliable"; 0 = "compound not detected"
- After setting NAs, re-normalize remaining peaks within each sample

### 1.3 Classification Correction Rules

Common systematic classification errors found in Py-GC-MS data:

| Pattern | Incorrect Classification | Correct Classification | Rationale |
|---------|--------------------------|------------------------|-----------|
| p-/o-Cymene (CAS 99-87-6, 527-84-4) | Alkenes | MAH | Aromatic ring with isopropyl + methyl substituents |
| 1-Nonene (CAS 124-11-8) | Alcohols/ketones/oxygenates | Alkenes | Straight-chain terminal alkene |
| Hexahydro-triazines (CAS 108-74-7) | N-MAH | Other_N | Saturated six-membered ring, NOT aromatic |
| Haloalkanes (Br/I-substituted) | Short/Long alkanes | Other_hydrocarbons | Halogen substitution changes chemical behavior |
| Benzonitrile (CAS 100-47-0) | Other_N | N-MAH | Cyano group attached to aromatic ring |
| Tetratetracontane (CAS 7098-22-8) | Other_hydrocarbons | Long_alkanes | C44 n-alkane |
| Compounds with invalid CAS (0-00-0) | Specific name | "Spectral type — specific name pending" | Downgrade annotation confidence |

**Workflow**: After automatic classification, run `verify_data.py --check-classification`
to detect these known patterns, then generate a corrections log before recomputing all
class composition statistics.

### 1.4 TMAH Thermochemolysis Artifact Removal

When samples were derivatized with **TMAH** (tetramethylammonium hydroxide,
20% in MeOH, pyrolyzed at 500°C), the reagent decomposes and its fragments
condense with formaldehyde (from carbohydrates) into N-methylated products
that are **reagent chemistry, not soil N**:

| TMAH artifact | NIST name | Mechanism |
|---------------|-----------|-----------|
| Trimethylamine | Methylamine, N,N-dimethyl- | TMAH Hofmann elimination |
| Tetramethylmethanediamine | Methanediamine, N,N,N',N'-tetramethyl- | TMA + formaldehyde |
| Trimethyltriazine | 1,3,5-Triazine, hexahydro-1,3,5-trimethyl- | TMA + formaldehyde cyclization |
| Triazine dione | 1,2,4-triazine, hexahydro-3,5-dione | TMA + formaldehyde cyclization |
| Methenamine | Methenamine (hexamethylenetetramine) | TMA + formaldehyde condensation |
| Dimethylaminoacetonitrile | Acetonitrile, (dimethylamino)- | Formaldehyde + cyanide |

**Why it matters**: `classify_compound()`'s keyword fallback sends compounds
containing `amine/amide/triazine/methanediamine` into **Other_N → Microbial**.
TMAH artifacts therefore inflate the microbial signal and R_MP ratio. Observed
impact (paddy soil biochar, 2026): Bulk/POC Other_N dropped 35-56 pp after
removal; MAOC was nearly unaffected (Δ < 9 pp).

**Workflow**:
1. `pipeline.py` removes these by default (`TMAH_ARTIFACTS` list).
2. `--keep_tmah` keeps them for a control run.
3. After removal, **renormalize** class composition to 100% of kept-peak total
   (default). Compare control vs. removed R_MP before interpreting any N-source signal.

**Always** run with TMAH data. **Report both** control and cleaned R_MP values;
verify the N-source narrative survives removal before writing Discussion.

### 1.5 Sensitivity Scenarios for Class Composition

When anomalous peaks are present, report class composition under multiple scenarios:

| Scenario | Description | Use Case |
|----------|-------------|----------|
| S0 | All peaks, raw | Reference only — shows the problem |
| S1 | SI ≥ 70 only | Remove low-confidence compound names |
| S2 | S1 + remove anomalous wide peaks | Basic cleaning |
| **S3** | **S2 + re-normalize per sample** | **Recommended main result (conservative)** |
| S4 | S3 + remove all integration-risk peaks (A/H > 10, window > 0.5 min) | Robustness check |

Report the **range** of key class percentages across S2–S4. If the main signal direction
is stable across scenarios, the result is robust.

---

## 2. Data Reconciliation & Closure Checks

### 2.1 Peak Count Reconciliation

Verify that peak counts are consistent across processing stages:

```
Raw EI spectra = Σ(per-sample [MS Spectrum] blocks)
SI<70 deleted = Σ(per-sample Hit#1 SI<70 counts)
Retained = Raw − Deleted
Consensus nodes = aligned feature count (NOT the same as retained peaks)
```

**Check**: `4×N₄ + 3×N₃ + 2×N₂ + 1×N₁ = total input spectra`
where Nₓ = nodes detected in x out of 4 treatments.

This verifies that network construction preserved all input spectra.

### 2.2 Network Metric Closure

Verify internal consistency of network metrics:

- **Mean degree**: k̄ = 2E/N (must match exactly, not approximately)
- **Annotation rates**: Total_nodes × rate% ≈ integer count of classified nodes
- **Propagation gain**: New_classified = After − Before (must match per-node count)
- **51-pair classification**: Σ(26 + 6 + 8 + 10 + 1) = 51 (all pairs accounted for)

### 2.3 Data Source Reconciliation

Track which data subset is used for each analysis:

| Analysis | Data Source | SI Filter | Anomalous Peaks |
|----------|------------|-----------|-----------------|
| Class composition (quantitative) | SI≥70 retained peaks | Yes | Areas set to NA |
| EI molecular network (nodes) | ALL EI spectra | No (names removed) | Spectra retained, area = NA |
| Network topology (edges) | EI cosine between consensus spectra | No | Included as nodes |
| Annotation propagation | Network + NIST names | SI≥70 for source | NA |

**Do NOT mix these**: The 916 SI≥70 peaks used for class composition are NOT the same
as the 1049 EI spectra used for network nodes. State this explicitly in Methods.

---

## 3. Cross-Fraction Spectral Matching (MAOC ↔ POC)

### 3.1 Matching Criteria

For identifying shared molecular features between carbon pools (e.g., MAOC vs POC):

```
Primary: EI cosine ≥ 0.90, mutual best match, matched ions ≥ 6
Secondary: Same CAS + same compound name + same chemical class
```

### 3.2 Classification of Matched Pairs

| Category | Criteria | Count in our data | Can be reported as |
|----------|----------|-------------------|-------------------|
| Same candidate compound | Same CAS, same name, same class, cosine ≥ 0.90 | 26 | "Same candidate molecule detected in both fractions" |
| Isomer / naming variant | Same formula, same class, different CAS/positional isomer, cosine ≥ 0.90 | 6 | "Structural isomers — merge as compound family" |
| Homolog / broad class | Same class, different formula/chain length, cosine ≥ 0.90 | 8 | "Homologous series member — same structural class" |
| Cross-class high similarity | Different classes, cosine ≥ 0.90 | 10 | "High spectral similarity but class conflict — requires manual verification" |
| Same class, ambiguous identity | Same broad class, specific name/CAS uncertain | 1 | "Same structural class, specific identity pending" |

### 3.3 Cross-Class Conflict Resolution

For the 10 cross-class pairs, check:
1. Are both SI values ≥ 80? If not, the low-SI assignment is suspicious.
2. Is there an RT/RI mismatch? Systems with different RT but same EI may be different compounds.
3. Does the EI spectrum contain characteristic ions? (e.g., m/z 31 for alcohols, m/z 77 for aromatics)
4. Could this be co-elution? Check A/H ratio and integration window.
5. Are both names from the same structural family? (e.g., cyclic alkene vs aromatic with same ring size)

---

## 4. Annotation Propagation Rules

### 4.1 Conservative Propagation Algorithm

For unannotated nodes (SI<70 or no NIST match) in the EI molecular network:

```
Step 1: Find all annotated neighbors (SI≥70) within cosine ≥ 0.70
Step 2: Weight each neighbor's class vote by its EI cosine
Step 3: If weighted majority class ≥ 70% AND highest neighbor cosine ≥ 0.80:
         → Assign that class with "putative" prefix
Step 4: If weighted majority class 60–70% OR highest neighbor cosine 0.70–0.80:
         → Assign "putative [class]-related spectral type" (lower confidence)
Step 5: If two classes have weight difference < 15%:
         → Mark as "class ambiguous" — do NOT assign
```

### 4.2 What Propagation CAN and CANNOT Do

| CAN do | CANNOT do |
|--------|-----------|
| Assign chemical class to unannotated node | Give a specific compound name |
| Increase overall annotation rate by 3–7 pp | Replace NIST library matching |
| Identify class-level patterns in unknown nodes | Confirm structural identity |
| Reduce "Unknown" category for class composition | Be reported as "identified compounds" |

**Writing rule**: Propagated classes must always carry "putative" or "related spectral type"
qualifier. Never write "identified as lignin compounds" for propagated nodes.

---

## 5. Parameter Sensitivity Framework

### 5.1 Network Parameters to Test

| Parameter | Main Value | Sensitivity Range |
|-----------|-----------|-------------------|
| Pair cosine threshold | 0.70 | 0.65, 0.70, 0.75, 0.80 |
| TopK | 10 | 5, 10, 15 |
| Minimum matched ions | 6 | 4, 6, 8 |
| Ion intensity filter | 1% relative | 1%, 2%, 5% |

### 5.2 Alignment Parameters to Test

| Parameter | Main Value | Sensitivity Range |
|-----------|-----------|-------------------|
| RT tolerance | 0.10 min | 0.05, 0.10, 0.20 min |
| Alignment cosine threshold | 0.85 | 0.80, 0.85, 0.90 |

### 5.3 Reporting Sensitivity Results

- Report the **direction** and **range** of key metrics across parameter space
- Use Jaccard index to quantify alignment stability (main vs alternative parameters)
- Present as: "Across cosine thresholds 0.65–0.80, MAOC consistently had fewer edges and
  a larger maximum molecular family than POC, indicating that the topological differences
  are robust to threshold choice."
- Do NOT write "results are completely stable" — metrics inevitably change with thresholds

---

## 6. Writing Guidelines for Papers

### 6.1 Expression Hierarchy by Data Quality

| Data Support Level | Allowed Expression | Example |
|-------------------|-------------------|---------|
| ≥3 biological replicates + significant test | "significantly increased" | "BC30 significantly increased PAH relative abundance (p < 0.05)" |
| No replicates, clear monotonic trend | "showed a consistent increase" | "Long-chain alkanes showed a consistent increase with biochar rate" |
| No replicates, non-monotonic | "varied without a clear dose-response pattern" | "MAOC class composition varied without a clear dose-response pattern" |
| Single anomalous peak | "was not included in quantitative analysis due to wide integration window" | — |
| Network-propagated class | "putative [class]-related spectral type" | — |

### 6.2 Key Distinctions to Maintain

- "Relative proportion" ≠ "absolute content" (no internal standard = only relative)
- "High EI spectral similarity" ≠ "same compound" (isomers, homologs can have similar EI)
- "Detected only in treatment X" ≠ "treatment-specific biomarker" (without replicates)
- "Putative class assignment" ≠ "identified compound" (for propagated annotations)

---

## 7. Deliverables Checklist

For a complete Py-GC-MS molecular network analysis, produce:

### Core data files
- [ ] Raw peak table with SI, CAS, formula, RT, area per sample
- [ ] SI≥70 filtered and anomalous-peak-cleaned quantitative table
- [ ] Classification corrections log
- [ ] Class composition under S0–S4 scenarios
- [ ] Cross-treatment aligned peak matrix (RT + EI cosine constrained)

### Molecular network files
- [ ] Consensus EI spectra (MGF format)
- [ ] Feature quantification matrix (CSV)
- [ ] Network node table (with class, SI level, QC flags)
- [ ] Network edge table (with cosine, matched ions)
- [ ] GraphML for Cytoscape
- [ ] Annotation propagation log (before/after counts, propagation rules applied)

### Verification & quality control
- [ ] Peak count reconciliation (raw → cleaned → network)
- [ ] Network metric closure check (degree, annotation rates, propagation)
- [ ] Cross-fraction spectral pair classification (26+6+8+10+1)
- [ ] Parameter sensitivity summary

### Figures
- [ ] Full molecular network (nodes colored by class, sized by abundance)
- [ ] Network topology comparison (MAOC vs POC)
- [ ] Cosine threshold sensitivity (edges and nodes vs threshold)
- [ ] Key compound class dose-response trends
- [ ] Annotation rate improvement (before/after propagation)

---

*Last updated: 2026-07-27. Based on paddy soil biochar Py-GC-MS MAOC/POC analysis pipeline.*
