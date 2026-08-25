# Shahriar 2026 Compound Classification Scheme

Source: Shahriar et al. (2026) *Environ. Sci. Technol.* 60, 9237-9249.
DOI: 10.1021/acs.est.5c16415

## 12 Compound Classes

| # | Class | Abbreviation | Source | Description |
|---|-------|-------------|--------|-------------|
| 1 | Monocyclic aromatic hydrocarbons | MAH | Mixed | Benzene, toluene, xylene, alkylbenzenes |
| 2 | N-containing MAH | N-MAH | Microbial | Pyrrole, pyridine, indole derivatives |
| 3 | Polycyclic aromatic hydrocarbons | PAH | Mixed | Naphthalene, fluorene, anthracene |
| 4 | Lignin | Lignin | Plant | Guaiacol, syringol, vanillin, eugenol |
| 5 | Phenols | Phenols | Mixed | Phenol, cresol, butylphenol |
| 6 | Degraded saccharides | Sugars | Mixed | Furfural, levoglucosan, furan derivatives |
| 7 | Fatty acids, alcohols and esters | Fatty_acids_lipids | Mixed | Long-chain acids, alcohols, esters |
| 8 | Alkenes | Alkenes | Plant | Unsaturated hydrocarbons C6-C30 |
| 9 | Short alkanes | Short_alkanes | Mixed | Saturated hydrocarbons C5-C18 |
| 10 | Long alkanes | Long_alkanes | Plant | Saturated hydrocarbons C19-C35 |
| 11 | Other hydrocarbons | Other_hydrocarbons | Mixed | Cyclic, branched, bicyclic, alkynes |
| 12 | Other N-containing compounds | Other_N | Microbial | Nitriles, amines, amides, N-heterocycles |

## Source Attribution

Only 5 of 12 classes reliably attributed (Section 3.2):
- Plant: Lignin, Long alkanes, Alkenes
- Microbial: N-MAH, Other N-containing
- Mixed: All others

R_MP = (N-MAH + Other_N) / (Lignin + Long_alkanes + Alkenes)

## Custom Library

1,678 compounds in `data/shahriar_library.json` (from SI Table S2).

## Key Caveats

1. Relative abundance based on peak heights, NOT carbon mass
2. Pyrolysis at 500C causes bond cleavage and secondary reactions
3. Single precursor can produce multiple products; different precursors can produce same product
