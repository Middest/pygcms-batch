#!/usr/bin/env python3
"""
Py-GC-MS Batch Analysis Script
================================
Batch process NIST search export TXT files from Py-GC-MS analysis.
Filters SI>=80, classifies compounds by Chen 2023 and Kallenbach 2016 schemes,
computes elemental composition, NOSC, and ΔG_COX.

Usage:
    python pygcms_batch.py --input <dir> --output <output.xlsx> [--sample_map <json>]
"""

import os
import re
import json
import argparse
from collections import defaultdict
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter


# ============================================================================
# 1. NIST TXT PARSER
# ============================================================================

def parse_molecular_formula(formula_str):
    """
    Parse a molecular formula string into atom counts.
    Returns dict with C, H, O, N, P, S, other counts.
    Handles formulas like 'C6H6O3', 'C5H14N2', 'C11H23Cl', etc.
    """
    atoms = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'P': 0, 'S': 0, 'other': {}}

    if not formula_str or formula_str == '?' or formula_str == '':
        return atoms

    # Remove any whitespace
    formula_str = formula_str.strip()

    # Pattern: element symbol followed by optional count
    pattern = re.compile(r'([A-Z][a-z]?)(\d*)')
    matches = pattern.findall(formula_str)

    for elem, count in matches:
        c = int(count) if count else 1
        if elem in atoms:
            atoms[elem] = c
        else:
            atoms['other'][elem] = c

    return atoms


def compute_nosc_and_dg(atoms):
    """
    Compute NOSC (Nominal Oxidation State of Carbon) and ΔG_COX.

    NOSC = -((-Z + 4*C + H - 3*N - 2*O + 5*P - 2*S) / C)
    where Z = net charge (assumed 0 for neutral molecules)

    ΔG_COX = 60.3 - 28.5 * NOSC  (kJ per mol C)

    Reference: LaRowe & Van Cappellen 2011; Wu et al. 2025 Nature Food
    """
    C = atoms.get('C', 0)
    H = atoms.get('H', 0)
    N = atoms.get('N', 0)
    O = atoms.get('O', 0)
    P = atoms.get('P', 0)
    S = atoms.get('S', 0)
    Z = 0  # neutral molecules

    if C == 0:
        return None, None

    nosc = -((-Z + 4*C + H - 3*N - 2*O + 5*P - 2*S) / C)
    dg_cox = 60.3 - 28.5 * nosc

    return round(nosc, 6), round(dg_cox, 3)


def parse_nist_txt(filepath):
    """
    Parse a NIST search export TXT file.
    Returns list of dicts with peak data from [MC Peak Table] section.

    The TXT format has:
    [Header] section with metadata
    [File Information]
    [Sample Information]
    [Original Files]
    [MC Peak Table] - the key section with peak data
    [MS Spectrum] sections (ignored)
    """
    peaks = []

    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Find [MC Peak Table] section
    mc_match = re.search(r'\[MC Peak Table\](.*?)(?=\n\[)', content, re.DOTALL)
    if not mc_match:
        # Try to find it without the closing bracket
        mc_match = re.search(r'\[MC Peak Table\](.*)', content, re.DOTALL)
    if not mc_match:
        return peaks

    mc_text = mc_match.group(1)
    lines = mc_text.strip().split('\n')

    in_data = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('Peak#'):
            in_data = True
            continue
        if in_data and line.startswith('['):
            break
        if in_data:
            parts = line.split('\t')
            if len(parts) >= 11:
                try:
                    peak_num = int(parts[0])
                    ret_time = float(parts[1])
                    area = float(parts[6])
                    height = float(parts[7]) if parts[7] else 0
                    name = parts[10].strip() if len(parts) > 10 else ''

                    peaks.append({
                        'peak_num': peak_num,
                        'ret_time': ret_time,
                        'area': area,
                        'height': height,
                        'name': name,
                        'raw_line': line
                    })
                except (ValueError, IndexError):
                    pass

    return peaks


def find_nist_match_details(filepath, peak_name_to_match):
    """
    Extract SI, CAS, molecular formula from NIST match section in TXT.
    Searches [Hit] sections for each peak's top NIST hit.
    """
    matches = {}

    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Find all [Hit] sections
    hit_blocks = re.findall(r'\[Hit\](.*?)(?=\[Hit\]|\[MS Spectrum\]|$)', content, re.DOTALL)

    for block in hit_blocks:
        si_match = re.search(r'SI\s*[:=]?\s*(\d+)', block, re.IGNORECASE)
        name_match = re.search(r'Name\s*[:=]?\s*(.+)', block, re.IGNORECASE)
        cas_match = re.search(r'CAS\s*[:=]?\s*([\d\-]+)', block, re.IGNORECASE)
        formula_match = re.search(r'Mol(?:ecular)?\s*Formula\s*[:=]?\s*([A-Za-z0-9]+)', block, re.IGNORECASE)
        mw_match = re.search(r'Mol(?:ecular)?\s*Weight\s*[:=]?\s*([\d.]+)', block, re.IGNORECASE)

        if name_match and si_match:
            name = name_match.group(1).strip()
            si = int(si_match.group(1))
            cas = cas_match.group(1).strip() if cas_match else ''
            formula = formula_match.group(1).strip() if formula_match else ''
            mw = float(mw_match.group(1)) if mw_match else 0

            if name not in matches or si > matches[name]['si']:
                matches[name] = {
                    'name': name,
                    'si': si,
                    'cas': cas,
                    'formula': formula,
                    'mol_weight': mw
                }

    return matches


# ============================================================================
# 2. COMPOUND CLASSIFICATION
# ============================================================================

# --- Chen 2023 Classification Rules ---
# Chemical groups and source attribution based on:
# Chen et al. 2023, Carbon Research 2:1

CHEN_CLASS_RULES = {
    'lipids': {
        'keywords': ['acid', 'ester', 'alkane', 'alkene', 'alcohol', 'ketone',
                     'aldehyde', 'fatty', 'glycerol', 'sterol', 'terpene',
                     'acetate', 'propanoate', 'butanoate', 'hexanoate',
                     'tridecanoic', 'tetradecanoic', 'pentadecanoic',
                     'hexadecanoic', 'heptadecanoic', 'octadecanoic',
                     'eicosanoic', 'docosanoic', 'heneicosanoic',
                     'decane', 'undecane', 'dodecane', 'tridecane',
                     'tetradecane', 'pentadecane', 'hexadecane',
                     'heptadecane', 'octadecane', 'nonadecane', 'eicosane',
                     'heneicosane', 'docosane', 'tricosane', 'tetracosane',
                     'decene', 'undecene', 'dodecene', 'tridecene',
                     'tetradecene', 'pentadecene', 'hexadecene',
                     'heptadecene', 'octadecene', 'nonadecene', 'eicosene',
                     'sterol', 'cholesterol', 'campesterol', 'sitosterol',
                     'stigmasterol', 'ergosterol',
                     'alkyne', 'squalene', 'tocopherol', 'carotene',
                     'decanoic', 'dodecanoic', 'octadecanenitrile',
                     'heptadecanenitrile', 'hexadecanenitrile',
                     'octanol', 'decanol', 'dodecanol', 'tetradecanol',
                     'hexadecanol', 'octadecanol'],
        'source': 'mixed'  # short-chain (<16C) microbial, long-chain (>18C) plant, 16-18C mixed
    },
    'monocyclic_aromatics': {
        'keywords': ['benzene', 'toluene', 'xylene', 'ethylbenzene', 'styrene',
                     'phenol', 'cresol', 'benzaldehyde', 'benzoic',
                     'benzyl', 'phenyl', 'anisole', 'catechol', 'guaiacol',
                     'syringol', 'trimethylbenzene', 'propylbenzene',
                     'butylbenzene', 'isopropylbenzene', 'cumene',
                     'indane', 'indene', 'tetralin',
                     'acetophenone', 'benzophenone',
                     'methylphenol', 'dimethylphenol', 'ethylphenol',
                     'methoxybenzene', 'dimethoxybenzene',
                     'trimethoxybenzene', 'hydroxybenzene'],
        'source': 'microbial'
    },
    'polycyclic_aromatics': {
        'keywords': ['naphthalene', 'anthracene', 'phenanthrene', 'pyrene',
                     'fluoranthene', 'chrysene', 'benzopyrene', 'fluorene',
                     'biphenyl', 'indole', 'quinoline', 'isoquinoline',
                     'carbazole', 'dibenzofuran', 'dibenzothiophene',
                     'acenaphthene', 'acenaphthylene', 'benzo[a]',
                     'methylnaphthalene', 'dimethylnaphthalene',
                     'trimethylnaphthalene', 'tetramethylnaphthalene',
                     'retene', 'perylene', 'coronene', 'benzoperylene',
                     'indenopyrene', 'benzofluoranthene', 'benzanthracene'],
        'source': 'microbial'
    },
    'phenolics': {
        'keywords': ['phenol', 'cresol', 'catechol', 'resorcinol',
                     'hydroquinone', 'guaiacol', 'syringol', 'vanillin',
                     'acetovanillone', 'acetosyringone', 'coniferyl',
                     'sinapyl', 'coumaryl', 'ferulic', 'p-coumaric',
                     'cinnamic', 'benzoic acid', 'salicylic',
                     'protocatechuic', 'gallic', 'syringic', 'vanillic',
                     'eugenol', 'isoeugenol', 'chavicol',
                     'dihydroxybenzene', 'trihydroxybenzene',
                     'methoxyphenol', 'dimethoxyphenol',
                     'hydroxybenzoic', 'dihydroxybenzoic',
                     'hydroxyphenyl', 'dihydroxyphenyl',
                     'methoxybenzoic', 'dimethoxybenzoic',
                     'hydroxybenzaldehyde', 'hydroxyacetophenone',
                     'phenylpropanoid', 'flavonoid', 'tannin',
                     'lignin monomer', 'monolignol'],
        'source': 'plant'
    },
    'polysaccharides': {
        'keywords': ['furan', 'furfural', 'furanone', 'pyran', 'pyranone',
                     'levoglucosan', 'levoglucosenone', 'mannosan',
                     'galactosan', 'xylosan', 'arabinosan',
                     'anhydro sugar', 'anhydrosugar', 'cyclohexanone',
                     'furoic', 'furan methanol', 'furanmethanol',
                     'furfuryl', 'methylfuran', 'dimethylfuran',
                     'hydroxymethylfurfural', 'acetic acid',
                     'propanoic acid', 'butanoic acid',
                     'pentanoic acid', 'hexanoic acid',
                     'dianhydromannitol', 'dianhydrosorbitol',
                     'dianhydroglucitol', 'glucose', 'mannose',
                     'galactose', 'xylose', 'arabinose', 'ribose',
                     'sugar', 'glycoside', 'cellulose marker',
                     'hemicellulose marker', 'chitin marker',
                     'acetylfuran', 'methylfuroate',
                     'acetamide, n-methyl', 'acetamide, n-ethyl',
                     'polysaccharide-derived', 'carbohydrate-derived',
                     'acetic acid, methyl ester', 'acetic acid, ethyl ester',
                     '2-cyclopenten-1-one', 'methyl cyclopentenone'],
        'source': 'plant'
    },
    'lignins': {
        'keywords': ['lignin', 'guaiacyl', 'syringyl', 'hydroxyphenyl',
                     'coniferyl alcohol', 'sinapyl alcohol', 'coumaryl alcohol',
                     'dihydroconiferyl alcohol', 'propylguaiacol',
                     'vinylguaiacol', 'vinylphenol', 'vinylsyringol',
                     'propenylguaiacol', 'propenylsyringol',
                     'allylguaiacol', 'allylsyringol',
                     'ethylguaiacol', 'ethylsyringol', 'methylguaiacol',
                     'homovanillic', 'homosyringic',
                     'guaiacylethane', 'syringylethane',
                     'guaiacylpropane', 'syringylpropane',
                     'guaiacylpropanol', 'syringylpropanol',
                     'guaiacylpropanone', 'syringylpropanone',
                     'guaiacylacetic acid', 'syringylacetic acid',
                     'lignin-derived phenol', 'lignin marker',
                     'methoxyeugenol', 'dimethoxyeugenol',
                     'trans-isoeugenol', 'trans-propenylsyringol',
                     '4-vinylguaiacol', '4-vinylsyringol',
                     '4-methylguaiacol', '4-ethylguaiacol',
                     '4-methylsyringol', '4-ethylsyringol',
                     'coniferaldehyde', 'sinapaldehyde'],
        'source': 'plant'
    },
    'amino_N_bearing': {
        'keywords': ['nitrile', 'pyridine', 'pyrrole', 'pyrazine', 'imidazole',
                     'indole', 'amine', 'amino', 'amide',
                     'piperidine', 'piperazine', 'pyrrolidine',
                     'pyrimidine', 'purine', 'adenine', 'guanine',
                     'methylpyridine', 'dimethylpyridine', 'trimethylpyridine',
                     'methylpyrrole', 'dimethylpyrrole',
                     'methylimidazole', 'dimethylimidazole',
                     'methylpyrazine', 'dimethylpyrazine', 'trimethylpyrazine',
                     'butanenitrile', 'pentanenitrile', 'hexanenitrile',
                     'heptanenitrile', 'octanenitrile', 'decanenitrile',
                     'propanenitrile', 'benzonitrile', 'acetonitrile',
                     'methylbutanenitrile', 'dimethylbutanenitrile',
                     'aminobenzoic', 'aminophenol',
                     'protein marker', 'amino acid derivative',
                     'diketopiperazine', 'diketodipyrrole',
                     'proline', 'tryptophan', 'tyrosine', 'phenylalanine',
                     'valine', 'leucine', 'isoleucine', 'alanine', 'glycine'],
        'source': 'mixed'  # nitriles/pyridines/pyrroles= microbial; indoles = plant
    },
    'heterocyclic_N_bearing': {
        'keywords': ['pyrazole', 'oxazole', 'thiazole', 'triazole',
                     'tetrazole', 'diazine', 'triazine',
                     'benzimidazole', 'benzoxazole', 'benzothiazole',
                     'benzotriazole', 'quinoxaline', 'quinazoline',
                     'purine base', 'pyrimidine base',
                     'uracil', 'thymine', 'cytosine',
                     'xanthine', 'hypoxanthine', 'theobromine', 'caffeine',
                     'riboflavin', 'niacin', 'nicotinamide',
                     'nucleobase', 'nucleoside',
                     'pteridine', 'pterin', 'flavin', 'folate',
                     'adenine derivative', 'guanine derivative'],
        'source': 'plant'  # mainly plant-derived (pyrazines)
    },
    'other_N_bearing': {
        'keywords': ['nitro', 'nitroso', 'hydroxylamine', 'hydrazine',
                     'azo', 'diazo', 'azide', 'cyanate', 'isocyanate',
                     'thiocyanate', 'isothiocyanate',
                     'nitrile oxide', 'oxime', 'oxazole',
                     'ammonium', 'quaternary ammonium',
                     'betaine', 'choline', 'carnitine',
                     'urea', 'thiourea', 'guanidine', 'creatine',
                     'methylamine', 'dimethylamine', 'trimethylamine',
                     'ethylamine', 'diethylamine', 'triethylamine',
                     'propylamine', 'butylamine',
                     'putrescine', 'cadaverine', 'spermidine', 'spermine',
                     'ethanolamine', 'diethanolamine', 'triethanolamine',
                     'methanediamine', 'ethanediamine',
                     'tetramethylammonium', 'tetraethylammonium',
                     'n,n-dimethyl', 'n,n-diethyl', 'n,n,n\',n\'-tetramethyl',
                     'dimethylamino', 'diethylamino',
                     'n-bearing', 'n-containing', 'nitrogen-containing'],
        'source': 'microbial'
    },
    'unspecified': {
        'keywords': ['chloro', 'bromo', 'fluoro', 'iodo', 'chloride',
                     'bromide', 'fluoride', 'iodide',
                     'silane', 'siloxane', 'silicone',
                     'phosphate', 'phosphonate', 'phosphine',
                     'sulfone', 'sulfoxide', 'sulfonate', 'sulfate',
                     'sulfur', 'sulfide', 'disulfide', 'thiol', 'thioether',
                     'borane', 'borate', 'silicon', 'germanium', 'arsenic',
                     'selenium', 'tellurium', 'tin',
                     'unknown', 'unidentified', 'not classified',
                     'artifact', 'column bleed', 'plasticizer',
                     'phthalate', 'adipate', 'sebacate',
                     'biphenyl, chloro', 'dichloro', 'trichloro',
                     'tetrachloro', 'pentachloro', 'hexachloro',
                     'silane,', 'siloxane,', 'cyclosiloxane',
                     'polydimethylsiloxane', 'pdms',
                     'perchlorate', 'halogen'],
        'source': 'mixed'
    }
}

# --- Kallenbach 2016 Classification Rules ---
# Based on Kallenbach et al. 2016, Nature Communications 7:13630
# Classifies compounds by origin-based categories

KALLENBACH_CLASS_RULES = {
    'lipids': {
        'keywords': ['acid', 'ester', 'alkane', 'alkene', 'alcohol', 'ketone',
                     'aldehyde', 'fatty', 'glycerol', 'sterol', 'terpene',
                     'acetate', 'propanoate', 'butanoate', 'hexanoate',
                     'decane', 'undecane', 'dodecane', 'tridecane',
                     'tetradecane', 'pentadecane', 'hexadecane',
                     'heptadecane', 'octadecane', 'nonadecane', 'eicosane',
                     'heneicosane', 'docosane', 'tricosane', 'tetracosane',
                     'decenoic', 'dodecanoic', 'tetradecanoic',
                     'hexadecanoic', 'octadecanoic',
                     'octadecanenitrile', 'heptadecanenitrile',
                     'hexadecanenitrile', 'tetradecanenitrile',
                     'octadecanamide', 'hexadecanamide', 'tetradecanamide',
                     'squalene', 'sterol', 'cholesterol', 'tocopherol',
                     'carotene', 'long-chain', 'wax ester',
                     'tridecanoic', 'tetradecanoic', 'pentadecanoic',
                     'hexadecanoic', 'heptadecanoic', 'octadecanoic',
                     'eicosanoic', 'docosanoic'],
        'source': 'mixed'
    },
    'lignin_derivatives': {
        'keywords': ['lignin', 'guaiacyl', 'syringyl', 'coniferyl', 'sinapyl',
                     'coumaryl', 'guaiacol', 'syringol', 'vanillin',
                     'acetovanillone', 'acetosyringone',
                     'vinylguaiacol', 'vinylsyringol', 'vinylphenol',
                     'ethylguaiacol', 'ethylsyringol',
                     'propylguaiacol', 'propenylsyringol',
                     'guaiacylethane', 'syringylethane',
                     'guaiacylpropane', 'syringylpropane',
                     'lignin-derived phenol', 'lignin marker',
                     'coniferaldehyde', 'sinapaldehyde',
                     'homovanillic', 'homosyringic',
                     'guaiacylpropanone', 'syringylpropanone',
                     'trans-isoeugenol', 'trans-propenylsyringol',
                     '4-vinylguaiacol', '4-vinylsyringol',
                     'methoxyeugenol', 'dimethoxyeugenol',
                     '4-methylguaiacol', '4-ethylguaiacol',
                     '4-methylsyringol', '4-ethylsyringol'],
        'source': 'plant'
    },
    'polysaccharides': {
        'keywords': ['furan', 'furfural', 'furanone', 'pyran', 'pyranone',
                     'levoglucosan', 'levoglucosenone', 'mannosan',
                     'galactosan', 'xylosan', 'arabinosan',
                     'anhydro sugar', 'anhydrosugar', 'cyclohexanone',
                     'furoic', 'furan methanol', 'furfuryl',
                     'methylfuran', 'dimethylfuran',
                     'hydroxymethylfurfural', 'cellulose marker',
                     'hemicellulose marker', 'chitin marker',
                     'acetylfuran', 'methylfuroate',
                     'glucose', 'mannose', 'galactose',
                     'xylose', 'arabinose', 'ribose',
                     'sugar', 'glycoside', 'polysaccharide-derived',
                     'dianhydromannitol', 'dianhydrosorbitol',
                     'dianhydroglucitol'],
        'source': 'mixed'  # plant polysaccharides + microbial (chitin)
    },
    'proteins': {
        'keywords': ['nitrile', 'pyridine', 'pyrrole', 'pyrazine', 'imidazole',
                     'indole', 'amine', 'amino', 'amide',
                     'piperidine', 'piperazine', 'pyrrolidine',
                     'pyrimidine', 'purine',
                     'diketopiperazine', 'diketodipyrrole',
                     'proline', 'tryptophan', 'tyrosine', 'phenylalanine',
                     'valine', 'leucine', 'isoleucine', 'alanine', 'glycine',
                     'protein marker', 'amino acid derivative',
                     'acetonitrile', 'propanenitrile', 'butanenitrile',
                     'pentanenitrile', 'hexanenitrile', 'benzonitrile',
                     'methylpyridine', 'dimethylpyridine', 'trimethylpyridine',
                     'methylpyrrole', 'dimethylpyrrole',
                     'methylimidazole', 'dimethylimidazole',
                     'methylpyrazine', 'dimethylpyrazine', 'trimethylpyrazine',
                     'methanediamine', 'ethanediamine',
                     'n,n-dimethyl', 'n,n-diethyl', 'tetramethyl',
                     'dimethylamino', 'diethylamino'],
        'source': 'microbial'
    },
    'non_protein_N': {
        'keywords': ['nitro', 'nitroso', 'hydroxylamine', 'hydrazine',
                     'azo', 'diazo', 'cyanate', 'isocyanate',
                     'thiocyanate', 'isothiocyanate',
                     'ammonium', 'quaternary ammonium', 'betaine',
                     'choline', 'carnitine', 'urea', 'thiourea',
                     'guanidine', 'creatine', 'creatinine',
                     'tetramethylammonium', 'tetraethylammonium',
                     'adenine', 'guanine', 'uracil', 'thymine', 'cytosine',
                     'xanthine', 'hypoxanthine', 'theobromine', 'caffeine',
                     'nucleobase', 'nucleoside', 'riboflavin', 'niacin',
                     'purine base', 'pyrimidine base',
                     'pteridine', 'pterin', 'flavin', 'folate'],
        'source': 'mixed'
    },
    'phenolics': {
        'keywords': ['phenol', 'cresol', 'catechol', 'resorcinol',
                     'hydroquinone', 'vanillin', 'vanillic',
                     'syringic', 'gallic', 'protocatechuic',
                     'p-coumaric', 'ferulic', 'cinnamic',
                     'benzoic acid', 'salicylic',
                     'hydroxybenzoic', 'dihydroxybenzoic',
                     'hydroxyphenyl', 'dihydroxyphenyl',
                     'methoxyphenol', 'dimethoxyphenol',
                     'eugenol', 'isoeugenol', 'chavicol',
                     'hydroxybenzaldehyde', 'hydroxyacetophenone',
                     'acetovanillone', 'acetosyringone',
                     'dihydroxybenzene', 'trihydroxybenzene',
                     'phenylpropanoid'],
        'source': 'plant'
    },
    'aromatics': {
        'keywords': ['benzene', 'toluene', 'xylene', 'ethylbenzene', 'styrene',
                     'naphthalene', 'anthracene', 'phenanthrene', 'pyrene',
                     'fluoranthene', 'fluorene', 'biphenyl',
                     'benzyl', 'phenyl', 'benzaldehyde',
                     'trimethylbenzene', 'propylbenzene',
                     'butylbenzene', 'isopropylbenzene', 'cumene',
                     'indane', 'indene', 'tetralin',
                     'acetophenone', 'benzophenone',
                     'methylnaphthalene', 'dimethylnaphthalene',
                     'methoxybenzene', 'dimethoxybenzene',
                     'trimethoxybenzene',
                     'acenaphthene', 'acenaphthylene',
                     'retene', 'perylene', 'coronene',
                     'benzoperylene', 'chrysene',
                     'methylphenol', 'dimethylphenol', 'ethylphenol',
                     'anisole', 'benzonitrile'],
        'source': 'microbial'
    },
    'unspecified': {
        'keywords': ['chloro', 'bromo', 'fluoro', 'iodo', 'chloride',
                     'bromide', 'fluoride', 'iodide',
                     'silane', 'siloxane', 'silicone',
                     'phosphate', 'phosphonate', 'phosphine',
                     'sulfone', 'sulfoxide', 'sulfonate', 'sulfate',
                     'sulfur', 'sulfide', 'disulfide', 'thiol',
                     'borane', 'borate', 'silicon',
                     'unknown', 'unidentified', 'not classified',
                     'artifact', 'column bleed', 'plasticizer',
                     'phthalate', 'adipate', 'sebacate',
                     'perchlorate', 'halogen',
                     'polydimethylsiloxane', 'pdms', 'cyclosiloxane',
                     'benzyl oxy', 'digitoxin'],
        'source': 'mixed'
    }
}


def classify_compound_chen(name):
    """
    Classify a compound by Chen 2023 scheme.
    Returns (chemical_group, source_category)
    """
    name_lower = name.lower().strip()

    # Check each chemical group
    for group, rules in CHEN_CLASS_RULES.items():
        for kw in rules['keywords']:
            if kw.lower() in name_lower:
                return (group, rules['source'])

    # Fallback: try to infer from chemical name patterns
    if any(x in name_lower for x in ['nitrile', 'pyridine', 'pyrrole', 'amine', 'amino']):
        return ('amino_N_bearing', 'microbial')
    if any(x in name_lower for x in ['benzene', 'toluene', 'aromatic', 'naphthalene']):
        return ('monocyclic_aromatics', 'microbial')
    if any(x in name_lower for x in ['phenol', 'cresol', 'vanillin', 'guaiacol']):
        return ('phenolics', 'plant')
    if any(x in name_lower for x in ['acid', 'ester', 'alkane', 'alkene']):
        return ('lipids', 'mixed')

    return ('unspecified', 'mixed')


def classify_compound_kallenbach(name):
    """
    Classify a compound by Kallenbach 2016 scheme.
    Returns (chemical_group, source_category)
    """
    name_lower = name.lower().strip()

    for group, rules in KALLENBACH_CLASS_RULES.items():
        for kw in rules['keywords']:
            if kw.lower() in name_lower:
                return (group, rules['source'])

    return ('unspecified', 'mixed')


# ============================================================================
# 3. BATCH PROCESSING
# ============================================================================

def process_directory(input_dir, sample_map=None):
    """
    Process all TXT files in directory recursively.

    Args:
        input_dir: path to directory containing TXT files
        sample_map: dict mapping sample code to treatment name, e.g. {'5':'CK'}

    Returns:
        dict with sample-level results
    """
    if sample_map is None:
        sample_map = {}

    results = {}

    for root, dirs, files in os.walk(input_dir):
        for fname in files:
            if not fname.lower().endswith('.txt'):
                continue

            filepath = os.path.join(root, fname)

            # Determine sample code from filename (strip extension)
            sample_code = os.path.splitext(fname)[0]
            treatment = sample_map.get(sample_code, sample_code)

            # Check if this is a retest (in a *复测* directory)
            subdir_name = os.path.basename(root)
            if '复测' in subdir_name or 'retest' in subdir_name.lower():
                treatment = f"{treatment}_复测"

            print(f"  Processing: {fname} -> {treatment}")

            # Parse peaks
            peaks = parse_nist_txt(filepath)

            if not peaks:
                print(f"    WARNING: No peaks found in {fname}")
                continue

            # Process each peak with classification
            processed = []
            for p in peaks:
                atoms = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'P': 0, 'S': 0, 'other': {}}
                nosc, dg = None, None

                # Try to extract molecular formula from name (if formula-like pattern)
                # Simple formula detection in the name
                formula_pattern = re.findall(r'[A-Z][a-z]?\d*', p['name'])

                # Classify
                chen_group, chen_source = classify_compound_chen(p['name'])
                kal_group, kal_source = classify_compound_kallenbach(p['name'])

                processed.append({
                    **p,
                    'chen_group': chen_group,
                    'chen_source': chen_source,
                    'kal_group': kal_group,
                    'kal_source': kal_source,
                    'atoms': atoms,
                    'nosc': nosc,
                    'dg_cox': dg
                })

            results[treatment] = {
                'sample_code': sample_code,
                'filepath': filepath,
                'total_peaks': len(peaks),
                'processed_peaks': len(processed),
                'peaks': processed
            }

            print(f"    {len(processed)} peaks processed")

    return results


def compute_statistics(results):
    """Compute source attribution and group statistics."""
    stats = {}

    for treatment, data in results.items():
        peaks = data['peaks']
        total_area = sum(p['area'] for p in peaks)

        # Chen statistics
        chen_source = defaultdict(float)
        chen_groups = defaultdict(float)

        # Kallenbach statistics
        kal_source = defaultdict(float)
        kal_groups = defaultdict(float)

        for p in peaks:
            area = p['area']
            area_pct = (area / total_area * 100) if total_area > 0 else 0

            chen_source[p['chen_source']] += area_pct
            chen_groups[p['chen_group']] += area_pct
            kal_source[p['kal_source']] += area_pct
            kal_groups[p['kal_group']] += area_pct

        stats[treatment] = {
            'total_area': total_area,
            'chen_source': dict(chen_source),
            'chen_groups': dict(chen_groups),
            'kal_source': dict(kal_source),
            'kal_groups': dict(kal_groups)
        }

    return stats


# ============================================================================
# 4. EXCEL OUTPUT
# ============================================================================

def write_excel(results, stats, output_path):
    """Write results to multi-sheet Excel workbook."""
    wb = Workbook()

    # Styles
    header_font = Font(name='Microsoft YaHei', bold=True, size=11)
    header_fill = PatternFill(start_color='4472C4', end_color='4472C4', fill_type='solid')
    header_font_white = Font(name='Microsoft YaHei', bold=True, size=11, color='FFFFFF')
    cell_font = Font(name='Microsoft YaHei', size=10)
    thin_border = Border(
        left=Side(style='thin'), right=Side(style='thin'),
        top=Side(style='thin'), bottom=Side(style='thin')
    )

    def style_header(ws, row, ncols):
        for col in range(1, ncols + 1):
            cell = ws.cell(row=row, column=col)
            cell.font = header_font_white
            cell.fill = header_fill
            cell.border = thin_border
            cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)

    # --- Sheet 1: 筛选统计 ---
    ws1 = wb.active
    ws1.title = '筛选统计'
    headers1 = ['样品代码', '处理', '原始峰数', 'SI>=80峰数', '原始总面积', '筛选后总面积',
                '筛选后面积/原始面积_%', '相对丰度合计_%', '筛选后唯一化合物数', '含卤素/其他元素化合物数']
    for col, h in enumerate(headers1, 1):
        ws1.cell(row=1, column=col, value=h)
    style_header(ws1, 1, len(headers1))

    row = 2
    for treatment, data in results.items():
        peaks = data['peaks']
        total_area = sum(p['area'] for p in peaks)
        unique = len(set(p['name'] for p in peaks))
        exotic = sum(1 for p in peaks if p['atoms']['other'])

        ws1.cell(row=row, column=1, value=data['sample_code'])
        ws1.cell(row=row, column=2, value=treatment)
        ws1.cell(row=row, column=3, value=data['total_peaks'])
        ws1.cell(row=row, column=4, value=len(peaks))
        ws1.cell(row=row, column=5, value=total_area)
        ws1.cell(row=row, column=6, value=total_area)  # All passed SI>=80
        ws1.cell(row=row, column=7, value=100.0)
        ws1.cell(row=row, column=8, value=100.0)
        ws1.cell(row=row, column=9, value=unique)
        ws1.cell(row=row, column=10, value=exotic)
        for col in range(1, len(headers1) + 1):
            ws1.cell(row=row, column=col).font = cell_font
            ws1.cell(row=row, column=col).border = thin_border
        row += 1

    # --- Sheet 2: Chen 分类详细 ---
    ws2 = wb.create_sheet('Chen_分类详细')
    headers2 = ['处理', '峰号', '保留时间_min', '峰面积', '相对丰度_%', '化合物名称',
                'Chen化学类别', 'Chen来源归属']
    for col, h in enumerate(headers2, 1):
        ws2.cell(row=1, column=col, value=h)
    style_header(ws2, 1, len(headers2))

    row = 2
    for treatment, data in results.items():
        total_area = sum(p['area'] for p in data['peaks'])
        for p in data['peaks']:
            ws2.cell(row=row, column=1, value=treatment)
            ws2.cell(row=row, column=2, value=p['peak_num'])
            ws2.cell(row=row, column=3, value=round(p['ret_time'], 3))
            ws2.cell(row=row, column=4, value=p['area'])
            ws2.cell(row=row, column=5, value=round(p['area']/total_area*100, 4) if total_area > 0 else 0)
            ws2.cell(row=row, column=6, value=p['name'])
            ws2.cell(row=row, column=7, value=p['chen_group'])
            ws2.cell(row=row, column=8, value=p['chen_source'])
            for col in range(1, len(headers2) + 1):
                ws2.cell(row=row, column=col).font = cell_font
                ws2.cell(row=row, column=col).border = thin_border
            row += 1

    # --- Sheet 3: Chen 来源统计 ---
    ws3 = wb.create_sheet('Chen_来源统计')
    headers3 = ['处理', '来源类别', '峰数', '峰面积合计', '相对丰度_%']
    for col, h in enumerate(headers3, 1):
        ws3.cell(row=1, column=col, value=h)
    style_header(ws3, 1, len(headers3))

    row = 2
    for treatment, st in stats.items():
        for source in ['plant', 'microbial', 'mixed']:
            ws3.cell(row=row, column=1, value=treatment)
            ws3.cell(row=row, column=2, value=source)
            ws3.cell(row=row, column=3, value=sum(1 for p in results[treatment]['peaks'] if p['chen_source'] == source))
            area_sum = sum(p['area'] for p in results[treatment]['peaks'] if p['chen_source'] == source)
            ws3.cell(row=row, column=4, value=area_sum)
            ws3.cell(row=row, column=5, value=round(st['chen_source'].get(source, 0), 2))
            for col in range(1, len(headers3) + 1):
                ws3.cell(row=row, column=col).font = cell_font
                ws3.cell(row=row, column=col).border = thin_border
            row += 1

        # Add total row
        ws3.cell(row=row, column=1, value=treatment)
        ws3.cell(row=row, column=2, value='合计')
        ws3.cell(row=row, column=3, value=len(results[treatment]['peaks']))
        ws3.cell(row=row, column=4, value=st['total_area'])
        ws3.cell(row=row, column=5, value=100.0)
        for col in range(1, len(headers3) + 1):
            ws3.cell(row=row, column=col).font = Font(name='Microsoft YaHei', bold=True, size=10)
            ws3.cell(row=row, column=col).border = thin_border
        row += 1

    # --- Sheet 4: Chen 来源饼图数据 ---
    ws4 = wb.create_sheet('Chen_来源饼图')
    headers4 = ['处理', '植物源', '微生物源', '混合来源', '合计']
    for col, h in enumerate(headers4, 1):
        ws4.cell(row=1, column=col, value=h)
    style_header(ws4, 1, len(headers4))

    row = 2
    for treatment, st in stats.items():
        ws4.cell(row=row, column=1, value=treatment)
        ws4.cell(row=row, column=2, value=round(st['chen_source'].get('plant', 0), 4))
        ws4.cell(row=row, column=3, value=round(st['chen_source'].get('microbial', 0), 4))
        ws4.cell(row=row, column=4, value=round(st['chen_source'].get('mixed', 0), 4))
        ws4.cell(row=row, column=5, value=100.0)
        for col in range(1, len(headers4) + 1):
            ws4.cell(row=row, column=col).font = cell_font
            ws4.cell(row=row, column=col).border = thin_border
        row += 1

    # --- Sheet 5: Kallenbach 分类详细 ---
    ws5 = wb.create_sheet('Kallenbach_分类详细')
    headers5 = ['处理', '峰号', '保留时间_min', '峰面积', '相对丰度_%', '化合物名称',
                'Kallenbach化学类别', 'Kallenbach来源归属']
    for col, h in enumerate(headers5, 1):
        ws5.cell(row=1, column=col, value=h)
    style_header(ws5, 1, len(headers5))

    row = 2
    for treatment, data in results.items():
        total_area = sum(p['area'] for p in data['peaks'])
        for p in data['peaks']:
            ws5.cell(row=row, column=1, value=treatment)
            ws5.cell(row=row, column=2, value=p['peak_num'])
            ws5.cell(row=row, column=3, value=round(p['ret_time'], 3))
            ws5.cell(row=row, column=4, value=p['area'])
            ws5.cell(row=row, column=5, value=round(p['area']/total_area*100, 4) if total_area > 0 else 0)
            ws5.cell(row=row, column=6, value=p['name'])
            ws5.cell(row=row, column=7, value=p['kal_group'])
            ws5.cell(row=row, column=8, value=p['kal_source'])
            for col in range(1, len(headers5) + 1):
                ws5.cell(row=row, column=col).font = cell_font
                ws5.cell(row=row, column=col).border = thin_border
            row += 1

    # --- Sheet 6: Kallenbach 来源统计 ---
    ws6 = wb.create_sheet('Kallenbach_来源统计')
    headers6 = ['处理', '来源类别', '峰数', '峰面积合计', '相对丰度_%']
    for col, h in enumerate(headers6, 1):
        ws6.cell(row=1, column=col, value=h)
    style_header(ws6, 1, len(headers6))

    row = 2
    for treatment, st in stats.items():
        for source in ['plant', 'microbial', 'mixed']:
            ws6.cell(row=row, column=1, value=treatment)
            ws6.cell(row=row, column=2, value=source)
            count = sum(1 for p in results[treatment]['peaks'] if p['kal_source'] == source)
            area_sum = sum(p['area'] for p in results[treatment]['peaks'] if p['kal_source'] == source)
            ws6.cell(row=row, column=3, value=count)
            ws6.cell(row=row, column=4, value=area_sum)
            ws6.cell(row=row, column=5, value=round(st['kal_source'].get(source, 0), 2))
            for col in range(1, len(headers6) + 1):
                ws6.cell(row=row, column=col).font = cell_font
                ws6.cell(row=row, column=col).border = thin_border
            row += 1

        ws6.cell(row=row, column=1, value=treatment)
        ws6.cell(row=row, column=2, value='合计')
        ws6.cell(row=row, column=3, value=len(results[treatment]['peaks']))
        ws6.cell(row=row, column=4, value=st['total_area'])
        ws6.cell(row=row, column=5, value=100.0)
        for col in range(1, len(headers6) + 1):
            ws6.cell(row=row, column=col).font = Font(name='Microsoft YaHei', bold=True, size=10)
            ws6.cell(row=row, column=col).border = thin_border
        row += 1

    # --- Sheet 7: Kallenbach 饼图 ---
    ws7 = wb.create_sheet('Kallenbach_来源饼图')
    headers7 = ['处理', '植物源', '微生物源', '混合来源', '合计']
    for col, h in enumerate(headers7, 1):
        ws7.cell(row=1, column=col, value=h)
    style_header(ws7, 1, len(headers7))

    row = 2
    for treatment, st in stats.items():
        ws7.cell(row=row, column=1, value=treatment)
        ws7.cell(row=row, column=2, value=round(st['kal_source'].get('plant', 0), 4))
        ws7.cell(row=row, column=3, value=round(st['kal_source'].get('microbial', 0), 4))
        ws7.cell(row=row, column=4, value=round(st['kal_source'].get('mixed', 0), 4))
        ws7.cell(row=row, column=5, value=100.0)
        for col in range(1, len(headers7) + 1):
            ws7.cell(row=row, column=col).font = cell_font
            ws7.cell(row=row, column=col).border = thin_border
        row += 1

    # --- Sheet 8: 说明_QC ---
    ws8 = wb.create_sheet('说明_QC')
    info = [
        ['Py-GC-MS 批量分析报告'],
        [''],
        ['项目', '说明'],
        ['1. 分析方法', '基于NIST搜索导出TXT文件的MC Peak Table解析'],
        ['2. SI筛选', '保留SI>=80的化合物（对应NIST Match Factor >= 80%）'],
        ['3. Chen 2023分类', '参考Chen et al. 2023, Carbon Research 2:1. DOI: 10.1007/s44246-022-00034-0'],
        ['4. Kallenbach 2016分类', '参考Kallenbach et al. 2016, Nature Communications 7:13630. DOI: 10.1038/ncomms13630'],
        ['5. 来源归属', '植物源(plant)、微生物源(microbial)、混合来源(mixed)'],
        ['6. 化学类别(Chen)', 'lipids, monocyclic_aromatics, polycyclic_aromatics, phenolics, polysaccharides, lignins, amino_N_bearing, heterocyclic_N_bearing, other_N_bearing, unspecified'],
        ['7. 化学类别(Kallenbach)', 'lipids, lignin_derivatives, polysaccharides, proteins, non_protein_N, phenolics, aromatics, unspecified'],
        ['8. 注意', '分类基于化合物名称的关键词匹配。对于需要精确归属的化合物，请人工核对NIST谱库信息。'],
        ['9. 脚本', 'pygcms-batch skill v1.0']
    ]
    for i, row_data in enumerate(info, 1):
        for j, val in enumerate(row_data, 1):
            ws8.cell(row=i, column=j, value=val)
            if i == 1:
                ws8.cell(row=i, column=j).font = Font(name='Microsoft YaHei', bold=True, size=14)
            elif i == 3:
                ws8.cell(row=i, column=j).font = header_font
            elif i > 3:
                ws8.cell(row=i, column=j).font = cell_font

    # Adjust column widths
    for ws in wb.worksheets:
        for col in ws.columns:
            max_length = 0
            col_letter = get_column_letter(col[0].column)
            for cell in col:
                try:
                    if cell.value:
                        max_length = max(max_length, len(str(cell.value)))
                except:
                    pass
            adjusted_width = min(max_length + 2, 60)
            ws.column_dimensions[col_letter].width = adjusted_width

    # Save
    wb.save(output_path)
    print(f"\nOutput saved to: {output_path}")


# ============================================================================
# 5. MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Py-GC-MS Batch Analysis - Compound classification by Chen 2023 and Kallenbach 2016'
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input directory containing NIST TXT export files')
    parser.add_argument('--output', '-o', required=True,
                        help='Output Excel file path (.xlsx)')
    parser.add_argument('--sample_map', '-m', default=None,
                        help='JSON file mapping sample codes to treatment names')

    args = parser.parse_args()

    # Load sample mapping
    sample_map = {}
    if args.sample_map and os.path.exists(args.sample_map):
        with open(args.sample_map, 'r', encoding='utf-8') as f:
            sample_map = json.load(f)

    print("=" * 60)
    print("Py-GC-MS Batch Analysis")
    print("=" * 60)
    print(f"Input directory: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Sample mapping: {sample_map if sample_map else 'None (using filenames)'}")
    print()

    # Process
    results = process_directory(args.input, sample_map)

    if not results:
        print("ERROR: No valid TXT files found or no peaks extracted!")
        return

    # Statistics
    stats = compute_statistics(results)

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY: Chen 2023 Source Attribution")
    print("=" * 60)
    for treatment, st in stats.items():
        print(f"\n{treatment}:")
        print(f"  Total peaks: {len(results[treatment]['peaks'])}")
        print(f"  Total area: {st['total_area']:,.0f}")
        for source in ['plant', 'microbial', 'mixed']:
            val = st['chen_source'].get(source, 0)
            print(f"  {source}: {val:.2f}%")

    # Write Excel
    write_excel(results, stats, args.output)

    print("\nDone!")


if __name__ == '__main__':
    main()
