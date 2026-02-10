#!/usr/bin/env python3
"""
genome_features_pipeline.py

Compute assembly metrics and annotation-derived candidate lists for one or more strains.

What it does:
- For each strain prefix provided (e.g. "E1" or "D4"):
  * Reads <PREFIX>_contig.fasta for assembly-level metrics:
      - genome size (bp), G+C content (%), contig N50 (bp), contig L50, number of contigs
  * Reads <PREFIX>_func_anno.xlsx (first sheet) for Prokka-style annotation:
      - discovers Prokka_ID / Feature / Product-like columns
      - counts unique loci with CDS, rRNA, tRNA
      - builds per-locus annotation text
  * Reads any CSV files matching <PREFIX>_Geno_Anno_*.csv and merges annotation text by Prokka_ID
  * Searches per-locus concatenated text for configurable keyword lists (adhesion, bacteriocin, virulence/defense)
      - reports unique Prokka_ID counts per category (deduplicated)
  * Saves outputs:
      - <PREFIX>_candidate_genes.csv (Prokka_ID, Feature, Product, matched categories, annotation text)
      - outputs combined summary for all strains to summary_table.csv
      - writes LaTeX table (table_E1_D4.tex) and Methods text (methods.txt)

Requirements:
- Python 3.8+
- pandas
- openpyxl (for reading .xlsx)
No network access required.

Usage example:
    python genome_features_pipeline.py --strains E1 D4 --indir /path/to/data --outdir ./results


Author: Alex
Date: 2026-02-10
"""

import argparse
import csv
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# -------------------------
# Utilities: fasta parsing
# -------------------------
def parse_fasta(fasta_path: Path) -> Tuple[List[str], List[int], int, int]:
    """
    Parse a FASTA and return sequence names list, lengths list, total GC count, total base count.
    Returns (names, lengths, gc_count, total_bases)
    """
    names = []
    lengths = []
    gc_count = 0
    total_bases = 0
    seq_parts = []
    header = None
    with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    seq = ''.join(seq_parts).upper()
                    l = len(seq)
                    names.append(header)
                    lengths.append(l)
                    gc_count += seq.count('G') + seq.count('C')
                    total_bases += l
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        # last
        if header is not None:
            seq = ''.join(seq_parts).upper()
            l = len(seq)
            names.append(header)
            lengths.append(l)
            gc_count += seq.count('G') + seq.count('C')
            total_bases += l
    return names, lengths, gc_count, total_bases


def assembly_metrics_from_lengths(lengths: List[int]) -> Tuple[int, Optional[float], Optional[int], Optional[int], int]:
    """
    Given a list of contig lengths, return:
    (genome_size, gc_percent placeholder None, N50, L50, n_contigs)
    Note: GC percent computed elsewhere when sequence GC available.
    """
    if not lengths:
        return 0, None, None, None, 0
    genome_size = sum(lengths)
    sorted_lengths = sorted(lengths, reverse=True)
    half = genome_size / 2.0
    cumsum = 0
    n50 = None
    l50 = None
    for i, L in enumerate(sorted_lengths):
        cumsum += L
        if cumsum >= half:
            n50 = L
            l50 = i + 1
            break
    return genome_size, None, n50, l50, len(lengths)


# -------------------------
# Annotation parsing utils
# -------------------------
def find_column_like(columns: List[str], candidates: List[str]) -> Optional[str]:
    """
    Find the first column in 'columns' that matches any candidate string (case-insensitive substring).
    Returns exact column name or None.
    """
    cols_lower = {c.lower(): c for c in columns}
    for cand in candidates:
        cand_low = cand.lower()
        # exact match first
        if cand_low in cols_lower:
            return cols_lower[cand_low]
    # substring match
    for c in columns:
        for cand in candidates:
            if cand.lower() in c.lower():
                return c
    return None


def load_func_anno_xlsx(path: Path) -> pd.DataFrame:
    """
    Load the first sheet of a functional annotation xlsx (Prokka-style).
    """
    # using pandas.read_excel (requires openpyxl)
    df = pd.read_excel(path, sheet_name=0)
    return df


def build_text_by_prokka(df_func: pd.DataFrame, additional_csvs: List[Path], prokka_col_candidates=None):
    """
    Return:
     - prokka_meta: dict[prokka_id] = {'Feature': ..., 'Product': ...}
     - text_by_prokka: dict[prokka_id] = concatenated text (lowercased)
    """
    if prokka_col_candidates is None:
        prokka_col_candidates = ['Prokka_ID', 'prokka_id', 'locus_tag', 'locus tag', 'ID', 'id']

    cols = list(df_func.columns)
    prokka_col = find_column_like(cols, prokka_col_candidates)
    feature_col = find_column_like(cols, ['Feature', 'feature', 'type', 'seqtype'])
    product_col = find_column_like(cols, ['Product', 'product', 'annotation', 'desc', 'description', 'note'])
    # Build metadata and text from df_func
    prokka_meta = dict()  # prokka_id -> dict(Feature=..., Product=...)
    text_by_prokka = defaultdict(list)
    if prokka_col:
        for _, row in df_func.iterrows():
            pid = str(row.get(prokka_col, '')).strip()
            if not pid or pid.lower() == 'nan':
                continue
            feat = str(row.get(feature_col, '')).strip() if feature_col else ''
            prod = str(row.get(product_col, '')).strip() if product_col else ''
            prokka_meta.setdefault(pid, {})['Feature'] = feat
            prokka_meta.setdefault(pid, {})['Product'] = prod
            # gather textual fields to help keyword searches
            texts = []
            for c in df_func.columns:
                if pd.api.types.is_string_dtype(df_func[c]) or 'product' in c.lower() or 'annot' in c.lower() or 'desc' in c.lower():
                    try:
                        val = row.get(c, '')
                        if pd.isna(val):
                            continue
                        texts.append(str(val))
                    except Exception:
                        continue
            if texts:
                text_by_prokka[pid].append(" || ".join(texts))
    else:
        # No Prokka_ID column; fallback to creating pseudo-IDs based on index (less reliable)
        for idx, row in df_func.iterrows():
            pid = f"row_{idx}"
            feat = str(row.get(feature_col, '')).strip() if feature_col else ''
            prod = str(row.get(product_col, '')).strip() if product_col else ''
            prokka_meta.setdefault(pid, {})['Feature'] = feat
            prokka_meta.setdefault(pid, {})['Product'] = prod
            texts = []
            for c in df_func.columns:
                if pd.api.types.is_string_dtype(df_func[c]) or 'product' in c.lower() or 'annot' in c.lower() or 'desc' in c.lower():
                    try:
                        val = row.get(c, '')
                        if pd.isna(val):
                            continue
                        texts.append(str(val))
                    except Exception:
                        continue
            if texts:
                text_by_prokka[pid].append(" || ".join(texts))

    # Merge text from additional CSVs where Prokka_ID-like column exists
    for csv_path in additional_csvs:
        try:
            df_csv = pd.read_csv(csv_path, low_memory=False)
        except Exception:
            # try with engine python
            df_csv = pd.read_csv(csv_path, engine='python')
        cols_csv = list(df_csv.columns)
        # try to detect prokka id column
        pcol = find_column_like(cols_csv, ['Prokka_ID', 'prokka_id', 'locus_tag', 'locus tag', 'ID', 'id', 'query id', 'query_id'])
        # description columns
        desc_cols = [c for c in cols_csv if any(x in c.lower() for x in ['product', 'desc', 'description', 'annotation', 'subject', 'function'])]
        if pcol:
            for _, row in df_csv.iterrows():
                pid = str(row.get(pcol, '')).strip()
                if not pid or pid.lower() == 'nan':
                    continue
                parts = []
                for c in desc_cols:
                    try:
                        v = row.get(c, '')
                        if pd.isna(v):
                            continue
                        parts.append(str(v))
                    except Exception:
                        continue
                if parts:
                    text_by_prokka[pid].append(" || ".join(parts))
                    # also conserve some metadata if missing
                    if pid not in prokka_meta:
                        # attempt to supply product/feature from CSV
                        prod = parts[0] if parts else ''
                        prokka_meta.setdefault(pid, {})['Feature'] = prokka_meta.get(pid, {}).get('Feature', '')
                        prokka_meta.setdefault(pid, {})['Product'] = prokka_meta.get(pid, {}).get('Product', prod)
        else:
            # If no prokka id column, skip merging (can't deduplicate)
            continue

    # Normalize text_by_prokka into strings
    text_by_prokka_str = {pid: " || ".join(texts).lower() for pid, texts in text_by_prokka.items()}
    return prokka_meta, text_by_prokka_str


def count_feature_types(prokka_meta: Dict[str, Dict[str, str]]) -> Tuple[int, int, int]:
    """
    Count unique CDS, rRNA, tRNA loci based on the prokka_meta Feature/Product fields.
    Returns (cds_count, rrna_count, trna_count)
    """
    cds_set = set()
    rrna_set = set()
    trna_set = set()
    for pid, meta in prokka_meta.items():
        feat = str(meta.get('Feature', '')).lower()
        prod = str(meta.get('Product', '')).lower()
        # heuristics for CDS
        if 'cds' in feat or feat.strip() == 'cds' or 'protein' in prod or 'hypothetical protein' in prod or 'coding' in feat:
            cds_set.add(pid)
        if 'rrna' in feat or 'rrna' in prod:
            rrna_set.add(pid)
        if 'trna' in feat or 'trna' in prod:
            trna_set.add(pid)
    return len(cds_set), len(rrna_set), len(trna_set)


# -------------------------
# Keyword search
# -------------------------
DEFAULT_KEYWORDS = {
    'adhesion': ['adhes', 'pili', 'fimbr', 'adhesin', 'adhesion'],
    'bacteriocin': ['bacteriocin', 'lantibiotic', 'lantipeptide', 'nisin', 'microcin', 'bovicin'],
    'virulence_defense': ['virulence', 'toxin', 'hemolysin', 'virulence factor', 'pathogen', 'resistance', 'efflux', 'capsule', 'disease']
}

def find_keyword_matches(text_by_prokka: Dict[str, str], keywords: Dict[str, List[str]]) -> Dict[str, set]:
    """
    For each category in keywords, return a set of Prokka_IDs matched (deduplicated).
    """
    matches = {cat: set() for cat in keywords}
    for pid, txt in text_by_prokka.items():
        if not txt:
            continue
        for cat, kwlist in keywords.items():
            for kw in kwlist:
                if kw.lower() in txt:
                    matches[cat].add(pid)
                    break
    return matches


# -------------------------
# AMRFinder parsing helper
# -------------------------
def count_amr_hits(amr_csv_path: Path) -> int:
    """
    Count unique AMR hits from an AMRFinder CSV. Try to detect a gene/element column.
    """
    try:
        df_amr = pd.read_csv(amr_csv_path, low_memory=False)
    except Exception:
        df_amr = pd.read_csv(amr_csv_path, engine='python')
    if df_amr.empty:
        return 0
    col = find_column_like(list(df_amr.columns), ['prokka_id', 'Prokka_ID', 'element', 'gene', 'gene_symbol', 'product', 'model'])
    if col:
        return df_amr[col].astype(str).dropna().unique().size
    else:
        return len(df_amr)


# -------------------------
# Output writers
# -------------------------
def write_candidates_csv(out_path: Path, candidate_rows: List[Dict]):
    df = pd.DataFrame(candidate_rows)
    df.to_csv(out_path, index=False)


def write_summary_csv(out_path: Path, summary_rows: List[Dict]):
    df = pd.DataFrame(summary_rows)
    df.to_csv(out_path, index=False)


def write_latex_table(out_path: Path, summary_rows: List[Dict], caption: str = "Genome assembly and annotation summary."):
    """
    summary_rows: list of dicts with keys: 'strain' plus metrics present in order
    We'll build a two-column table (feature x strains).
    """
    if not summary_rows:
        return
    # Determine order of keys for rows (stable)
    example = summary_rows[0]
    keys_order = [
        'Genome size (bp)', 'G+C content (%)', 'Contig N50 (bp)', 'Contig L50',
        'Number of contigs', 'Number of protein-coding sequences (CDSs)',
        'Number of rRNAs', 'Number of tRNAs',
        'Genes related to adhesion (candidate loci)',
        'Genes related to bacteriocin (candidate loci)',
        'Genes related to virulence, disease and defense (candidate loci)',
        'AMR hits (AMRFinder)'
    ]
    # Map summary_rows into dict: strain -> metrics
    strain_map = {r['strain']: r for r in summary_rows}
    strains = list(strain_map.keys())
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write("% Generated by genome_features_pipeline.py\n")
        fh.write("\\begin{table*}[htbp]\n")
        fh.write("\\centering\n")
        fh.write(f"\\caption{{{caption}}}\n")
        fh.write("\\label{tab:genome-features}\n")
        fh.write("\\small\n")
        # header
        fh.write("\\begin{tabular}{l" + "r" * len(strains) + "}\n")
        fh.write("\\toprule\n")
        fh.write("\\textbf{Feature} ")
        for s in strains:
            fh.write(f" & \\textbf{{{s}}}")
        fh.write(" \\\\\n")
        fh.write("\\midrule\n")
        # rows
        for k in keys_order:
            fh.write(k.replace('_', ' ') )
            for s in strains:
                val = strain_map[s].get(k, "N/A")
                fh.write(f" & {val}")
            fh.write(" \\\\\n")
        fh.write("\\bottomrule\n")
        fh.write("\\end{tabular}\n")
        fh.write("\\end{table*}\n")


def write_methods_txt(out_path: Path):
    methods_text = """Methods (automatically generated)
Assembly statistics (genome size, G+C content, contig N50, contig L50, and contig count) were computed directly from each strain's contig FASTA. Genome size is the sum of contig lengths; G+C content is the proportion of G and C nucleotides in the total assembly. N50 and L50 were obtained from the cumulative distribution of sorted contig lengths (N50 = contig length at which cumulative length >= 50%% of the assembly; L50 = number of contigs to reach that threshold).

Feature counts (CDS, rRNA, tRNA) were extracted from Prokka-style functional annotation spreadsheets. We detect the Prokka_ID/locus identifier and count unique loci with Feature or Product fields indicating CDS, rRNA or tRNA. Candidate functional categories (adhesion, bacteriocin, virulence/defense) were identified by searching per-locus concatenated annotation text across spreadsheets and CSVs for curated keyword lists. Matches are deduplicated by Prokka_ID and reported as candidate loci counts (these are annotation-based candidates and require manual curation/validation).

AMR gene counts are taken from AMRFinder CSV output (unique reported hits).

This script preserves reproducibility by saving per-strain candidate CSVs and the summary table. Use the original Prokka GFF for definitive locus-type counts when possible, and specialized tools (BAGEL/antiSMASH, VFDB) to validate bacteriocin/virulence candidates.
"""
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(methods_text)


# -------------------------
# Main
# -------------------------
def process_strain(prefix: str, indir: Path, outdir: Path, keywords=None) -> Dict:
    """
    Process one strain with prefix (e.g. 'E1' or 'D4').
    Returns a summary dict with metrics and 'strain' key.
    Also writes per-strain candidate CSV to outdir.
    """
    if keywords is None:
        keywords = DEFAULT_KEYWORDS

    prefix = str(prefix)
    contig_fa = indir / f"{prefix}_contig.fasta"
    func_xlsx = indir / f"{prefix}_func_anno.xlsx"
    # gather CSVs
    csv_glob = list(indir.glob(f"{prefix}_Geno_Anno_*.csv"))
    amr_csv = indir / f"{prefix}_Geno_Anno_AMRFinder.csv"
    # check files
    if not contig_fa.exists():
        raise FileNotFoundError(f"Contig FASTA not found for {prefix}: {contig_fa}")
    # parse fasta
    names, lengths, gc_count, total_bases = parse_fasta(contig_fa)
    genome_size = sum(lengths)
    gc_pct = round((gc_count / total_bases * 100) if total_bases > 0 else 0, 2)
    # compute N50/L50
    _, _, n50, l50, n_contigs = assembly_metrics_from_lengths(lengths)
    # parse annotation
    if not func_xlsx.exists():
        # proceed but warning
        df_func = pd.DataFrame()
    else:
        df_func = load_func_anno_xlsx(func_xlsx)
    prokka_meta, text_by_prokka = build_text_by_prokka(df_func, csv_glob)
    cds_c, rrna_c, trna_c = count_feature_types(prokka_meta)
    matches = find_keyword_matches(text_by_prokka, keywords)
    # count matches per category (deduplicated)
    adhesion_count = len(matches.get('adhesion', set()))
    bacteriocin_count = len(matches.get('bacteriocin', set()))
    virulence_count = len(matches.get('virulence_defense', set()))
    # AMR hits
    amr_count = 0
    if amr_csv.exists():
        amr_count = count_amr_hits(amr_csv)
    # produce candidate list (rows)
    candidate_rows = []
    for cat, pidset in matches.items():
        for pid in pidset:
            meta = prokka_meta.get(pid, {})
            candidate_rows.append({
                'Prokka_ID': pid,
                'Feature': meta.get('Feature', ''),
                'Product': meta.get('Product', ''),
                'Matched_category': cat,
                'Annotation_text': text_by_prokka.get(pid, '')
            })
    # deduplicate candidate rows by Prokka_ID, merge matched categories
    df_cand = pd.DataFrame(candidate_rows)
    if not df_cand.empty:
        df_cand = df_cand.groupby('Prokka_ID').agg({
            'Feature': 'first',
            'Product': 'first',
            'Matched_category': lambda s: ";".join(sorted(set(s))),
            'Annotation_text': 'first'
        }).reset_index()
    # write candidate CSV
    outdir.mkdir(parents=True, exist_ok=True)
    cand_path = outdir / f"{prefix}_candidate_genes.csv"
    df_cand.to_csv(cand_path, index=False)
    # summary dict
    summary = {
        'strain': prefix,
        'Genome size (bp)': int(genome_size),
        'G+C content (%)': float(gc_pct),
        'Contig N50 (bp)': int(n50) if n50 is not None else 'N/A',
        'Contig L50': int(l50) if l50 is not None else 'N/A',
        'Number of contigs': int(n_contigs),
        'Number of protein-coding sequences (CDSs)': int(cds_c),
        'Number of rRNAs': int(rrna_c),
        'Number of tRNAs': int(trna_c),
        'Genes related to adhesion (candidate loci)': int(adhesion_count),
        'Genes related to bacteriocin (candidate loci)': int(bacteriocin_count),
        'Genes related to virulence, disease and defense (candidate loci)': int(virulence_count),
        'AMR hits (AMRFinder)': int(amr_count)
    }
    return summary


def main():
    parser = argparse.ArgumentParser(description="Genome assembly + annotation summary pipeline")
    parser.add_argument('--strains', nargs='+', required=True,
                        help="strain prefixes to process (e.g. E1 D4)")
    parser.add_argument('--indir', type=str, default='.',
                        help="input directory where <PREFIX>_contig.fasta and annotation files live")
    parser.add_argument('--outdir', type=str, default='./results',
                        help="output directory for candidate CSVs, summary and LaTeX")
    parser.add_argument('--latex-caption', type=str, default="Summary of assembly and annotation features.",
                        help="caption to use in generated LaTeX table")
    parser.add_argument('--keywords-file', type=str, default=None,
                        help="optional path to a JSON file with keyword dict to override defaults")
    args = parser.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    # optionally read keywords from JSON
    keywords = DEFAULT_KEYWORDS
    if args.keywords_file:
        import json
        with open(args.keywords_file, 'r', encoding='utf-8') as fh:
            keywords = json.load(fh)

    summaries = []
    for s in args.strains:
        print(f"Processing strain {s} ...")
        try:
            summary = process_strain(s, indir, outdir, keywords=keywords)
            summaries.append(summary)
            print(f"  - done: {s}")
        except Exception as e:
            print(f"  - ERROR processing {s}: {e}")

    # write combined summary and latex + methods
    summary_csv = outdir / "summary_table.csv"
    write_summary_csv(summary_csv, summaries)
    latex_path = outdir / "table_genome_features.tex"
    write_latex_table(latex_path, summaries, caption=args.latex_caption)
    methods_path = outdir / "methods.txt"
    write_methods_txt(methods_path)
    print(f"\nOutputs written to {outdir.resolve()}:")
    print(f" - summary CSV: {summary_csv}")
    print(f" - LaTeX table: {latex_path}")
    print(f" - Methods text: {methods_path}")
    print(" - Per-strain candidate CSVs: <PREFIX>_candidate_genes.csv")


if __name__ == "__main__":
    main()
