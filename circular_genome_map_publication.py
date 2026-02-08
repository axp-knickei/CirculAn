#!/usr/bin/env python3
"""
Whole-genome circular map - PUBLICATION VERSION
Enhanced styling (Colors, Fonts, High-DPI) while preserving 
the working geometry of the original script.
"""

import os
import argparse
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# --- Utilities ---

def read_fasta_to_dict(path):
    seqs = {}
    with open(path, 'r') as fh:
        name = None
        parts = []
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(parts).upper()
                name = line[1:].strip()
                parts = []
            else:
                parts.append(line)
        if name:
            seqs[name] = ''.join(parts).upper()
    return seqs

def file_md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()

def normalize_name(s):
    return str(s).lower().replace('|','_').replace('-','_').replace(' ','_')

def build_concatenated_sequence(fasta_dict, spacer=100, order=None):
    if order is None:
        order = list(fasta_dict.keys())
    offsets = []
    parts = []
    offset = 0
    for h in order:
        seq = fasta_dict[h]
        start = offset
        parts.append(seq)
        offset += len(seq)
        end = offset - 1
        offsets.append((h, len(seq), start, end))
        if spacer > 0:
            parts.append('N' * spacer)
            offset += spacer
    return ''.join(parts), offsets

def write_offsets_tsv(offsets, path):
    with open(path,'w') as fh:
        fh.write('scaffold_name\tscaffold_length\toffset_start\toffset_end\n')
        for h,l,s,e in offsets:
            fh.write(f"{h}\t{l}\t{s}\t{e}\n")

def remap_annotations(df, offsets):
    norm_to_idx = {normalize_name(h): i for i,(h,_,_,_) in enumerate(offsets)}
    rows = []
    for _, r in df.iterrows():
        unitig = r.get('Unitig','')
        if pd.isna(unitig) or str(unitig).strip()=='':
            continue
        norm = normalize_name(unitig)
        idx = None
        if norm in norm_to_idx:
            idx = norm_to_idx[norm]
        else:
            for k,v in norm_to_idx.items():
                if k.find(norm) != -1 or norm.find(k) != -1:
                    idx = v
                    break
        if idx is None:
            continue
        h,l,sc_start,sc_end = offsets[idx]
        try:
            s = int(r['Start']); e = int(r['End'])
        except Exception:
            continue
        gstart = sc_start + (s-1)
        gend = sc_start + (e-1)
        rows.append({'scaffold_idx': idx, 'feature': r.get('Feature',''), 'strand': str(r.get('Strand','')).strip(), 'gstart': gstart, 'gend': gend})
    return pd.DataFrame(rows)

def extract_cds(remap):
    if remap.empty:
        return [], []
    fvals = remap['feature'].astype(str).str.upper()
    mask = fvals.str.contains('CDS')
    if mask.sum() == 0:
        mask = fvals.str.contains('GENE|CODING')
    cds = remap[mask]
    fwd = []
    rev = []
    for _, r in cds.iterrows():
        s = int(r['gstart']); e = int(r['gend']); st = r['strand']
        if st in ('+','1','+1'):
            fwd.append((s,e))
        elif st in ('-','-1'):
            rev.append((s,e))
        else:
            fwd.append((s,e))
    return fwd, rev

def merge_intervals(intervals, gap=20):
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    cur_s, cur_e = intervals[0]
    for s,e in intervals[1:]:
        if s <= cur_e + gap:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged

def sliding_gc_and_skew(seq, window=1000, step=None):
    if step is None:
        step = max(1, window//4)
    L = len(seq)
    pos = []
    gc = []
    skew = []
    i = 0
    while i < L:
        wseq = seq[i:i+window]
        if len(wseq) == 0:
            break
        g = wseq.count('G')
        c = wseq.count('C')
        gc_val = (g + c) / len(wseq)
        skew_val = (g - c) / (g + c) if (g + c) > 0 else 0.0
        pos.append(i + len(wseq)//2)
        gc.append(gc_val)
        skew.append(skew_val)
        i += step
    return np.array(pos), np.array(gc), np.array(skew)

def N50_calc(lengths):
    x = sorted(lengths, reverse=True)
    tot = sum(x)
    cum = 0
    for L in x:
        cum += L
        if cum >= tot/2:
            return L
    return None

# --- Plotting ---

def plot_map(seq, fwd, rev, pos_gc, gc_vals, skew_vals, out_png, out_svg, window_size=1000):
    L = len(seq)
    
    # 1. High-Tier Styling (Fonts)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
    
    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_axis_off()

    # 2. Publication-Grade Color Palette
    COLOR_FWD = '#377eb8'  # Solid Blue
    COLOR_REV = '#e41a1c'  # Solid Red
    COLOR_GC  = '#252525'  # Near Black
    COLOR_SKEW_POS = '#4daf4a' # Soft Green
    COLOR_SKEW_NEG = '#984ea3' # Soft Purple
    COLOR_SKEW_LINE = '#542788' # Dark Purple Outline

    def pos2ang(p):
        return 2*np.pi*(p/L)

    # 3. PRESERVED GEOMETRY (Working Radii)
    outer = 1.0
    h = 0.06
    inner = outer - (h + 0.03)
    gc_r   = inner - 0.09
    skew_r = gc_r - 0.18

    # 4. CDS Tracks
    fwd_m = merge_intervals(fwd)
    rev_m = merge_intervals(rev)
    min_ang = 2*np.pi/L*3

    for s,e in fwd_m:
        ang = pos2ang(s)
        w = pos2ang(max(1, e-s))
        if w < min_ang: w = min_ang
        ax.bar(ang, h, width=w, bottom=outer, align='edge', color=COLOR_FWD, linewidth=0)

    for s,e in rev_m:
        ang = pos2ang(s)
        w = pos2ang(max(1, e-s))
        if w < min_ang: w = min_ang
        ax.bar(ang, h, width=w, bottom=inner, align='edge', color=COLOR_REV, linewidth=0)

    # 5. GC and Skew Tracks
    ang_gc = pos2ang(pos_gc)
    gc_mean = gc_vals.mean()
    gc_line = gc_r + (gc_vals - gc_mean) * 0.35
    skew_line = skew_r + (skew_vals) * 0.9

    # GC Content Filled Area
    ax.fill_between(ang_gc, gc_line, gc_r, where=(gc_line >= gc_r), color=COLOR_GC, alpha=0.5)
    ax.fill_between(ang_gc, gc_line, gc_r, where=(gc_line < gc_r), color=COLOR_GC, alpha=0.5)
    ax.plot(ang_gc, gc_line, color=COLOR_GC, linewidth=0.8)

    # GC Skew Filled Area
    ax.fill_between(ang_gc, skew_line, skew_r, where=(skew_vals >= 0), color=COLOR_SKEW_POS, alpha=0.6)
    ax.fill_between(ang_gc, skew_line, skew_r, where=(skew_vals < 0), color=COLOR_SKEW_NEG, alpha=0.6)
    ax.plot(ang_gc, skew_line, color=COLOR_SKEW_LINE, linewidth=0.8)

    # 6. Professional Legend (Preserved location, refined text)
    legend_elems = [
        Patch(facecolor=COLOR_FWD, label='Forward CDS'),
        Patch(facecolor=COLOR_REV, label='Reverse CDS'),
        Patch(facecolor=COLOR_GC, alpha=0.5, label=f'GC Content (Mean: {gc_mean*100:.1f}%)'),
        Patch(facecolor=COLOR_SKEW_POS, alpha=0.6, label='GC Skew Positive (G > C)'),
        Patch(facecolor=COLOR_SKEW_NEG, alpha=0.6, label='GC Skew Negative (C > G)')
    ]
    ax.legend(
        handles=legend_elems, 
        loc='lower center', 
        bbox_to_anchor=(0.5, -0.15), 
        fontsize=10, 
        frameon=False,
        title=f"Genomic Features (Window: {window_size} bp)",
        title_fontproperties={'weight': 'bold', 'size': 11}
    )

    # 7. Title
    ax.set_title(f"Whole-Genome Circular Map\nLength: {L:,} bp", 
                 y=1.08, fontsize=15, fontweight='bold')

    plt.savefig(out_png, dpi=600, bbox_inches='tight')
    plt.savefig(out_svg, bbox_inches='tight')
    plt.close(fig)

# --- Main ---

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--fasta', required=True)
    p.add_argument('--anno', required=True)
    p.add_argument('--window', type=int, default=1000)
    p.add_argument('--step', type=int, default=None)
    p.add_argument('--spacer', type=int, default=100)
    p.add_argument('--outdir', default='.')
    p.add_argument('--order_by', choices=['fasta','size','annotation'], default='fasta')
    return p.parse_args()

def main():
    args = parse_args()
    print('Reading FASTA...')
    seqs = read_fasta_to_dict(args.fasta)
    
    headers = list(seqs.keys())
    if args.order_by == 'size':
        headers = sorted(headers, key=lambda h: len(seqs[h]), reverse=True)
    
    whole_seq, offsets = build_concatenated_sequence(seqs, spacer=args.spacer, order=headers)
    os.makedirs(args.outdir, exist_ok=True)
    
    write_offsets_tsv(offsets, os.path.join(args.outdir, 'scaffold_offsets_pub.tsv'))
    
    print('Processing annotations...')
    ann = pd.read_excel(args.anno)
    remapped = remap_annotations(ann, offsets)
    fwd, rev = extract_cds(remapped)
    
    print(f'Computing GC & skew (window {args.window})...')
    pos_gc, gc_vals, skew_vals = sliding_gc_and_skew(whole_seq, window=args.window, step=args.step)
    
    out_png = os.path.join(args.outdir, 'circular_map_publication.png')
    out_svg = os.path.join(args.outdir, 'circular_map_publication.svg')

    plot_map(whole_seq, fwd, rev, pos_gc, gc_vals, skew_vals, out_png, out_svg, window_size=args.window)
    print(f'Success! Figures saved to {args.outdir}')

if __name__ == '__main__':
    main()