#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import matplotlib.colors as mcolors

def generate_nature_style_plot(csv_path, output_prefix="figure_eggNOG_nature_style"):
    # ===============================
    #   Load & Normalize Input Data
    # ===============================
    df = pd.read_csv(csv_path)

    # Normalize column names
    df.columns = [c.strip().lower() for c in df.columns]

    # Flexible column detection
    col_category = next((c for c in df.columns if c.startswith("category")), "category")
    col_count = next((c for c in df.columns if c in ("count","counts","n")), "count")
    col_percent = next((c for c in df.columns if c in ("percentage","percent","pct")), None)

    # Compute percent if missing
    if col_percent is None:
        df["percent"] = df[col_count] / df[col_count].sum() * 100
        col_percent = "percent"

    # Sort categories by size
    df = df.sort_values(col_count, ascending=False).reset_index(drop=True)

    # ============================================================
    #   Optional: Map COG letters → full descriptions if present
    # ============================================================
    COG_DESCRIPTIONS = {
        'J': 'Translation, ribosomal structure and biogenesis',
        'A': 'RNA processing and modification',
        'K': 'Transcription',
        'L': 'Replication, recombination and repair',
        'B': 'Chromatin structure and dynamics',
        'D': 'Cell cycle control, cell division, chromosome partitioning',
        'Y': 'Nuclear structure',
        'V': 'Defense mechanisms',
        'T': 'Signal transduction mechanisms',
        'M': 'Cell wall/membrane/envelope biogenesis',
        'N': 'Cell motility',
        'Z': 'Cytoskeleton',
        'W': 'Extracellular structures',
        'U': 'Intracellular trafficking, secretion, and vesicular transport',
        'O': 'Posttranslational modification, protein turnover, chaperones',
        'C': 'Energy production and conversion',
        'G': 'Carbohydrate transport and metabolism',
        'E': 'Amino acid transport and metabolism',
        'F': 'Nucleotide transport and metabolism',
        'H': 'Coenzyme transport and metabolism',
        'I': 'Lipid transport and metabolism',
        'P': 'Inorganic ion transport and metabolism',
        'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
        'R': 'General function prediction only',
        'S': 'Function unknown'
    }

    df["description"] = df[col_category].map(COG_DESCRIPTIONS).fillna(df[col_category])

    # ========================================
    #   Group tiny categories into "Other"
    # ========================================
    threshold = 1.0  # percent
    small = df[df[col_percent] < threshold]

    if len(df) > 15 and not small.empty:
        other_count = small[col_count].sum()
        other_percent = small[col_percent].sum()
        df = df[df[col_percent] >= threshold].copy()

        df = pd.concat([
            df,
            pd.DataFrame({
                col_category: ["Other"],
                col_count: [other_count],
                col_percent: [other_percent],
                "description": ["Other functional categories"]
            })
        ], ignore_index=True)

    # ============================================
    #   Colorblind-friendly palette (tab20)
    # ============================================
    cmap = plt.get_cmap("tab20")
    df["color"] = [cmap(i % cmap.N) for i in range(len(df))]

    # ==========================
    #   Plot Styling Defaults
    # ==========================
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "svg.fonttype": "none",     # keep text selectable in Illustrator
        "pdf.fonttype": 42,
        "ps.fonttype": 42
    })

    fig, ax = plt.subplots(figsize=(14, 10), dpi=600)

    wedgeprops = {
        "width": 0.42,
        "edgecolor": "white",
        "linewidth": 1.2,
        "antialiased": True
    }

    # ========================
    #   Draw Donut Pie Chart
    # ========================
    wedges, texts, autotexts = ax.pie(
        df[col_count],
        colors=df["color"],
        startangle=90,
        autopct=lambda p: f"{p:.1f}%" if p > 2.5 else "",
        pctdistance=0.75,
        wedgeprops=wedgeprops,
        textprops={"fontsize": 11}
    )

    # Adaptive label color (white on dark, black on light)
    def luminance(color):
        r, g, b = mcolors.to_rgb(color)
        return (0.2126*r + 0.7152*g + 0.0722*b)

    for t, col in zip(autotexts, df["color"]):
        t.set_color("white" if luminance(col) < 0.45 else "#222222")
        t.set_fontweight("bold")

    # ========================
    #   Legend Construction
    # ========================
    handles = [Patch(facecolor=c, edgecolor='none') for c in df["color"]]
    labels = [
        f"[{row[col_category]}] {row['description']} ({row[col_percent]:.1f}%)"
        for _, row in df.iterrows()
    ]

    ncols = 1 if len(df) <= 12 else 2

    ax.legend(
        handles,
        labels,
        title="Functional Categories",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=10,
        title_fontsize=12,
        ncol=ncols,
        labelspacing=0.8
    )

    # ===========================
    #   Center text of donut
    # ===========================
    total = df[col_count].sum()
    center_text = f"Total\n{int(total):,}\ngenes"

    ax.text(
        0, 0, center_text,
        ha='center', va='center',
        fontsize=13, fontweight="bold",
        color="#222222"
    )

    # ===================
    #   Title
    # ===================
    plt.title(
        "EggNOG Functional Annotation Distribution",
        fontsize=18, fontweight="bold", pad=25
    )

    plt.tight_layout()

    # ============
    #   Save
    # ============
    plt.savefig(f"{output_prefix}.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"{output_prefix}.svg", bbox_inches="tight")

    print(f"✓ Generated: {output_prefix}.png and .svg")


if __name__ == "__main__":
    import argparse
    import os
    parser = argparse.ArgumentParser(description='Generate Nature-style eggNOG plot.')
    parser.add_argument('--input', default='functional_category_distribution.csv', help='Input CSV file')
    parser.add_argument('--output_prefix', default='figure_eggNOG_nature_style', help='Prefix for output figures')
    args = parser.parse_args()

    if os.path.exists(args.input):
        generate_nature_style_plot(args.input, args.output_prefix)
    else:
        print(f"Error: {args.input} not found.")