# Circular Genome Mapping & Functional Annotation Pipeline

A professional suite of Python tools for bacterial genome visualization and functional annotation analysis, optimized for high-tier journal publications (Nature, ISME, etc.).

## 🛠 Features

- **Circular Genome Mapping**: High-resolution maps with GC content/skew fills and publication-grade aesthetics.
- **Nature-Style Visualization**: High-tier donut charts for COG/eggNOG functional distributions.
- **Data Processing**: Robust scripts for protein extraction and EggNOG-mapper results parsing.

## 📁 Repository Structure

- `circular_genome_map_publication.py`: The primary tool for creating circular maps.
- `plot_nature_figure.py`: Generates publication-quality functional distribution figures.
- `generate_scientific_figure.py`: Alternative COG distribution script using Seaborn styling.
- `create_functional_table.py`: Parses eggNOG-mapper tabular outputs into summarized distributions.
- `create_fasta.py`: Extracts validated protein sequences from XLSX annotation files.

## 🚀 Quick Start

Ensure you have [uv](https://github.com/astral-sh/uv) installed for dependency management.

### 1. Extract Proteins for Annotation
```bash
uv run python create_fasta.py --input your_annotation.xlsx --output proteins.faa
```

### 2. Generate Circular Map
```bash
uv run python circular_genome_map_publication.py \
   --fasta your_genome.fasta \
   --anno your_annotation.xlsx \
   --window 9000 \
   --outdir results
```

### 3. Create Functional Distribution Figure
```bash
# First, summarize your EggNOG results
uv run python create_functional_table.py --input emapper_out.annotations --output summary.csv

# Then, generate the Nature-style donut chart
uv run python plot_nature_figure.py --input summary.csv --output_prefix functional_plot
```

## 📜 License
MIT

---
**Note**: Post-analysis for **Genomics NGS Analysis Team (Taiwan © 2017-2018 Genomics)**.
