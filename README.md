# Circular Genome Mapping & Functional Annotation Pipeline

A comprehensive and user-friendly toolkit designed specifically for researchers working with bacterial genomes. This pipeline creates publication-quality circular genome visualizations and functional annotation analyses suitable for high-impact journals such as *Nature*, *ISME Journal*, *Science*, and *Cell*.

## 🎯 Who This Tool Is For

This toolkit is specifically designed for:
- **Microbiologists and environmental microbiologists** analyzing bacterial genomes
- **Genomics researchers** who need professional visualizations for publications
- **Bioinformatics beginners** who may not have extensive programming experience
- **Research groups** working with bacterial isolates from environmental or clinical samples

## 🌟 Key Features

### 🧬 Circular Genome Visualization
- **High-resolution circular maps** with customizable color schemes
- **GC content and GC skew analysis** with automatic calculation and visualization
- **Gene annotation display** with color-coded functional categories
- **Publication-ready output** in multiple formats (PNG, SVG, PDF)
- **Automatic scaling** optimized for different genome sizes

### 📊 Functional Analysis Tools
- **Nature-style donut charts** for COG/eggNOG functional category distributions
- **Professional color palettes** designed for scientific publications
- **Statistical summaries** of functional annotations
- **Multiple visualization options** for different presentation needs

### 🛠 Data Processing Utilities
- **Protein sequence extraction** from Excel annotation files
- **EggNOG-mapper results parsing** with comprehensive summarization
- **Data validation** and error checking to ensure reliable results
- **Flexible input/output formats** compatible with common bioinformatics tools

---

## 📋 Prerequisites & System Requirements

### What You Need Before Starting

**Basic Requirements:**
- **Computer**: Windows 10+, macOS 10.15+, or Linux (Ubuntu 18.04+)
- **Python**: Version 3.8 or higher
- **Memory**: Minimum 4GB RAM (8GB+ recommended for large genomes)
- **Storage**: At least 1GB free space for installation and output files

**No Programming Experience Required**: All tools can be run using simple copy-and-paste commands.

### Installing the Required Tools

This toolkit uses `uv` for easy dependency management (similar to how app stores work on your phone):

1. **Install uv** (one-time setup):
   ```bash
   # For macOS/Linux
   curl -LsSf https://astral.sh/uv/install.sh | sh
   
   # For Windows
   powershell -c "irm https://astral.sh/uv/install.sh | iex"
   ```

2. **Restart your terminal** after installation.

---

## 🚀 Step-by-Step Guide

This section walks you through the entire process from raw data to publication-ready figures.

### Step 1: Prepare Your Input Files

**You will need:**

1. **Genome Sequence File** (`genome.fasta`)
   - Contains your bacterial genome DNA sequence
   - Must be in FASTA format (standard format used by most sequencing tools)
   - Can be assembled from your NGS data or downloaded from databases

2. **Gene Annotation File** (`annotation.xlsx`)
   - Contains information about genes, their locations, and functions
   - Must be in Excel format (.xlsx)
   - Should include columns for gene name, start position, end position, strand, and functional category

**Example of a properly formatted annotation file:**

| Gene | Start | End | Strand | Function | COG_Category |
|------|-------|-----|--------|----------|--------------|
| gyrA | 100 | 2500 | + | DNA gyrase subunit A | L |
| recA | 2600 | 4000 | - | Recombinase A | L |
| rpoB | 4100 | 7000 | + | RNA polymerase beta | K |

### Step 2: Extract Protein Sequences (Optional but Recommended)

If you plan to perform functional annotation using EggNOG-mapper:

```bash
uv run python create_fasta.py --input your_annotation.xlsx --output proteins.faa
```

**What this does:**
- Reads your Excel annotation file
- Extracts protein sequences for each gene
- Creates a FASTA file compatible with functional annotation tools

### Step 3: Create Circular Genome Map

**Basic Usage:**
```bash
uv run python circular_genome_map_publication.py \
   --fasta your_genome.fasta \
   --anno your_annotation.xlsx \
   --window 9000 \
   --outdir results
```

**Understanding the Parameters:**
- `--fasta`: Your genome sequence file
- `--anno`: Your gene annotation Excel file
- `--window`: Window size for GC content calculation (9000 works well for most bacteria)
- `--outdir`: Folder where results will be saved

**Output Files:**
- `circular_map.png`: High-resolution image for presentations
- `circular_map.svg`: Vector format for publications (infinite resolution)
- `circular_map.pdf`: PDF format for journal submissions
- `gc_content.txt`: Tabular data for further analysis

### Step 4: Functional Annotation Analysis (Optional)

**Part A: Run EggNOG-mapper** (External tool, not included)
- Upload your `proteins.faa` file to the [EggNOG-mapper web server](http://eggnog-mapper.embl.de/)
- Download the results as a tab-separated file

**Part B: Summarize Functional Categories:**
```bash
uv run python create_functional_table.py \
   --input emapper_out.annotations \
   --output functional_summary.csv
```

**Part C: Create Publication-Quality Figures:**

**Option 1: Nature-Style Donut Chart**
```bash
uv run python plot_nature_figure.py \
   --input functional_summary.csv \
   --output_prefix functional_analysis
```

**Option 2: Scientific Bar Chart**
```bash
uv run python generate_scientific_figure.py \
   --input functional_summary.csv \
   --output_prefix functional_bars
```

---

## 📈 Understanding Your Results

### Interpreting Circular Maps

- **Outer Ring**: Gene annotations color-coded by functional category
- **Middle Ring**: GC content (higher values in red, lower in blue)
- **Inner Ring**: GC skew (indicative of replication origin and terminus)
- **Scale Marks**: Genome position in kilobases (kb)

### Functional Category Colors

| Category | Description | Color |
|----------|-------------|-------|
| J | Translation, ribosomal structure | Blue |
| K | Transcription | Green |
| L | Replication, recombination and repair | Red |
| C | Energy production and conversion | Orange |
| E | Amino acid transport and metabolism | Purple |
| F | Nucleotide transport and metabolism | Brown |
| G | Carbohydrate transport and metabolism | Pink |
| H | Coenzyme transport and metabolism | Gray |
| I | Lipid transport and metabolism | Yellow |
| P | Inorganic ion transport and metabolism | Cyan |
| Q | Secondary metabolites biosynthesis | Magenta |
| R | General function prediction only | Light Blue |
| S | Function unknown | Light Gray |

---

## 🛠 Advanced Customization

### Changing Color Schemes

You can modify colors by editing the script files:
- Open the relevant Python script in a text editor
- Look for color definitions (usually hex codes like `#FF0000`)
- Replace with your preferred colors
- Save and re-run the script

### Adjusting Figure Resolution

For extremely high-resolution outputs:
```bash
uv run python circular_genome_map_publication.py \
   --fasta your_genome.fasta \
   --anno your_annotation.xlsx \
   --window 9000 \
   --outdir results \
   --dpi 600
```

### Batch Processing Multiple Genomes

Create a simple text file with your genome paths and use:
```bash
while IFS= read -r genome; do
    uv run python circular_genome_map_publication.py \
       --fasta "$genome" \
       --anno "${genome%.fasta}_annotation.xlsx" \
       --window 9000 \
       --outdir "results_$(basename "$genome" .fasta)"
done < genome_list.txt
```

---

## 🆘 Troubleshooting & FAQ

### Common Issues and Solutions

**Q: "Python command not found" error**
A: Install Python 3.8+ from [python.org](https://www.python.org/downloads/) or use your system's package manager.

**Q: "uv command not found" error**
A: Make sure you installed uv correctly and restarted your terminal. Try running `uv --version` to verify installation.

**Q: Memory errors with large genomes**
A: Increase the `--window` parameter to smaller values (e.g., 5000) or use a computer with more RAM.

**Q: Empty or corrupted output files**
A: Check that your input files are properly formatted and not corrupted. Try opening them in their respective applications first.

**Q: Colors don't match my journal's requirements**
A: Use the SVG output format and import it into vector graphics software (Adobe Illustrator, Inkscape) for final adjustments.

### Getting Help

**For Technical Issues:**
- Check the error messages carefully for specific file or parameter issues
- Verify your input files match the expected format
- Ensure all required software is properly installed

**For Scientific Questions:**
- Consult with your bioinformatics core facility
- Reference the original EggNOG-mapper publications for methodology
- Consider collaborating with computational biologists for complex analyses

---

## 📚 Scientific Background

### What Are Circular Genome Maps?

Circular genome maps are visual representations of bacterial chromosomes that display:
- **Gene organization and density**
- **Functional categorization of genes**
- **Nucleotide composition patterns (GC content/skew)**
- **Evolutionary relationships between genomic regions**

These visualizations help researchers identify:
- **Replication origin and terminus regions**
- **Genomic islands and horizontal gene transfer events**
- **Metabolic capabilities and lifestyle adaptations**
- **Comparative genomic features between strains**

### Functional Annotation Categories

The COG (Clusters of Orthologous Groups) classification system groups genes based on:
- **Evolutionary relationships** across different species
- **Functional similarities** and biochemical pathways
- **Cellular processes** and molecular functions

This standardized system enables meaningful comparisons between different bacterial genomes and facilitates meta-analyses across multiple studies.

---

## 📄 Citation Information

If you use this toolkit in your research, please cite:

```
Circular Genome Mapping & Functional Annotation Pipeline
Genomics NGS Analysis Team, Taiwan © 2017-2018
```

**Additional citations to consider:**
- EggNOG-mapper methodology paper (Huerta-Cepas et al., 2017, 2019)
- Original COG database publications (Tatusov et al., 1997, 2003)
- Specific journals' guidelines for figure preparation

---

## 📜 License

This project is released under the MIT License, allowing for both academic and commercial use with proper attribution.

---

## 👥 About This Project

**Developed by**: Genomics NGS Analysis Team, Taiwan  
**Development Period**: 2017-2018  
**Purpose**: To provide accessible, high-quality genome visualization tools for the microbiology research community  
**Target Users**: Researchers needing publication-ready figures without extensive bioinformatics expertise  

**Contributions and feedback are welcome!** This toolkit continues to evolve based on user needs and emerging visualization standards in scientific publishing.

---

## 🔄 Version History

**Version 1.0** (2017-2018): Initial release with core circular mapping functionality  
**Current Version**: Enhanced with additional visualization options and improved user experience  

---

*Last updated: February 2026*  
*Compatible with Python 3.8+ and modern operating systems*