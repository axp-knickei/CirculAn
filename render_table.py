import matplotlib.pyplot as plt
import os

# Create a figure and axis
fig, ax = plt.subplots(figsize=(8, 6))

# Define the table content based on E1_D4.tex
table_data = [
    ["Feature", "E1", "D4"],
    ["Genome size (bp)", "3,114,003", "3,156,065"],
    ["G+C content (%)", "46.19", "46.15"],
    ["Contig N50 (bp)", "86,629", "78,029"],
    ["Contig L50", "11", "12"],
    ["Number of contigs", "166", "195"],
    ["Number of protein-coding sequences (CDSs)", "2,934", "2,972"],
    ["Number of rRNAs", "12", "12"],
    ["Number of tRNAs", "111", "111"],
    ["Genes related to adhesion (candidate loci)^a", "18", "18"],
    ["Bacteriocin-related candidate loci^a", "13", "14"],
    ["Virulence / disease / defense candidate loci^a", "52", "54"],
    ["AMR hits (AMRFinder)^b", "1", "1"]
]

# Add the table to the axis
table = ax.table(cellText=table_data, loc='center', cellLoc='left')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.2)

# Remove axes
ax.axis('off')

# Add caption/notes below the table
notes = ("Notes. Numbers are produced from strain-specific files (E1: E1_contig.fasta, E1_func_anno.xlsx, E1_Geno_Anno_*.csv; "
         "D4: D4_contig.fasta, D4_func_anno.xlsx, D4_Geno_Anno_*.csv). Candidate loci counts (a) were identified by keyword search "
         "of the annotation text and are reported as unique Prokka_ID matches (deduplicated across annotation files). "
         "AMR hits (b) were taken from the AMRFinder CSV output.")
plt.figtext(0.1, 0.05, notes, wrap=True, horizontalalignment='left', fontsize=8)

# Save as PNG and PDF
plt.savefig("E1_D4.png", dpi=300, bbox_inches='tight')
plt.savefig("E1_D4.pdf", bbox_inches='tight')

print("Created E1_D4.png and E1_D4.pdf")
