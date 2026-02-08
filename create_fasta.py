import pandas as pd
import re
import argparse

def create_protein_fasta(input_xlsx, output_faa):
    try:
        df = pd.read_excel(input_xlsx)
    except Exception as e:
        print(f"Error loading {input_xlsx}: {e}")
        return

    seq_col = "Prot seq"
    if seq_col not in df.columns:
        print(f"Error: Column '{seq_col}' not found in the file.")
        return

    id_col = "Prokka_ID"
    has_id_col = id_col in df.columns
    aa_regex = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY\*X]+$')

    count = 0
    with open(output_faa, 'w') as f:
        for index, row in df.iterrows():
            raw_seq = row[seq_col]
            if pd.isna(raw_seq):
                continue
                
            seq = str(raw_seq).strip().upper()
            if not seq or seq == 'NAN' or not aa_regex.match(seq):
                continue

            if has_id_col and pd.notna(row[id_col]) and str(row[id_col]).strip() != "":
                gene_id = str(row[id_col]).strip()
            else:
                gene_id = f"gene{index}"

            f.write(f">{gene_id}\n{seq}\n")
            count += 1

    print(f"Success: Wrote {count} protein sequences to {output_faa}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract proteins from XLSX to FASTA.')
    parser.add_argument('--input', required=True, help='Input XLSX file (e.g. data/E1_func_anno.xlsx)')
    parser.add_argument('--output', required=True, help='Output FASTA filename (e.g. E1_proteins.faa)')
    args = parser.parse_args()
    
    create_protein_fasta(args.input, args.output)
