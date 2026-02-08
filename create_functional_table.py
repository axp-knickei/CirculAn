import pandas as pd
import sys

def create_table(input_file):
    try:
        # eggNOG-mapper files often have a lot of header lines starting with #
        # We read with comment='#' to skip them.
        df = pd.read_csv(input_file, sep='	', comment='#', header=None)
        
        # Read the header line to get column names if possible
        # Typically the last line starting with # is the header
        col_names = []
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('#query'):
                    col_names = line.lstrip('#').strip().split('	')
                    break
        
        if col_names:
            df.columns = col_names[:len(df.columns)]
        
        # Identify column
        target_col = None
        if 'eggNOG_functional_category' in df.columns:
            target_col = 'eggNOG_functional_category'
        elif 'COG_category' in df.columns:
            target_col = 'COG_category'
        else:
            # Fallback to column index if names are missing
            # In emapper v2, COG category is usually column 21 (0-indexed)
            if len(df.columns) > 21:
                target_col = df.columns[21]
            else:
                raise ValueError("Could not find functional category column.")

        # Clean and count
        # Some rows have multiple categories (e.g., 'S/V' or 'SV')
        all_cats = []
        for val in df[target_col].dropna():
            val = str(val).strip()
            if val and val != '-':
                # Handle concatenated categories
                for char in val:
                    if char.isalpha():
                        all_cats.append(char)
        
        counts = pd.Series(all_cats).value_counts().reset_index()
        counts.columns = ['Category', 'Count']
        
        total = counts['Count'].sum()
        counts['Percentage'] = (counts['Count'] / total) * 100
        
        counts.to_csv('functional_category_distribution.csv', index=False)
        print("Successfully created functional_category_distribution.csv")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process eggNOG annotations.')
    parser.add_argument('--input', default='data/Galaxy3-[eggNOG_dataset 1_ anno].tabular', help='Input tabular file')
    parser.add_argument('--output', default='functional_category_distribution.csv', help='Output CSV file')
    args = parser.parse_args()
    
    # We need to pass the output filename to the create_table function
    # I will slightly modify the function call to accept it.
    def create_table_with_output(input_file, output_file):
        try:
            df = pd.read_csv(input_file, sep='\t', comment='#', header=None)
            col_names = []
            with open(input_file, 'r') as f:
                for line in f:
                    if line.startswith('#query'):
                        col_names = line.lstrip('#').strip().split('\t')
                        break
            if col_names:
                df.columns = col_names[:len(df.columns)]
            
            target_col = None
            if 'eggNOG_functional_category' in df.columns:
                target_col = 'eggNOG_functional_category'
            elif 'COG_category' in df.columns:
                target_col = 'COG_category'
            else:
                if len(df.columns) > 21:
                    target_col = df.columns[21]
                else:
                    raise ValueError("Could not find functional category column.")

            all_cats = []
            for val in df[target_col].dropna():
                val = str(val).strip()
                if val and val != '-':
                    for char in val:
                        if char.isalpha():
                            all_cats.append(char)
            
            counts = pd.Series(all_cats).value_counts().reset_index()
            counts.columns = ['Category', 'Count']
            total = counts['Count'].sum()
            counts['Percentage'] = (counts['Count'] / total) * 100
            
            counts.to_csv(output_file, index=False)
            print(f"Successfully created {output_file}")
        except Exception as e:
            print(f"Error: {e}")
            import sys
            sys.exit(1)

    create_table_with_output(args.input, args.output)
