import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='Generate reference matrix from each reference compound')
parser.add_argument('--input-dir', '-i', required=True, help='Input directory path containing files')
parser.add_argument('--ref-list', '-m', required=True, help='File containing compound list, one per line')
parser.add_argument('--output', '-o', default='ref.txt', help='Output file path')
parser.add_argument('--nrows', '-n', type=int, default=1181, help='Number of rows in matrix')

args = parser.parse_args()

with open(args.ref_list, 'r') as f:
    compounds = [line.strip() for line in f if line.strip()]
print(f"Read {len(compounds)} compounds")

ref = pd.DataFrame(index=range(args.nrows), columns=compounds)

for compound in compounds:
    file_path = os.path.join(args.input_dir, f'{compound}.txt')
    auc = pd.read_csv(file_path, header=None, sep='\t')
    
    if len(auc) < args.nrows:
        values = auc.iloc[:, 1].values
        values = np.pad(values, (0, args.nrows - len(values)), constant_values=np.nan)
        ref[compound] = values
    else:
        ref[compound] = auc.iloc[:args.nrows, 1].values
    
    print(f"Read {compound} data")

ref.to_csv(args.output, sep='\t', index=False, header=False)
print(f"Processing completed. Output saved to {args.output} with dimensions {ref.shape}")
