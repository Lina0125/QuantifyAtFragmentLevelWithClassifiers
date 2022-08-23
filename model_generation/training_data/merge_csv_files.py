"""
INTRO:
    This script take a directory contains all csv files transformed from openswath output as input, then output a merged csv file
USAGE:
    python merge_csv_file.py -d [directory] -o [path to output file]
"""

#!/usr/bin/env python3

import argparse
import sys, os
import pandas as pd

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description = print('''
Merge all csv files.
'''),usage='%(prog)s -d -o [-h for help]')

parser.add_argument('-d', '--directory', type=str, help='Path to directory that contains csv files', metavar='')
parser.add_argument('-o', '--output', type=str, help='Path to merged csv file', metavar='')

# If enough parameters?
if len(sys.argv[1:]) == 4:
    args = parser.parse_args()
else:
    parser.print_help()
    sys.exit()


# 
csv_dir = args.directory
merged_file_path = args.output

# Function that merge csv data to df
def merge_data(file_dir):
    files = os.listdir(file_dir)
    data = pd.read_csv(os.path.join(file_dir, files[0]), sep=',')
    data['FILE_ID'] = [files[0].split('.')[0]] * len(data.index.tolist())
    for file in files[1:]:
        tmp = pd.read_csv(os.path.join(file_dir, file), sep=',')
        tmp['FILE_ID'] = [file.split('.')[0]] * len(tmp.index.tolist())
        data = pd.concat([data,tmp], join='inner')
    return data
# 
merged_df = merge_data(csv_dir)
merged_df.to_csv(merged_file_path, index=False)