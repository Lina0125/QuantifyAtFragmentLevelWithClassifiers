#!/usr/bin/env python3

"""
This script excute q value calculation for snakemake pipeline
"""

from main import score
import argparse
import sys

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description = print(''),usage='%(prog)s -i -c -s -o [-h for help]')

parser.add_argument('-i','--tsv_file',type=str,help='path to input.tsv',metavar='')
parser.add_argument('-c','--classifier',type=str,help='path to classifier',metavar='')
parser.add_argument('-s','--scaler',type=str,help='path to standard scaler',metavar='')
parser.add_argument('-o','--output',type=str,help='path to output.tsv',metavar='')

if len(sys.argv[1:]) == 8:
    args = parser.parse_args()
    score(fp=args.tsv_file, Classifier=args.classifier, Scaler=args.scaler, out_path=args.output)
else:
    parser.print_help()
    sys.exit('Please pass enough parameters')

