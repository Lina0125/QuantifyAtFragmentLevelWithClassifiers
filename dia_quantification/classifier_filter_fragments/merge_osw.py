#!/usr/bin/env python3

"""
This script excute merge_osw_tables for snakemake pipeline
"""

from preprocessing import merge_osw_tables
import argparse, sys

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description = print(''),usage='%(prog)s -i -o [-h for help]')

parser.add_argument('-i','--input',type=str,required=True,help='path to .osw file',metavar='')
parser.add_argument('-o','--output',type=str,required=True, help='path to .tsv file',metavar='')

if len(sys.argv[1:]) == 4:
    args = parser.parse_args()
    merge_osw_tables(args.input, args.output)
else:
    parser.print_help()
    sys.exit('\nPlease pass enough parameters')