#!/usr/bin/env python3

"""
This script takes a directory contains files with fragment ions q values as input, 
Process each run in the following steps:
1) Align matrix back to global model to control FDR at proteins and peptides level;
2) Filter unreliable fragment ions and calculate MS1 intensity using different combination of remaining fragment ions;
3) Merge every run into a matrix on precursors;
4) Align runs to calculate RT and output quantitative matrix.
"""

def filter_workflow(in_dir, template, out_path, method='top3', max_transition_qvalue=0.01, max_fragment_qvalue=0.01):

    from main import merge_runs, align_runs
    # Conduct 1)~3)
    merged = merge_runs(in_dir=in_dir, template=template, method=method,max_transition_qvalue=max_transition_qvalue, max_fragment_qvalue=max_fragment_qvalue)
    # Conduct 4)
    align_runs(merged=merged, out_path=out_path)

import argparse, sys

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description = print(''),usage='%(prog)s -i -t -o[-h for help]')

parser.add_argument('-i','--in_dir',type=str,required=True,help='directory to *.tsv',metavar='')
parser.add_argument('-t','--template',type=str,required=True,help='model_global.osw',metavar='')
parser.add_argument('-o','--output',type=str,required=True, help='path to quantitative matrix',metavar='')
parser.add_argument('-m','--method',type=str,required=False, default='top3', help='method to choose optimal fragment(topN, mean, median, sum)',metavar='')
parser.add_argument('-q1','--transition_qvalue',type=float,required=False,default=0.01,help='max transition qvalue',metavar='')
parser.add_argument('-q2','--fragment_qvalue',type=float,required=False,default=0.01,help='max fragment qvalue',metavar='')

if len(sys.argv[1:]) >= 6:
    args = parser.parse_args()
    filter_workflow(
        in_dir = args.in_dir, 
        template =  args.template,
        out_path = args.output, 
        method = args.method, 
        max_transition_qvalue=args.transition_qvalue, 
        max_fragment_qvalue=args.fragment_qvalue
        )

else:
    parser.print_help()
    sys.exit('At least pass -i, -t, -o')