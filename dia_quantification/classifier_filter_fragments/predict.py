#!/usr/bin/env python3

def score(fp, out_path, Classifier, Scaler):
    '''
    -fp merged tsv file from pyprophet scored osw tables
    -out_path output tsv with scored qvalue
    '''
    from preprocessing import data_cleaning
    from metrics import calculate_probability
    from gscore.fdr import ScoreDistribution
    import os, joblib
    import pandas as pd

    os.path.expanduser(os.getcwd())
    
    # Loading data
    print('Info: Loading classifier', Classifier, sep='\t', end='\n')
    classifier = joblib.load(Classifier)
    scaler = joblib.load(Scaler)

    # Testing data preprocessing
    print('Info: Preparing data for classifier', end='\n')
    # Here we use fp instead of df, because data_cleaning also inspect data, it'll report if sth wrong with this file
    feature, target = data_cleaning(fp)
    feature_standard = scaler.transform(feature)

    # Predicting and get decision scores
    print('Info: Predicting and calculating decision scores', end='\n')
    y_scores = calculate_probability(classifier, feature_standard, target)

    # Calculate q values
    print('Info: Calculting q values', end='\n')
    score_distribution = ScoreDistribution()
    score_distribution.fit(
        y_scores.ravel(), # Decision function scores from the classifier
        target.ravel() # TARGET labels
        )

    q_values = score_distribution.calculate_q_values(y_scores)

    # Reread and merge qvalue to file
    osw_df = pd.read_csv(fp, sep='\t')
    osw_df['FRAGMENT_QVALUE'] = q_values

    # Output scored tsv
    osw_df.to_csv(out_path, sep='\t', index=False)


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
    sys.exit('\nPlease pass enough parameters')

