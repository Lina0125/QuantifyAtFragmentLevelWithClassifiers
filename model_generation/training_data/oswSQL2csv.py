# This script take a openswath osw file as input, then output a csv file that well merged according to columns in different tables

#!/usr/bin/env python3
import sqlite3
import pandas as pd

def ows_tables_2_csv(osw_file, csv_file):
    osw_tables = sqlite3.connect(osw_file, uri=True)
    db_df = pd.read_sql_query("""
    SELECT *
        FROM FEATURE_TRANSITION
        INNER JOIN(
        SELECT
            ID FEATURE_ID,
            PRECURSOR_ID PRECUSOR_ID,
            EXP_RT RT_APEX,
            LEFT_WIDTH RT_START,
            RIGHT_WIDTH RT_END
        from FEATURE
    ) FEATURE ON FEATURE_TRANSITION.FEATURE_ID =  FEATURE.FEATURE_ID
    INNER JOIN (
        SELECT
            ID,
            CHARGE,
            PRECURSOR_MZ MZ,
            DECOY
        FROM PRECURSOR
    ) PRECURSOR ON FEATURE.PRECUSOR_ID = PRECURSOR.ID
    ORDER BY PRECURSOR.ID ASC,
         FEATURE.RT_APEX ASC""", osw_tables)
    db_df.to_csv(csv_file, index=False)

import argparse
import sys
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description = print(''),usage='%(prog)s -i -o [-h for help]')

parser.add_argument('-i','--input',type=str,help='path to input.osw',metavar='')
parser.add_argument('-o','--output',type=str,help='path to output.csv',metavar='')

if len(sys.argv[1:]) >= 1:
    args = parser.parse_args()
else:
    parser.print_help()
    sys.exit('\nPlease pass enough parameters')
    
ows_tables_2_csv(args.input, args.output)
print(args.input.split('.')[0]+' is merged!')
