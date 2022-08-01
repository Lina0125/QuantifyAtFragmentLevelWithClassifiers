#!/usr/bin/env python3

"""
INTRO:
    1)merge_osw_tables is used to merge osw tables to a dataframe
    2)merge_global_tables is used to merge pyProphet global osw tables to a dataframe
    3)data_cleaning is used to generate feature and label nparray that needed to be feeded to the classifier
    3)calculate_gscore is used to calculate q values
"""

def merge_osw_tables(osw_file, out_path=None):

    import pandas as pd
    import sqlite3

    print("Info: Loading", osw_file, sep=' ')

    osw_tables = sqlite3.connect(osw_file, uri=True)
    print("Info: Merging ows tables")
    db_df = pd.read_sql_query("""
    SELECT FEATURE.FEATURE_ID, FEATURE_TRANSITION.TRANSITION_ID,PRECURSOR.PRECURSOR_ID, 
    PEPTIDE.PEPTIDE_ID, PROTEIN.PROTEIN_ID,AREA_INTENSITY, TOTAL_AREA_INTENSITY, APEX_INTENSITY, 
    TOTAL_MI,VAR_INTENSITY_SCORE, VAR_INTENSITY_RATIO_SCORE, VAR_LOG_INTENSITY, VAR_XCORR_COELUTION,
    VAR_XCORR_SHAPE, VAR_LOG_SN_SCORE, VAR_MASSDEV_SCORE, VAR_MI_SCORE, VAR_MI_RATIO_SCORE,
    VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE, RT_APEX, RT_START, RT_END,
    PEAKGROUP_SCORE,PEAKGROUP_QVALUE, PEAKGROUP_RANK, TRANSITION_ANNOTATION, PRECURSOR_CHARGE,
    PRECURSOR_MZ,DECOY,UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, PROTEIN_ACCESSION, 
    PEPTIDE_QVALUE, PROTEIN_QVALUE

    FROM FEATURE_TRANSITION
        INNER JOIN(
        SELECT
            ID FEATURE_ID,
            PRECURSOR_ID,
            EXP_RT RT_APEX,
            LEFT_WIDTH RT_START,
            RIGHT_WIDTH RT_END
        from FEATURE
    ) FEATURE ON FEATURE_TRANSITION.FEATURE_ID =  FEATURE.FEATURE_ID
        INNER JOIN (
        SELECT
            FEATURE_ID FEATURE_ID,
            SCORE PEAKGROUP_SCORE,
            QVALUE PEAKGROUP_QVALUE,
            RANK PEAKGROUP_RANK
        FROM SCORE_MS2
    ) SCORE_MS2 ON FEATURE_TRANSITION.FEATURE_ID = SCORE_MS2.FEATURE_ID
        INNER JOIN (
        SELECT
            ID TRANSITION_ID,
            ANNOTATION TRANSITION_ANNOTATION
        FROM TRANSITION
    ) TRANSITION ON FEATURE_TRANSITION.TRANSITION_ID = TRANSITION.TRANSITION_ID
        INNER JOIN (
        SELECT
            ID PRECURSOR_ID,
            CHARGE PRECURSOR_CHARGE,
            PRECURSOR_MZ PRECURSOR_MZ,
            DECOY
        FROM PRECURSOR
    ) PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.PRECURSOR_ID
        INNER JOIN (
        SELECT
            PRECURSOR_ID PRECURSOR_ID,
            PEPTIDE_ID
        FROM PRECURSOR_PEPTIDE_MAPPING
    ) PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.PRECURSOR_ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID
        INNER JOIN (
        SELECT
            PEPTIDE_ID PEPTIDE_ID,
            PROTEIN_ID
        FROM PEPTIDE_PROTEIN_MAPPING
    ) PEPTIDE_PROTEIN_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID
        INNER JOIN (
        SELECT
            PEPTIDE_ID PEPTIDE_ID,
            QVALUE PEPTIDE_QVALUE
        FROM SCORE_PEPTIDE
    ) SCORE_PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID
        INNER JOIN (
        SELECT
            ID PEPTIDE_ID,
            UNMODIFIED_SEQUENCE,
            MODIFIED_SEQUENCE
        FROM PEPTIDE
    ) PEPTIDE ON PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID = PEPTIDE.PEPTIDE_ID
        INNER JOIN (
        SELECT
            ID PROTEIN_ID,
            PROTEIN_ACCESSION
        FROM PROTEIN
    ) PROTEIN ON PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID = PROTEIN.PROTEIN_ID
        INNER JOIN (
                SELECT
                    PROTEIN_ID,
                    QVALUE PROTEIN_QVALUE
                FROM SCORE_PROTEIN
            ) SCORE_PROTEIN ON SCORE_PROTEIN.PROTEIN_ID = PROTEIN.PROTEIN_ID
    ORDER BY PRECURSOR.PRECURSOR_ID ASC,
         FEATURE.RT_APEX ASC""", osw_tables)

    # TODO 0504 add this to avoid value error
    #print('Info: NaN is droped to ensure qscore calculation do not raise error')
    #db_df = db_df.dropna()
    #db_df = db_df.reset_index().drop(columns='index')
    
    if out_path:
        db_df.to_csv(out_path, sep='\t', index=False)
    else:
        return db_df

def merge_global_tables(global_osw_file, out_path=None):

    import pandas as pd
    import sqlite3

    print("Info: Loading", global_osw_file, "to get global protein list", sep=' ')

    global_osw_tables = sqlite3.connect(global_osw_file, uri=True)
    db_df = pd.read_sql_query("""
    SELECT ID PROTEIN_ID, PROTEIN_ACCESSION, PROTEIN_QVALUE,
       SCORE_PEPTIDE.PEPTIDE_ID, MODIFIED_SEQUENCE,
       PEPTIDE_QVALUE, PEPTIDE.DECOY, PRECURSOR_ID
    FROM PROTEIN
        INNER JOIN(
        SELECT
            PROTEIN_ID,
            QVALUE PROTEIN_QVALUE
        from SCORE_PROTEIN
    ) SCORE_PROTEIN ON SCORE_PROTEIN.PROTEIN_ID =  PROTEIN.ID
        INNER JOIN (
        SELECT
            PEPTIDE_ID,
            PROTEIN_ID
        FROM PEPTIDE_PROTEIN_MAPPING
    ) PEPTIDE_PROTEIN_MAPPING ON PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID = PROTEIN.ID
        INNER JOIN (
        SELECT
            PEPTIDE_ID,
            QVALUE PEPTIDE_QVALUE
        FROM SCORE_PEPTIDE
    ) SCORE_PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID
        INNER JOIN (
        SELECT
            ID PEPTIDE_ID,
            MODIFIED_SEQUENCE,
            DECOY
        FROM PEPTIDE
    ) PEPTIDE ON PEPTIDE.PEPTIDE_ID = SCORE_PEPTIDE.PEPTIDE_ID
        INNER JOIN (
        SELECT
            PRECURSOR_ID,
            PEPTIDE_ID
        FROM PRECURSOR_PEPTIDE_MAPPING
    ) PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = SCORE_PEPTIDE.PEPTIDE_ID
    ORDER BY PROTEIN.ID ASC,
         PEPTIDE.PEPTIDE_ID ASC""", global_osw_tables)
    
    if out_path:
        db_df.to_csv(out_path, sep='\t', index=False)
    else:
        return db_df


def data_cleaning(fp):

    import numpy as np,pandas as pd

    print(f"Info: Loading and preprocessing {fp}")
    df = pd.read_csv(fp, sep='\t')

    # Use TARGET instead of DECOY 
    df["TARGET"] = np.abs(df["DECOY"] - 1)
    # Calculate RT_LENGTH
    df.loc[:, 'RT_LENGTH'] = df.loc[:, 'RT_END'] - df.loc[:, 'RT_START']

    # Only keep sub socres, RT_LENGTH and TARGET
    all_cols = [x for x in df.columns.tolist() if x[:4] == 'VAR_']
    all_cols.append('RT_LENGTH')
    all_cols.append('TARGET')
    df = df[all_cols]
    
    # Drop columns that all 0
    all_0_cols = df.loc[:, (df == 0).all(axis=0)].columns.tolist()

    if len(all_0_cols) != 0:
        print(f"Warning: {','.join(list(all_0_cols))} columns are all 0 in {fp}")
        #df = df.loc[:, (df != 0).any(axis=0)]

    feature = df[all_cols[:-1]]
    feature = feature.values
    
    target = df[['TARGET']]
    target = target.values
    target = target[:, 0]
    
    if len(all_0_cols) != 0:
        print('Warning: send data to the classifier, but the score may not be so accuracy because of the all 0 column feature')
    
    return feature, target

def calculate_gscore(feature, target, classifier, scaler):
    from metrics import calculate_probability
    from gscore.fdr import ScoreDistribution
    
    feature_standard = scaler.transform(feature)
    # Scoring
    print('Info: Predicting and calculating decision scores', end='\n')
    y_scores = calculate_probability(classifier, feature_standard, target)

    print('Info: Calculting q values', end='\n')
    score_distribution = ScoreDistribution()
    score_distribution.fit(
        y_scores.ravel(), # Decision function scores from the classifier
        target.ravel() # TARGET labels
        )

    q_values = score_distribution.calculate_q_values(y_scores)

    return q_values