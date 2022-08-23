#!/usr/bin/env python3

def global_align(df, template):

    """
    This function control FDR on protein and peptide level by aligning scored data 
    back to global model.
    """

    from preprocessing import merge_global_tables
    import pandas as pd
    import os
    
    os.path.expanduser(os.getcwd())

    global_df = merge_global_tables(template)
    global_data = global_df[global_df['DECOY'] == 0]
    global_data = global_data[global_data['PROTEIN_QVALUE'] < 0.01]
    global_data = global_data[global_data['PEPTIDE_QVALUE'] < 0.01]
    global_data['global_index'] = global_data['PROTEIN_ACCESSION'] + global_data['MODIFIED_SEQUENCE']
    global_data = global_data[['global_index']]
    global_data = global_data.drop_duplicates()

    df['global_index'] = df['PROTEIN_ACCESSION'] + df['MODIFIED_SEQUENCE']
    df = pd.merge(global_data, df, on='global_index', how='left')

    return df

def filter_fragments(fp, template, method='mean', max_transition_qvalue=0.01, max_fragment_qvalue=0.01):

    """
    This function align data back to global model to filter out unconfident proteins and peptides;
    Then it filters fragment ions according to its q value, and calculate precursor intensity using optional methods;
    The rest proteins and peptides that don't have any confident fragment ions are assigned NAN on RT and intensity
    """

    import pandas as pd, numpy as np
    import os
    from collections import Counter
    import time
    import warnings

    from pandas.core.common import SettingWithCopyWarning
    warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

    start_cpu = time.process_time()
    os.path.expanduser(os.getcwd())

    print(f'Info: Loading {fp}, quantify precursor using {method} method.')
    df = pd.read_csv(fp, sep='\t')
    sample = os.path.split(fp)[-1].split('.')[0]

    # Align this sample to global model to get the data frame of global protieins and peptides in each run
    print(f'Info: Align {fp} to global model.')
    df = global_align(df, template)
    
    precursor_list = df[['PRECURSOR_ID', 'PROTEIN_ACCESSION', 'MODIFIED_SEQUENCE', 'PRECURSOR_CHARGE']]
    precursor_list = precursor_list.drop_duplicates()

    # Filtering unconfident precusrsors and fragments
    data = df[df['PEAKGROUP_QVALUE'] < max_transition_qvalue]
    data = data[data['FRAGMENT_QVALUE'] < max_fragment_qvalue]

    # Precursor intensity calculation metheds and it's process
    if method[:3] == 'top':
        topn = method
        if topn == 'top1':
            precursor_top1 = data[["PRECURSOR_ID", "APEX_INTENSITY", 'FRAGMENT_QVALUE', 'RT_APEX']]
            precursor_top1 = precursor_top1.sort_values('FRAGMENT_QVALUE', ascending=True).groupby(['PRECURSOR_ID']).first().reset_index()
            filtered = precursor_top1[["PRECURSOR_ID", "APEX_INTENSITY", 'RT_APEX']]
        else:
            # Get precursor that fragment number is lower than TopN
            count_number = {'PRECURSOR_ID': Counter(data['PRECURSOR_ID']).keys(), 
                            'FRAGMENT_COUNT': Counter(data['PRECURSOR_ID']).values()}
            count_number = pd.DataFrame(count_number)
            under = count_number[count_number['FRAGMENT_COUNT'] < int(topn[3:])]['PRECURSOR_ID'].tolist()

            # For these precursor that has more than N confident fragment ions, remove extra, only keep N fragment ions
            sub1 = data[data.PRECURSOR_ID.isin(under)]
            sub2 = data[~data.PRECURSOR_ID.isin(under)]
            sub2 = sub2.sort_values('FRAGMENT_QVALUE', ascending=True).groupby(['PRECURSOR_ID']).head(int(topn[3:]))
            data = pd.concat([sub1, sub2], join='inner')

            precursor_intensity = data.groupby('PRECURSOR_ID', as_index=False)['APEX_INTENSITY'].sum()
            precursor_rt = data.groupby('PRECURSOR_ID', as_index=False)['RT_APEX'].agg(np.nanmedian)

            filtered = pd.merge(precursor_intensity, precursor_rt, on='PRECURSOR_ID', how='outer')

    elif method == 'mean':
        precursor_intensity = data.groupby('PRECURSOR_ID', as_index=False)['APEX_INTENSITY'].agg(np.nanmean)
        precursor_rt = data.groupby('PRECURSOR_ID', as_index=False)['RT_APEX'].agg(np.nanmedian)
        # Merge to get RT
        filtered = pd.merge(precursor_intensity, precursor_rt, on='PRECURSOR_ID',how='left')

    elif method == 'median':
        precursor_intensity = data.groupby('PRECURSOR_ID', as_index=False)['APEX_INTENSITY'].agg(np.nanmedian)
        precursor_rt = data.groupby('PRECURSOR_ID', as_index=False)['RT_APEX'].agg(np.nanmedian)
        # Merge to get RT
        filtered = pd.merge(precursor_intensity, precursor_rt, on='PRECURSOR_ID',how='left')

    elif method == 'sum':
        precursor_intensity = data.groupby('PRECURSOR_ID', as_index=False)['APEX_INTENSITY'].sum()
        precursor_rt = data.groupby('PRECURSOR_ID', as_index=False)['RT_APEX'].agg(np.nanmedian)
        # Merge to get RT
        filtered = pd.merge(precursor_intensity, precursor_rt, on='PRECURSOR_ID',how='left')


    # Merge back to global model protein list
    filtered = pd.merge(precursor_list, filtered, on='PRECURSOR_ID', how='left')
    filtered = filtered.rename(columns={'APEX_INTENSITY':sample, 'RT_APEX':f'RT_{sample}'})
    
    print(f'Info: {sample} quantified by {method} method, takes {round(time.process_time()-start_cpu, 2)}s.')
    
    return filtered

def merge_runs(in_dir, template, method='mean', max_transition_qvalue=0.01, max_fragment_qvalue=0.01):

    """
    This function take a directory contains all scored runs dataframe as input, process every run separately,
    and merge them based on precursors. The output contains precursor information, and the precursor RT and 
    intensity of each sample.
    """

    import pandas as pd
    import os, warnings

    from pandas.core.common import SettingWithCopyWarning
    warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

    os.path.expanduser(os.getcwd())

    files = [os.path.join(x,z) for x,y,zs in os.walk(in_dir) for z in zs]

    # Process a scored matrix, filter out fragments and calculate MS1 intensity
    filtered = filter_fragments(fp=files[0], template=template, method=method, max_transition_qvalue=max_transition_qvalue, max_fragment_qvalue=max_fragment_qvalue)
    
    # Separate data into precursor information part, RT and intensity part, merge them separately
    info_cols = ['PRECURSOR_ID', 'PROTEIN_ACCESSION', 'MODIFIED_SEQUENCE', 'PRECURSOR_CHARGE']
    info = filtered[info_cols]
    intensity_cols = ['PRECURSOR_ID'] + [x for x in filtered.columns.tolist() if x not in info_cols]
    intensity = filtered[intensity_cols]
     
    print(f'Info: Merging {files[0]}')
    for file in files[1:]:
        filtered = filter_fragments(fp=file, template=template, method=method, max_transition_qvalue=max_transition_qvalue, max_fragment_qvalue=max_fragment_qvalue)
        tmp_info = filtered[info_cols]
        intensity_cols = ['PRECURSOR_ID'] + [x for x in filtered.columns.tolist() if x not in info_cols]
        tmp_intensity = filtered[intensity_cols]

        intensity = pd.merge(intensity, tmp_intensity, on='PRECURSOR_ID', how='outer')
        info = pd.concat((info, tmp_info), join='inner')
        print(f'Info: Merging {file}')

    info = info.drop_duplicates()
    merged = pd.merge(intensity, info, on='PRECURSOR_ID', how='left')
    
    return merged

def align_runs(merged, out_path):

    """
    This function calculate median RT for merged quantitative matrix
    """

    import pandas as pd, numpy as np

    # Separate data for calculating  median RT
    rts = merged[[x for x in merged.columns if x[:2]=='RT']]
    rest = merged[[x for x in merged.columns if x[:2]!='RT']]
    rts_median = np.nanmedian(np.array(rts.values), axis=1)
    data = rest.copy()
    data.insert(0, 'RT', rts_median)

    # If precursor intensity are nan in all samples, drop this precursor
    info_cols = ["PRECURSOR_ID", "PROTEIN_ACCESSION", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE", "RT"]
    info = data[info_cols]
    intensity = data[[x for x in data.columns.tolist() if x not in info_cols]]
    intensity_cols = intensity.columns.tolist()
    info = info.values
    intensity = intensity.values
    new_intensity = intensity[~np.isnan(intensity).all(axis=1)]
    new_info = info[~np.isnan(intensity).all(axis=1)]
    data = np.hstack((new_info,new_intensity))
    data = pd.DataFrame(data, columns=info_cols+intensity_cols)

    # Put fragment peptide charge protein columns at the first columns
    mid = data.PRECURSOR_CHARGE.astype("int64")
    data = data.drop('PRECURSOR_CHARGE', axis=1)
    data.insert(0, 'Charge', mid)

    mid = data.PRECURSOR_ID.astype("int64")
    data = data.drop('PRECURSOR_ID', axis=1)
    data.insert(0, 'Precursor', mid)

    mid = data.MODIFIED_SEQUENCE
    data = data.drop('MODIFIED_SEQUENCE', axis=1)
    data.insert(0, 'Peptide', mid)
    
    mid = data.PROTEIN_ACCESSION
    data = data.drop('PROTEIN_ACCESSION', axis=1)
    data.insert(0, 'Protein', mid)

    data = data.sort_values(by=['Protein', 'Peptide', 'Precursor', 'Charge'])

    # Save quantitative matrix to file
    data.to_csv(out_path, sep='\t', index=True)
    print(f'Info: {out_path} is saved.')
