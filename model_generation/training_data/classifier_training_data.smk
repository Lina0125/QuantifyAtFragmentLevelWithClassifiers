#!/usr/bin/env python3
SAMPLES, = glob_wildcards(f"{config['base_file_paths']['mzml']}/{{sample}}.mzML")

rule all:
    input:
        f"{config['training_data']['merged_csv_path']}"
'''TODO before 0326
rule openswath:
    input:
        pqp = config['libraries']['spectral_library'],
        mzml = f"{config['base_file_paths']['mzml']}/{{sample}}.mzML",
        rt = config['libraries']['rt_library'],
        swath_windows_file = config['swath_windows_file']
    output:
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw",
        chrom = f"{config['base_file_paths']['osw']}/{{sample}}.sqMass" 
    params:
        min_upper_edge_dist = 1,
        mz_extraction_window = 30,
        mz_extraction_window_unit = 'ppm',
        mz_extraction_window_ms1 = 20,
        mz_extraction_window_ms1_unit = 'ppm',
        irt_mz_extraction_window = 50,
        irt_mz_extraction_window_unit = 'ppm',
        rt_extraction_window = 600,
        rt_normalization_alignment_method = 'lowess',
        scoring_stop_report_after_feature = 5
    threads: 50
    container:
        "/srv/data1/li7186lu/data/container/openms.sif"
    shell:
        (
            "OpenSwathWorkflow "
            "-in {input.mzml} "
            "-tr {input.pqp} "
            "-tr_irt {input.rt} "
            "-out_osw {output.osw} "
            "-out_chrom {output.chrom} "
            "-enable_ms1 true "
            "-enable_ipf true "
            "-threads {threads} "
            "-Scoring:Scores:use_uis_scores "
            "-Scoring:Scores:use_total_mi_score "
            "-min_upper_edge_dist {params.min_upper_edge_dist} "
            "-mz_extraction_window {params.mz_extraction_window} "
            "-mz_extraction_window_unit {params.mz_extraction_window_unit} "
            "-mz_extraction_window_ms1 {params.mz_extraction_window_ms1} "
            "-mz_extraction_window_ms1_unit {params.mz_extraction_window_ms1_unit} "
            "-Library:retentionTimeInterpretation minutes "
            "-RTNormalization:alignmentMethod {params.rt_normalization_alignment_method} "
            "-Scoring:stop_report_after_feature {params.scoring_stop_report_after_feature}"
        )
'''

rule openswath:
    input:
        pqp = config['libraries']['spectral_library'],
        mzml = f"{config['base_file_paths']['mzml']}/{{sample}}.mzML",
        rt = config['libraries']['rt_library']
    output:
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw",
        chrom = f"{config['base_file_paths']['osw']}/{{sample}}.sqMass" 
    params:
        swath_windows_file = config['swath_windows_file'],
        min_upper_edge_dist = 1,
        mz_extraction_window = 50,
        mz_extraction_window_unit = 'ppm',
        mz_extraction_window_ms1 = 20,
        mz_extraction_window_ms1_unit = 'ppm',
        irt_mz_extraction_window = 50,
        irt_mz_extraction_window_unit = 'ppm',
        rt_extraction_window = 600,
        rt_normalization_alignment_method = 'lowess',
        scoring_stop_report_after_feature = 5,
        batch_size = 1000, 
        min_rsq = 0.95
    threads: 50
    container:
        "/srv/data1/li7186lu/data/container/openms.sif"
    shell:
        (
        "OpenSwathWorkflow "
        "-in {input.mzml} "
        "-tr {input.pqp} "
        "-out_osw {output.osw} "
        "-swath_windows_file {params.swath_windows_file} "
        "-tr_irt {input.rt} "
        "-out_chrom {output.chrom} "
        "-threads {threads} "
        "-enable_ms1 true "
        "-enable_ipf true "
        "-Scoring:Scores:use_uis_scores "
        "-Scoring:Scores:use_total_mi_score "
        "-Scoring:stop_report_after_feature {params.scoring_stop_report_after_feature} "
        "-Scoring:TransitionGroupPicker:compute_peak_quality "
        "-min_rsq {params.min_rsq} "
        "-batchSize {params.batch_size}"
        )

rule osw_to_csv:
    input:
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw",
        chrom = f"{config['base_file_paths']['osw']}/{{sample}}.sqMass"
    output:
        f"{config['training_data']['csvs']}/{{sample}}.csv"
    shell:
        (
    	"python "
        "/srv/data1/li7186lu/scripts/classifier/training_data/oswSQL2csv.py "
        "-i {input.osw} "
        "-o {output}"
    	)


rule merge_csvs:
    params:
        input_dir = f"{config['training_data']['csvs']}"
    input:
        expand(f"{config['training_data']['csvs']}/{{sample}}.csv", sample=SAMPLES)
    output:
        f"{config['training_data']['merged_csv_path']}"
    shell:
        (
    	"python "
        "/srv/data1/li7186lu/scripts/classifier/training_data/merge_csv_files.py "
        "-d {params.input_dir} "
        "-o {output}"
    	)

