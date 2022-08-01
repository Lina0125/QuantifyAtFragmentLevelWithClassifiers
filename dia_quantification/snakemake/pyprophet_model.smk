#!/usr/bin/env python3

SAMPLES, = glob_wildcards(f"{config['base_file_paths']['mzml']}/{{sample}}.mzML")


rule all:
    input:
        model = f"{config['base_file_paths']['pyprophet']}/models/scoring_model.osw"

rule openswath:
    input:
        pqp = config['libraries']['spectral_library'],
        mzml = expand(f"{config['base_file_paths']['mzml']}/{{sample}}.mzML", sample=SAMPLES),
        rt = config['libraries']['rt_library']
    output:
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw",
        chrom = f"{config['base_file_paths']['osw']}/{{sample}}.sqMass" 
    params:
        min_upper_edge_dist = 1,
        mz_extraction_window = 30,
        mz_extraction_window_unit = 'ppm',
        mz_extraction_window_ms1 = 20,
        mz_extraction_window_ms1_unit = 'ppm',
        irt_mz_extraction_window = 100,
        irt_mz_extraction_window_unit = 'ppm',
        rt_extraction_window = 600,
        rt_normalization_alignment_method = 'lowess',
        scoring_stop_report_after_feature = 5,
        min_rsq = 0.7,
        batch_size = 10000,
        swath_windows_file = config['swath_windows_file']
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
            "-swath_windows_file {params.swath_windows_file} "
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


rule pyprophet_subsample:
    input:
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw"
    output:
        subsampled = f"{config['base_file_paths']['pyprophet']}/subsampled/{{sample}}.osws"
    #container:
        #config['containers']['pyprophet']
    threads: 1
    params:
        subsample_ratio = 0.1
    shell:
        (
            "pyprophet subsample "
            "--subsample_ratio {params.subsample_ratio} "
            "--in {input.osw} "
            "--out {output.subsampled}"
        )
rule pyprophet_merge:
    input:
        template = config['libraries']['spectral_library'],
        osws = expand(f"{config['base_file_paths']['pyprophet']}/subsampled/{{sample}}.osws", sample=SAMPLES)
    output:
        merged = f"{config['base_file_paths']['pyprophet']}/training_data/merged.osw"
    wildcard_constraints:
        loop="[0-9]+"
    #container:
        #config['containers']['pyprophet']
    threads: 1
    shell:
        (
            "pyprophet merge "
            "--template {input.template} "
            "--out {output.merged} "
            "-- {input.osws}"
        )
rule pyprophet_learn:
    input:
        rules.pyprophet_merge.output.merged
    output:
        model = f"{config['base_file_paths']['pyprophet']}/models/scoring_model.osw"
    wildcard_constraints:
        loop="[0-9]+"
    #container:
        #config['containers']['pyprophet']
    threads: 50
    params:
        level = 'ms1ms2',
        xevel_num_iter = 3,
        ss_initial_fdr=0.05,
        ss_iteration_fdr=0.01,
    shell:
        (
            "pyprophet score "
            "--in {input} "
            "--out {output.model} "
            "--level {params.level} "
            "--xeval_num_iter {params.xevel_num_iter} "
            "--ss_initial_fdr {params.ss_initial_fdr} "
            "--ss_iteration_fdr {params.ss_iteration_fdr} "
            "--threads {threads}"
        )
        
''''
rule pyprophet_apply:
    input:
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw",
        model = f"{config['base_file_paths']['pyprophet']}/models/scoring_model.osw"
    output:
        score_placeholder = f"{config['base_file_paths']['pyprophet']}/scored/{{sample}}"
    #container:
        #config['containers']['pyprophet']
    params:
        group_id = 'feature_id',
        level = 'ms1ms2'
    threads: 1
    shell:
        (
            "pyprophet score "
            "--in {input.osw} "
            "--group_id {params.group_id} "
            "--apply_weights {input.model} "
            "--level {params.level} && "
            "touch {output.score_placeholder}"
        )

rule pyprophet_reduce:
    input:
        score_placeholder = f"{config['base_file_paths']['pyprophet']}/scored/{{sample}}",
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw"
    output:
        oswr = f"{config['base_file_paths']['pyprophet']}/models/global/{{sample}}.oswr"
    #container:
        #config['containers']['pyprophet']
    threads: 1
    shell:
        (  
            "pyprophet reduce "
            "--in {input.osw} "
            "--out {output.oswr}"
        )
rule pyprophet_merge:
    input:
        oswr = expand(f"{config['base_file_paths']['pyprophet']}/models/global/{{sample}}.oswr", sample=SAMPLES)
    output:
        global_model = f"{config['base_file_paths']['pyprophet']}/models/global/global_model.oswr",
        merged = f"{config['base_file_paths']['pyprophet']}/models/global/merged"
    #container:
        #config['containers']['pyprophet']
    params:
        context = 'global',
        template = config['libraries']['spectral_library']
    shell:
        (
            "pyprophet merge "
            "--template {params.template} "
            "--out {output.global_model} "
            "-- {input.oswr} && "
            "touch {output.merged}"
        )
rule pyprophet_global:
    input:
        merged = f"{config['base_file_paths']['pyprophet']}/models/global/merged",
        global_model = f"{config['base_file_paths']['pyprophet']}/models/global/global_model.oswr",
    output:
        modeled = f"{config['base_file_paths']['pyprophet']}/models/global/modeled"
    #container:
        #config['containers']['pyprophet']
    params:
        context = 'global'
    shell:
        (
            "pyprophet peptide "
            "--context {params.context} "
            "--in {input.global_model} && "
            "pyprophet protein "
            "--context {params.context} "
            "--in {input.global_model} && "
            "touch {output.modeled}"
        )
rule pyprophet_backpropagate:
    input:
        modeled = f"{config['base_file_paths']['pyprophet']}/models/global/modeled",
        merged = f"{config['base_file_paths']['pyprophet']}/models/global/merged",
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw",
    output:
        backpropagate_placeholder = f"{config['base_file_paths']['pyprophet']}/scored/{{sample}}_global"
    #container:
        #config['containers']['pyprophet']
    params:
        global_model = f"{config['base_file_paths']['pyprophet']}/models/global/global_model.oswr"
    threads: 1
    shell:
        (            
            "pyprophet backpropagate "
            "--apply_scores {params.global_model} "
            "--in {input.osw} && "
            "touch {output.backpropagate_placeholder}"
        )
rule pyprophet_export:
    input:
        backpropagate_placeholder = f"{config['base_file_paths']['pyprophet']}/scored/{{sample}}_global",
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw"
    output:
        tsv = f"{config['base_file_paths']['pyprophet']}/tric/{{sample}}.tsv"
    wildcard_constraints:
        loop="[0-9]+"
    #container:
        #config['containers']['pyprophet']
    params:
        format = 'legacy_merged',
        max_rs_peakgroup_qvalue = 0.01,
        max_global_peptide_qvalue = 0.01,
        max_global_protein_qvalue = 0.01
    shell:
        (
            "pyprophet export "
            "--in {input.osw} "
            "--out {output.tsv} "
            "--format {params.format} "
            "--max_rs_peakgroup_qvalue {params.max_rs_peakgroup_qvalue} "
            "--max_global_peptide_qvalue {params.max_global_peptide_qvalue} "
            "--max_global_protein_qvalue {params.max_global_protein_qvalue}"
        )
rule tric_feature_alignment:
    input:
        tsv = expand(f"{config['base_file_paths']['pyprophet']}/tric/{{sample}}.tsv", sample=SAMPLES)
    output:
        long = f"{config['base_file_paths']['pyprophet']}/tric/feature_alignment.tsv",
        matrix = f"{config['base_file_paths']['pyprophet']}/tric/feature_alignment_matrix.tsv"
    #container:
        #config['containers']['msproteomicstools']
    params:
        max_rt_diff = 30,
        mst_use_rt_correction = "False",
        mst_stdev_multiplier = 3.0,
        fdr_cutoff = 0.01,
        max_fdr_quality = 0.01
    shell:
        (
            "python /srv/data1/li7186lu/bin/masterproject/lib/python3.8/site-packages/msproteomicstools/feature_alignment.py "
            "--in {input} "
            "--out {output.long} "
            "--out_matrix {output.matrix} "
            "--max_rt_diff {params.max_rt_diff} "
            "--mst:useRTCorrection {params.mst_use_rt_correction} "
            "--fdr_cutoff {params.fdr_cutoff} "
            "--max_fdr_quality {params.max_fdr_quality}"
        )
'''