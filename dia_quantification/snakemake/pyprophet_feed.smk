#!/usr/bin/env python3

SAMPLES, = glob_wildcards(f"{config['base_file_paths']['mzml']}/{{sample}}.mzML")


rule all:
    input:
        long = f"{config['base_file_paths']['pyprophet']}/tric/feature_alignment.tsv",
        matrix = f"{config['base_file_paths']['pyprophet']}/tric/feature_alignment_matrix.tsv"

rule pyprophet_apply:
    input:
        osw = f"{config['base_file_paths']['osw']}/{{sample}}.osw",
    output:
        score_placeholder = f"{config['base_file_paths']['pyprophet']}/scored/{{sample}}.osw"
    #container:
        #config['containers']['pyprophet']
    params:
        group_id = 'feature_id',
        level = 'ms1ms2',
        scoring_model = f"{config['base_file_paths']['pyprophet']}/models/scoring_model.osw"
    threads: 1
    shell:
        (
            "pyprophet score "
            "--in {input.osw} "
            "--out {output.score_placeholder} "
            "--group_id {params.group_id} "
            "--apply_weights {params.scoring_model} "
            "--level {params.level}"
        )

rule pyprophet_reduce:
    input:
        score_placeholder = f"{config['base_file_paths']['pyprophet']}/scored/{{sample}}.osw"
    output:
        oswr = f"{config['base_file_paths']['pyprophet']}/reduced/{{sample}}.oswr"
    shell:
        (  
            "pyprophet reduce "
            "--in {input.score_placeholder} "
            "--out {output.oswr}"
        )

rule pyprophet_global:
    input:
        reduced = expand(f"{config['base_file_paths']['pyprophet']}/reduced/{{sample}}.oswr", sample=SAMPLES)
    output:
        global_model = f"{config['base_file_paths']['pyprophet']}/models/model_global.osw"
    params:
        context = 'global',
        template = f"{config['base_file_paths']['pyprophet']}/models/scoring_model.osw"
    shell:
        (
            "pyprophet merge "
            "--template {params.template} "
            "--out {output.global_model} "
            "{input.reduced} && "
            "pyprophet peptide "
            "--context {params.context} "
            "--in {output.global_model} && "
            "pyprophet protein "
            "--context {params.context} "
            "--in {output.global_model}"
        )

rule pyprophet_backpropagate:
    input:
        score_placeholder = f"{config['base_file_paths']['pyprophet']}/scored/{{sample}}.osw",
        global_model = f"{config['base_file_paths']['pyprophet']}/models/model_global.osw"
    output:
        backpropagate_placeholder = f"{config['base_file_paths']['pyprophet']}/backpropagated/{{sample}}.osw"
    shell:
        (            
            "pyprophet backpropagate "
            "--apply_scores {input.global_model} "
            "--in {input.score_placeholder} "
            "--out {output.backpropagate_placeholder}"
        )

rule pyprophet_export:
    input:
        backpropagate_placeholder = f"{config['base_file_paths']['pyprophet']}/backpropagated/{{sample}}.osw"
    output:
        tsv = f"{config['base_file_paths']['pyprophet']}/export/{{sample}}.tsv"
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
            "--in {input.backpropagate_placeholder} "
            "--out {output.tsv} "
            "--format {params.format} "
            "--max_rs_peakgroup_qvalue {params.max_rs_peakgroup_qvalue} "
            "--max_global_peptide_qvalue {params.max_global_peptide_qvalue} "
            "--max_global_protein_qvalue {params.max_global_protein_qvalue}"
        )
        
rule tric_feature_alignment:
    input:
        tsv = expand(f"{config['base_file_paths']['pyprophet']}/export/{{sample}}.tsv", sample=SAMPLES)
    output:
        long = f"{config['base_file_paths']['pyprophet']}/tric/feature_alignment.tsv",
        matrix = f"{config['base_file_paths']['pyprophet']}/tric/feature_alignment_matrix.tsv"
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