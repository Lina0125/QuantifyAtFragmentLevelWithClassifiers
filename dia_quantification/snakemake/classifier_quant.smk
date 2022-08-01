#!/usr/bin/env python3

SAMPLES, = glob_wildcards(f"{config['base_file_paths']['osw']}/{{sample}}.osw")


rule all:
    input:
        f"{config['classifier_results']['merged_matirx']}"
        
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

rule preprocess_merge_osw:
    input:
        f"{config['base_file_paths']['pyprophet']}/backpropagated/{{sample}}.osw"
    output:
        f"{config['classifier_results']['feeding_tsv']}/{{sample}}.tsv"
    shell:
        (
            "python "
            "/srv/data1/li7186lu/scripts/quantification/classifier_filter_fragments/merge_osw.py "
            "-i {input} "
            "-o {output}"
    	)

rule fragment_qvalue:
    params:
        classifier = f"{config['classifier_model']['classifier_file']}",
        scaler = f"{config['classifier_model']['scaler_file']}"
    input:
        f"{config['classifier_results']['feeding_tsv']}/{{sample}}.tsv"
    output:
        f"{config['classifier_results']['scored']}/{{sample}}.tsv"
    shell:
        (
            "python "
            "/srv/data1/li7186lu/scripts/quantification/classifier_filter_fragments/predict.py "
            "-i {input} "
            "-c {params.classifier} "
            "-s {params.scaler} "
            "-o {output}"
    	)

rule merge:
    params:
        in_dir = f"{config['classifier_results']['scored']}",
        method = f"{config['classifier_results']['global_method']}",
        max_transition_qvalue = 0.01,
        max_fragment_qvalue = 0.01,
        global_model = f"{config['base_file_paths']['pyprophet']}/models/model_global.osw"
    input:
        expand(f"{config['classifier_results']['scored']}/{{sample}}.tsv", sample=SAMPLES)
    output:
        f"{config['classifier_results']['merged_matirx']}"
    threads: 50
    shell:
        (
            "python "
            "/srv/data1/li7186lu/scripts/quantification/classifier_filter_fragments/global_quantify.py "
            "-i {params.in_dir} "
            "-t {params.global_model} "
            "-o {output} "
            "-m {params.method} "
            "-q1 {params.max_transition_qvalue} "
            "-q2 {params.max_fragment_qvalue}"
    	)