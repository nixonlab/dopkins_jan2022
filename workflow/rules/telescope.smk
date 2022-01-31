#! /usr/bin/env python
# -*- coding utf-8 -*-

################################### TELESCOPE ###################################

rule telescope:
    conda: "../envs/telescope.yaml"
    output:
        "results/telescope/{s}/{s}-telescope_report.tsv"
    input:
        bam = "results/star_alignment/{s}/{s}_GDC38.Aligned.out.bam",
        annotation = rules.telescope_annotation.output
    benchmark: "benchmarks/telescope/{s}_telescope.tsv"
    log:
        "results/telescope/{s}/telescope.log"
    threads: config['telescope_threads']
    params:
        outdir = "results/telescope/{s}",
        exp_tag = "{s}"
    shell:
        """
        telescope assign\
         --exp_tag {params.exp_tag}\
         --theta_prior 200000\
         --max_iter 200\
         --updated_sam\
         --outdir {params.outdir}\
         --stranded_mode RF\
         {input[0]}\
         {input[1]}\
         2>&1 | tee {log[0]}
        chmod 660 {output[0]}
        """

rule sample_complete:
    input:
        rules.telescope.output
    output:
        touch("results/completed/{s}_completed.txt")
    params:
        star_out = "results/star_alignment/{s}/",
        telescope_out = "results/telescope/{s}/",
        dir_to_move_telescope = config['efs_dir_telescope'],  # "/efs/projects/DLBCL/telescope/"
        dir_to_move_star = config['efs_dir_star'] # "/efs/projects/DLBCL/star_algn/"
    shell:
        """
        mkdir -p {params.dir_to_move_telescope}
        mkdir -p {params.dir_to_move_star}
        cp -R {params.telescope_out} {params.dir_to_move_telescope}
        cp -R {params.star_out} {params.dir_to_move_star}
        rm -r {params.telescope_out}
        rm -r {params.star_out}
        """
