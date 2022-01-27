#! /usr/bin/env python
# -*- coding utf-8 -*-

################################### FASTQC ###################################

rule fastqc:
    conda: "../envs/fastqc.yaml"
    output:
        R1 = "results/fastqc/{s}/{s}_R1_001_fastqc.html",
        R2 = "results/fastqc/{s}/{s}_R2_001_fastqc.html"
    input:
        R1 = config['R1'],
        R2 = config['R2']
    benchmark:
        "benchmarks/fastqc/{s}_fastqc.tsv"
    params:
        output_prefix = "results/fastqc/{s}"
    log:
        "results/fastqc/{s}/{s}_fastqc.log"
    threads: config['fastqc_threads']
    shell:
        """
        fastqc -t {threads} -o {params.output_prefix} {input} 2> {log}
        """
