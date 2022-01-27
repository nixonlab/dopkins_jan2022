#! /usr/bin/env python
# -*- coding utf-8 -*-

# Add quant mode for gene TE_counts
# Make sure that star runs with stringtie
################################ STAR ALIGNMENT ################################

rule star_alignment:
    conda:
        "../envs/star.yaml"
    input:
        R1 = config['R1'],
        R2 = config['R2'],
        genome = config['indexes']['star']
    output:
        aligned_bam = "results/star_alignment/{s}/{s}_GDC38.Aligned.out.bam"
    params:
        out_prefix="results/star_alignment/{s}/{s}_GDC38.",
        outFilterMultimapNmax=config['outFilterMultimapNmax'],
        winAnchorMultimapNmax=config['winAnchorMultimapNmax']
    threads: config['star_alignment_threads']
    benchmark: "benchmarks/star_alignment/{s}_star_alignment.tsv"
    shell:
        '''
        STAR\
            --runThreadN {threads}\
            --genomeDir {input.genome}\
            --readFilesIn {input.R1} {input.R2}\
            --outSAMattributes NH HI NM MD AS XS\
            --outSAMtype BAM Unsorted\
            --outFileNamePrefix {params.out_prefix}\
            --quantMode GeneCounts\
            --outSAMstrandField intronMotif\
            --outFilterMultimapNmax {params.outFilterMultimapNmax}\
            --winAnchorMultimapNmax {params.winAnchorMultimapNmax}\
            --outSAMunmapped Within KeepPairs
        '''
