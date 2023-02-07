#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Align reads with metaphlan to estimate microbial composition
process metaphlan_align {
    // Docker/Singularity container used to run the process
    container "${params.container__metaphlan}"

    // Resources used
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"
    
    input:
    // Input from a pair of FASTQ files
    tuple val(sample_name), path(R1), path(R2)
    // Reference Database Files
    path "db/"

    output:
    // Capture just the aligned reads
    tuple val(sample_name), path("${sample_name}.bowtie2.bz2"), emit: alignment

"""#!/bin/bash

set -e

echo Processing sample : '${sample_name}'

metaphlan \
    --input_type fastq \
    --bowtie2db db \
    --index ${params.db.replaceAll(".*/", "")} \
    ${R1},${R2} \
    -o ${sample_name}.metaphlan \
    --bowtie2out ${sample_name}.bowtie2.bz2 \
    --sample_id_key "${sample_name}" \
    --sample_id "${sample_name}" \
    --nproc ${task.cpus}
"""

}

// Estimate microbial composition from pre-aligned reads
process metaphlan_call {
    // Docker/Singularity container used to run the process
    container "${params.container__metaphlan}"
    // Write output files to the output directory
    publishDir "${params.output}/bz2/", pattern: "*.bz2", mode: "copy", overwrite: true
    publishDir "${params.output}/mpl/", pattern: "*.metaphlan", mode: "copy", overwrite: true
    publishDir "${params.output}/biom/", pattern: "*.biom", mode: "copy", overwrite: true
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"
    
    input:
    // Input from a single file of alignments
    tuple val(sample_name), path(bowtie2out)
    // Reference Database Files
    path "db/"

    output:
    // Capture all output files
    tuple val(sample_name), path("${sample_name}.bowtie2.bz2"), emit: alignment
    path "${sample_name}.metaphlan", emit: metaphlan
    path "${sample_name}.biom"

"""#!/bin/bash

set -e

echo Processing sample : '${sample_name}'

metaphlan \
    ${bowtie2out} \
    --nproc ${task.cpus} \
    --input_type bowtie2out \
    -o ${sample_name}.metaphlan

#    --biom ${sample_name}.biom \
#    -t rel_ab_w_read_stats \
#    --unclassified_estimation \

"""

}

process combine {
    // Docker/Singularity container used to run the process
    container "${params.container__pandas}"
    // Write output files to the output directory
    publishDir "${params.output}", mode: "copy", overwrite: true
    
    input:
    path "*"

    output:
    // Capture all output files
    path "metaphlan.*.csv.gz", emit: all
    path "metaphlan.long.csv.gz", emit: long

"""#!/bin/bash

set -e

combine.py
"""

}

process merge {
    container "${params.container__metaphlan}"
    // Write output files to the output directory
    publishDir "${params.output}", mode: "copy", overwrite: true
    
    input:
    path "inputs/*"

    output:
    // Capture all output files
    path "merged_abundance_table.txt"

"""
merge_metaphlan_tables.py inputs/* > merged_abundance_table.tsv
"""

}

process report {
    // Docker/Singularity container used to run the process
    container "${params.container__pandas}"
    // Write output files to the output directory
    publishDir "${params.output}", mode: "copy", overwrite: true
    
    input:
    path "metaphlan.long.csv.gz"
    path "template.jinja"

    output:
    path "metaphlan_report.html"

"""#!/bin/bash

set -e

prep.py
"""
}

process concat {
    container "${params.container__pandas}"
    
    input:
    tuple val(sample_name), path("inputs/*.bz2")

    output:
    tuple val(sample_name), path("${sample_name}.bz2")

"""
cat inputs/*.bz2 > "${sample_name}.bz2"
"""
}