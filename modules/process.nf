#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Syntax for defining a process
process metaphlan {
    // Docker/Singularity container used to run the process
    container "${params.container__metaphlan}"
    // Write output files to the output directory
    publishDir "${params.output}", mode: "copy", overwrite: true
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"
    
    input:
    // Input from a single file
    tuple val(sample_name), path(R1), path(R2)
    // Reference Database Files
    path "db/"

    output:
    // Capture all output files
    path "${sample_name}.metaphlan"

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
    --nproc ${task.cpus}
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