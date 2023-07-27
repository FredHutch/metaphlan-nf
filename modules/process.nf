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
    // Input from a FASTQ file
    tuple val(sample_name), path(fastq)
    // Reference Database Files
    path "db/"

    output:
    // Capture just the aligned reads
    tuple val(sample_name), path("${sample_name}.bowtie2.bz2"), emit: bwt
    tuple val(sample_name), path("${sample_name}.sam.bz2"), emit: sam

"""#!/bin/bash

set -e

echo Processing sample : '${sample_name}'

metaphlan \
    --input_type fastq \
    --bowtie2db db \
    --index ${params.db.replaceAll(".*/", "")} \
    ${fastq} \
    -o ${sample_name}.metaphlan \
    -s ${sample_name}.sam.bz2 \
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
    path "${sample_name}.metaphlan", emit: metaphlan
    path "${sample_name}.biom"

"""#!/bin/bash

set -e

echo Processing sample : '${sample_name}'

metaphlan \
    --nproc ${task.cpus} \
    --input_type bowtie2out \
    --biom ${sample_name}.biom \
    -t rel_ab_w_read_stats \
    --unclassified_estimation \
    --bowtie2db db \
    --index ${params.db.replaceAll(".*/", "")} \
    --sample_id_key "${sample_name}" \
    --sample_id "${sample_name}" \
    ${bowtie2out} \
    ${sample_name}.metaphlan

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
    path "merged_abundance_table.tsv"

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

process concat_sam {
    container "${params.container__samtools}"
    publishDir "${params.output}/sam/", mode: "copy", overwrite: true
    
    input:
    tuple val(sample_name), path("inputs/*.sam.bz2")

    output:
    tuple val(sample_name), path("${sample_name}.sam.bz2")

"""#!/bin/bash
set -e
echo STARTING
find .
if (( \$(find inputs -name '*.sam.bz2' | wc -l) == 1 )); then
    echo "Only a single SAM file found"
    cp inputs/*.sam.bz2 ./
else
    echo "Multiple SAM files found"
    echo "Decompressing SAM files"
    bunzip2 inputs/*.bz2
    echo "Merging SAM files"
    samtools merge "${sample_name}.sam" inputs/*.sam
    find .
    echo "Compressing merged SAM file"
    bzip2 "${sample_name}.sam"
fi
echo DONE
find .
"""
}

process concat_bwt {
    container "${params.container__pandas}"
    publishDir "${params.output}/bowtie2/", mode: "copy", overwrite: true
    
    input:
    tuple val(sample_name), path("inputs/*.bowtie2.bz2")

    output:
    tuple val(sample_name), path("${sample_name}.bowtie2.bz2"), emit: bowtie

"""#!/bin/bash
set -e
concat_alignments.py "${sample_name}.bowtie2.bz2"
"""
}

process sample2markers {
    container "${params.container__metaphlan}"
    publishDir "${params.output}", mode: "copy", overwrite: true
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"
    
    input:
    path "db/"
    path "sams/"

    output:
    path "consensus_markers/*", optional: true

"""#!/bin/bash
set -e

# Only produce outputs if the inputs contain reads to start with
if (( \$(bunzip2 -c sams/*.sam.bz2 | grep -v '^@' | wc -l) > 1 )); then

    mkdir -p consensus_markers
    sample2markers.py -d db/ -i sams/*.sam.bz2 -o consensus_markers -n ${task.cpus}

fi
"""
}

process extract_markers {
    container "${params.container__metaphlan}"
    publishDir "${params.output}", mode: "copy", overwrite: true
    
    input:
    path "db/"
    val clade

    output:
    tuple val(clade), path("db_markers/*")

"""#!/bin/bash
set -e
mkdir -p db_markers
extract_markers.py -d db/*.pkl -c ${clade} -o db_markers/
"""
}

process strainphlan {
    container "${params.container__metaphlan}"
    publishDir "${params.output}", mode: "copy", overwrite: true
    
    input:
    path "db/"
    path "consensus_markers/"
    tuple val(clade), path(db_markers)
    path "reference_genomes/"

    output:
    path "${clade}/*", emit: all
    tuple val(clade), path("${clade}/RAxML_bestTree.*.tre"), emit: tre

"""#!/bin/bash
set -e
mkdir -p ${clade}
if [ -d reference_genomes ] && (( \$(find reference_genomes | wc -l) > 0 )); then
    FLAGS="-r reference_genomes/*"
else
    FLAGS=""
fi

strainphlan \
    -s consensus_markers/*.pkl \
    -d db/*.pkl \
    -m ${db_markers} \
    -o ${clade} \
    -n ${task.cpus} \
    -c ${clade} \
    --mutation_rates \
    \$FLAGS \
"""
}

process get_metadata {
    container "${params.container__pandas}"
    
    input:
    path "metadata.tsv"

    output:
    path "metadata_categories.txt"

    script:
"""
set -e
get_metadata_categories.py
"""
}

process add_metadata {
    container "${params.container__metaphlan}"
    publishDir "${params.output}/${clade}/${metadata_category}/", mode: "copy", overwrite: true, pattern: "*.tre.metadata"
    
    input:
    tuple val(clade), path(tre), val(metadata_category), path("metadata.tsv")

    output:
    tuple val(clade), path("${tre}.*"), val(metadata_category)

"""
echo "Adding metadata for ${metadata_category} to ${tre}"
add_metadata_tree.py \
    -t ${tre} \
    -f metadata.tsv \
    -m "${metadata_category}"

plot_tree_graphlan.py \
    -t ${tre}.metadata \
    -m "${metadata_category}" || true
"""
}

process plot_metadata {
    container "${params.container__graphlan}"
    publishDir "${params.output}/${clade}/${metadata_category}", mode: "copy", overwrite: true
    
    input:
    tuple val(clade), path("*"), val(metadata_category)

    output:
    path "${clade}.${metadata_category}*"

"""#!/bin/bash
set -e
echo "Making a plot for ${metadata_category} from ${clade}"

graphlan_annotate.py \
    --annot *.annot \
    *.graphlantree

for suffix in pdf png svg; do
    graphlan.py \
        *.graphlantree \
        "${clade}.${metadata_category}.\$suffix"
done
"""
}