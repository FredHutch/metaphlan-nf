#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// All of the default parameters are being set in `nextflow.config`

// Import the process
include {
    metaphlan_paired;
    metaphlan_single;
    metaphlan_call;
    combine;
    report;
    concat_bwt;
    concat_sam;
    merge
} from './modules/process'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/metaphlan-nf <ARGUMENTS>

Arguments:

  Input Data:
  --input_folder    Folder containing all input files (FASTQ pairs)
  --file_suffix     File ending for all input FASTQ files (default: ${params.file_suffix})
  or
  --samplesheet     CSV listing input files with header sample,fastq_1,fastq_2

  Reference Database
  --db              Path to reference database
                    (must include the prefix shared across reference files)

  Output Location:
  --output          Folder for output files

    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.output == false || params.db == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    // Raise an error if neither input type is specified
    if ( params.samplesheet == false && params.input_folder == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    // Raise an error if BOTH input types are specified
    if ( params.samplesheet && params.input_folder ){
        log.info"""
        You may only specify one input type -- samplesheet OR input_folder
        """.stripIndent()
        exit 1
    }

    // Make a channel with the reference database files
    Channel
        .fromPath("${params.db}*")
        .ifEmpty { error "No database files found at ${params.db}*" }
        .toSortedList()
        .set {db_ch}

    // If the samplesheet input was specified
    if ( params.samplesheet ){
        samplesheet = Channel
            .fromPath(
                "${params.samplesheet}",
                checkIfExists: true,
                glob: false
            )
            .splitCsv(
                header: true
            )
            .branch {
                single: it.fastq_2 == null || it.fastq_2 == ""
                paired: true
            }

        // Paired-end alignment
        metaphlan_paired(
            samplesheet
                .paired
                .map {
                    row -> [
                        row.sample,
                        file(row.fastq_1, checkIfExists: true),
                        file(row.fastq_2, checkIfExists: true)
                    ]
                },
            db_ch
        )

        // Single-end
        metaphlan_single(
            samplesheet
                .single
                .map {
                    row -> [
                        row.sample,
                        file(row.fastq_1, checkIfExists: true)
                    ]
                },
            db_ch
        )

        bwt_ch = metaphlan_paired.out.bwt
            .mix(metaphlan_single.out.bwt)

        sam_ch = metaphlan_paired.out.sam
            .mix(metaphlan_single.out.sam)

    } else {

        // Make a channel with all of the files from the --input_folder
        inputs = Channel
            .fromFilePairs([
                "${params.input_folder}/*${params.file_spacer}{1,2}${params.file_suffix}"
            ])
            .ifEmpty { error "No file pairs found at ${params.input_folder}/*${params.file_spacer}{1,2}${params.file_suffix}" }

        // Run paired-end alignment using metaphlan
        metaphlan_paired(
            inputs.map {
                it -> [it[0], it[1][0], it[1][1]]
            },
            db_ch
        )

        bwt_ch = metaphlan_paired.out.bwt
        sam_ch = metaphlan_paired.out.sam

    }

    // Transform the output to yield a tuple of sample_name, list(alignments.bz2)
    // Then split the channel based on whether there is just one alignment file
    // or multiple for a given sample.
    bwt_ch
        .groupTuple()
        .branch {
            single: it[1].size() == 1
            multiple: it[1].size() > 1
        }
        .set {
            bwt_ch
        }

    sam_ch
        .groupTuple()
        .branch {
            single: it[1].size() == 1
            multiple: it[1].size() > 1
        }
        .set {
            sam_ch
        }

    // If there are any samples with multiple sets of read pairs,
    // concat those alignments into a single file
    concat_bwt(bwt_ch.multiple)
    concat_sam(sam_ch.multiple)

    // Run the metaphlan community profiling algorithm on the combined
    // set of (1) samples which only had a single pair of reads, and
    // (2) the merged alignments from samples with multiple pairs of reads
    metaphlan_call(
        concat_bwt.out.mix(bwt_ch.single),
        db_ch
    )

    // Combine the results
    combine(
        metaphlan_call.out.metaphlan.toSortedList()
    )

    // Merge tables using the metaphlan utility
    merge(
        metaphlan_call.out.metaphlan.toSortedList()
    )

    // Make a summary report
    report(
        combine.out.long,
        file(
            "$projectDir/bin/template.jinja",
            checkIfExists: true,
            glob: false
        )
    )
}