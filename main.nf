#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// All of the default parameters are being set in `nextflow.config`

// Import the process
include { metaphlan; combine; report } from './modules/process'


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

    // If the samplesheet input was specified
    if ( params.samplesheet ){
        Channel
            .fromPath(
                "${params.samplesheet}",
                checkIfExists: true,
                glob: false
            )
            .splitCsv(
                header: true
            )
            .map {
                row -> [
                    row.sample,
                    file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true)
                ]
            }
            .set { input_ch }
    } else {

        // Make a channel with all of the files from the --input_folder
        Channel
            .fromFilePairs([
                "${params.input_folder}/*${params.file_spacer}{1,2}.${params.file_suffix}"
            ])
            .ifEmpty { error "No file pairs found at ${params.input_folder}/*${params.file_spacer}{1,2}.${params.file_suffix}" }
            .map {it -> [it[0], it[1][0], it[1][1]]}
            .set { input_ch }

    }

    // Make a channel with the reference database files
    Channel
        .fromPath("${params.db}*")
        .ifEmpty { error "No database files found at ${params.db}*" }
        .toSortedList()
        .set {db_ch}

    // Run the process on the data
    metaphlan(input_ch, db_ch)

    // Combine the results
    combine(
        metaphlan
            .out
            .toSortedList()
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