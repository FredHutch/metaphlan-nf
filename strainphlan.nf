#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// All of the default parameters are being set in `nextflow.config`

// Import the process
include {
    sample2markers;
    extract_markers;
    strainphlan;
    get_metadata;
    add_metadata;
    plot_metadata
} from './modules/process'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/metaphlan-nf <ARGUMENTS>

Arguments:

  Input Data:
  --input_folder    Folder containing all input files (*.sam.bz2)
  --clades          Comma-separated list of clades (e.g. t__SGB1877)

  Reference Database
  --db              Path to reference database
                    (must include the prefix shared across reference files)

  Optional Reference Genomes
  --reference_genomes   Comma-delimited list of reference genomes to include
                        (accepts wildcard globs)

  Output Location:
  --output          Folder for output files

    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.output == false || params.input_folder == false || params.clades == false || params.db == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }


    // Make a channel with the reference database files
    Channel
        .fromPath("${params.db}*")
        .ifEmpty { error "No database files found at ${params.db}*" }
        .toSortedList()
        .set {db_ch}

    // Make a channel with all of the .sam.bz2 files from the --input_folder
    inputs = Channel
        .fromPath([
            "${params.input_folder}/*.sam.bz2"
        ])
        .ifEmpty { error "No files found in ${params.input_folder}/*.sam.bz2" }

    // Reconstruct strains present in each sample
    sample2markers(
        db_ch,
        inputs
    )

    // Extract markers from the database for each of the organisms provided by the user
    extract_markers(
        db_ch,
        Channel
            .of(
                "${params.clades}".split(",").toList()
            )
            .flatten()
    )

    // Collect any optionally-provided reference genomes
    if (params.reference_genomes != ""){
        reference_genomes = Channel
            .fromPath(
                "${params.reference_genomes}".split(",").toList(),
                checkIfExists: true
            )
    }else{
        reference_genomes = Channel.empty()
    }

    // Get the phylophlan configuration
    phylophlan_config = file(
        "${projectDir}/lib/phylophlan.config",
        checkIfExists: true
    )

    strainphlan(
        db_ch,
        sample2markers.out.toSortedList(),
        extract_markers.out,
        reference_genomes.toSortedList(),
        phylophlan_config
    )

    if (params.metadata){

        metadata_tsv = file(
            params.metadata,
            checkIfExists: true,
            glob: false
        )

        get_metadata(metadata_tsv)

        strainphlan
                .out
                .tre
                .combine(
                    get_metadata
                        .out
                        .splitText()
                        .map { it.strip("\n") }
                )
                .combine( Channel.of(metadata_tsv) )
                .set { metadata_ch }
        
        add_metadata(metadata_ch) | plot_metadata
    }

}