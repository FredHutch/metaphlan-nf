# metaphlan-nf
Quantify microbial abundance from whole-genome shotgun sequencing data

## References

- MetaPhlAn-3.0 - [https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0)
- Running Nextflow at Fred Hutch - [SciWiki - Nextflow Run Script](https://sciwiki.fredhutch.org/hdc/workflows/running/run_script/)
- Configuring Singularity for Nextflow - [SciWiki - Singularity / Gizmo](https://sciwiki.fredhutch.org/hdc/workflows/running/on_gizmo/)

## Running the Workflow

```
Usage:

nextflow run FredHutch/metaphlan-nf <ARGUMENTS>

Arguments:

  Input Data:
  --input_folder    Folder containing all input files (FASTQ pairs)
  --file_suffix     File ending for all input FASTQ files (default: ${params.file_suffix})

  Reference Database
  --db              Path to reference database
                    (must include the prefix shared across reference files)

  Output Location:
  --output          Folder for output files

```