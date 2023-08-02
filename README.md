# metaphlan-nf
Quantify microbial abundance from whole-genome shotgun sequencing data

## References

- MetaPhlAn-3.0 - [https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0)
- Running Nextflow at Fred Hutch - [SciWiki - Nextflow Run Script](https://sciwiki.fredhutch.org/hdc/workflows/running/run_script/)
- Configuring Singularity for Nextflow - [SciWiki - Singularity / Gizmo](https://sciwiki.fredhutch.org/hdc/workflows/running/on_gizmo/)

## Running MetaPhlAn

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

## Running StrainPhlAn

To run the StrainPhlan tool, use the `strainphlan.nf` entrypoint
as the `-main-script` argument for the Nextflow CLI.

```
Usage:

nextflow run FredHutch/metaphlan-nf -main-script strainphlan.nf <ARGUMENTS>

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

  Optional Metadata Table
  --metadata        Path to TSV containing sample metadata
                    First column must have the header: sampleID

  StrainPhlAn Parameters
  --breadth_threshold   Threshold defining the minimum breadth of coverage
                        for the markers as a percentage value (default: 80)
  --min_reads_aligning  Minimum number of reads aligning to a clade
                        per sample (default: 8)
  --min_base_coverage   Minimum threshold for the number of reads aligning
                        to a given position on the marker sequence (default: 1)
  --min_base_quality    Minimum base quality for alignment (default: 30)
  --min_mapping_quality Minimum mapping quality score (default: 10)
  --marker_in_n_samples Minimum number of samples required to keep a marker
                        (default: 2, note this is an absolute number)
  --sample_with_n_markers
                        Minimum percentage of markers required to retain
                        a sample (default: 80, note this is a percentage)

  Output Location:
  --output          Folder for output files
  ```
