#!/usr/bin/env python3

import logging
from pathlib import Path
import pandas as pd
import sys
from living_figures.bio.fom.widgets.microbiome.main import MicrobiomeExplorer

logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

logging.info("Saving HTML for the MicrobiomeExplorer")

logging.info("Instantiating the object")
microbiome_explorer = MicrobiomeExplorer()

# Turn on a few more plots, one of each type
plots = microbiome_explorer._get_child("plots")
# This is the list of bools which control whether the plot is being shown
plot_flags = plots.get_value()
# Iterate over each type of plot which has been set up
for plot_i, plot_type in enumerate(plots.children[0].options):
    # Show a plot in the list
    plot_flags[plot_i] = True
    # Set that plot to this type
    plots.children[plot_i].set_value(plot_type.label)
# Apply the change which shows all of these plots
plots.set_value(plot_flags)

# Add the microbial abundances
logging.info("Adding microbial abundances")
microbiome_explorer._get_child(
    "data",
    "abund"
).parse_files(
    Path("merged_abundance_table.tsv")
)

# Check to see if there are any samplesheets provided
for samplesheet in Path("samplesheet").rglob("*.csv"):

    # Read the table
    logging.info(f"Reading {samplesheet.name}")
    annots = pd.read_csv(
        samplesheet
    )
    logging.info(annots.to_csv(index=None))

    # Drop any columns with reads
    for cname in ["fastq_1", "fastq_2"]:
        if cname in annots.columns.values:
            logging.info(f"Dropping column: {cname}")
            annots = annots.drop(columns=[cname])
            logging.info(annots.to_csv(index=None))

    # Drop any duplicate rows
    logging.info("Dropping any duplicate rows")
    annots = annots.groupby('sample').head(1)

    # Set the index
    logging.info("Setting the index")
    annots = annots.set_index("sample")
    logging.info(annots.to_csv())

    # Add the samplesheet information to the explorer
    logging.info("Adding annotations to the explorer")
    microbiome_explorer._get_child(
        "data",
        "annots"
    ).value = annots

    # Don't add more than one samplesheet
    break

# Update the submenu options
logging.info("Updating submenus")
microbiome_explorer.update_options()

fp = "MicrobiomeExplorer.html"
logging.info(f"Saving to {fp}")
microbiome_explorer.to_html(Path(fp))
