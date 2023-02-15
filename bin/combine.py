#!/usr/bin/env python3
"""Combine a collection of metaphlan output files"""

import pandas as pd
import os

def combine_metaphlan(fp_dict):
    
    # Join all of the results together
    df = pd.concat([
        read_metaphlan(fp).assign(name=name)
        for name, fp in fp_dict.items()
    ])

    # Write them out in long format
    df.to_csv("metaphlan.long.csv.gz", index=None)
    print("Wrote out long")

    # Write out each taxonomic level in wide format
    for level, level_df in df.groupby("level"):

        level_df.pivot_table(
            index="org_name",
            columns="name",
            values="relative_abundance"
        ).fillna(
            0
        ).to_csv(
            f"metaphlan.wide.{level}.csv.gz"
        )
        print(f"Wrote out {level}")

    print("Done")


def get_org_name(clade_name: str):

    if clade_name == "UNCLASSIFIED":
        return

    else:

        # Get the last name in the list
        org = clade_name.rsplit("|", 1)[-1]

        # The name should always start with X__
        msg = f"Clade name does not fit the expected pattern: {clade_name}"
        assert len(org) > 3, msg

        return org[3:]


def get_tax_level(clade_name: str):

    if clade_name.upper() == "UNCLASSIFIED":
        return

    else:

        # Get the last name in the list
        org = clade_name.rsplit("|", 1)[-1]

        # The name should always start with X__
        msg = f"Clade name does not fit the expected pattern: {clade_name}"
        assert len(org) > 3, msg

        tax_level = dict(
            k="kingdom",
            p="phylum",
            c="class",
            o="order",
            f="family",
            g="genus",
            s="species",
            t="strain"
        ).get(
            org[0]
        )

        assert tax_level is not None, msg

        return tax_level


def read_metaphlan(fp):

    # Read the table
    df = pd.read_csv(
        fp,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "clade_name",
            "clade_taxid",
            "relative_abundance",
            "coverage",
            "estimated_number_of_reads_from_the_clade"
        ]
    )

    # Add the taxonomic level
    df = df.assign(
        level=df["clade_name"].apply(
            get_tax_level
        ),
        org_name=df["clade_name"].apply(
            get_org_name
        )
    )

    return df

if __name__ == "__main__":
    combine_metaphlan(
        {
            fp.replace(".metaphlan", ""): fp
            for fp in os.listdir(".")
            if fp.endswith(".metaphlan")
        }
    )