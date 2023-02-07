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


def read_metaphlan(fp):

    # Read the table
    df = pd.read_csv(
        fp,
        sep="\t",
        comment="#",
        header=None,
        names=["clade_name", "NCBI_tax_id", "relative_abundance", "additional_species"]
    )

    # Add the taxonomic level
    df = df.assign(
        level=df["clade_name"].apply(
            lambda n: dict(
                k="kingdom",
                p="phylum",
                c="class",
                o="order",
                f="family",
                g="genus",
                s="species",
                t="strain"
            ).get(
                n.rsplit("|")[-1][0],
                'unknown'
            )
        ),
        org_name=df["clade_name"].apply(
            lambda n: n.rsplit("|")[-1][3:] if len(n.rsplit("|")[-1]) > 3 else n
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