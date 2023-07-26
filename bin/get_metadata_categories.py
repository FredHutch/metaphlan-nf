#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

metadata = Path("metadata.tsv")
assert metadata.exists(), "File does not exist: metadata.tsv"

print("Reading from metadata.tsv")
metadata = pd.read_csv(metadata, sep="\t")
print(metadata.to_csv(index=None))

msg = "The first column must be 'sampleID'"
assert metadata.columns.values[0] == "sampleID", msg

msg = "There must be more than one column"
assert metadata.shape[1] > 1, msg

# Write out the other column names to a file
with open("metadata_categories.txt", "w") as handle:
    for cname in metadata.columns.values[1:]:
        print(cname)
        handle.write(cname + "\n")

print("DONE")
