#!/usr/bin/env python3

import pandas as pd
import jinja2

csv_fp = "metaphlan.long.csv.gz"
print(f"Reading in {csv_fp}")
df = pd.read_csv(csv_fp)
print("Removing unused columns")
df = df.drop(
    columns=["clade_name", "clade_taxid"]
)

print("Loading the template")
templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "template.jinja"
template = templateEnv.get_template(TEMPLATE_FILE)

print("Adding the data")
report = template.render(
    data=df.to_csv(
        index=None
    )
)

print("Saving the report")
with open("metaphlan_report.html", "w") as handle:
    handle.write(report)
