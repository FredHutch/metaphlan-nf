#!/usr/bin/env python3

import bz2
from pathlib import Path
import sys

def combine_mpl_alignments(
    input_folder="inputs",
    output_file="metaphlan.bz2"
):

    # Keep track of the running totals for reads and avg_len
    running_nreads = 0
    running_avg_len = 0

    with bz2.open(output_file, "wt") as output_handle:

        for fp in Path(input_folder).iterdir():

            if fp.name.endswith(".bz2"):

                nreads, avg_len = add_reads(fp, output_handle)

                running_nreads += nreads
                running_avg_len += (nreads * avg_len)

        add_stats(
            running_nreads,
            running_avg_len / running_nreads,
            output_handle
        )


def add_stats(nreads, avg_len, output_handle):

    output_handle.write(f"nreads\t{nreads}\n")
    output_handle.write(f"avg_read_length\t{avg_len}\n")


def add_reads(input_fp, output_handle):

    nreads = None
    avg_len = None

    with bz2.open(input_fp, "rt") as input_handle:

        for line in input_handle:

            if line.startswith("#"):

                tag, val = line[1:].rstrip("\n").split("\t", 1)

                try:
                    val = float(val)
                except Exception as e:
                    msg = f"Could not convert value to float: {val}\n{str(e)}"
                    raise Exception(msg)

                if tag == "nreads":
                    nreads = val
                elif tag == "avg_read_length":
                    avg_len = val
                else:
                    raise Exception(f"Unrecognized tag: {tag}")

        else:
            output_handle.write(line)

    if nreads is None or avg_len is None:

        raise Exception("Did not find tags for nreads and avg_read_length")

    return int(nreads), avg_len


if __name__ == "__main__":

    # Get the output path from sys argv
    combine_mpl_alignments(
        output_file=sys.argv[1]
    )
