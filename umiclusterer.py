#!/usr/bin/env python3

import click

from src import UMIclusterer


@click.command()
@click.argument(
    "bam",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--threads",
    "-j",
    type=click.INT,
    help="Threads to use to speed up clustering process.",
    required=False,
    default=1,
)
@click.option(
    "--threshold",
    "-t",
    type=click.INT,
    help="Hamming distance threshold between UMIs.",
    required=False,
    default=1,
)
@click.option(
    "--window",
    "-w",
    type=click.INT,
    help="Allowed window around genomic start and end to consider two reads the same.",
    required=False,
    default=5,
)
@click.option(
    "--debug",
    "-d",
    is_flag=True,
    help="Enables debug mode.",
    required=False,
)
def main(bam, threads, threshold, window, debug):
    UMIclusterer(bam, threads, threshold, window, debug)


if __name__ == "__main__":
    main()
