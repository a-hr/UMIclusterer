#!/usr/bin/env python3

import logging
import logging.config
from time import perf_counter as time
from typing import TypeVar

import click
import pandas as pd

from clusterer import Clusterer
from consensus import Consensus
from utils import LogMessages


ConsensusRead = TypeVar("ConsensusRead")


@click.command()
@click.argument(
    "bam",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--target_regions",
    "-t",
    type=click.Path(exists=True),
    help="Path to the file containing the target regions.",
    required=True,
)
@click.option(
    "--saf",
    "-s",
    is_flag=True,
    help="Target regions are in SAF format.",
    required=False,
)
@click.option(
    "--debug",
    "-d",
    is_flag=True,
    help="Enables debug mode.",
    required=False,
)
def main(bam, target_regions, saf, debug):
    """
    Takes the path to a BAM file and its associated FASTQ, finds the reads in
    the specified targets and clusters them based on their umis' similarity.
    """

    # ------------------ Logging ------------------
    logging.config.dictConfig(LogMessages.get_config(debug=debug))
    logger = logging.getLogger(__name__)

    logger.info(LogMessages.init_log(bam, target_regions))
    ti = time()
    # ------------------ Parsing Inputs ------------------
    if not saf:
        target_regions = pd.read_csv(target_regions, sep=";")
    else:
        target_regions = pd.read_csv(target_regions, sep="\t")

    target_regions = [
        (
            row["Chr"] if "chr" in str(row["Chr"]) else f"chr{row['Chr']}",
            row["Start"],
            row["End"],
        )
        for _, row in target_regions.iterrows()
    ]

    # ------------------ START ------------------
    uc = Clusterer(bam, target_regions)
    bam_reads = uc.read_bam()

    # ------------------ Clustering ------------------
    logger.info("Starting clustering...")
    t = time()

    clusters = dict()  # key: target_n, val: dict of {clusterN: [bam_reads]}

    # TODO: use multiprocessing to cluster in parallel
    for target_n, read_objs in bam_reads.items():
        logger.info(f"Clustering target {target_n + 1}/{len(bam_reads)}")
        UMIs = uc.get_UMIs(read_objs)

        clusters[target_n] = uc.compute_clusters(UMIs, read_objs)

    logger.info(f"Clustering completed in {(time() - t):2f}s.\n{'-' * 50}")

    # ------------------ Consensus sequences ------------------
    logger.info("Starting consensus sequence computing...")

    t = time()

    consensus_reads = list()

    for target_n, target_clusters in clusters.items():
        logger.info(
            f"Computing consensus sequences for target {target_n + 1}/{len(clusters)}"
        )

        for cluster_n, reads in enumerate(target_clusters.values()):
            logger.debug(
                f"Computing consensus for cluster {cluster_n + 1}/{len(target_clusters)}"
            )
            cs = Consensus(reads)
            consensus_reads.append(cs.compute_consensus())

    logger.info(
        f"{len(consensus_reads)} consensus sequences computed in {(time() - t):2f}s.\n{'-' * 50}"
    )

    # ------------------ Exporting ------------------
    logger.info("Exporting consensus sequences to STDOUT...")

    for read in consensus_reads:
        print(read)

    # ------------------ END ------------------
    logger.info(f"Execution competed in {(time() - ti):2f}s.")


if __name__ == "__main__":
    main()
