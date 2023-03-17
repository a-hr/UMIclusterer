import logging
import logging.config
from time import time

import click
import pandas as pd

from UMIClusterer import UMIClusterer
from utils import LogMessages


@click.command()
@click.argument(
    "bam",
    type=click.Path(exists=True),
    required=True,
)
@click.argument(
    "fastq",
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
    "--outdir",
    "-o",
    type=click.Path(exists=False),
    help="Path to the output directory.",
    default="output/",
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
def main(bam, fastq, target_regions, outdir, saf, debug):
    """
    Takes the path to a BAM file and its associated FASTQ, finds the reads in
    the specified targets and clusters them based on their umis' similarity.
    """

    logging.config.dictConfig(LogMessages.get_config(debug=debug))
    logger = logging.getLogger(__name__)
        
    logger.info(LogMessages.init_log(bam, target_regions, outdir))

    if not saf:
        target_regions = pd.read_csv(target_regions, sep=";")
    else:
        target_regions = pd.read_csv(target_regions, sep="\t")
    
    target_regions = [
        (row["chr"], row["start"], row["end"]) for _, row in target_regions.iterrows()
    ]

    uc = UMIClusterer(bam, target_regions, outdir)
    uc.read_bam()

    logger.info("Starting clustering...")
    t = time()

    all_clusters = uc.compute_clusters()

    logger.info(f"Clustering completed in {time() - t}s.\n{'-' * 50}")

    logger.info("Writing fastqs...")
    t = time()

    for i, clusters in enumerate(all_clusters, start=1):
        # todo: use asyncio to write fastqs in parallel
        for cluster, reads in clusters.items():
            uc.gen_fastq(reads, fastq, f"target{i}_{cluster}.fastq.gz")

    logger.info(f"Fastqs written in {time() - t}s.\n{'-' * 50}\n")
    logger.info("Execution competed.")


if __name__ == "__main__":
    main()
