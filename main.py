import logging
import logging.config
from time import time

from UMI_clusterer import UMI_clusterer
from utils import LogMessages

logging.config.dictConfig(
    LogMessages.get_config()
)
logger = logging.getLogger(__name__)


if __name__ == "__main__":
    # ------------ INPUTS ------------
    bam = "data/forward_lib1_2838_KI_STR_nanopore.bam"
    fastq = "data/forward_lib1_2838_KI_STR.fastq.gz"
    target_regions = [
        ("chr5", 35069758, 35069874),
        ("chr5", 34920043, 34920162),
        ("chr5", 34920505, 34920620),
    ]
    outdir = "output"

    # ------------ EXECUTION ------------
    logger.info(LogMessages.init_log(bam, target_regions, outdir))

    uc = UMI_clusterer(bam, target_regions, outdir)
    uc.read_bam()

    logger.info("Starting clustering...")
    t = time()

    all_clusters = uc.compute_clusters()
    
    logger.info(f"Clustering completed in {time() - t}s.\n{'-' * 50}")

    logger.info("Writing fastqs...")
    t = time()

    for i, clusters in enumerate(all_clusters, start=1):
        # use asyncio to write fastqs in parallel
        for cluster, reads in clusters.items():
            uc.gen_fastq(reads, fastq, f"target{i}_{cluster}.fastq.gz")

    logger.info(f"Fastqs written in {time() - t}s.\n{'-' * 50}\n")
    logger.info("Execution competed.")