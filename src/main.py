#!/usr/bin/env python3

import logging
import logging.config
import multiprocessing
from time import perf_counter as time
from typing import TypeVar, List

from pysam import AlignedSegment

from .clusterer import Clusterer
from .consensus import Consensus
from .utils import LogMessages, EmptyClusterError, PickableRead, group_reads

ConsensusRead = TypeVar("ConsensusRead")

def main(bam: str, threads: int, threshold: int, window: int, debug: bool):
    """
    Takes the path to a UMI tagged BAM file, parses the reads
    and clusters them based on their UMI's similarity and genomic
    coordinates.
    """

    # ------------------ Logging ------------------
    logging.config.dictConfig(LogMessages.get_config(debug=debug))
    logger = logging.getLogger(__name__)

    logger.info(LogMessages.init_log(bam, ""))
    ti = time()

    # ------------------ START ------------------
    if min(threads, multiprocessing.cpu_count()) != threads:
        threads = multiprocessing.cpu_count()
        logger.warning(
            f"Invalid number of threads. Setting threads to {threads}."
        )

    uc = Clusterer(bam)
    bam_reads: List[List[AlignedSegment]] = uc.read_bam(threads=threads)

    if threads > 1:
        total_reads = sum([len(contig) for contig in bam_reads])
    else:
        total_reads = len(bam_reads)

    # ------------------ Clustering ------------------
    logger.info("Starting clustering...")
    t = time()

    # create an object to hold relevant data from each read and still be pickable
    if threads > 1:
        pk_reads = list()
        for contig in bam_reads:
            pk_reads.append([PickableRead(read) for read in contig])

    # cluster the reads
    if threads > 1:
        with multiprocessing.Pool(processes=threads) as pool:
            # Apply the process_element function to each element in the list
            pk_clustered_reads = pool.map(uc.cluster, pk_reads)
    else:
        clustered_reads: List[List[AlignedSegment]] = uc.cluster(bam_reads)

    # trace back the original reads if using threading
    if threads > 1:
        clustered_reads = group_reads(pk_clustered_reads, bam_reads)
        del pk_clustered_reads
        del pk_reads
        del bam_reads
    else:
        del bam_reads

    logger.info(f"Clustering completed in {(time() - t):2f}s.\n{'-' * 50}")
    logger.info(f"{len(clustered_reads)} clusters found out of {total_reads} reads.")

    # check integrity of the clustering
    def count_objects(nested_list, obj_type):
        count = 0
        if isinstance(nested_list, obj_type):
            count += 1
        elif isinstance(nested_list, list):
            for sublist in nested_list:
                count += count_objects(sublist, obj_type)
        return count

    logger.info("Checking integrity of the clustering...")

    try:
        assert total_reads == count_objects(clustered_reads, AlignedSegment)
    except AssertionError:
        logger.error("Some reads were duplicated during the clustering process.")
        logger.error(f"Total reads: {total_reads}")
        logger.error(f"Clustered reads: {count_objects(clustered_reads, AlignedSegment)}")
        raise AssertionError
    logger.info(f"Clustering integrity check passed.\n{'-' * 50}")

    # ------------------ Consensus sequences ------------------
    logger.info("Starting consensus sequence computing...")

    t = time()

    consensus_reads = list()
    errors = 0
    for reads in clustered_reads:
        try:
            cs = Consensus(reads)
        except EmptyClusterError:
            errors += 1
            continue
        consensus_reads.append(cs.compute_consensus())

    logger.info(f"Consensus sequences computed in {(time() - t):2f}s.\n{'-' * 50}")
    logger.info(f"{errors} clusters were empty.")

    # ------------------ Exporting ------------------
    logger.info("Exporting consensus sequences to STDOUT...")
    print(*consensus_reads, sep="\n")

    # ------------------ END ------------------
    logger.info(f"Execution competed in {(time() - ti):2f}s.")


if __name__ == "__main__":
    main()
