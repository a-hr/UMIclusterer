import logging

from pathlib import Path
from typing import List, Union, TypeVar

import numpy as np

from pysam import AlignmentFile
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, hamming

logger = logging.getLogger(__name__)

AlignedSegment = TypeVar("AlignedSegment")


class Clusterer:
    def __init__(self, bam: str) -> None:
        self.bam = Path(bam)

        if not self.bam.exists():
            logger.error(f"File {self.bam} not found.")
            raise FileNotFoundError

    def read_bam(self, threads: int) -> Union[List[AlignedSegment], List[List[AlignedSegment]]]:
        """Reads the bam file and returns the reads for each target region.

        Returns:
            Dict[str, List[AlignedSegment]]: A dictionary with a list of reads for each target region.
        """
        logger.info("Reading bam file...")

        # Check if the bam file is paired-end or single-end
        with AlignmentFile(self.bam, "rb") as bam:
            sample = list(bam.fetch())[:1000]
            counts = [read.is_paired for read in sample if not read.is_unmapped]

            if sum(counts) > 0:
                logger.error(f"{self.bam} BAM file is paired-end, but only single-end is supported.")
                raise ValueError(f"{self.bam} BAM file is paired-end, but only single-end is supported.")

        # Read the bam file
        logger.info("Fetching reads...")
        if threads > 1:
            reads = self._split_bam(self.bam)
        else:
            with AlignmentFile(self.bam, "rb") as bam:
                reads = [read for read in bam.fetch() if not read.is_unmapped]
        logger.info("Bam file parsed.")
        
        if not any(reads):
            logger.error("No reads found.")
            raise ValueError("No reads found.")

        return reads

    def cluster(self, reads: List[AlignedSegment], threshold: int = 1, window: int = 5) -> List[List[AlignedSegment]]:
        """Clusters the reads by their UMI and coordinates.

        Args:
            reads (List[AlignedSegment]): list of reads to cluster.
            threshold (int, optional): maximum UMI distance to consider same cluster. Defaults to 1.
            window (int, optional): window size to consider for the coordinates. Defaults to 5.

        Returns:
            List[List[AlignedSegment]]: list of lists containing the reads for each cluster.
        """
        # Compute the distance matrix based on the UMI distance and genomic coordinates
        if not reads:
            logger.warning("No reads found. Skipping clustering.")
            return []
        if len(reads) == 1:
            logger.warning(f"Skipping clustering. Only one read found in contig {reads[0].reference_name}.")
            return reads
        
        distance_matrix = self._generate_matrix(reads, threshold, window)

        if not distance_matrix.any():
            logger.warning("Empty distance matrix. Skipping clustering.")
            return reads

        # Perform hierarchical clustering with complete linkage and threshold = max umi distance + allowed window
        linkage_matrix = linkage(distance_matrix, method='complete')
        clusters = fcluster(linkage_matrix, t=threshold + window, criterion='distance')

        if not clusters.any():
            logger.warning("No clusters found. Skipping clustering.")
            return []

        # Fetch the reads associated with each cluster
        grouped_reads = []
        for c in np.unique(clusters):
            assoc_reads: List = self._fetch_by_cluster_idx(reads, clusters, c)
            grouped_reads.append(assoc_reads)

            if not assoc_reads:
                logger.warning(f"No reads found for cluster {c}.")

        return grouped_reads

    def _generate_matrix(self, reads: List[AlignedSegment], threshold: int = 1, window: int = 5) -> np.ndarray:
        """Generates a similarity matrix for the reads in the target region.

        Args:
            reads (List[AlignedSegment]): list of reads in the target region.

        Returns:
            np.array: similarity matrix.
        """
        read_keys = [
            (
                read.query_name.split("_")[-1],
                read.reference_name,
                read.reference_start,
                read.reference_end
            ) for read in reads
        ]
        return pdist(read_keys, lambda x, y: self._metric(x, y, threshold, window))

    @staticmethod
    def _fetch_by_cluster_idx(reads: List, clusters: np.array, idx: int) -> List[AlignedSegment]:
        reads = np.array(reads)
        return reads[clusters == idx].tolist()

    @staticmethod
    def _metric(x: tuple, y: tuple, threshold: int = 1, window: int = 5) -> int:
        """Custom metric to compute the similarity between two reads.
        - If the UMI distance is greater than the threshold, the reads are not considered to be from the same origin.
        - The same applies to coordinate distance: if dist > window, the similarity score is 999.

        Args:
            x (tuple): tuple containing the UMI, chr, start and end coordinate of first read.
            y (tuple): tuple containing the UMI, chr, start and end coordinate of second read.
            threshold (int, optional): maximum UMI distance to consider same cluster. Defaults to 1.
            window (int, optional): window size to consider for the coordinates. Defaults to 5.

        Returns:
            int: similarity score.
        """
        umi_dist = hamming(list(x[0]), list(y[0])) * len(x[0])  # Hamming distance (0-1) * length of UMI

        if umi_dist > threshold:
            return 999

        if x[1] != y[1]:  # Different chromosome
            return 999

        cdist = (abs(int(x[2]) - int(y[2])) + abs(int(x[3]) - int(y[3])))/2
        coord_dist = cdist if cdist <= window else 999

        return umi_dist + coord_dist

    @staticmethod
    def _split_bam(bam: str) -> List[List[AlignedSegment]]:
        """Fetches the reads from each chromosome in the bam file to a separate list.
        This allows for multiprocessing without running the risk of splitting a cluster into multiple processes,
        thus artificially increasing the number of clusters."""

        with AlignmentFile(bam, "rb") as bamfile:
            contigs = bamfile.references
            logger.info(f"Found {len(contigs)} contigs.")
            logger.info(f"Found {bamfile.mapped} mapped reads in {bam}.")
            logger.info("Splitting bam file by contig...")

            bam_reads = []
            for contig in contigs:
                reads = [read for read in bamfile.fetch(contig=contig) if not read.is_unmapped]
                if reads:
                    bam_reads.append(reads)
            logger.info(f"Found {len(bam_reads)} mapped contigs.")                
        return bam_reads