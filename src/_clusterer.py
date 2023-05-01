import logging
from pathlib import Path
from typing import List, Dict

# import Levenshtein
import numpy as np
from pysam import AlignmentFile, AlignedSegment
from scipy.cluster import hierarchy

logger = logging.getLogger(__name__)


class Clusterer:
    def __init__(self, bam: str, regions: List) -> None:
        self.bam = Path(bam)

        self.regions = regions  # list of tuples (chr, start, end)

        self.UMIs = list()  # will hold an array of umis for each region
        self.reads = list()  # will hold an array of reads for each region

        if not self.bam.exists():
            logger.error(f"File {self.bam} not found.")
            quit(1)

    def read_bam(self) -> Dict[str, List[AlignedSegment]]:
        """Reads the bam file and returns the reads for each target region.

        Returns:
            Dict[str, List[AlignedSegment]]: A dictionary with a list of reads for each target region.
        """
        reads = dict()

        logger.info("Reading bam file...")

        # Check if the bam file is paired-end or single-end
        with AlignmentFile(self.bam, "rb") as bam:
            sample = list(bam.fetch())[:1000]
            counts = [read.is_paired for read in sample if not read.is_unmapped]

            if sum(counts) > 0:
                logger.error(f"{self.bam} BAM file is paired-end, but only single-end is supported.")
                quit(1)
        
        # Read the bam file
        with AlignmentFile(self.bam, "rb") as bamfile:
            for i, region in enumerate(self.regions):
                logger.info(f"Fetching reads from {region[0]}:{region[1]}-{region[2]}")
                reads[i] = [
                    read for read in bamfile.fetch(*region) if not read.is_unmapped
                ]

        logger.info("Bam file parsed.")

        if not reads[0]:
            logger.error("No reads found.")
            quit(1)

        return reads

    def compute_clusters(self, UMIs: np.array, reads: List, threshold: int = 1) -> Dict:
        """Computes the UMI humming-distance-based clusters for a single target region.

        Args:
            UMIs (np.array[str]): numpy array containing the UMIs for each read.
            reads (List[AlignedSegment]): numpy array containing the AlignedSegment objects.
            threshold (int, optional): humming-distance threshold to consider same origin. Defaults to 1.

        Returns:
            dict[str, list]: dict containing the AlignedSegments for each cluster.
        """
        clusters: Dict[str, List[AlignedSegment]] = dict()  # key: cluster id

        # create clusters of read indexes based on their UMI similarity
        _clusters = self._run_clustering(UMIs, threshold)

        # fetch the reads associated with each UMI by index
        for c in np.unique(_clusters):
            assoc_reads: List = self._fetch_by_cluster_idx(reads, _clusters, c)
            clusters[f"cluster_{c}"] = assoc_reads

            if not assoc_reads:
                logger.warning(f"No reads found for cluster {c}.")

        return clusters

    @staticmethod
    def get_UMIs(reads: List[AlignedSegment]) -> np.array:
        return np.array([read.query_name.split("_")[-1] for read in reads])

    @staticmethod
    def get_read_ids(reads: List[AlignedSegment]) -> np.array:
        return np.array([read.query_name for read in reads])

    # @staticmethod
    # def _get_distances(umis: np.array) -> np.array:
    #     return np.array(
    #         [
    #             Levenshtein.hamming(umis[i], umis[j])
    #             for i in range(len(umis))
    #             for j in range(len(umis))
    #         ]
    #     ).reshape(len(umis), len(umis))

    @staticmethod
    def _fetch_by_cluster_idx(reads: List, clusters: np.array, idx: int) -> List:
        reads = np.array(reads)
        return reads[clusters == idx].tolist()

    def _run_clustering(self, umis: np.array, threshold: int = 1) -> np.array:
        # Compute the Levenshtein distance matrix
        dist_matrix = self._get_distances(umis)
        if not dist_matrix.any():
            logger.warning("No distances computed. Skipping clustering.")
            return None

        # Convert distance matrix to condensed form
        distance_condensed = hierarchy.distance.squareform(dist_matrix)

        # Perform hierarchical clustering with complete linkage
        z = hierarchy.linkage(distance_condensed, method="complete")

        # Extract clusters from the dendrogram
        clusters = hierarchy.fcluster(z, t=threshold, criterion="distance")
        if not clusters.any():
            logger.warning("No clusters found. Skipping clustering.")
            return None

        return clusters
