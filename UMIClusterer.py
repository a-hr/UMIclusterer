import gzip
import logging
from pathlib import Path

import Levenshtein
import numpy as np
from pysam import AlignmentFile, FastxFile
from scipy.cluster import hierarchy

logger = logging.getLogger(__name__)


class UMIClusterer:
    def __init__(self, bam: str, regions: list, outdir: str) -> None:
        self.bam = Path(bam)
        self.outdir = Path(outdir)

        self.regions = regions  # list of tuples (chr, start, end)

        self.UMIs = list()  # will hold an array of umis for each region
        self.reads = list()  # will hold an array of reads for each region

        if not self.bam.exists():
            logger.error(f"File {self.bam} not found.")
            quit(1)

    def read_bam(self) -> None:
        self.reads = list()
        self.UMIs = list()

        logger.info("Reading bam file...")
        with AlignmentFile(self.bam, "rb") as bamfile:
            for region in self.regions:
                logger.info(
                    f"Fetching reads and umis from {region[0]}:{region[1]}-{region[2]}"
                )
                self.reads.append(
                    np.array([read.query_name for read in bamfile.fetch(*region)])
                )
                self.UMIs.append(
                    np.array(
                        [
                            read.query_name.split("_")[-1]
                            for read in bamfile.fetch(*region)
                        ]
                    )
                )
        logger.info("Bam file parsed.")

        if not self.UMIs and not self.reads:
            logger.error("No umis or reads found.")
            quit(1)

    def compute_clusters(self, threshold: int = 1) -> list:
        if not self.UMIs and not self.reads:
            raise ValueError("umis and reads not initialized. Run read_bam() first.")

        clusters_per_region = list()

        # TODO: run _run_clustering for each region in parallel
        target = 1
        for umis, reads in zip(self.UMIs, self.reads):
            logger.info(f"Clustering region {target}/{len(self.regions)}")
            region_cluster_reads = (
                {}
            )  # key: cluster index, value: list of associated reads
            _clusters = self._run_clustering(
                umis, threshold
            )  # all the clusters on this region

            # for each cluster, fetch the associated reads and store them
            for c in np.unique(_clusters):
                assoc_reads = self._fetch_by_cluster_idx(reads, _clusters, c)
                region_cluster_reads[f"cluster_{c}"] = assoc_reads

                if not assoc_reads:
                    logger.warning(f"No reads found for cluster {c}.")

            clusters_per_region.append(region_cluster_reads)
            target += 1

        return clusters_per_region

    def gen_fastq(self, reads: list, in_f: str, out_f: str) -> None:
        if not self.outdir.exists():
            logger.info(f"Creating output directory {str(self.outdir)}")
            self.outdir.mkdir(parents=True, exist_ok=True)

        out_f = self.outdir / out_f
        logger.debug(f"Generating {out_f}...")

        with FastxFile(in_f, persist=False) as fq, gzip.open(out_f, "wb") as of:
            for read in fq:
                if read.name in reads:
                    _r = bytes(str(read) + "\n", "utf-8")
                    of.write(_r)
        logger.debug(f"{out_f} generated.")

    @staticmethod
    def _get_distances(umis: np.array) -> np.array:
        return np.array(
            [
                Levenshtein.hamming(umis[i], umis[j])
                for i in range(len(umis))
                for j in range(len(umis))
            ]
        ).reshape(len(umis), len(umis))

    @staticmethod
    def _fetch_by_cluster_idx(umis, clusters, idx):
        return umis[clusters == idx].tolist()

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

    def _plot_clusters(self):
        # import matplotlib.pyplot as plt

        # plt.figure(figsize=(25, 10))
        # plt.title("Hierarchical Clustering Dendrogram")
        # plt.xlabel("sample index")
        # plt.ylabel("distance")
        # hierarchy.dendrogram(
        #     Z,
        #     leaf_rotation=90.0,  # rotates the x-axis labels
        #     leaf_font_size=8.0,  # font size for the x-axis labels
        # )
        # plt.savefig("dendrogram.png")
        pass
