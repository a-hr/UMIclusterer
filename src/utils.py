import os
import sys
import logging
from typing import List, Tuple, TypeVar

import pysam

AlignedSegment = TypeVar("AlignedSegment")
logger = logging.getLogger(__name__)


class LogMessages:
    @staticmethod
    def init_log(bam: str, regions: list):
        r = f"Initializing UMIclusterer.\n"
        r += f"Python version: {sys.version}\n"
        r += f"Pysam version: {pysam.__version__}\n\n"
        r += f"Input bam: {bam}\n"
        r += f"Target regions: {regions}\n"
        r += f"{'-' * 50}\n"
        return r

    @staticmethod
    def get_config(debug: bool):
        return {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {
                "default": {
                    "format": "%(asctime)s - %(levelname)s(%(name)s): %(message)s",
                    "datefmt": "%d/%m/%Y %H:%M:%S",
                }
            },
            "handlers": {
                "file": {
                    "class": "logging.FileHandler",
                    "filename": os.path.join(os.getcwd(), "UMIClusterer.log"),
                    "mode": "w",
                    "encoding": "utf-8",
                    "formatter": "default",
                }
            },
            "loggers": {
                "": {
                    "level": "DEBUG" if debug else "INFO",
                    "handlers": ["file"],
                }
            },
        }


class EmptyClusterError(Exception):
    pass


class CustomAlignedSegment:
    def __init__(self, read: AlignedSegment):
        if not isinstance(read, pysam.libcalignedsegment.AlignedSegment):
            raise EmptyClusterError("Read is not an AlignedSegment object.")
        
        self.seq: str = read.query_sequence
        self.int_qual: List[int] = list(read.query_qualities)  # Q scores as integers, not ASCII values (base 33)

        self.id: str = read.query_name
        self.cigar: Tuple[Tuple[int, int]] = read.cigartuples

        self.__seq_padded = False
        self.__qual_padded = False

    def seq_pad(self) -> None:
        """
        Takes a read and its cigar and formats the read accordingly:
            - Matches and mismatches are left as is.
            - Insertions are converted to lowecase.
            - Deletions are padded with 'p'.
        """
        if not self.__seq_padded:
            self.__seq_padded = True
            pos = 0
            _seq = ""

            for op, length in self.cigar:
                if op == 0:  # match or mismatch
                    _seq += self.seq[pos : pos + length]
                    pos += length
                elif op == 1:  # insertion
                    _seq += self.seq[pos : pos + length].lower()
                    pos += length
                elif op == 2:  # deletion
                    _seq += "p" * length
                elif op == 3:  # skipped region
                    continue
                elif op == 4:  # soft clipping
                    continue
                elif op == 5:  # hard clipping
                    continue
            self.seq = _seq

    def qual_pad(self) -> None:
        """
        Takes a quality string and a read and pads the quality string accordingly:
            - Matches and mismatches are left as is.
            - Insertions left as is.
            - Deletions are padded with 'p'.
        """
        if not self.__seq_padded:
            raise Exception("Sequence must be padded before quality.")

        if not self.__qual_padded:
            self.__qual_padded = True
            pos = 0
            _int_qual = []

            for char in self.seq:
                if char == "p":
                    _int_qual.append("p")
                else:
                    _int_qual.append(self.int_qual[pos])
                    pos += 1
            self.int_qual = _int_qual

    def __str__(self):
        return f"{self.id=} {self.seq=} {self.int_qual=}"


class ConsensusRead:
    def __init__(self, seq: str, qual: List[int], _id: str) -> None:
        self.seq: str = seq
        self.q_score: str = "".join([str(q) for q in qual])
        self.ascii_qual: str = "".join([chr(score + 33) for score in qual])
        self.id: str = _id

    def __str__(self):
        return f"@{self.id}\n{self.seq}\n+{self.q_score}\n{self.ascii_qual}"


class PickableRead:
    def __init__(self, read: AlignedSegment) -> None:
        self.query_name: str = read.query_name
        self.reference_name: str = read.reference_name
        self.reference_start: int = read.reference_start
        self.reference_end: int = read.reference_end
        self._id = (self.query_name, self.reference_name, self.reference_start, self.reference_end)


def group_reads(ordered_reads: List[List[PickableRead]], original_reads: List[AlignedSegment]) -> List[List[AlignedSegment]]:
    """
    Groups the original reads into the same groups as the ordered reads, using their query_name as a key.
    It also flattens the list of lists of lists of reads into a list of lists of reads, removing the 
    groping by contig, that is no longer needed.
    """
    # flatten the first level of ordererd lists
    _reads = []
    for contig_reads in ordered_reads:
        _reads.extend(contig_reads)
    ordered_reads = _reads
    del _reads

    # flatten the first level of original lists
    _reads = []
    for contig_reads in original_reads:
        _reads.extend(contig_reads)
    original_reads = _reads

    #####################################
    infer_id = lambda read: (read.query_name, read.reference_name, read.reference_start, read.reference_end)
    grouped_reads = []
    for cluster in ordered_reads:

        # get the ids of the reads in the cluster
        try:
            cluster_read_ids = {read.query_name for read in cluster}
        except TypeError:
            # cluster contains a single read
            cluster_read_ids = {cluster.query_name}
            cluster = [cluster]

        # create a new cluster with matching query_names (multimappers will be duplicated)
        new_cluster = [read for read in original_reads if read.query_name in cluster_read_ids]

        # check the presence of multimappers, and ask the object for more info to resolve them
        if isinstance(cluster, list) and len(cluster) != len(new_cluster):
            # compare with ids (slower but definitive
            cluster_read_ids = [read._id for read in cluster]
            new_cluster = [read for read in new_cluster if infer_id(read) in cluster_read_ids]

        grouped_reads.append(new_cluster)

    return grouped_reads
