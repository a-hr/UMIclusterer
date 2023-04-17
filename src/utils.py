import os
import sys
from typing import List, Tuple, TypeVar, Union

import pysam

AlignedSegment = TypeVar("AlignedSegment")


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
        self.quality: List[Union[str, int]] = list(read.query_qualities)  # no 33 offset
        self.id: str = read.query_name
        self.cigar: Tuple[Tuple[int, int]] = read.cigartuples

        self.pad_seq: str = ""
        self.pad_qual: List[Union[str, int]] = []

        self.__ascii = False
        self.__seq_padded = False
        self.__qual_padded = False

    def seq_pad(self) -> str:
        """
        Takes a read and its cigar and formats the read accordingly:
            - Matches and mismatches are left as is.
            - Insertions are converted to lowecase.
            - Deletions are padded with 'p'.
        """
        if not self.__seq_padded:
            self.__seq_padded = True
            self.pad_seq = ""
            pos = 0

            for op, length in self.cigar:
                if op == 0:  # match or mismatch
                    self.pad_seq += self.seq[pos : pos + length]
                    pos += length
                elif op == 1:  # insertion
                    self.pad_seq += self.seq[pos : pos + length].lower()
                    pos += length
                elif op == 2:  # deletion
                    self.pad_seq += "p" * length
                elif op == 3:  # skipped region
                    continue
                elif op == 4:  # soft clipping
                    continue
                elif op == 5:  # hard clipping
                    continue
            return self.pad_seq

        return self.pad_seq

    def qual_pad(self) -> List[Union[str, int]]:
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
            self.pad_qual = []

            for char in self.pad_seq:
                if char == "p":
                    self.pad_qual.append("p")
                else:
                    self.pad_qual.append(self.quality[pos])
                    pos += 1
            return self.pad_qual

        return self.pad_qual

    def to_phred(self) -> List[int]:
        """Converts an ASCII-encoded quality string to a list of Phred quality scores.
        - Note: this class works withouth the 33 offset. It is only added when exporting.
        """
        if not self.__ascii:
            return self.quality
        else:
            self.__ascii = False
            self.quality = [char if char == "p" else ord(char) for char in self.quality]

        return self.quality

    def to_ascii(self) -> List[str]:
        """Converts a list of Phred quality scores to an ASCII-encoded quality string.
        - Note: this class works withouth the 33 offset. It is only added when exporting.
        """
        if self.__ascii:
            return self.quality
        else:
            self.__ascii = True
            self.quality = [
                score if score == "p" else chr(score) for score in self.quality
            ]

        return self.quality


class ConsensusRead:
    def __init__(self, seq: str, qual: List[int], _id: str) -> None:
        self.seq: str = seq
        self.qual: str = self.__unphred(qual)
        self.id: str = _id

    @staticmethod
    def __unphred(qual: List[int]) -> str:
        return "".join([chr(score + 33) for score in qual])

    def __str__(self):
        return f"@{self.id}\n{self.seq}\n+\n{self.qual}"
