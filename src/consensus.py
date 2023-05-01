import logging
from typing import List, Tuple, TypeVar

import numpy as np

from .utils import CustomAlignedSegment, ConsensusRead, EmptyClusterError


logger = logging.getLogger(__name__)
AlignedSegment = TypeVar("AlignedSegment")


class Consensus:
    def __init__(self, reads: List[AlignedSegment]) -> None:
        # convert the reads to CustomAlignedSegment objects and sort them by their length
        self.reads = [CustomAlignedSegment(read) for read in reads]

    def compute_consensus(self) -> ConsensusRead:
        """Computes the consensus read for a cluster.

        Returns:
            CustomAlignedSegment: An object containing the consensus sequence, quality string and read_id.
        """
        if len(self.reads) == 1:
            return ConsensusRead(self.reads[0].seq, self.reads[0].int_qual, self.reads[0].id)

        # pad the sequences to account for indels
        [read.seq_pad() for read in self.reads]

        # sort the reads by their sequence length
        self.reads = sorted(self.reads, key=lambda x: len(x.seq), reverse=True)
        try:
            seq, qual = self._consensus()
        except Exception as e:
            logger.error(f"{e}\n\t{[str(read) for read in self.reads]}")
            raise Exception(e)

        return ConsensusRead(seq, qual, self.reads[0].id)

    def _consensus(self) -> Tuple[str, List[int]]:
        """
        Takes a list of AlignedSegment objects and computes their consensus sequence.

        Returns:
            Tuple[str, List[int]]: A tuple containing the consensus sequence and quality string.
        """

        # Align the reads with the longest read using the Needleman-Wunsch algorithm to account for insertions
        self.reads = self._align_strings(self.reads)

        # Format the quality strings to match their sequences
        [read.qual_pad() for read in self.reads]

        # Find the consensus sequence
        consensus = ""
        quality = []

        for bases, qualities in zip(zip(*[r.seq for r in self.reads]), zip(*[r.int_qual for r in self.reads])):
            # for each A, T, G, C, and p in the column, count the number of times it appears
            counts = {base: bases.count(base) / len(bases) for base in set(bases)}
            n_scores = {key: value * 10 for key, value in counts.items()}

            # for each A, T, G, C, and p in the column, get the mean quality score for this position
            q_scores = {base: 0 for base in set(bases)}
            for i, base in enumerate(bases):
                if base != "p":
                    q_scores[base] += qualities[i]

            get_q_score = lambda base: q_scores[base] / (counts[base] * len(bases))
            qual_per_base: dict = {b: get_q_score(b) for b in set(bases) if b != "p" and counts[b] != 0}

            # !penalty: in case of tie, real base is taken. remove penalty if needed
            try:
                # deletions quality is the average quality of the other bases minos a penalty
                qual_per_base["p"] = sum(qual_per_base.values()) / len(qual_per_base) - 5
            except ZeroDivisionError:
                # all bases are deletions, skip this position
                continue

            # assign a score to each base based on its quality
            q_scores = self._assign_q_score(qual_per_base)

            # combine the scores from the quantity and quality
            scores = {base: 0.5 * n_scores[base] + 0.5 * q_scores[base] for base in n_scores.keys()}

            # pick the base with the highest score
            base = max(scores, key=scores.get)
            consensus += base if base != "p" else ""

            if base != "p":
                # asign the quality of the selected base to the consensus
                quality.append(int(qual_per_base[base]))

        if len(consensus) != len(quality):
            logger.error("Consensus and quality strings are not the same length.")
            raise ValueError("Consensus and quality strings are not the same length.")
        
        return consensus, quality

    @staticmethod
    def _align_strings(reads: List[CustomAlignedSegment]) -> List[CustomAlignedSegment]:
        # Set the gap penalty and the substitution matrix
        gap_penalty = -1
        substitution_matrix = np.zeros((256, 256), dtype=int)
        # 'p' is a special character that represents a deletion and is treated as wildcard
        for i in range(256):
            substitution_matrix[i, i] = 1
            substitution_matrix[i, ord("p")] = 1
            substitution_matrix[ord("p"), i] = 1

        # Find the length of the longest string
        max_length = max(len(r.seq) for r in reads)

        # Align each string with the longest string using the Needleman-Wunsch algorithm
        for i in range(len(reads)):
            # Initialize the scoring matrix and the traceback matrix
            rows = len(reads[i].seq) + 1
            cols = max_length + 1
            score_matrix = np.zeros((rows, cols), dtype=int)
            traceback_matrix = np.zeros((rows, cols), dtype=int)

            # Fill in the scoring matrix
            for r in range(1, rows):
                for c in range(1, cols):
                    if reads[i].seq[r - 1] == "p" or reads[0].seq[c - 1] == "p":
                        match_score = 1
                    else:
                        match_score = substitution_matrix[
                            ord(reads[i].seq[r - 1]), ord(reads[0].seq[c - 1])
                        ]
                    diagonal_score = score_matrix[r - 1, c - 1] + match_score
                    up_score = score_matrix[r - 1, c] + gap_penalty
                    left_score = score_matrix[r, c - 1] + gap_penalty
                    score = max(diagonal_score, up_score, left_score)
                    score_matrix[r, c] = score

                    if score == diagonal_score:
                        traceback_matrix[r, c] = 1
                    elif score == up_score:
                        traceback_matrix[r, c] = 2
                    else:
                        traceback_matrix[r, c] = 3

            # Trace back the alignment from the bottom-right corner of the matrix
            aligned_i = ""
            aligned_j = ""
            r = rows - 1
            c = cols - 1

            while r > 0 or c > 0:
                if traceback_matrix[r, c] == 1:
                    aligned_i = reads[i].seq[r - 1] + aligned_i
                    aligned_j = reads[0].seq[c - 1] + aligned_j
                    r -= 1
                    c -= 1
                elif traceback_matrix[r, c] == 2:
                    aligned_i = reads[i].seq[r - 1] + aligned_i
                    aligned_j = "p" + aligned_j
                    r -= 1
                else:
                    aligned_i = "p" + aligned_i
                    aligned_j = reads[0].seq[c - 1] + aligned_j
                    c -= 1

            # Add the aligned string to the output
            reads[i].seq = aligned_i

        return reads

    @staticmethod
    def _assign_q_score(q_score: dict) -> dict:
        r = {}
        for key, value in q_score.items():
            if value >= 30:
                r[key] = 8
            elif value >= 20:
                r[key] = 6
            elif value >= 15:
                r[key] = 4
            else:
                r[key] = 2
        return r
