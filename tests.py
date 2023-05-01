#!/usr/bin/env python3

import os
from time import perf_counter


def fast_test() -> None:
    """
    Run a quick test to check if the program is working.
    """
    cmd = "python umiclusterer.py test_files/sample.bam -j 4 -t 1 -w 5 --debug > output/test_sample.fastq"
    print("Starting short test run in debug mode...")
    print(f"Command: {cmd}", flush=True)

    start = perf_counter()
    os.system(cmd)
    end = perf_counter()

    print(f"Time elapsed: {end - start:.2f} seconds", flush=True)


def full_test() -> None:
    """
    Run a full test to check if the program is working.
    """
    bam = "test_files/full.bam"
    cmd = f"python umiclusterer.py {bam} -j 10 -t 1 -w 5 --debug | gzip > output/test_full.fastq.gz"

    print("Starting full test run in debug mode...")
    print(f"Command: {cmd}", flush=True)

    start = perf_counter()
    os.system(cmd)
    end = perf_counter()

    print(f"Time elapsed: {end - start:.2f} seconds", flush=True)


def profiling() -> None:
    """
    Run a profiling test to check bottlenecks.
    """
    ...


if __name__ == "__main__":
    fast_test()
    full_test()
    print("Test runs completed. Check the log file for more information.")
