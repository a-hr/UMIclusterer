#!/usr/bin/env python3

import os
from time import perf_counter

if __name__ == "__main__":
    bam = "test_files/mousev1_183a-WT-1a_Ctx.sortedByCoord.out.bam"
    target_regions = "test_files/targets.csv"
    target_regions_saf = "test_files/targets.saf"

    if not os.path.exists("output"):
        os.mkdir("output")

    cmd_1 = f"python umiclusterer.py {bam} -t {target_regions} | gzip > output/test_output.fastq.gz"
    cmd_2 = f"python umiclusterer.py {bam} -t {target_regions_saf} -s | gzip > output/test_output_saf.fastq.gz"

    print("Starting test run in debug mode...")
    print(f"Command 1: {cmd_1}", flush=True)
    print(f"Command 2: {cmd_2}", flush=True)

    start = perf_counter()
    os.system(cmd_1)
    os.system(cmd_2)
    os.system("rm -rf output/")
    # os.system("clear")

    end = perf_counter()
    print(f"Time elapsed: {end - start:.2f} seconds", flush=True)
    print("Test runs completed. Check the log file for more information.")
    