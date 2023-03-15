import os

if __name__ == "__main__":
    bam = "test_files/forward_lib1_2838_KI_STR_nanopore.bam"
    fastq = "test_files/forward_lib1_2838_KI_STR.fastq.gz"
    target_regions = "test_files/targets.csv"
    outdir = "test_output/"

    cmd = f"python main.py {bam} {fastq} -t {target_regions} -o {outdir} -d"

    print("Starting test run in debug mode...")
    print(f"Commad: {cmd}", flush=True)

    os.system(cmd)
    os.system("clear")

    print("Test run completed. Check the log file for more information.")
    