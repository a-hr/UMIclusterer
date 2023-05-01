# Introduction

UMIclusterer is a tool for clustering reads by UMI similarity in Nanopore single-end sequencing data. It uses the `pysam` library to parse `bam` files, clusters the reads using ``numpy`` and ``scipy`` and finally outputs a new *consensus FASTQ* to the STDOUT.

The aim of this project is to provide a tool for seamlessly deduplicating and correcting the error-prone Nanopore sequences, clustering reads into groups by UMI and genomic coordinates, and then creating a consensus read after evaluating the abundance and quality of each possible base.

After generating the consensu file, it can then be realigned with standard aligning tools such as [STAR](https://github.com/alexdobin/STAR).

# Why UMIclusterer?

UMIclusterer was developed to perform consensus-basd error-correction on Nanopore single-end sequencing data, which means having to deal with the following issues:

* **Error-prone reads**: Nanopore reads are known to be error-prone, with a high rate of indels and substitutions.
* **UMIs with sequencing errors**: UMIs are short sequences that are added to the 5' end of the reads to identify PCR duplicates. However, they are also prone to sequencing errors,a nd can lead to false positives when clustering reads.
* **Coordinate discordance**: due to the high error rate of Nanopore reads, (specially in low-complexity regions), reads from the same cluster can have different start and end coordinates, which can lead to false negatives when clustering reads. 

Other deduplication tools, such as UMI-tools or Gencore, either:
* select the highest quality duplicate, discarding the rest without correcting the errors
* only allow paired-end reads
* don't allow for coordinate discordance

UMIclusterer aims to solve these issues by clustering reads by UMI and genomic coordinates, and then creating a consensus read after evaluating the abundance and quality of each possible base.

# Installation

The easiest way to install UMIclusterer is to simply clone this repo and install the dependencies in a virtual environment.
Using conda:

```bash
# clone the repo
git clone https://github.com/a-hr/UMIclusterer.git
# install dependencies in a virtual environment
conda create -n umiclusterer -c anaconda -c conda-forge -c bioconda python=3.8 pysam numpy scipy click
source ~/.bashrc
# add the repo to your path
echo "export PATH=$PATH:/path/to/UMIclusterer" >> ~/.bashrc
```

A docker image is also available on [Docker Hub](https://hub.docker.com/r/ahr1/umiclusterer).

## Dependencies

UMIclusterer requires the following dependencies:

* Python 3.8 or higher
* pysam
* numpy
* scipy
* click

# Usage

To use UMIclusterer, simply run the script with the path to a BAM file:

```bash
python umiclusterer.py [path-to-bam-file] [--options] | gzip > [output-file].fastq.gz
```

> Note that the BAM file must be sorted by read coordinates and indexed, and only contain single-end reads. 

## Input

UMIclusterer requires a BAM file as input. The BAM file must be sorted by read coordinates and indexed, and only contain single-end reads. It can handle multimapping reads without any issues.

## Options

The following options are available:
* **-j | --threads**: Number of threads to use. Defaults to 1.
* **-t | --threshold**: Threshold for the Hamming distance between UMIs. Defaults to 1.
* **-w | --window**: Window size for the genomic coordinates. Creates a *safe-zone* of *-w [INT]* bases around both start and end coordinates, inside of which reads are considered to be elegible to be clustered. Defaults to 5, change based on the expected size of your reads.
* **-d | --debug**: Enables debug mode, which includes additional information in the log file.

## Output

UMIclusterer outputs a ``FASTQ`` file to the STDOUT containing either the consensus read or the original read, depending on whether a cluster was found or not. All the mapped reads in the BAM file are processed and passed to the output.

# How it works

1. The BAM file is parsed and the reads are grouped by genomic coordinates. This enables the use of multiprocessing without having to worry about reads from the same cluster being processed in different threads.
2. Then, the paiwise distance between all the reads in the contig is calculated, where the distance is the sum of the Hamming distance between the UMIs and the genomic distance between the reads.
3. The reads are then grouped using hierarchical clustering, based on the threshold and window given by the user.
4. The clusters are then passed to the consensus module, that aligns the reads in the cluster between them (global alignment using Needleman-Wunsch), and analyses the sequences on a per-base basis, giving each possible base (A, C, G, T or None) a score based on their abundance rate in the cluster and their basecall quality.

# TODO

* Add support for specific target regions, speeding up the process
* Translate bottleneck functions to C

# Credits
UMIclusterer was developed by Álvaro Herrero as part of his final thesis carried out at the Biodonostia HRI.

# License
UMIclusterer is licensed under the GPLv3 license. See the LICENSE file for more information.