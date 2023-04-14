# Introduction

UMIclusterer is a tool for clustering reads by UMI similarity in Nanopore single-end sequencing data. It uses the `pysam` library to parse `bam` files, clusters the reads using ``numpy`` and ``scipy`` and finally outputs a new *consensus FASTQ* to the STDOUT.

The aim of this project is to provide a tool for seamlessly deduplicating and correcting Nanopore sequences, searching for target mapped regions and clustering reads into groups by UMI.
UMIs are grouped by their Humming distance (default 1) and the grouped reads are converted to  `FASTQ` format. These files can then be realigned with standard aligning tools such as [STAR](https://github.com/alexdobin/STAR).

# Installation

The easiest way to install UMIclusterer is to simply clone this repo and install the dependencies in a virtual environment.
Using conda:

```bash
# clone the repo
git clone https://github.com/a-hr/UMIclusterer.git
# install dependencies in a virtual environment
conda create -n umiclusterer python=3.8 pysam numpy pandas scipy Levenshtein click
source ~/.bashrc
# add the repo to your path
echo "export PATH=$PATH:/path/to/UMIclusterer" >> ~/.bashrc
```

A docker image is also available on [Docker Hub](https://hub.docker.com/r/ahr1/umiclusterer). TODO: upload image

## Dependencies

UMIclusterer requires the following dependencies:

* Python 3.8 or higher
* pysam
* pandas
* numpy
* scipy
* Levenshtein
* click

# Usage

To use UMIclusterer, run the script with the path to a BAM file as well as the path to the file containing the target regions:

```bash
python umiclusterer.py [path-to-bam-file] --target_regions [path-to-target-regions-file] | gzip > [output-file].fastq.gz
```

> Note that the BAM file must be sorted by read coordinates and indexed, and only contain single-end reads. 

The target regions file must be a table-like file with the following columns:

* `Chr`: The chromosome number. It can either have the format `chr1` or `1`. `<str>|<int>`
* `Start`: The start position of the target region. `<int>`
* `End`: The end position of the target region. `<int>`

It can contain more columns, that will be ignored. Additionally, the file can be either in **tsv** or **csv** (``;``) format. The default expected format is the **csv**. Optionally, you can include the following flags:

--saf or -s: Indicates that the target regions are in **tsv** format. The file used in featureCounts to generate the target regions file is in this *SAF* format.
--debug or -d: Enables debug mode, which includes additional information in the log file.

# Output

UMIclusterer outputs a ``FASTQ`` file to the STDOUT containing either the consensus read or the original read, depending on whether a cluster was found or not. All the mapped reads in the BAM file are processed and passed to the output.

# How it works

1. For each target region specified in the input file, UMIclusterer will extract all the reads that map to that region.
2. Then, a Humming-distance matrix is calculated for all the reads in the region.
3. The reads are then clustered using hierarchical clustering, with a threshold of 1 Humming distance to be considered in the same cluster.
4. The clusters are then passed to the consensus module, that aligns the reads in the cluster between them, and analyses the sequences a a per-base basis, giving each possible base (A, C, G, T or None) a score based on their occurrence rate in the cluster and their basecall quality.

# Credits
UMIclusterer was developed by Álvaro Herrero as part of his final project carried out at the Biodonostia HRI.

# License
UMIclusterer is licensed under the GPLv3 license. See the LICENSE file for more information.