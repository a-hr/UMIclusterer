# Introduction

UMIclusterer is a tool for clustering reads by UMI in Nanopore sequencing data. It uses the `pysam` library to parse `bam` files and `fastq` files. The clustering is performed using ``numpy`` and ``scipy``.

The aim of this project is to provide a tool for seamlessly parsing ``fastq`` and ``bam`` files, searching for target mapped regions and clustering reads into groups by UMI.
UMIs are grouped by their Humming distance (default 1) and the grouped reads are written to separate output `fastq` files. These files can then be feeded to consensus calling tools such as `SPOA`.

# Installation

The easiest way to install UMIclusterer is to simply clone this repo:

```bash
git clone https://github.com/a-hr/UMIclusterer.git
```

A docker image is also available on [Docker Hub](https://hub.docker.com/r/ahr1/umiclusterer). TODO: upload image

## Dependencies

UMIclusterer requires the following dependencies:

* Python 3.8 or higher
* pysam
* numpy
* scipy
* Levenshtein

# Usage

TODO: write usage