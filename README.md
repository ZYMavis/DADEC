# DADEC

DADEC is a hybrid error correction tool specifically designed for long-read sequencing data. It combines multiple correction strategies to effectively reduce errors while preserving biological variations.

## Installation

### Install from source code

To compile DADEC from source, follow these steps:

#### Clone the repository
- `git clone https://github.com/ZYMavis/DADEC.git`
- `cd DADEC`

#### Initialize and update submodules
- `git submodule update --init --recursive`

#### Create and activate the conda environment
- `conda env create -f environment.yaml`
- `conda activate DADEC`

#### Build the executable
- `make all`

Note: The conda environment will install all required dependencies. If you prefer to manage dependencies manually, please ensure all libraries listed in environment.yaml are available.


## Usage

### Running

#### Quick start
Run DADEC with default parameters:

##### bash
 `DADEC -s short_reads.fa -l long_reads.fa -t 16`

-s : short reads file (FASTA format)

-l : long reads file (FASTA format)

-t : number of threads for alignment (default: 1)

Quickstart: `DADEC -s short_reads.fa -l long_reads.fa -t 16`

See Parameters, the option DADEC --help and the subsections below for more information and options

#### Demo

Test the basic functionality using the provided sample data:

- `cd dataDemo`
- `sh demo.sh`

The demo script will run DADEC on a small dataset and output the corrected reads.

### Parameters

- `s` input short reads Format .fa, uncompressed or gzipped
- `l` input long reads Format .fa, uncompressed or gzipped
All parameters below are optional.

- `t` number of aligner threads. The program also uses two IO threads in addition to these. (default: 1)
- `o` output file name. Format .fa (default: DADEC.fa)
- `S` number of splits for long-read files (default: 5)
- `r` Haplotype filtering threshold for the step two correction (default: 0.08)
- `k` k-mer size for the step one correction (default: 39)
- `K` k-mer size for the step three correction (default: 39)
- `a` Abundance threshold for the step one correction (default: 2)
- `A` Abundance threshold for the step three correction (default: 1)