# DADEC

### Installation
Install via [bioconda](https://bioconda.github.io/):

### Compilation
If you however want to compile GraphAligner yourself, run these:


- `git clone https://github.com/YCMavis/DADEC.git`
- `cd DADEC`
- `git submodule update --init --recursive`
- `conda env create -f neta.yaml`
- `source activate DADEC`
- `make all`

### Running

Quickstart: `DADEC -s short_reads.fa -l long_reads.fa -t 16`

See Parameters, the option DADEC --help and the subsections below for more information and options

### Parameters
- `s` input short reads Format .fa, uncompressed or gzipped
- `l` input long reads Format .fa, uncompressed or gzipped
All parameters below are optional.

- `t` number of aligner threads. The program also uses two IO threads in addition to these.
- `o` output file name. Format .fa
- `S` number of splits for long-read files
- `r` Haplotype filtering threshold for the step two correction
- `k` k-mer size for the step one correction
- `K` k-mer size for the step three correction
- `a` Abundance threshold for the step one correction
- `A` Abundance threshold for the step three correction