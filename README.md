# DADEC

### Compilation
If you however want to compile DADEC yourself, run these:

- `git clone https://github.com/YCMavis/DADEC.git`
- `cd DADEC`
- `git submodule update --init --recursive`
- `conda env create -f environment.yaml`
- `conda activate DADEC`
- `make all`

### Running

Quickstart: `DADEC -s short_reads.fa -l long_reads.fa -t 16`

See Parameters, the option DADEC --help and the subsections below for more information and options

### Demo

Run basic functionality with sample data

- `cd dataDemo`
- `sh demo.sh`

### Parameters
- `s` input short reads Format .fa, uncompressed or gzipped
- `l` input long reads Format .fa, uncompressed or gzipped
All parameters below are optional.

- `t` number of aligner threads. The program also uses two IO threads in addition to these. (default: 1)
- `o` output file name. Format .fa (default: DADEC.fa)
- `S` number of splits for long-read files (default: DADEC.fa)
- `r` Haplotype filtering threshold for the step two correction (default: 5)
- `k` k-mer size for the step one correction (default: 39)
- `K` k-mer size for the step three correction (default: 39)
- `a` Abundance threshold for the step one correction (default: 2)
- `A` Abundance threshold for the step three correction (default: 1)