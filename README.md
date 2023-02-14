# atac-barcode-matcher

For checking the structure of MAT-seq reads. Quantifies the number of reads with:
- adapters
- spacers
- adapters and spacers
- adapters, spacers, and valid barcodes
Does so for both forward and reverse reads.

## Usage:
```
./atac-barcode-matcher input.fq.gz output.fq.gz barcodes.csv
```
- `input.fq.gz` is the input fastq file
- `output.fq.gz` is currently unused, will be used to collect reads with valid structure
- `barcodes.csv` is a newline-separated list of barcodes to search for

## Installation:
```
git clone https://github.com/jakob-schuster/atac-barcode-matcher.git
cd atac-barcode-matcher
make
```
