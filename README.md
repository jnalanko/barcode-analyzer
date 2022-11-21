# Barcode analyzer

A fast tool to search for barcode sequences inside fasta/fastq data. Can also remove the reads that have a barcode. Internally, uses the [Aho-Corasick](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm) multiple string matching algorithm implementation from [here](https://github.com/cjgdev/aho_corasick).

## Compiling

Compile by running `make`. Compiling requires a C++ compiler with support for the C++17 standard. Tested to work with g++ version 9. 

## Quick start

There is some example data provided in the repository. Running barcode analysis on reads at `example_data/reads.fastq` with barcodes at `example_data/barcodes.txt`:

```
./barcode_analyzer analyze -i example_data/reads.fastq -b example_data/barcodes.txt
```

Filtering reads with barcodes with the same inputs, writing output to example_data/filtered.fastq:

```
./barcode_analyzer filter -i example_data/reads.fastq -b example_data/barcodes.txt -o example_data/filtered.fastq
```

## Usage

There are two commands:

```
Available commands: 
   ./barcode_analyzer analyze
   ./barcode_analyzer filter
Running a command without arguments prints the usage instructions for the command.
```

### Analyze

```
Search for barcode sequences inside a fasta/fastq file.
Usage:
  analyze [OPTION...]

  -i arg         The sequence file in fasta or fastq format.
  -o arg         Output file. If not given, prints to stdout.
  -b arg         A file containing the barcodes, one per line. Do not give 
                 reverse complements.
  -v, --verbose  Verbose output.
  -h, --help     Print usage
```

### Filter

```
Remove all sequences from a fasta/fastq file that have a barcode sequence.
Usage:
  filter [OPTION...]

  -i arg      The sequence file in fasta or fastq format.
  -o arg      Output file.
  -b arg      A file containing the barcodes, one per line. Do not give 
              reverse complements.
  -h, --help  Print usage
```

