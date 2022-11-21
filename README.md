# Barcode analyzer

A fast tool to search for barcode sequences inside fasta/fastq data. Can also remove the reads that have a barcode. Internally, uses the [Aho-Corasick](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm) multiple string matching algorithm implementation from [here](https://github.com/cjgdev/aho_corasick).

## Compiling

Compile by running `make`. Compiling requires a C++ compiler with support for the C++17 standard. Tested to work with g++ version 9. 

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
  -b arg         A file containing the barcodes, one per line.
  -v, --verbose  Verbose output.
  -h, --help     Print usage
```

### Filter

```
Remove all sequence from a fasta/fastq file that have a barcode sequence.
Usage:
  filter [OPTION...]

  -i arg      The sequence file in fasta or fastq format.
  -o arg      Output file.
  -b arg      A file containing the barcodes, one per line.
  -h, --help  Print usage
```

