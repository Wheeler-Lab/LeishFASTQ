# LeishFASTQ

Find sequencing barcodes in gzipped FASTQ files.

## Mode of operation

This tool counts barcodes of a given length, flanked by two pre-defined sequences:
```
5' - GTGTATCGGATGTCAGTTGCATGCGTACGTCAGTACGTATAATGCAGACCTGCTGC - 3'
     \__________________/\______________/\__________________/
         Left Flank          Barcode         Right Flank
```
LeishFASTQ peruses the given FILENAME.fastq.gz file and finds if a sequence matches the structure above, with the flanking sequences enclosing a barcode sequence of the required lengths. While matching, it allows for up to the given number of mismatches or deletions (using `--allow-mismatches`) in either flanking sequence. If a match is positive, the tool then adds the barcode sequence to its database. The number of times each barcode is seen is then output as a CSV file with the name FILENAME.csv. If the name of a file containing a list of known barcode sequences is given, LeishFASTQ will also allow single substitutions of the given sequences if `--allow-mismatches` is non-zero. Allowing more than one substition might lead to ambiguities in most barcode databases, therefore the tool disallows this.

If the file name contains the string "__R2__" LeishFASTQ assumes it has been given a reverse read file and matches the reverse complement of the scheme above.

## Usage

```
usage: leishfastq [-h] [--allow-mismatches MISMATCHES]
                  [--known-barcodes BARCODE_LIST]
                  [--left-flanking-sequence LEFT_FLANKING_SEQUENCE]
                  [--right-flanking-sequence RIGHT_FLANKING_SEQUENCE]
                  [--barcode-length BARCODE_LENGTH] [--nr-workers NR_WORKERS]
                  fastq-gz-file [fastq-gz-file ...]

Find barcodes in gzipped fastq files.

positional arguments:
  fastq-gz-file         List of FASTQ files to process

options:
  -h, --help            show this help message and exit
  --allow-mismatches MISMATCHES
                        Allow up to the given number of mismatches in each
                        flank and the barcode, if a list of known barcodes is
                        given (default 0)
  --known-barcodes BARCODE_LIST
                        File containing a list of known barcodes to match
                        against. This helps if mismatches are allowed.
  --left-flanking-sequence LEFT_FLANKING_SEQUENCE
                        Left barcode flanking sequence (default:
                        GTGTATCGGATGTCAGTTGC)
  --right-flanking-sequence RIGHT_FLANKING_SEQUENCE
                        Right barcode flanking sequence (default:
                        GTATAATGCAGACCTGCTGC)
  --barcode-length BARCODE_LENGTH
                        Barcode length in nucleotides (default: 17)
  --nr-workers NR_WORKERS
                        Number of parallel processes to use (default: same as
                        the number of cores)
```

## Example

```
leishfastq --allow-mismatches 2 --left-flanking-sequence GTGTATCGGATGTCAGTTGC --right-flanking-sequence GTATAATGCAGACCTGCTGC --barcode-length 17 --nr-workers 12
```
