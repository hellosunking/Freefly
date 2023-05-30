## An example pipeline for size and end motif analysis of cfDNA using reference-based approach

The main script is `conventional.pipeline`; call it without parameters to show the usage:
```
Usage: ./conventional.pipeline [options] <-x bowtie2.index> <-g genome.fa> <-o output.prefix> <-1 read1.fq> <-2 read2.fq>

Options:
  -x path    Set the path to bowtie2 index
  -g path    Set the path to the fasta sequence of the genome
  -s size    Set minimum read size. Default: 36
  -t thread  Set running threads. Default: 8
  -k kit     Set kit for trimming adaptors. Default: illumina
  -m kmer    Set k-mer for motif analysis. Default: 4
  -q score   Mininum mapping score for analysis. Default: 30

```
Note that you need to prepare bowtie2 index files as well as genome sequences (in FASTA format) first.

