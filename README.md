# SHIMS3 assembly pipeline

### Description

Pipeline for assembling large-insert clones (e.g. BACs, Fosmids) from Nanopore and Illumina data.

Software and documentation written by Daniel W. Bellott.

### Prerequisites

Install the following tools used by the pipeline:

- [gzip](https://www.gnu.org/software/gzip/) (≥1.6; tested with 1.6 and 1.10)
- [Minimap2](https://github.com/lh3/minimap2#install) (≥2.17; tested with 2.17-r941)
- [Racon](https://github.com/lbcb-sci/racon) (≥1.4.17; tested with 1.4.17)
- [SAMtools](http://www.htslib.org/download/) (≥1.11; tested with 1.11)
- [Gap5](http://staden.sourceforge.net) (≥1.2.14-r; tested with 1.2.14-r)
- (Optional) [Consed](http://www.phrap.org/consed/consed.html#howToGet) (≥29.0; tested with 29.0)

Your computer should already have perl installed

All the required perl modules are packaged in the /vendor/cache directory for your convenience.

Tested in Ubuntu 18.04.6 LTS and 20.04.3 LTS.

### Installation

Navigate to the directory where you want to install the pipeline and type:

```
git clone
https://github.com/dwbellott/shims3_assembly_pipeline.git
```

You will be prompted for your username and password.

After the download completes, type:

```
cd shims3_assembly_pipeline
vendor/bin/carton install --cached --deployment
```

To install the cached perl modules.

Installation should take less than a minute on a typical desktop computer.

### Usage

```
$ ./shims3.pl -h
shims3: automate the process of finding full-length nanopore reads, creating a
consensus, polishing with other data, and generating an index for gap5.
(Version: 1.0.0)


USAGE: ./shims3.pl -b vector.fa -n nanopore_reads.fq

  Basic Arguments

  -b  <fasta file>         fasta file of vector backbone sequence

  -n  <fastq file(s)>      fasta or fastq files of reads; may be supplied
                           multiple times or once as a comma separated list

  Optional Arguments

  -o  <output name>        base name for output files (default: shims3)
  -d  <outpur directory>   directory for output  (default: current directory)

  Illumina Reads

  -f  <fastq files(s)>     fasta or fastq files of reads; may be supplied
  -r  <fastq file(s)>      multiple times or once as a comma separated list
  -s  <fastq file(s)>

  PacBio Reads

  -p  <fastq file(s)>      fasta or fastq files of reads; may be supplied
                           multiple times or once as a comma separated list
  Assembly Editor

  -g, --gap5                Generate a gap5 database (default)
  --no-g, --no-gap5         Do not generate a gap5 database

  -c, --consed              Generate an ace file for consed
  --no-c, --no-consed       Do not generate an ace file for consed (default)

  Other:

  -v, --version            print version number and exit

  -h, --help, -u, --usage  print this helpful screen and exit
```

### Demo

The `demo` directory contains real data from clone CH250-214O8 (GenBank accession: AC280424.1), from the rhesus macaque X chromosome. To test your install, run:

```bash
shims3.pl -b demo/input/backbone.fa -n demo/input/nanopore.fq.gz -f demo/input/illumina_1.fq.gz -r demo/input/illumina_2.fq.gz -d demo/output
```

The resulting file at `demo/output/shims3.polish.5.fa` should be identical to the one at `demo/expected_result/shims3.polish.5.fa`

The demo should take 5-10 minutes to run on a typical desktop computer.

### Additional scripts

This package includes additional scripts that automate a subset of the pipeline.
They are occasionally useful in the finishing process

- `identify_full_length.pl` provides a list of full-length nanopore reads
- `iterate_racon.pl` automates the process of polishing a consensus
- `minimap2_and_gap5.pl` automates the process of building a gap5 database
- `minimap2_and_consed.pl` as above, but for consed
- `rotate_fasta.pl` for rotating full-length reads with the vector in the middle
