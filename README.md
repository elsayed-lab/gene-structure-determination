RNA-Seq ORF detector
====================

Overview
--------

sequences.

The script takes as input a reference genome and annotation, a bedtools
coverage map for a collection of one or more mapped samples and, optionally, a
set of detected 5' and 3'UTR boundaries.

This script was designed primarily for use with trypanasome RNA-Seq data, for
which most genes contain only a single exon, and genes exists in the genome in
long tracts on the same strand.

Requirements
------------

- [BioPython](http://biopython.org/wiki/Main_Page)
- [Numpy](http://www.numpy.org/)
- [Matplotlib](http://matplotlib.org/)
- [bcbio-gff](https://github.com/chapmanb/bcbb/tree/master/gff)
- [pybedtools](https://pythonhosted.org/pybedtools/)


Running
-------

### Input files

This script requires the following files for input:

1. Reference FASTA
2. Reference GFF
3. RNA-Seq coverage map
4. Trans-splice acceptor site GFF
5. Polyadenylation site GFF

### Construct a genome coverage map

A genome coverage map is a file which contains counts of reads at each position
in a genome, e.g.:

```
...
TcChr1-S	98	0
TcChr1-S	99	0
TcChr1-S	100	0
TcChr1-S	101	1
TcChr1-S	102	3
TcChr1-S	103	3
TcChr1-S	104	3
TcChr1-S	105	6
...
```

You can use the [bedtools
genomecov](http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html)
function to create such a map.

First, of you have multiple samples you wish you wish you use, first combine them
using [samtools view](http://www.htslib.org/doc/samtools.html):

```sh
samtools merge combined.bam */accepted_hits.bam
```

Next, call `genomecov` with the `-d` (single-nt resolution) option:

```sh
bedtools genomecov -d -ibam combined.bam | gzip > combined.coverage.gz
```

Usage Example
-------------

Example call to `rnaseq_orf_detector.py`:

```sh
./rnaseq_orf_detector.py \
    -c /path/to/combined.coverage.gz \
    -g /path/to/TriTrypDB-27_TcruziCLBrenerEsmeraldo-like.gff \
    -f /path/to/TriTrypDB-27_TcruziCLBrenerEsmeraldo-like_Genome.fasta \
    output.gff
```

Diagnostic Images
-----------------

TODO...


