#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Gene structure and novel ORF detection using RNA-Seq data
Keith Hughitt (khughitt@umd.edu)
2016/03/29

The goal of this script is to attempt to use RNA-Seq information from RNA-Seq
reads to define the UTR boundaries for genes, and to detect putative novel ORFs
with transcriptional and other types of evidence.

Currently, bedtools does not support multiple input BAM files when using the
genomecov function to compute read coverage across the genome. As such, if you
wish to use multiple bam files, you must first combine them into a single bam
file.

For example, using samtools, you can do:

    samtools merge combined.bam */accepted_hits.bam
    samtools sort combined.bam > combined_sorted.bam

To generate a genome coverage maps:

    bedtools genomecov -d -ibam combined_sorted.bam > combined_sorted.coverage

"""
import warnings
from genestructure.analyzer import GeneStructureAnalyzer
from genestructure.io import parse_args
from Bio import BiopythonWarning

def main():
    """Main"""
    # Ignore Biopython warnings
    warnings.simplefilter('ignore', BiopythonWarning)

    # Get command-line arguments and instanstiate GeneStructureAnalyzer
    args = parse_args()
    analyzer = GeneStructureAnalyzer(args['fasta'], args['gff'],
                                     args['coverage'], args['sl_gff'],
                                     args['polya_gff'],
                                     args['min_protein_length'], args['outdir'])

if __name__ == "__main__":
    main()

