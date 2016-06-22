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
import os
import pandas
import warnings 
from Bio import Seq, SeqIO, BiopythonWarning
from Bio import  
from BCBio import GFF
from genestructure.genome import get_cds_regions, get_inter_cds_regions, find_primary_sites
from genestructure.plots import plot_genome_features
from genestructure.io import load_gff, parse_args, write_output_gff

def main():
    """Main"""
    # Ignore Biopython warnings
    warnings.simplefilter('ignore', BiopythonWarning) 

    # Get command-line arguments
    args = parse_args()

    print(args)

    # Load genome sequence
    genome_sequence = {x.id:x for x in SeqIO.parse(args['fasta'], 'fasta')}

    # Load genome annotations
    genome_annotations = load_gff(args['gff'])

    # Load RNA-Seq coverage map
    rnaseq_coverage = {}

    df = pandas.read_table(args['coverage'], names=['chr_id', 'loc', 'coverage'])

    for chr_id in genome_sequence:
        rnaseq_coverage[chr_id] = df[df.chr_id == chr_id].coverage

    # Load SL and Poly(A) sites
    sl_sites = load_gff(args['sl_gff'])
    polya_sites = load_gff(args['polya_gff'])

    # Load RNA-Seq reads
    # from pybedtools import BedTool
    # rnaseq_reads = BedTool(bam)
    # Generate coverage maps for each chromosome
    # rnaseq_vec = rnaseq_reads.genome_coverage(d=True)

    # Determine CDS and inter-CDS boundaries
    cds_regions = get_cds_regions(genome_annotations)
    inter_cds_regions = get_inter_cds_regions(genome_annotations)

    # Optimize SL / polyadenylation primary site selection for each
    # annotated CDS, assigning novel ORFs where appropriate
    find_primary_sites(genome_sequence, genome_annotations, sl_sites,
                       polya_sites, args['min_protein_length'])

    # generate genome coverage and novel ORF plots
    plot_genome_features(genome_sequence, cds_regions, rnaseq_coverage,
                         args['plot_type'])

    # write_output_gff()

if __name__ == "__main__":
    main()

