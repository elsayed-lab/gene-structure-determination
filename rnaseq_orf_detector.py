#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Novel ORF detection using RNA-Seq data
Keith Hughitt (khughitt@umd.edu)
2016/03/29

The goal of this tool is to provide a way to search for unannotated ORFs which
fall between existing annotated ORFs with RNA-Seq read support.

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
from Bio import Seq, SeqIO
from BCBio import GFF
from orfdetector.orfs import search_genome_for_orfs
from orfdetector.genome import get_cds_regions, get_inter_cds_regions
from orfdetector.plots import plot_genome_features
from orfdetector.util import max_filter
from orfdetector.io import parse_args, write_output_gff

def main():
    """Main"""
    # Get command-line arguments
    args = parse_args()

    # Store parameters
    # bam = args['bam']
    fasta = args['fasta']
    gff = args['gff']
    output = args['output']
    min_protein_length = args['min_protein_length']

    # Load genome sequence
    genome_sequence = {x.id:x for x in SeqIO.parse(fasta, 'fasta')}

    # Load genome annotations
    annotations = {}

    for entry in GFF.parse(open(gff)):
        if len(entry.features) > 0 and entry.features[0].type == 'chromosome':
            annotations[entry.id] = entry

    # Load RNA-Seq coverage map
    rnaseq_coverage = {}

    df = pandas.read_table(args['coverage'], names=['chr_id', 'loc', 'coverage'])
    for chr_id in genome_sequence:
        rnaseq_coverage[chr_id] = df[df.chr_id == chr_id].coverage

    # Load SL and Poly(A) sites
    # TODO
    # sl_sites = {}
    # polya_sites = {}

    # for entry in GFF.parse(open(sl_gff)):
        # if len(entry.features) > 0 and entry.features[0].type == 'chromosome':
            # sl_sites[entry.id] = entry

    # Load RNA-Seq reads
    # from pybedtools import BedTool
    # rnaseq_reads = BedTool(bam)
    # Generate coverage maps for each chromosome
    # rnaseq_vec = rnaseq_reads.genome_coverage(d=True)

    # Determine CDS and inter-CDS boundaries
    cds_regions = get_cds_regions()
    inter_cds_regions = get_inter_cds_regions()

    # Detect ORFs in inter-CDS regions
    orfs = search_genome_for_orfs(genome_sequence, inter_cds_regions,
                                  min_protein_length)

    # Optimize SL / polyadenylation primary site selection for each
    # annotated CDS, assigning novel ORFs where appropriate
    find_primary_sites()

    # discrete or continuous plot
    # TODO: refactor this logic / make an option
    plot_type = 'discrete'
    # _plot_type = 'continuous'

    # generate genome coverage and novel ORF plots
    plot_genome_features(genome_sequence, cds_regions, rnaseq_coverage, plot_type)

    # write_output_gff()

if __name__ == "__main__":
    main()

