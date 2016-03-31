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
import csv
import sys
import matplotlib
import numpy as np
import pandas
from argparse import ArgumentParser
from matplotlib.colors import ListedColormap
from pybedtools import BedTool
from Bio import Seq, SeqIO
from BCBio import GFF

matplotlib.use('Agg')
from matplotlib import pyplot as plt

class ORFDetector(object):
    """ORFDetector class definition"""
    def __init__(self):
        """Initialize a new ORFDetector instance"""
        # Get command-line arguments
        args = self._get_args()

        # Store parameters
        # self.bam = args['bam']
        self.fasta = args['fasta']
        self.gff = args['gff']
        self.output = args['output']

        # Load genome sequence
        self.sequence = {x.id:x for x in SeqIO.parse(self.fasta, 'fasta')}

        # Load genome annotations
        self.annotations = {}

        for entry in GFF.parse(open(self.gff)):
            if len(entry.features) > 0 and entry.features[0].type == 'chromosome':
                self.annotations[entry.id] = entry

        # Load RNA-Seq coverage map
        self.coverage = {}

        df = pandas.read_table(args['coverage'], 
                               names=['chr_id', 'loc', 'coverage'])
        for chr_id in self.sequence:
            self.coverage[chr_id] = df[df.chr_id == chr_id].coverage

        # Load RNA-Seq reads
        # self.rnaseq_reads = BedTool(self.bam)
        # Generate coverage maps for each chromosome
        # rnaseq_vec = self.rnaseq_reads.genome_coverage(d=True)

        self.inter_cds_regions = self.get_inter_cds_regions()
        self.search_genome_for_orfs()

        # generate genome coverage and novel ORF plots
        self._plot_features()
        # self.write_output_gff()


    def _create_cds_location_vector(self, chr_id):
        """Create a single-nt resolution ternary vector indicating locations
           covered by known CDSs on either strand"""
        START_IDX = 0
        END_IDX = 1

        cds_vec = np.zeros(len(self.sequence[chr_id]))

        # negative strand
        for loc in self.inter_cds_regions[chr_id][-1]:
            cds_vec[loc[START_IDX]:loc[END_IDX]] = 1 

        # positive strand
        for loc in self.inter_cds_regions[chr_id][1]:
            cds_vec[loc[START_IDX]:loc[END_IDX]] = 2

        return cds_vec

    def _get_args(self):
        """Parses input and returns arguments"""
        parser = ArgumentParser(description='Detect novel ORFs using RNA-Seq data.')

        # parser.add_argument('-b', '--bam', help='Input bam file', required=True)
        parser.add_argument('-c', '--coverage', 
                             help='Single nucleotide resolution genome coverage map',
                             required=True)
        parser.add_argument('-f', '--fasta', help='Input FASTA file',
                            required=True)
        parser.add_argument('-g', '--gff', help='Input GFF file', required=True)
        parser.add_argument('output', help='Location to save output GFF to',
                            metavar='OUTFILE')

        # parse arguments into a dict
        args = vars(parser.parse_args())

        # @TODO
        # check to make sure valid filepaths specified
        # if not os.path.exists(fasta):
            # print("Incorrect genome filepath specified")
            # sys.exit(0)
        # if not os.path.exists(gff):
            # print("Incorrect annotations filepath specified")
            # sys.exit(0)

        return args

    def get_inter_cds_regions(self):
        """Returns a chromosome-indexed dict of inter-CDS regions for a
        specified GFF.
        """
        # Determine locations of inter-CDS regions for each chromosome
        inter_cds_regions = {}

        for chr_id, chromosome in self.annotations.items():
            # get chromosome dimensions (the first feature represents the
            # chromosome itself)
            ch_end = int(chromosome.features[0].location.end)

            # filter out everything except for genes
            genes = [x for x in chromosome.features if x.type == 'gene']

            # order by position along chromosome
            genes.sort(key=lambda x: x.location.start)

            # add range before first gene
            start = 0

            # keep track of strand of polycistronic transcriptional unit
            strand = None

            inter_cds_regions[chr_id] = {
                -1: [],
                +1: []
            }

            # iterate through genes and store the ranges between them;
            # for TriTrypDB files, the gene boundaries are generally the same
            # as the CDS boundaries.
            for gene in genes:
                # Determine location for the region up to start of the gene
                end = int(gene.location.start)

                # Skip over snoRNAs, etc. that are nested inside of other genes
                # For example: TcCLB TcChr22-2 179,000:180,000
                if end <= start:
                    next

                # Add CDS to relevant list based on strand
                if strand is None:
                    # Left-most gene
                    inter_cds_regions[chr_id][gene.location.strand].append((start, end))
                elif strand != gene.location.strand:
                    # Add ORFs in both directions at transcription switch sites (TSSs)
                    inter_cds_regions[chr_id][+1].append((start, end))
                    inter_cds_regions[chr_id][-1].append((start, end))
                else:
                    # Within PTU; look for ORFs on same strand
                    inter_cds_regions[chr_id][strand].append((start, end))

                # update start counter and strand
                start = int(gene.location.end)
                strand = gene.location.strand

            # add region after last gene
            inter_cds_regions[chr_id][strand].append((start, ch_end))

        return inter_cds_regions

    def search_genome_for_orfs(self, min_protein_length=30):
        # Create a list to store ORF entries
        self.orfs = []

        # Iterate through inter-CDS regions
        for i, chr_id in enumerate(self.inter_cds_regions, 1):
            print("Processing %s (%d/%d)" % (chr_id, i, len(self.sequence)))
            for strand in [-1, 1]:
                for j, region in enumerate(self.inter_cds_regions[chr_id][strand]):
                    # get sequence record for the range
                    record = self.sequence[chr_id][region[0]:region[1]]

                    # Find ORFs in each of the six possible reading frames
                    # that are at least the specified length
                    inter_cds_orfs = self.find_orfs(record.seq,
                                                    min_protein_length, strand)

                    # Write GFF entries for each match
                    for k, orf in enumerate(inter_cds_orfs, 1):
                        # Convert coordinates to chromosomal position
                        #if orf[2] == "+":
                        start = orf[0] + region[0] + 1
                        stop = orf[1] + region[0]
                        str_strand = orf[2]

                        self.orfs.append({
                            'id': "%s.ORF.%d.%d" % (chr_id, j, k),
                            'chr': chr_id,
                            'start': start,
                            'stop': stop,
                            'strand': str_strand
                        })

    def find_orfs(self, seq, min_protein_length, strand=1, trans_table=1,
                  ignore_ambiguous_orfs=True):
        """
        Finds ORFs of a specified minimum protein length in a SeqRecord.

        Based on: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec360
        """
        answer = []
        seq_len = len(seq)

        # Get sequence associated with the specified location and strand
        if strand == 1:
            dna_seq = seq
        else:
            dna_seq = seq.reverse_complement()

        for frame in range(3):
            trans = str(dna_seq[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0

            # Iterate through ORFS in reading frame
            while aa_start < trans_len:
                # Set end counter to position of next stop codon
                aa_start = trans.find("M", aa_start)
                aa_end = trans.find("*", aa_start)

                # If no start or stop codons found, stop here
                if aa_start == -1 or aa_end == -1:
                    break

                if (aa_end < aa_start):
                    raise Exception('wtf')

                # Compute coordinates of ORF
                if strand == 1:
                    start = frame + aa_start * 3
                    end = min(seq_len, frame + aa_end * 3 + 3)
                else:
                    start = seq_len - frame - aa_end * 3 - 3
                    end = seq_len - frame - aa_start * 3

                # Add to output
                str_strand = "+" if strand == 1 else '-'

                # Check to make sure ORF doesn't contain a bunch of N's
                if ignore_ambiguous_orfs:
                    num_unknown = trans[aa_start:aa_end].count('X')
                    if (num_unknown / (aa_end - aa_start)) > 0.25:
                        aa_start = aa_end + 1
                        continue

                # increment start counter
                aa_start = aa_end + 1

                # Add ORF coordinates and continue
                answer.append((start, end, str_strand))

        # Sort results
        answer.sort()

        return answer

    def write_output_gff(self):
        # Iterate over inter-CDS regions and find ORFs of at least the specified
        # length in any of the six possible reading frames and output as GFF entries
        fp = open(self.output, 'w')

        # Write csv header
        fp.write("##gff-version\t3\n")
        fp.write("##feature-ontology\tsofa.obo\n")
        fp.write("##attribute-ontology\tgff3_attributes.obo\n")

        # Write header to output
        writer = csv.writer(fp, delimiter='\t')

        # Iterate through inter_cds regions
        for orf in self.orfs:
            # Write entry
            gff_attrs = "ID=%s;Name=%s" % (orf['id'], orf['id'])
            writer.writerow([orf['chr'],
                            "ElSayedLab",
                            'ORF',
                            orf['start'],
                            orf['stop'],
                            '.',
                            orf['strand'],
                            '.',
                            gff_attrs])

        # clean up
        fp.close()

    def _plot_features(self):
        """Creates image plots of specified features (CDSs, RNA-Seq coverage,
           or novel ORFs)"""
        # create a separate plot for each chromosome
        #
        #  0 empty
        #  1 CDS on negative strand
        #  2 CDS on positive strand
        #  3 RNA-Seq reads / novel ORF
        #
        for chr_id in self.sequence:
            # create single-nt vectors indicating location for each feature
            # of interest
            cds_vec = self._create_cds_location_vector(chr_id)

            # rna-seq coverage map
            rnaseq_vec = self.coverage[chr_id]

            # for now, convert to binary
            rnaseq_vec[rnaseq_vec > 0] = 3

            # output filepath
            filename = '%s_orig.png' % chr_id
            outfile = os.path.join('output', 'image', 'raw', filename) 

            self.plot_genome_image(outfile, cds_vec, rnaseq_vec)

    def plot_genome_image(self, outfile, cds_vec, rnaseq_vec):
        """Creates an image plot for all chromosomes in a genome with known
        CDS's shown with one color, and regions of coverage shown in another"""
        # convert input to numpy arrays
        cds_vec = np.array(cds_vec)
        rnaseq_vec = np.array(rnaseq_vec)

        # combine cds and rna-seq arrays
        rnaseq_vec[cds_vec != 0] = cds_vec[cds_vec != 0]

        # convert vector to a zero-filled square matrix
        mat_dim = int(np.ceil(np.sqrt(len(rnaseq_vec))))
        fill = np.zeros(mat_dim**2 - len(rnaseq_vec))

        mat = np.concatenate((rnaseq_vec, fill)).reshape(mat_dim, mat_dim)

        # flip so that 0,0 is top-left
        mat = mat[::-1]

        # colormap
        colors = ['#000000', '#a2e803', '#0219e8', '#a302e8']
        cmap = ListedColormap(colors, name='genome_cmap')

        plt.figure()
        plt.matshow(mat, cmap=cmap)

        # create output directory and save
        outdir = os.path.dirname(outfile)

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(11.25, 11.25)
        plt.savefig(outfile, dpi=96)

        # flip to display from top-left corner
        # mat = t(apply(mat, 2, rev))
        # image(z=mat, col=c("#333333", "red", "green", "blue", "yellow", "cyan",
                            # "magenta", "orange"))

def _colormap():
    """Custom colormap to use for genome plots"""
    

if __name__ == "__main__":
    orf_detector = ORFDetector()
