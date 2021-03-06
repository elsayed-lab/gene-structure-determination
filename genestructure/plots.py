"""
Plotting related functionality
"""
import os
import numpy as np
import matplotlib
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

matplotlib.use('Agg')
from matplotlib import pyplot as plt

def get_cds_regions(annotations):
    """Returns a chromosome-indexed dict of CDS locations for a
    specified GFF.

    Arguments
    ---------
    annotations: dict
        A dictionary of chromomsome Seq instances as parsed from a GFF file.

    Returns
    -------
    cds_regions: dict
        A dictionary of CDS coordinates across the genome.
    """
    # Determine locations of CDS regions for each chromosome
    cds_regions = {}

    for chr_id, chromosome in annotations.items():
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

        cds_regions[chr_id] = {
            -1: [],
            +1: []
        }

        # iterate over genes and add CDS coordinates
        for gene in genes:
            coords = (int(gene.location.start), int(gene.location.end))
            cds_regions[chr_id][gene.location.strand].append(coords)

    return cds_regions

def get_inter_cds_regions(annotations):
    """Returns a chromosome-indexed dict of inter-CDS regions for a
    specified GFF.

    Arguments
    ---------
    annotations: dict
        A dictionary of chromomsome Seq instances as parsed from a GFF file.

    Returns
    -------
    cds_regions: dict
        A dictionary of CDS coordinates across the genome.
    """
    # Determine locations of inter-CDS regions for each chromosome
    inter_cds_regions = {}

    for chr_id, chromosome in annotations.items():
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
                continue

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


def plot_genome_features(genome_sequence, cds_regions, rnaseq_coverage, 
                         plot_type='discrete'):
    """Creates image plots of specified features (CDSs, RNA-Seq coverage,
        or novel ORFs)"""
    # create a separate plot for each chromosome
    #
    #
    for chr_id in genome_sequence:
        # create single-nt vectors indicating location for each feature
        # of interest
        # inter_cds_vec = _create_feature_vector(inter_cds_regions[chr_id],
                                                    # len(genome_sequence[chr_id]))
        cds_vec = create_feature_vector(cds_regions[chr_id],
                                        len(genome_sequence[chr_id]))

        # rna-seq coverage map
        rnaseq_vec = rnaseq_coverage[chr_id]

        # threshold
        # TODO: make a command-line option; may not ultimately need to do
        # any thresholding if density overy each ORF is used for selection,
        # but could still be useful for visualizing potential novel ORFs
        if plot_type == 'discrete':
                threshold = 500
                rnaseq_vec[rnaseq_vec <= threshold] = 0
                rnaseq_vec[rnaseq_vec > threshold] = 3

        # output filepath
        filename = '%s_orig.png' % chr_id
        outfile = os.path.join('output', 'image', 'raw', filename) 

        plot_genome_image(outfile, cds_vec, rnaseq_vec)

def plot_genome_image(outfile, cds_vec, rnaseq_vec, plot_type='discrete'):
    """Creates an image plot for all chromosomes in a genome with known
    CDS's shown with one color, and regions of coverage shown in another"""
    # convert input to numpy arrays
    cds_vec = np.array(cds_vec)
    rnaseq_vec = np.array(rnaseq_vec)

    # combine cds and rna-seq arrays
    if plot_type == 'discrete':
        # for discrete plot, show each feature using a single color
        rnaseq_vec[cds_vec != 0] = cds_vec[cds_vec != 0]
    else:
        # for the continuous plot, coverage is plotted across a range
        # of colors, while CDSs are set to a single color (max)
        rnaseq_vec = np.log2(rnaseq_vec)
        rnaseq_vec[cds_vec != 0] = max(rnaseq_vec)

    # convert vector to a zero-filled square matrix
    mat_dim = int(np.ceil(np.sqrt(len(rnaseq_vec))))
    fill = np.zeros(mat_dim**2 - len(rnaseq_vec))

    mat = np.concatenate((rnaseq_vec, fill)).reshape(mat_dim, mat_dim)

    # flip so that 0,0 is top-left
    # mat = mat[::-1]

    #
    # discrete colormap
    #
    #  0 black      Empty
    #  1 green      CDS on negative strand
    #  2 blue       CDS on positive strand
    #  3 magenta    RNA-Seq reads / novel ORF
    #
    if plot_type == 'discrete':
        # default colors
        # colors = ['#000000', '#a2e803', '#0219e8', '#a302e8']

        # dream magnet (CL)
        colors = ['#343838', '#00B4CC', '#005F6B', '#008C9E']
        cmap = ListedColormap(colors, name='genome_cmap')
    else:
        # for the continuous plot, use the hot colormap, but add an
        # extra color at the end to use for CDS's
        hot_cmap = matplotlib.cm.get_cmap("hot", 1024)
        cds_color = list(matplotlib.colors.hex2color('#0219e8')) + [1.0]
        hot_vals = hot_cmap(np.arange(1024))
        hot_vals[-1] = cds_color
        cmap = LinearSegmentedColormap.from_list("genome_cmap", hot_vals) 

    plt.figure()
    plt.matshow(mat, cmap=cmap)

    # create output directory and save
    outdir = os.path.dirname(outfile)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(40, 20)
    plt.savefig(outfile, dpi=96)
    # fig.set_size_inches(60, 20)
    # plt.savefig(outfile, dpi=300)
    plt.clf()

def create_feature_vector(features, length):
    """Create a single-nt resolution ternary vector representation of
    a given set of features across a chromosome"""
    START_IDX = 0
    END_IDX = 1

    output_vector = np.zeros(length)

    # negative strand
    for loc in features[-1]:
        output_vector[loc[START_IDX]:loc[END_IDX]] = 1 

    # positive strand
    for loc in features[1]:
        output_vector[loc[START_IDX]:loc[END_IDX]] = 2

    return output_vector

