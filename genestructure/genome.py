"""
Functions for retrieving features across a set of genome annotations.
"""
import numpy as np
from .orfs import find_orfs
from .util import max_filter

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

def find_primary_sites(genome_sequence, genome_annotations, sl_sites,
                       polya_sites, min_protein_length=30):
    """
    Simultaneously chooses optimal SL and Poly(A) sites for a set of two
    neighboring CDS's and detecting novel ORFs where appropriate.

    The goal of this function is to attempt to determine the most likely
    primary SL/Poly(A) sites for the set of known CDS's for which we have
    sites detected for and, in cases where evidence exists for possible
    novel ORFs (e.g. additional SL and Poly(A) sites surrounding an ORF in
    between two CDS's), to extend the current genome annotation to include
    these ORFs. This also serves to provide us with better estimates about
    the true 5' and 3'UTR lengths for parasite genes.

    Depending on the amount of detected SL and Poly(A) sites between a set
    of known CDS's, there may be multiple possible configurations of
    feature assignments to both the known CDS's and various possible novel
    ORFs.

    In order to attempt to select an "optimal" combination of novel ORF and
    SL / Poly(A) site assignments, we will use a scoring function which
    takes into account different types of support.

    For each possible configuration, we will compute a score defined by:

        10E6 * # Assigned features +
        10E3 * sum(reads mapped to each feature) +
        10E1 * sum(read density spanning the ORF) +
        1    / sum(distance from each site to associated feature)

    This optimization criterion will favor the selection of configurations
    where the maximal number of CDS/ORFs are assigned SL and Poly(A) sites,
    and will attempt to choose the sites with the most support (largest
    number of reads supporting those sites). Finally, if there are multiple
    ORFs which can be chosen for a given SL/Poly(A) combination, the ORF
    which captures the most read density (usually the longest ORF), will be
    selected.

    Arguments
    ---------
    genome_sequence: Seq
        Genome sequence as a BioPython Seq instance
    genome_annotations: dict
        Genome annotations as a dictionary of chromosome entries
    sl_sites: dict
        Detected spliced-leader sites as a chromosome-indexed dict
    polya_sites: dict
        Detected polyadenylation sites as a chromosome-indexed dict
    min_protein_length: int
        Minimum length in amino acids to allow for novel ORFs
    """
    # Iterate over chromosomes
    for chr_id, chromosome in genome_annotations.items():
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

        # iterate over genes
        # for gene in genes:
        for right_gene in genes:
            # Determine location for the region up to start of the gene
            end = int(right_gene.location.start)

            # Skip over snoRNAs, etc. that are nested inside of other genes
            # For example: TcCLB TcChr22-2 179,000:180,000
            if end <= start:
                continue

            # Get sequence record for the range
            record = genome_sequence[chr_id][start:end]

            # Find ORFs in each of the six possible reading frames
            # that are at least the specified length
            inter_cds_orfs = find_orfs(record.seq, min_protein_length, strand)

            # Retrieve SL and Poly(A) sites in region
            inter_cds_sl_sites = get_features_by_range(sl_sites[chr_id], start, end)
            inter_cds_polya_sites = get_features_by_range(polya_sites[chr_id], start, end)

            # If no sites found, skip
            if len(inter_cds_sl_sites) == 0 and len(inter_cds_polya_sites) == 0:
                continue

            # Filter sites down to ~3 optimal SL's / Poly(A)'s
            filtered_sl_sites = filter_features(inter_cds_sl_sites, start, end)
            filtered_polya_sites = filter_features(inter_cds_polya_sites, start, end)

            # Add CDS to relevant list based on strand
            if strand is None:
                # Left-most gene
                if right_gene.location.strand == 1:
                    if len(inter_cds_sl_sites) == 0:
                        continue
                    sl = get_max_feature(filtered_sl_sites,
                                         location_preference='right')
                else:
                    polya = get_max_feature(filtered_polya_sites,
                                            location_preference='left')
                pass
            elif strand != right_gene.location.strand:
                # Transcription Switch Sites (TSSs)
                # TODO: assign features at TSS
                pass
            else:
                # Within PTU; look for ORFs on same strand
                # TODO: assign within PTU features 
                pass

            # update gene, start counter and strand
            left_gene = right_gene
            start = int(left_gene.location.end)
            strand = left_gene.location.strand


def get_features_by_range(annotations, start, stop):
    """Returns all annotation features which fall within the specified range"""
    features = []

    for feature in annotations.features:
        if feature.location.start > start and feature.location.end < stop:
            features.append(feature)

    return features

def get_max_feature(features, location_preference='right'):
    """Takes a list of features and finds the one with the highest coverage
    support. In the case of ties it will choose the site closest to the
    direction specified."""
    # get scores for each feature
    scores = np.array([int(site.qualifiers['score'][0]) for site in features])

    # maximum coverage score
    max_score = max(scores)

    # number of sites with the maximum coverage score
    # num_matches = len([x for x in scores if x == max_score])
    max_sites = np.array(features)[scores == max_score]

    # If only one site found with maximum coverage, return it
    if len(max_sites) == 1:
        return max_sites[0]

    # Otherwise find the site closest to the edge specified
    positions = np.array([site.location.start for site in features])

    if location_preference == 'left':
        return np.array(features)[positions == min(positions)][0]
    else:
        return np.array(features)[positions == max(positions)][0]

def features_to_1d_array(features, start, end):
    """Takes a list of GFF features such as SL sites and creates a 1d array
    representation of them with the value at each position in the array
    equaling the score of the feature, or 0 if no features exist at that
    location.

    Arguments
    ---------
    features: list
        A list of GFF features
    start: int
        Start coordinate for feature range
    end: int
        End coordinate for feature range

    Returns
    -------
    np.array: A 1-dimensional array representation of the features.
    """
    arr = np.zeros(end - start)

    for feature in features:
        arr[feature.location.start - start] = int(feature.qualifiers.get('score')[0])

    return arr

def filter_features(features, start, end, window_size=10, max_features=3):
    """Filters a set of features (e.g. SL sites) to remove any low-support or
    closeby sites.
    
    Parameters
    ----------
    features: list
        A list of SeqFeature instances representing the locations of putative
        SL accepted sites and polyadenylation sites in the genome.
    start: int
        Start coordinate along the chromosome for the region currently scanning
    end: int
        End coordinate along the chromosome for the region currently scanning
    window_size: int
        Sliding window size to use when smoothing feature data
    max_features: int
        Maximum number of features to return - will find at most this number of
        features, based on coverage support.
    """
    # Convert features to a 1d array representation and filter
    feature_arr = features_to_1d_array(features, start, end)

    # Smooth out everything except for local peaks
    feature_arr = max_filter(feature_arr, width=window_size)

    # Rank by score and keep only the top N highest-scoring features
    rank_cutoff = len(feature_arr) - max_features
    feature_arr[feature_arr.argsort().argsort() < rank_cutoff] = 0

    # Get the indices of coverage peaks
    indices = np.arange(len(feature_arr))[feature_arr != 0] + start

    return [x for x in features if x.location.start in indices]

