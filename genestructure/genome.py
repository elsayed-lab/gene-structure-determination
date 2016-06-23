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
                       polya_sites, rnaseq_coverage, min_protein_length=30):
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

        10E6  * # Assigned features +
        10E0  * sum(reads mapped to each feature) +
        10E-6 * sum(read density spanning the ORF) +
       -10E6  * if ORF included

    This optimization criterion will favor the selection of configurations
    where the maximal number of CDS/ORFs are assigned SL and Poly(A) sites,
    and will attempt to choose the sites with the most support (largest
    number of reads supporting those sites). If there are multiple ORFs which
    can be chosen for a given SL/Poly(A) combination, the ORF which captures
    the most read density (usually the longest ORF), will be selected.
    Finally, if two configurations score equally well, but one uses only the
    existing annotations, and another assigns features to a novel ORF, the
    configuration using only the existing annotations will be favored.

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
    rnaseq_coverage: dict
        Dictionary of 1d chromosome RNA-Seq read coverage vectors
    min_protein_length: int
        Minimum length in amino acids to allow for novel ORFs
    """
    # GFF entries
    gff_entries = []

    # counters
    total_counter    = 0
    assigned_counter = 0

    # Iterate over chromosomes
    for chr_id, chromosome in genome_annotations.items():
        print("=======" + chr_id + "=======")

        # get chromosome length (the first feature represents the
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
        left_gene = None

        # iterate over genes
        # for gene in genes:
        for right_gene in genes:
            # Determine location for the region up to start of the gene
            end = int(right_gene.location.start)

            total_counter = total_counter + 1

            # Gene ids (used for logging)
            left_gene_id = left_gene.id if left_gene != None else 'None'
            right_gene_id = right_gene.id
            gene_ids = "%s - %s" % (left_gene_id, right_gene_id)

            # Skip over snoRNAs, etc. that are nested inside of other genes
            # For example: TcCLB TcChr22-2 179,000:180,000
            if end <= start:
                print(" [SKIPPING] %s: Nested genes" % gene_ids)
                left_gene = right_gene
                start = int(left_gene.location.end)
                strand = left_gene.location.strand
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
                print(" [SKIPPING] %s: No features detected" % gene_ids)
                left_gene = right_gene
                start = int(left_gene.location.end)
                strand = left_gene.location.strand
                continue

            # Filter sites down to ~3 optimal SL's / Poly(A)'s
            filtered_sl_sites = filter_features(inter_cds_sl_sites, start, end)
            filtered_polya_sites = filter_features(inter_cds_polya_sites, start, end)

            # TODO: deal with cases where only sl site or polya sites were
            # found, but not found.
            if len(filtered_sl_sites) == 0:
                print(" [SKIPPING] %s: No SL sites detected" % gene_ids)
                left_gene = right_gene
                start = int(left_gene.location.end)
                strand = left_gene.location.strand
                continue
            elif len(filtered_polya_sites) == 0:
                print(" [SKIPPING] %s: No Poly(A) sites detected" % gene_ids)
                left_gene = right_gene
                start = int(left_gene.location.end)
                strand = left_gene.location.strand
                continue

            # Add CDS to relevant list based on strand
            if strand is None:
                # Left-most gene; not looking for undetected ORFs in this
                # region at the moment.
                if right_gene.location.strand == 1:
                    if len(inter_cds_sl_sites) == 0:
                        print(" [SKIPPING] %s: No SL sites detected" % gene_ids)
                        left_gene = right_gene
                        start = int(left_gene.location.end)
                        strand = left_gene.location.strand
                        continue
                    sl = get_max_feature(filtered_sl_sites,
                                         location_preference='right')
                    entries = [build_gff_utr_entry(sl, right_gene, chr_id)]
                    gff_entries = gff_entries + entries
                else:
                    polya = get_max_feature(filtered_polya_sites,
                                            location_preference='left')
                    entries = [build_gff_utr_entry(polya, right_gene, chr_id)]
                    gff_entries = gff_entries + entries
                print(" [ASSIGNED] %s: Annotated features only" % gene_ids)
            elif strand != right_gene.location.strand:
                # Transcription Switch Sites (TSSs)
                # TODO: assign features at TSS
                print(" [SKIPPING] %s: Inter-PTU region" % gene_ids)
                left_gene = right_gene
                start = int(left_gene.location.end)
                strand = left_gene.location.strand
                continue
            else:
                # Assign within PTU features 
                # In this case, we will evaluate each possible combination of
                # SL/Poly(A) sites to choose the configuration which is most
                # appropriate for the neighboring pair of genes.
                entries = assign_inter_ptu_sites(filtered_sl_sites, 
                                                 filtered_polya_sites, 
                                                 inter_cds_orfs, 
                                                 rnaseq_coverage[chr_id],
                                                 start, end, strand, gene_ids)
                gff_entries = gff_entries + entries

            # update gene, start counter and strand
            left_gene = right_gene
            start = int(left_gene.location.end)
            strand = left_gene.location.strand

            if len(entries) > 0:
                assigned_counter = assigned_counter + 1

    print("Total inter-CDS regions scanned: %d (%d assigned)" % (total_counter,
                                                                 assigned_counter))

    return gff_entries

def assign_inter_ptu_sites(sl_sites, polya_sites, orfs, coverage,
                           region_start, region_end, strand, gene_ids):
    """
    Attempts to choose the most like primary SL and Poly(A) sites for a pair of
    adjacenct genes on the same polycistronic transcriptional unit (PTU),
    allowing for the possibility of novel ORFs.

    For a specific set of SL sites, Poly(A) sites, and ORFs in a inter-CDS
    region, each possible assignment of features is considered, including a
    configuration with no novel ORFs.

    The highest scoring combination of assignments is then used to assign
    primary sites to each feature.
    """
    max_score = 0
    sl = None
    polya = None

    # Iterate over any novel ORFs
    for orf in orfs:
        # ORF boundaries relevant to chromosome start
        orf_start = region_start + orf[0]
        orf_end   = region_start + orf[1]

        # compute density of reads covering ORF
        orf_coverage = sum(coverage[orf_start:orf_end])
        
        # compute score for region to the left of ORF
        lhs_sl_sites    = [sl for sl in sl_sites if sl.location.start < orf_start]
        lhs_polya_sites = [polya for polya in polya_sites if polya.location.start < orf_start]
        lhs = select_optimal_features(lhs_sl_sites, lhs_polya_sites, strand)

        # compute score for region to the right of ORF
        rhs_sl_sites    = [sl for sl in sl_sites if sl.location.start > orf_end]
        rhs_polya_sites = [polya for polya in polya_sites if polya.location.start > orf_end]
        rhs = select_optimal_features(rhs_sl_sites, rhs_polya_sites, strand)

        # total number of features assigned for configuration
        total_assigned = bool(lhs['sl']) + bool(lhs['polya']) + bool(rhs['sl']) + bool(rhs['polya'])

        # compute score
        # one substracted from total_assigned to penalize inclusion of ORFs
        score = (10E6  * (total_assigned - 1) +
                 10E0  * lhs['read_support'] + rhs['read_support'] +
                 10E-6 * orf_coverage)

        if  score > max_score:
            max_score = score

            # for now, just keep track of feature assignments for annotated
            # genes and don't worry about the ORFs themselves...
            if strand == 1:
                sl = rhs['sl']
                polya = lhs['polya']
            else:
                sl = lhs['sl']
                polya = rhs['polya']

    # Also compute the score for the best arrangement without any novel ORFs
    no_orfs = select_optimal_features(sl_sites, polya_sites, strand)

    # total number of features assigned for configuration
    total_assigned = bool(no_orfs['sl']) + bool(no_orfs['polya'])

    # compute score
    score = (10E6  * total_assigned +
             10E0  * no_orfs['read_support'])

    if  score > max_score:
        print(" [ASSIGNED] %s: Annotated features only" % gene_ids) 
        max_score = score

        # for now, just keep track of feature assignments for annotated
        # genes and don't worry about the ORFs themselves...
        sl = no_orfs['sl']
        polya = no_orfs['polya']
    else:
        if max_score > 0:
            print(" [ASSIGNED] %s: Putative novel ORF detected in range: %d - %d" % (gene_ids, 
                                                                                     region_start,
                                                                                     region_end)) 
    # If no valid configurations were encountered (e.g. SL/Poly(A) sites
    # entirely in wrong orientation) then don't add any GFF entries
    if max_score == 0:
        return []

    # Create GFF entries
    if strand == 1:
        sl_gff_entry = build_gff_utr_entry(sl, right_gene, chr_id)
        polya_gff_entry = build_gff_utr_entry(polya, left_gene, chr_id)
    else:
        sl_gff_entry = build_gff_utr_entry(sl, left_gene, chr_id)
        polya_gff_entry = build_gff_utr_entry(polya, right_gene, chr_id)

    # Return GFF entries for each feature
    return [sl_gff_entry, polya_gff_entry]

def select_optimal_features(sl_sites, polya_sites, strand):
    """
    Finds the optimal SL/Poly(A) primary site assignments for a given region,
    without considering the effect of novel ORFs.
    """
    max_read_support = 0
    sl = None
    polya = None

    # Iterate over all combinations of filtered SL/Poly(A) sites
    # TODO: handle case where only one type of feature exists...
    for sl_site in sl_sites:
        for polya_site in polya_sites:
            # TODO: Decide how to assign features when this is the only
            # configuration.

            #  If SL/Poly(A) are in wrong in wrong orientation,
            #  which would lead to negative-sized intergenic
            #  regions, skip.
            if ((strand ==  1 and polya_site.location.start >= sl_site.location.start) or 
                (strand == -1 and polya_site.location.start <= sl_site.location.start)):
                continue

            # compute read support for features
            read_support = (int(sl_site.qualifiers['score'][0]) + 
                            int(polya_site.qualifiers['score'][0]))

            # keep feature pair if highest scoring encountered
            if read_support > max_read_support:
                max_read_support = read_support
                sl = sl_site
                polya = polya_site

    return {'sl': sl, 'polya': polya, 'read_support': max_read_support}

def build_gff_utr_entry(feature, gene, chr_id):
    """Takes a single SeqFeature corresponding to either a primary SL or
    polyadenylation site, as well as its associated gene, and return a GFF 5'
    or 3'UTR entry for the pair.

    Parameters
    ----------
    feature: SeqFeature
        A single-coordinate SeqFeature representing an SL or Poly(A) site.
    gene: SeqFeature
        Gene for which the feature belongs to.
    chr_id: str
        Chromosome identifier.
        
    Return
    ------
    str: A GFF entry of type five_prime_UTR or three_prime_UTR.
    """
    # 3'UTR
    if feature.type == 'polyA_site':
        entry_type = 'three_prime_UTR'

        # Positive strand
        if gene.strand == 1:
            start = gene.location.end + 1
            end = feature.location.start
        else:
            # Negative strand
            start = feature.location.end + 1
            end = gene.location.start
    else:
        # 5'UTR
        entry_type = 'five_prime_UTR'

        # Positive strand
        if gene.strand == 1:
            start = feature.location.end + 1
            end = gene.location.start
        else:
            # Negative strand
            start = gene.location.end + 1
            end = feature.location.start

    # Description
    # ID=TcCLB.511911.98_5utr;Name=TcCLB.511911.98;description=hypothetical+protein,+conserved
    short_type = '5utr' if entry_type == 'five_prime_UTR' else '3utr'

    desc = 'ID=%s_%s;Name=%s;description=%s' % (gene.id, short_type, gene.id,
                                                gene.qualifiers['description'][0])

    # GFF parts
    strand = '+' if gene.strand == 1 else '-'
    score = feature.qualifiers['score'][0]

    parts = [chr_id, 'El-Sayed', entry_type, str(start), str(end), score, strand, '.', desc]

    return "\t".join(parts) 

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

def filter_features(features, start, end, window_size=10, max_features=3,
                    min_ratio_alt_to_prim=0.01):
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
    min_ratio_alt_to_prim: float
        The minimum ratio of alternative site to primary site reads required
        for an alternative site to be considered valid. For example, if the
        primary site was found to have 1000 reads mapped to it, and the ratio
        is 0.01, then any secondary sites must have at least 0.1 * 1000 = 10
        reads to avoid being filtered out.
    """
    # Convert features to a 1d array representation and filter
    feature_arr = features_to_1d_array(features, start, end)

    # Smooth out everything except for local peaks
    feature_arr = max_filter(feature_arr, width=window_size)

    # Rank by score and keep only the top N highest-scoring features
    rank_cutoff = len(feature_arr) - max_features
    feature_arr[feature_arr.argsort().argsort() < rank_cutoff] = 0

    # Drop any non-primary sites with much lower read coverage compared
    # to the primary site
    min_num_reads = max(feature_arr) * min_ratio_alt_to_prim
    feature_arr[feature_arr < min_num_reads] = 0

    # Get the indices of coverage peaks
    indices = np.arange(len(feature_arr))[feature_arr != 0] + start

    return [x for x in features if x.location.start in indices]

