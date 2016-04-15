"""
Functions for retrieving features across a set of genome annotations.
"""
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

def find_primary_sites(annotations):
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
        10E1 * sum(read density spanning the ORF)

    This optimization criterion will favor the selection of configurations
    where the maximal number of CDS/ORFs are assigned SL and Poly(A) sites,
    and will attempt to choose the sites with the most support (largest
    number of reads supporting those sites). Finally, if there are multiple
    ORFs which can be chosen for a given SL/Poly(A) combination, the ORF
    which captures the most read density (usually the longest ORF), will be
    selected.
    """
    # Iterate over chromosomes
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

        # iterate over genes
        # for gene in genes:
        for right_gene in genes:
            # Determine location for the region up to start of the gene
            end = int(right_gene.location.start)

            # Skip over snoRNAs, etc. that are nested inside of other genes
            # For example: TcCLB TcChr22-2 179,000:180,000
            if end <= start:
                next

            # TODO: retrieve SL and Poly(A) sites in region
            
            # TODO: if no sites found, skip

            # TODO: filter sites down to ~3 optimal SL's / Poly(A)'s

            # Add CDS to relevant list based on strand
            if strand is None:
                # Left-most gene
                # TODO: assign SL/Poly(A) to left-most gene
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

