"""
GeneStructureAnalyzer class definition
Keith Hughitt <khughitt@umd.edu>
July, 2016
"""
import os
import sys
import logging
import pandas
import re
import numpy as np
from Bio import Seq, SeqIO
from BCBio import GFF
from .orfs import find_orfs
from .util import max_filter
from .plots import plot_genome_features
from .io import build_gff_utr_entry, create_extended_gff, create_summary_csv_files, load_gff, create_alt_site_csv_files, create_polypyrimidine_tract_csv, create_intercds_csv

class GeneStructureAnalyzer(object):
    """
    GeneStructureAnalyzer class definition

    Class for orchestrating the analysis of structural features in
    Trypanosomatid genomes.

    This class combines multiple sources of data and annotations, and attempt
    to predict primary and alternative trans-splicing acceptors sites, primary
    and alternative polyadenylation sites, and the corresponding 5' and 3'UTR
    boundaries where evidence is available.

    Additional features such as the splice site acceptor dinucleotide and
    polypyrimidine tract are also identified.

    In order to optimally predict UTR boundaries, predicted trans-splicing and
    polyadenylation are considered in the inter-CDS regions between pairs of
    annotated genes. This ensures that the UTR boundary of one gene will not
    overlap that of its neighbor.

    For simplicity, this script currently does not analyzes UTRs at the ends
    of chromosomes, or in transcript switch sites (TSSs).
    """
    def __init__(self, fasta, gff, coverage, sl_gff, polya_gff,
                 min_protein_length, max_intercds_length,
                 polypyrimidine_window, outdir):
        """Creates a new GeneStructureAnalyzer instance

        Arguments
        ---------
        fasta: str
            Genome sequence filepath
        gff: str
            Genome annotations filepath
        coverage: str
            RNA-Seq read coverage filepath
        sl_gff: str
            filepath to GFF containing predicted SL sites
        polya_gff: str
            filepath to GFF containing predicted Poly(A) sites
        min_protein_length: int
            Minimum length in amino acids to allow for novel ORFs
        max_intercds_length: int
            Maximum distance between two genes considered to be adjacent.
        poly_pyrimidine_window: int
            Window size to scan for polypyrimidine tracts 
        outdir: str
            Location to save outputs to
        """
        self.min_protein_length = min_protein_length
        self.max_intercds_length = max_intercds_length
        self.polypyrimidine_window = polypyrimidine_window
        self.outdir = outdir

        # Setup loggers
        self.setup_logger()

        # Load genome sequence
        logging.info("Loading Genome sequence...")
        self.genome_sequence = {x.id:x for x in SeqIO.parse(fasta, 'fasta')}

        # Load genome annotations
        logging.info("Loading Genome annotations...")
        self.genome_annotations = load_gff(gff)

        # Load RNA-Seq coverage map
        logging.info("Loading RNA-Seq coverage data...")
        df = pandas.read_table(coverage, names=['chr_id', 'loc', 'coverage'])

        self.rnaseq_coverage = {}
        for chr_id in self.genome_sequence:
            #pylint: disable=maybe-no-member
            self.rnaseq_coverage[chr_id] = df[df.chr_id == chr_id].coverage

        # Load SL and Poly(A) sites
        logging.info("Loading SL/Poly(A) site data...")
        self.sl_sites = load_gff(sl_gff)
        self.polya_sites = load_gff(polya_gff)

        # Variables for keeping track of our position and count as we move
        # across the chromosomes in the genome

        # counters
        self.total_counter = 0
        self.assigned_counter = 0

        # position
        self.chr_id = None
        self.strand = None
        self.start = 0
        self.end = None

        # genes
        self.gene_ids = None
        self.left_gene = None
        self.right_gene = None

        # inter-CDS regions
        self.inter_cds_regions = []

        # Optimize SL / polyadenylation primary and alternative site selection for
        # each annotated CDS, assigning novel ORFs where appropriate
        logging.info("Detecting UTR boundaries...")
        self.find_utr_boundaries()

        # generate genome coverage and novel ORF plots
        # Determine CDS and inter-CDS boundaries
        # cds_regions = get_cds_regions(self.genome_annotations)
        # inter_cds_regions = get_inter_cds_regions(self.genome_annotations)
        # plot_genome_features(genome_sequence, cds_regions, rnaseq_coverage,
        #                      plot_type)

        # Get a list of GFF/CSV entries to be written to output
        gff_entries = []

        intercds_csv_entries = []

        utr5_csv_entries = []
        utr3_csv_entries = []

        utr5_csv_all_sites = []
        utr3_csv_all_sites = []

        polypyrimidine_tract_entries = []

        # Iterate over inter-CDS regions and generate GFF/CSV output entries
        for region in self.inter_cds_regions:
            # GFF entry
            gff_entries = gff_entries + region.to_gff()

            # 5'UTR CSV entries
            utr5 = region.get_utr5()
            if utr5 is not None:
                utr5_csv_entries.append(utr5.to_primary_utr_csv())
                utr5_csv_all_sites = utr5_csv_all_sites + utr5.all_sites_csv()

            # 3'UTR CSV entries
            utr3 = region.get_utr3()
            if utr3 is not None:
                utr3_csv_entries.append(utr3.to_primary_utr_csv())
                utr3_csv_all_sites = utr3_csv_all_sites + utr3.all_sites_csv()

            # Inter-CDS CSV entry
            if utr5 is not None and utr3 is not None:
                intercds_csv_entries.append(region.to_csv())

            # Polypyrimidine tract entries
            if region.polypyrimidine_tract is not None:
                polypyrimidine_tract_entries.append(region.polypyrimidine_tract.to_csv())

        # Drop any empty entries
        utr5_csv_entries = [x for x in utr5_csv_entries if x is not None]
        utr3_csv_entries = [x for x in utr3_csv_entries if x is not None]
        utr5_csv_all_sites = [x for x in utr5_csv_all_sites if x is not None]
        utr3_csv_all_sites = [x for x in utr3_csv_all_sites if x is not None]

        # Write output
        logging.info("Saving results...")

        out_dir = os.path.join('output', outdir)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir, mode=0o755)

        create_extended_gff(out_dir, gff, gff_entries)
        create_intercds_csv(out_dir, intercds_csv_entries)
        create_summary_csv_files(out_dir, utr5_csv_entries, utr3_csv_entries)
        create_alt_site_csv_files(out_dir, utr5_csv_all_sites,
                                  utr3_csv_all_sites)
        create_polypyrimidine_tract_csv(out_dir, polypyrimidine_tract_entries)

    def find_utr_boundaries(self):
        """
        Simultaneously detects primray and alternative  SL and Poly(A) sites for a
        set of two neighboring CDS's and detecting novel ORFs where appropriate.

        The goal of this function is to attempt to determine the most likely
        primary SL/Poly(A) sites for the set of known CDS's for which we have
        sites detected for and, in cases where evidence exists for possible
        novel ORFs (e.g. additional SL and Poly(A) sites surrounding an ORF in
        between two CDS's), to extend the current genome annotation to include
        these ORFs. This also serves to provide us with better estimates about
        the true 5' and 3'UTR lengths for parasite genes.

        Further, this function also keeps track of alternative (non-primary)
        trans-splicing and polyadenylation sites.

        Depending on the amount of detected SL and Poly(A) sites between a set
        of known CDS's, there may be multiple possible configurations of
        feature assignments to both the known CDS's and various possible novel
        ORFs.

        In order to attempt to select an "optimal" combination of novel ORF and
        SL / Poly(A) site assignments, we will use a scoring function which
        takes into account different types of support.

        For each possible configuration, we will compute a score defined by:

            10E6   * # Assigned features +
            10E0   * sum(reads mapped to each feature) +
            10E-6  * sum(read density spanning the ORF) +
           -10E6   * if ORF included

        This optimization criterion will favor the selection of configurations
        where the maximal number of CDS/ORFs are assigned SL and Poly(A) sites,
        and will attempt to choose the sites with the most support (largest
        number of reads supporting those sites). If there are multiple ORFs which
        can be chosen for a given SL/Poly(A) combination, the ORF which captures
        the most read density (usually the longest ORF), will be selected.
        Finally, if two configurations score equally well, but one uses only the
        existing annotations, and another assigns features to a novel ORF, the
        configuration using only the existing annotations will be favored.
        """
        # Iterate over chromosomes
        for chr_id, chromosome in self.genome_annotations.items():
            logging.info("=======" + chr_id + "=======")
            self.chr_id = chr_id

            # filter out everything except for genes and order by position
            # along chromosome
            genes = [x for x in chromosome.features if x.type == 'gene']
            genes.sort(key=lambda x: x.location.start)

            # iterate over genes
            for right_gene in genes:
                self.right_gene = right_gene

                # Determine location for the region up to start of the gene
                self.end = int(right_gene.location.start)
                self.total_counter = self.total_counter + 1

                # Update gene range description (used for logging)
                left_gid = self.left_gene.id if self.left_gene != None else 'None'
                right_gid = self.right_gene.id
                self.gene_ids = "%s - %s" % (left_gid, right_gid)

                # attempt to detect UTR boundaries in the region
                self.process_inter_cds_region()

        # Print message when finished
        txt = "Total inter-CDS regions scanned: %d (%d assigned)"
        logging.info(txt % (self.total_counter, self.assigned_counter))

    def process_inter_cds_region(self):
        """Analyzes a single inter-CDS region"""
        # Skip over rare snoRNAs, etc. that are nested inside of other genes
        # For example: TcCLB TcChr22-2 179,000:180,000
        if self.end <= self.start:
            self.next("[SKIPPING] %s: Nested genes" % self.gene_ids)
            return

        # Skip very large inter-CDS regions (usually these correspond to
        # inter-PTU regions for PTUs on the same strand)
        if self.end - self.start > self.max_intercds_length:
            self.next("[SKIPPING] %s: region > 30,000nt wide." % self.gene_ids)
            return

        # Get sequence record for the range
        seq_record = self.genome_sequence[self.chr_id][self.start:self.end]

        # Retrieve SL and Poly(A) sites in region
        inter_cds_sl_sites    = self.get_features_by_range(self.sl_sites[self.chr_id])
        inter_cds_polya_sites = self.get_features_by_range(self.polya_sites[self.chr_id])

        # If no sites found, skip
        if len(inter_cds_sl_sites) == 0 and len(inter_cds_polya_sites) == 0:
            self.next("[SKIPPING] %s: No features detected" % self.gene_ids)
            return

        # Add CDS to relevant list based on strand
        if self.strand is None:
            # Left-most gene
            if self.right_gene.location.strand == 1 and len(inter_cds_sl_sites) > 0:
                # Positive strand with sites detected
                logging.info("[ASSIGNED] %s: Start of chromosome" % self.gene_ids)
                inter_cds = self.get_sl_only_inter_cds(inter_cds_sl_sites)
            elif self.right_gene.location.strand == -1 and len(inter_cds_polya_sites) > 0:
                # Negative strand with sites detected
                logging.info("[ASSIGNED] %s: Start of chromosome" % self.gene_ids)
                inter_cds = self.get_polya_only_inter_cds(inter_cds_polya_sites)
            else:
                # No sites detected for left-most gene
                self.next("[SKIPPING] %s: No features detected" % self.gene_ids)
                return
        elif self.strand != self.right_gene.location.strand:
            # Transcription Switch Sites (TSSs)
            self.next("[SKIPPING] %s: Inter-PTU region" % self.gene_ids)
            return
        elif len(inter_cds_sl_sites) == 0:
            # Only Poly(A) sites detected for region
            logging.info("[ASSIGNED] %s: 3'UTR only" % self.gene_ids)
            inter_cds = self.get_polya_only_inter_cds(inter_cds_polya_sites)
        elif len(inter_cds_polya_sites) == 0:
            # Only SL sites detected for region
            logging.info("[ASSIGNED] %s: 5'UTR only" % self.gene_ids)
            inter_cds = self.get_sl_only_inter_cds(inter_cds_sl_sites)
        else:
            # Within PTU features
            inter_cds = self.create_inter_cds_region(inter_cds_sl_sites,
                                                     inter_cds_polya_sites,
                                                     seq_record)

        # if features were assigned, add to list and continue
        if inter_cds is not None:
            self.inter_cds_regions.append(inter_cds)
            self.assigned_counter = self.assigned_counter + 1

        # update gene, start counter and strand
        self.next()

    def create_inter_cds_region(self, sl_sites, polya_sites, seq_record):
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
        # Filter sites down to ~3 optimal SL's / Poly(A)'s
        filtered_sl_sites    = self.filter_features(sl_sites)
        filtered_polya_sites = self.filter_features(polya_sites)

        # Find the optimal ORF-containing configuration, if ORFs can be
        # detected in the inter-CDS region
        orf_result = self.select_optimal_features_with_orf(filtered_sl_sites,
                                                           filtered_polya_sites,
                                                           seq_record)
        max_score = orf_result['max_score']
        sl = orf_result['sl']
        polya = orf_result['polya']

        # Also compute the score for the best arrangement without any novel ORFs
        no_orfs = self.select_optimal_features(filtered_sl_sites,
                                               filtered_polya_sites)

        # total number of features assigned for configuration
        total_assigned = bool(no_orfs['sl']) + bool(no_orfs['polya'])

        # compute score
        no_orf_score = (10E6  * total_assigned +
                        10E0  * no_orfs['read_support'])

        # check to see if configuration without a novel ORF performs best
        if no_orf_score > max_score:
            logging.info("[ASSIGNED] %s: Annotated features only" % self.gene_ids)
            max_score = no_orf_score
            sl = no_orfs['sl']
            polya = no_orfs['polya']
        elif max_score > 0:
            msg = "[ASSIGNED] %s: Putative novel ORF detected in range: %d - %d"
            logging.info(msg % (self.gene_ids, self.start, self.end))

        # If no valid configurations were encountered (e.g. SL/Poly(A) sites
        # entirely in wrong orientation) then don't add any GFF entries
        if max_score == 0:
            return None

        # Get alternative trans-splicing and polyadenylation sites
        alt_sl_sites = [x for x in sl_sites if x.id != sl.id]
        alt_polya_sites = [x for x in polya_sites if x.id != polya.id]

        # Create InterCDSRegion instances
        if self.strand == 1:
            # positive strand (3'UTR before 5'UTR)
            return InterCDSRegion(
                ThreePrimeUTR(self.left_gene, polya, alt_polya_sites, self.strand,
                              self.chr_id, self.genome_sequence),
                FivePrimeUTR(self.right_gene, sl, alt_sl_sites, self.strand,
                             self.chr_id, self.genome_sequence),
                self.chr_id, self.genome_sequence, self.polypyrimidine_window
            )
        else:
            # negative strand (5'UTR before 3'UTR)
            return InterCDSRegion(
                FivePrimeUTR(self.left_gene, sl, alt_sl_sites, self.strand,
                             self.chr_id, self.genome_sequence),
                ThreePrimeUTR(self.right_gene, polya, alt_polya_sites, self.strand,
                              self.chr_id, self.genome_sequence),
                self.chr_id, self.genome_sequence, self.polypyrimidine_window
            )

    def get_sl_only_inter_cds(self, sl_sites):
        """Returns an InterCDS instance for a region with only SL sites"""
        if self.strand == 1:
            # positive strand
            primary_sl = self.get_max_feature(sl_sites, location_preference='right')
            alt_sl_sites = [x for x in sl_sites if x.id != primary_sl.id]

            return InterCDSRegion(
                None,
                FivePrimeUTR(self.right_gene, primary_sl, alt_sl_sites,
                             self.strand, self.chr_id, self.genome_sequence),
                self.chr_id, self.genome_sequence, self.polypyrimidine_window
            )
        else:
            # negative strand
            primary_sl = self.get_max_feature(sl_sites, location_preference='left')
            alt_sl_sites = [x for x in sl_sites if x.id != primary_sl.id]

            return InterCDSRegion(
                FivePrimeUTR(self.left_gene, primary_sl, alt_sl_sites,
                             self.strand, self.chr_id, self.genome_sequence),
                None,
                self.chr_id, self.genome_sequence, self.polypyrimidine_window
            )

    def get_polya_only_inter_cds(self, polya_sites):
        """Returns an InterCDS instance for a region with only polya sites"""
        if self.strand == 1:
            # positive strand
            primary_polya = self.get_max_feature(polya_sites, 
                                                 location_preference='left')
            alt_polya_sites = [x for x in polya_sites if x.id != primary_polya.id]

            return InterCDSRegion(
                ThreePrimeUTR(self.left_gene, primary_polya, alt_polya_sites,
                              self.strand, self.chr_id, self.genome_sequence),
                None,
                self.chr_id, self.genome_sequence, self.polypyrimidine_window
            )
        else:
            # negative strand
            primary_polya = self.get_max_feature(polya_sites,
                                                 location_preference='right')
            alt_polya_sites = [x for x in polya_sites if x.id != primary_polya.id]

            return InterCDSRegion(
                None,
                ThreePrimeUTR(self.right_gene, primary_polya, alt_polya_sites,
                              self.strand, self.chr_id, self.genome_sequence),
                self.chr_id, self.genome_sequence, self.polypyrimidine_window
            )


    def next(self, msg=''):
        """Go to next inter-CDS region"""
        if len(msg) > 0:
            logging.info(msg)
        self.left_gene = self.right_gene
        self.start = int(self.left_gene.location.end)
        self.strand = self.left_gene.location.strand

    def select_optimal_features_with_orf(self, sl_sites, polya_sites, seq_record):
        """
        Find the ORF-containing inter-CDS region configuration which scores the
        best.

        Args
        ----
        sl_sites: list
            List of SL sites to use when considering presence of ORFs
        polya_sites: list
            List of polyadenylation sites to use when considering presence of ORFs
        seq_record: Bio.SeqRecord.SeqRecord
            Sequence of the inter-CDS region.

        Returns
        -------
        out: tuple
            Ordered pair containing <F9><F9>
        """
        # Find ORFs in each of the six possible reading frames
        # that are at least the specified length
        orfs = find_orfs(seq_record.seq, self.min_protein_length, self.strand)

        coverage = self.rnaseq_coverage[self.chr_id]

        max_score = 0
        max_score_orf = None
        sl = None
        polya = None

        # Iterate over any novel ORFs and score configurations
        for orf in orfs:
            # ORF boundaries relevant to chromosome start
            orf_start = self.start + orf[0]
            orf_end   = self.start + orf[1]

            # compute density of reads covering ORF
            orf_coverage = sum(coverage[orf_start:orf_end])

            # compute score for region to the left of ORF
            lhs_sl_sites    = [sl for sl in sl_sites if sl.location.end < orf_start]
            lhs_polya_sites = [polya for polya in polya_sites if polya.location.end < orf_start]
            lhs = self.select_optimal_features(lhs_sl_sites, lhs_polya_sites)

            # compute score for region to the right of ORF
            rhs_sl_sites    = [sl for sl in sl_sites if sl.location.end > orf_end]
            rhs_polya_sites = [polya for polya in polya_sites if polya.location.end > orf_end]
            rhs = self.select_optimal_features(rhs_sl_sites, rhs_polya_sites)

            # total number of features assigned for configuration
            total_assigned = bool(lhs['sl']) + bool(lhs['polya']) + bool(rhs['sl']) + bool(rhs['polya'])

            # compute score
            # one substracted from total_assigned to penalize inclusion of ORFs
            score = (10E6  * (total_assigned - 1) +
                     10E0  * lhs['read_support'] + rhs['read_support'] +
                     10E-6 * orf_coverage)

            # if including this ORF improves score, keep configuration
            if  score > max_score:
                max_score = score
                max_score_orf = orf

                if self.strand == 1:
                    sl = rhs['sl']
                    polya = lhs['polya']
                else:
                    sl = lhs['sl']
                    polya = rhs['polya']

        return {'sl': sl, 'polya': polya, 'max_score': max_score,
                'orf': max_score_orf}

    def select_optimal_features(self, sl_sites, polya_sites):
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
                if ((self.strand ==  1 and polya_site.location.end >= sl_site.location.end) or
                    (self.strand == -1 and polya_site.location.end <= sl_site.location.end)):
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

    def get_features_by_range(self, annotations):
        """Returns all annotation features which fall within the specified range"""
        features = []

        # retrieve all features that fall in the range of interest
        for feature in annotations.features:
            # skip chromosome entries
            if feature.type == 'chromosome':
                continue

            # check to see if SL/Poly(A) site falls in inter-CDS range
            if feature.location.end > self.start and feature.location.end < self.end:
                # get putative UTR sequence corresponding to feature and check
                # to make sure there are not a large number of N's
                if ((feature.type == 'trans_splice_site' and feature.strand == 1) or
                    (feature.type == 'polyA_site' and feature.strand == -1)):
                    # SL / negative strand Poly(A)
                    seq = self.genome_sequence[self.chr_id][feature.location.end:self.end + 1]
                else:
                    # Poly(A) / negative strand SL
                    seq = self.genome_sequence[self.chr_id][self.start:feature.location.end]

                # Determine number of ambiguous sequence positions
                num_ambiguous = seq.seq.count('N')

                # Include all sites with less than 100 N's
                if num_ambiguous < 100:
                    features.append(feature)

        return features

    def get_max_feature(self, features, location_preference='right'):
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
        positions = np.array([site.location.end for site in features])

        if location_preference == 'left':
            return np.array(features)[positions == min(positions)][0]
        else:
            return np.array(features)[positions == max(positions)][0]

    def features_to_1d_array(self, features, start, end):
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
            arr[feature.location.end - start] = int(feature.qualifiers.get('score')[0])

        return arr

    def filter_features(self, features, window_size=10, max_features=3,
                        min_ratio_alt_to_prim=0.01):
        """Filters a set of features (e.g. SL sites) to remove any low-support or
        closeby sites.

        Parameters
        ----------
        features: list
            A list of SeqFeature instances representing the locations of putative
            SL accepted sites and polyadenylation sites in the genome.
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
        feature_arr = self.features_to_1d_array(features, self.start, self.end)

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
        indices = np.arange(len(feature_arr))[feature_arr != 0] + self.start

        return [x for x in features if x.location.end in indices]
    
    def setup_logger(self):
        """Sets up global logger for gene structure analysis"""
        # get master logger and set log-level to DEBUG
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)

        # add a stream handler to have log output to STDOUT
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.DEBUG)

        # log format
        formatter = logging.Formatter('%(asctime)s (%(levelname)s) %(message)s')
        handler.setFormatter(formatter)

        root.addHandler(handler)

class InterCDSRegion(object):
    """InterCDSRegion class definition"""
    def __init__(self, left_feature, right_feature, chr_id, genome_seq,
                 polypyrimidine_window_size):
        self.left_feature = left_feature
        self.right_feature = right_feature

        self.intergenic_seq = None
        self.polypyrimidine_tract = None

        # Get intergenic sequence and detect polypyrimidine tract
        if left_feature is not None and right_feature is not None:
            start = left_feature.end
            end = right_feature.start - 1
            self.intergenic_seq = str(genome_seq[chr_id].seq[start:end])

            # detect polypyrimidine tract in region
            self.detect_polypyrimidine_tract(polypyrimidine_window_size)

    def to_gff(self):
        """Returns a GFF representation of the UTR boundaries"""
        entries = []

        for feature in [self.left_feature, self.right_feature, 
                        self.polypyrimidine_tract]:
            if feature is not None:
                entries.append(feature.to_gff())
        return entries

    def to_csv(self):
        """Returns a CSV representation of the inter-CDS region"""
        intercds_len = self.right_feature.end - self.left_feature.start + 1
        intergenic_len = self.right_feature.start - self.left_feature.end - 1

        return [self.left_feature.gene.id, self.right_feature.gene.id,
                str(self.left_feature.gene.strand), str(intercds_len), 
                str(intergenic_len), str(self.intergenic_seq)]

    def get_utr5(self):
        """Returns the 5'UTR"""
        if (self.left_feature is not None) and (self.left_feature.entry_type == 'five_prime_UTR'):
            return self.left_feature
        elif (self.right_feature is not None) and (self.right_feature.entry_type == 'five_prime_UTR'):
            return self.right_feature

    def get_utr3(self):
        """Returns the 3'UTR"""
        if (self.left_feature is not None) and (self.left_feature.entry_type == 'three_prime_UTR'):
            return self.left_feature
        elif (self.right_feature is not None) and (self.right_feature.entry_type == 'three_prime_UTR'):
            return self.right_feature

    def detect_polypyrimidine_tract(self, window_size):
        """Detects the longest upstream polypyrimidine tract for a 5'UTR"""
        # limit search window to specified size
        search_seq = self.intergenic_seq
        offset = 0

        if window_size < len(search_seq):
            if self.left_feature.strand == 1:
                # positive strand (right-side of intergenic region)
                offset = len(search_seq) - window_size
                search_seq = search_seq[-window_size:]
            else:
                # negative strand (left-side of intergenic region)
                search_seq = search_seq[:window_size]

        # split string by stretches of 3 or more purine and return longest match
        seq = max(re.split('[AGN]{3,}', search_seq), key=len)

        # strip any purines or N's at edge of match
        seq = seq.strip('AGN')

        # For very short intergenic regions, there may be no tracts detected
        if seq == '':
            self.polypyrimidine_tract = None
            return

        # get coordinates of match
        m = re.search(seq, search_seq)

        start =  self.left_feature.end + m.start() - 1 + offset
        end =  self.left_feature.end + m.end() - 1 + offset

        # distance to SL site
        # add one to account second base in AG dinucleotide
        utr5 = self.get_utr5()
        sl_dist = min(abs(start - utr5.primary_site_loc), 
                      abs(end - utr5.primary_site_loc)) + 1

        # distance to Poly(A) site
        utr3 = self.get_utr3()
        polya_dist = min(abs(start - utr3.primary_site_loc), 
                         abs(end - utr3.primary_site_loc))

        # store as a tuple with the start loc, end loc, and sequence
        self.polypyrimidine_tract = PolypyrimidineTract(utr5, start, end,
                                                        sl_dist, polya_dist, seq)

class UntranslatedRegion(object):
    def __init__(self, gene, primary_site, alternative_sites, strand, chr_id, genome_seq):
        self.gene = gene
        self.primary_site = primary_site
        self.alternative_sites = alternative_sites
        self.chr_id = chr_id
        self.strand = strand
        self.entry_type = None

        self.score = primary_site.qualifiers['score'][0]

        # set UTR boundaries and primary splice/polya site location
        start, end, site_loc = self.get_site_locations(primary_site)

        self.start = start
        self.end = end
        self.primary_site_loc = site_loc

        # Get UTR sequence
        self.seq = str(genome_seq[chr_id].seq[self.start - 1:self.end])

    def to_gff(self):
        """Returns a GFF representation of UTR"""
        # Description
        desc = 'ID=%s;Name=%s;description=%s;Parent=%s' % (self.id,
                                                           self.entry_type_short, 
                                                           self.entry_type_short, 
                                                           self.gene.id)
        # GFF parts
        strand = '+' if self.strand == 1 else '-'

        return "\t".join([self.chr_id, 'El-Sayed', self.entry_type, 
                          str(self.start), str(self.end), self.score, 
                          strand, '.', desc])

    def to_primary_utr_csv(self):
        """Returns a CSV representation of primary UTR boundaries"""
        # compute utr length
        utr_length = len(self.seq)

        # TODO: look into edge case where SL site appears to be directly
        # adjacenct to the CDS (ex: TcCLB.509233.50)
        if utr_length == 0:
            return

        # Get GC- and CT-richness
        gc_richness = round((self.seq.count('G') + self.seq.count('C')) /
                            len(self.seq), 3)
        ct_richness = round((self.seq.count('C') + self.seq.count('T')) /
                            len(self.seq), 3)

        return [self.gene.id, self.primary_site.id, utr_length, self.score, 
                gc_richness, ct_richness, self.seq]

    def all_sites_csv(self):
        """Returns a CSV representation of all trans-splicing /
        polyadenylation sites detected for a given feature."""
        csv_entries = []

        # compute utr length
        utr_length = len(self.seq)

        # read_support
        num_reads = self.primary_site.qualifiers['score'][0]

        csv_entries.append([self.gene.id, self.primary_site.id, 'primary',
                            num_reads, self.primary_site_loc, self.start, self.end, 
                            utr_length])

        # add secondary sites
        for site in self.alternative_sites:
            start, end, site_loc = self.get_site_locations(site)
            utr_length = end - start + 1
            num_reads = site.qualifiers['score'][0]

            csv_entries.append([self.gene.id, site.id, 'alternative',
                                num_reads, site_loc, start, end, utr_length])

        return csv_entries

class FivePrimeUTR(UntranslatedRegion):
    """FivePrimeUTR class definition"""
    def __init__(self, gene, primary_sl, alt_sites, strand, chr_id, genome_seq):
        super().__init__(gene, primary_sl, alt_sites, strand, chr_id, genome_seq)

        self.entry_type = 'five_prime_UTR'
        self.entry_type_short = '5utr'

        # Set ID
        self.id = "%s_%s" % (self.entry_type_short, self.gene.id)

    def get_site_locations(self, site):
        """
        Returns the UTR start and end locations, and the feature location.
       
        Note
        ----
        When parsing GFF files with BCBio, the single locations used for
        start and end coordinates in the GFF files get converted to a start
        position that is one less than what is in the GFF, and and end position
        equal to the GFF value.  Thus we always want to use the "site.end"
        position when determining where the splice site is actually located.
        """
        site_loc = site.location.end

        if self.gene.strand == 1:
            # Positive strand
            utr_start = site_loc + 1
            utr_end = self.gene.location.start
        else:
            # Negative strand
            utr_start = self.gene.location.end + 1
            utr_end = site_loc - 1

        return([utr_start, utr_end, site_loc])

class ThreePrimeUTR(UntranslatedRegion):
    """ThreePrimeUTR class definition"""
    def __init__(self, gene, primary_polya, alt_sites, strand, chr_id, genome_seq):
        super().__init__(gene, primary_polya, alt_sites, strand, chr_id, genome_seq)

        self.entry_type = 'three_prime_UTR'
        self.entry_type_short = '3utr'

        # Set ID
        self.id = "%s_%s" % (self.entry_type_short, self.gene.id)

    def get_site_locations(self, site):
        """Returns the UTR start and end locations, and the feature location"""
        site_loc = site.location.end

        if self.gene.strand == 1:
            # Positive strand
            utr_start = self.gene.location.end + 1
            utr_end = site_loc
        else:
            # Negative strand
            utr_start = site_loc + 1
            utr_end = self.gene.location.start

        return([utr_start, utr_end, site_loc])

class PolypyrimidineTract(object):
    """PolypyrimidineTract class definition"""
    def __init__(self, utr5, start, end, sl_dist, polya_dist, seq):
        """Creates a new PolypyrimidineTract instance"""
        self.utr5 = utr5
        self.start = start
        self.end = end
        self.sl_dist = sl_dist
        self.polya_dist = polya_dist
        self.seq = seq

    def to_csv(self):
        """Returns a list representation of the tract ready for csv output"""
        return [self.utr5.primary_site.id, str(self.start), str(self.end),
                str(self.sl_dist), str(self.polya_dist), str(self.seq)]

    def to_gff(self):
        """Returns a GFF representation of polypyrimidine tract"""
        # Description
        desc_template = 'ID=polypyrimidine_tract_%s;Name=polypyrimidine_tract;description=polypyrimidine_tract;size=%d;Parent=%s'
        desc = desc_template % (self.utr5.gene.id, self.end - self.start, self.utr5.id)

        # GFF parts
        strand = '+' if self.utr5.strand == 1 else '-'

        return "\t".join([self.utr5.chr_id, 'El-Sayed', 'polypyrimidine_tract',
                          str(self.start), str(self.end), '.',
                          strand, '.', desc])

