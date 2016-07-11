"""
GeneStructureAnalyzer class definition
Keith Hughitt <khughitt@umd.edu>
July, 2016
"""
import os
import pandas
import re
import numpy as np
from Bio import Seq, SeqIO, BiopythonWarning
from BCBio import GFF
from .orfs import find_orfs
from .util import max_filter
from .plots import plot_genome_features
from .io import build_gff_utr_entry, create_extended_gff, create_summary_csv_files, load_gff

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
                 min_protein_length, outdir):
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
        outdir: str
            Location to save outputs to
        """
        self.min_protein_length = min_protein_length
        self.outdir = outdir
        
        # Load genome sequence
        logging.info("- Loading Genome sequence...")
        self.genome_sequence = {x.id:x for x in SeqIO.parse(fasta, 'fasta')}

        # Load genome annotations
        logging.info("- Loading Genome annotations...")
        self.genome_annotations = load_gff(gff)

        # Load RNA-Seq coverage map
        logging.info("- Loading RNA-Seq coverage data...")
        df = pandas.read_table(coverage, names=['chr_id', 'loc', 'coverage'])

        self.rnaseq_coverage = {}
        for chr_id in self.genome_sequence:
            self.rnaseq_coverage[chr_id] = df[df.chr_id == chr_id].coverage

        # Load SL and Poly(A) sites
        logging.info("- Loading SL/Poly(A) site data...")
        self.sl_sites = load_gff(sl_gff)
        self.polya_sites = load_gff(polya_gff)

        # Variables for keeping track of our position and count as we move
        # across the chromosomes in the genome

        # inter-CDS regions
        self.inter_cds_regions = []

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
        logging.info("- Detecting UTR boundaries...")
        self.find_utr_boundaries()

        # generate genome coverage and novel ORF plots
        # Determine CDS and inter-CDS boundaries
        # cds_regions = get_cds_regions(self.genome_annotations)
        # inter_cds_regions = get_inter_cds_regions(self.genome_annotations)
        # plot_genome_features(genome_sequence, cds_regions, rnaseq_coverage,
        #                      plot_type)

        # Write output
        logging.info("- Saving results...")
        out_dir = os.path.join('output', outdir])
        if not os.path.exists(out_dir):
            os.makedirs(out_dir, mode=0o755)

        create_extended_gff(out_dir, gff, gff_entries)
        create_summary_csv_files(out_dir, utr5_csv_entries, utr3_csv_entries)

    def find_utr_boundaries():
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
                left_gid = left_gene.id if self.left_gene != None else 'None'
                right_gid = self.right_gene.id
                self.gene_ids = "%s - %s" % (left_gid, right_gid)

                # attempt to detect UTR boundaries in the region
                self.process_intercds_region()

        # Print message when finished
        txt = "# Total inter-CDS regions scanned: %d (%d assigned)"
        logging.info(txt % (self.total_counter, self.assigned_counter))

    def process_inter_cds_region(self):
        """Analyzes a single inter-CDS region"""
        # Skip over snoRNAs, etc. that are nested inside of other genes
        # For example: TcCLB TcChr22-2 179,000:180,000
        if self.end <= self.start:
            self.next(" [SKIPPING] %s: Nested genes" % self.gene_ids)
            return

        # Get sequence record for the range
        seq_record = genome_sequence[self.chr_id][self.start:self.end]

        # Retrieve SL and Poly(A) sites in region
        inter_cds_sl_sites    = self.get_features_by_range(self.sl_sites[self.chr_id])
        inter_cds_polya_sites = self.get_features_by_range(self.polya_sites[self.chr_id])

        # If no sites found, skip
        if len(inter_cds_sl_sites) == 0 and len(inter_cds_polya_sites) == 0:
            self.next(" [SKIPPING] %s: No features detected" % self.gene_ids)
            return

        # DEBUGGING
        # if right_gene.id == 'TcCLB.508727.60':
        #     import pdb; pdb.set_trace();

        # TODO: deal with cases where only sl site or polya sites were found 
        # if len(filtered_sl_sites) == 0:
        #     self.next(" [SKIPPING] %s: No SL sites detected" % self.gene_ids)
        #     return
        # elif len(filtered_polya_sites) == 0:
        #     self.next(" [SKIPPING] %s: No Poly(A) sites detected" % self.gene_ids)
        #     return

        # Add CDS to relevant list based on strand
        if self.strand is None:
            # Left-most gene
            self.next(" [SKIPPING] %s: Start of chromosome" % self.gene_ids)
            return
        elif self.strand != self.right_gene.location.strand:
            # Transcription Switch Sites (TSSs)
            self.next(" [SKIPPING] %s: Inter-PTU region" % self.gene_ids)
            return
        else:
            # Within PTU features 
            inter_cds = self.create_inter_cds_region(inter_cds_sl_sites,
                                                     inter_cds_polya_sites,
                                                     seq_record) 

        # if features were assigned, add to list and continue
        if inter_cds is not None:
            self.inter_cds_regions.append(inter_cds)
            self.assigned_counter = assigned_counter + 1

        # update gene, start counter and strand
        self.next()

    def create_inter_cds_region(sl_sites, polya_sites, seq_record)
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
        filtered_sl_sites    = filter_features(sl_sites)
        filtered_polya_sites = filter_features(polya_sites)

        # Find the optimal ORF-containing configuration, if ORFs can be
        # detected in the inter-CDS region
        orf_result = select_optimal_features_with_orf(filtered_sl_sites,
                                                      filtered_polya_sites, 
                                                      seq_record)
        max_score = orf_result['max_score']
        sl = orf_result['sl']
        polya = orf_result['polya']

        # Also compute the score for the best arrangement without any novel ORFs
        no_orfs = select_optimal_features(filtered_sl_sites, filtered_polya_sites, strand)

        # total number of features assigned for configuration
        total_assigned = bool(no_orfs['sl']) + bool(no_orfs['polya'])

        # compute score
        no_orf_score = (10E6  * total_assigned +
                        10E0  * no_orfs['read_support'])

        # check to see if configuration without a novel ORF performs best
        if no_orf_score > max_score:
            logging.info(" [ASSIGNED] %s: Annotated features only" % self.gene_ids) 
            max_score = no_orf_score
            sl = no_orfs['sl']
            polya = no_orfs['polya']
        elif max_score > 0:
            logging.info(" [ASSIGNED] %s: Putative novel ORF detected in range: %d - %d" % (self.gene_ids, 
                                                                                        self.start,
        # If no valid configurations were encountered (e.g. SL/Poly(A) sites
        # entirely in wrong orientation) then don't add any GFF entries
        if max_score == 0:
            return None

        # Create InterCDSRegion instances
        if strand == 1:
            # positive strand (3'UTR before 5'UTR)
            return InterCDSRegion(
                ThreePrimeUTR(left_gene, polya_sites, polya, self.chr_id, strand),
                FivePrimeUTR(right_gene, sl_sites, sl, self.chr_id, strand)
            )
        else:
            # negative strand (5'UTR before 3'UTR)
            return InterCDSRegion(
                FivePrimeUTR(left_gene, sl_sites, sl, self.chr_id, strand),
                ThreePrimeUTR(right_gene, polya_sites, polya, self.chr_id, strand)
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
            lhs_sl_sites    = [sl for sl in sl_sites if sl.location.start < orf_start]
            lhs_polya_sites = [polya for polya in polya_sites if polya.location.start < orf_start]
            lhs = select_optimal_features(lhs_sl_sites, lhs_polya_sites)

            # compute score for region to the right of ORF
            rhs_sl_sites    = [sl for sl in sl_sites if sl.location.start > orf_end]
            rhs_polya_sites = [polya for polya in polya_sites if polya.location.start > orf_end]
            rhs = select_optimal_features(rhs_sl_sites, rhs_polya_sites)

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

                if strand == 1:
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
                if ((self.strand ==  1 and polya_site.location.start >= sl_site.location.start) or 
                    (self.strand == -1 and polya_site.location.start <= sl_site.location.start)):
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

    def get_features_by_range(annotations):
        """Returns all annotation features which fall within the specified range"""
        features = []

        for feature in annotations.features:
            if feature.location.start > self.start and feature.location.end < self.end:
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
        feature_arr = features_to_1d_array(features, self.start, self.end)

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

class Gene(object):
    """Base Gene class"""
    def __init__(strand):
        pass

class UpstreamGene(Gene):
    """Upstream gene"""
    def __init__(strand):
        pass

class DownstreamGene(Gene):
    """Downstream gene"""
    def __init__(strand):
        pass

class PrimarySite(object):
    def __init__(self, gene, chr_id):
        """Base constructor"""
        self.gene = gene
        self.chr_id = chr_id

class UntranslatedRegion(object):
    def __init__(self, gene, features, primary_feature, chr_id, strand):
        self.gene = gene
        self.features = features
        self.primary_feature = primary_feature
        self.chr_id = chr_id
        self.strand = strand

class FivePrimeUTR(UntranslatedRegion):
    """FivePrimeUTR class definition"""
    def __init__(self, gene, sl_sites, primary_sl, chr_id, strand):
        super().__init__(gene, sl_sites, primary_sl, chr_id, strand)

class ThreePrimeUTR(UntranslatedRegion):
    """ThreePrimeUTR class definition"""
    def __init__(self):
        super().__init__(gene, polya_sites, primary_polya, chr_id, strand)

class InterCDSRegion(object):
    """InterCDSRegion class definition"""
    def __init__(self, left_feature, right_feature):
        self.left_feature = left_feature
        self.right_feature = right_feature
    def to_gff(self):
        """Returns a GFF representation of the UTR boundaries"""
        # Indices in GFF row list
        TYPE_IDX  = 2
        START_IDX = 3
        END_IDX   = 4
        SCORE_IDX = 5
        DESC_IDX  = 8

        # Generate row in CSV output
        for utr in [self.left_feature, self.right_feature]:
            utr_start = utr[START_IDX]
            utr_end   = utr[END_IDX]

            # TODO: look into edge case where SL site appears to be directly
            # adjacenct to the CDS (ex: TcCLB.509233.50)
            if utr_end - utr_start == 0:
                return

            gene_id = re.search('Name=([^;]*)', utr[DESC_IDX]).groups()[0]

            # Get GC- and CT-richness
            seq = str(genome_sequence[self.chr_id].seq[utr_start:utr_end])
            gc_richness = round((seq.count('G') + seq.count('C')) / len(seq), 3)
            ct_richness = round((seq.count('C') + seq.count('T')) / len(seq), 3)

            if utr[TYPE_IDX] == 'five_prime_UTR':
                utr5_csv_entries.append([gene_id, utr_end - utr_start,
                                        utr[SCORE_IDX], gc_richness, 
                                        ct_richness])
            else:
                utr3_csv_entries.append([gene_id, utr_end - utr_start,
                                        utr[SCORE_IDX], gc_richness, 
                                        ct_richness])

class PrimaryTransSplicingSite(PrimarySite):
    def __init__(self, gene, chr_id):
        """Creates a new PrimaryTransSplicingSite instance"""
        super().__init__(gene, chr_id)

        # 5'UTR
        entry_type = 'five_prime_UTR'

        # Positive strand
        if gene.strand == 1:
            start = feature.location.end
            end = gene.location.start
        else:
            # Negative strand
            start = gene.location.end
            end = feature.location.start

class PrimaryPolyAdenylationSite(PrimarySite):
    def __init__(self, gene, chr_id):
        """Creates a new PrimaryPolyAdenylationSite instance"""
        super().__init__(gene, chr_id)

        self.entry_type = 'three_prime_UTR'

        # Positive strand
        if gene.strand == 1:
            start = gene.location.end
            end = feature.location.start
        else:
            # Negative strand
            start = feature.location.end
            end = gene.location.start

