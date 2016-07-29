"""
I/O functions
"""
import os
import csv
import shutil
from argparse import ArgumentParser
from BCBio import GFF

def parse_args():
    """Parses input and returns arguments"""
    parser = ArgumentParser(description=('Determine gene structure and detect '
                                         'novel ORF using RNA-Seq data.'))

    # parser.add_argument('-b', '--bam', help='Input bam file', required=True)
    parser.add_argument('-c', '--coverage', 
                            help='Single nucleotide resolution genome coverage map',
                            required=True)
    parser.add_argument('-f', '--fasta', help='Input genome FASTA file',
                        required=True)
    parser.add_argument('-g', '--gff', help='Input genome GFF file', required=True)
    parser.add_argument('-s', '--sl-gff', help='Spliced leader site GFF')
    parser.add_argument('-p', '--polya-gff', help='Polyadenylation site GFF')
    parser.add_argument('-l', '--min-protein-length', type=int, default=30,
                        help=('Minimum size in amino acids allowed for '
                                'novel ORFs. (default=30)'))
    parser.add_argument('-t', '--plot-type', default='discrete',
                        help=('Type of colormap to use when generating '
                              'coverage plots. [discrete|continuous]'))
    parser.add_argument('outdir', help='Location to save results to',
                        metavar='OUTDIR')

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

def load_gff(gff):
    """Parses a single GFF file and returns a chromosome-indexed dict for
       that file.

    Arguments
    ---------
    gff: str
        Filepath to GFF

    Returns
    -------
    dict: A dictionary representation of the GFF entries, indexed by
            chromosome ID
    """
    annotations = {}

    if gff.endswith('.gz'):
        import gzip
        from io import TextIOWrapper
        fp = TextIOWrapper(gzip.open(gff))
    else:
        fp = open(gff)

    for entry in GFF.parse(fp):
        if len(entry.features) > 0 and entry.features[0].type == 'chromosome':
            annotations[entry.id] = entry
    fp.close()

    return annotations

def create_extended_gff(out_dir, gff, entries):
    """
    Extends an existing GFF file by adding 5' and 3'UTR entries.

    This function takes an existing GFF file, e.g. a TriTrypDB genome
    annotation GFF containing entries for mRNAs, genes, etc. and adds UTR
    entries predicted in this program.

    Parameters
    ----------
    out_dir: str
        Location where modified GFF should be written to.
    gff: str
        Filepath to original GFF.
    entries: list
        List of strings representing single line GFF entries
    """
    # TODO: Add UTR entries as children of mRNA GFF entries
    # http://biopython.org/wiki/GFF_Parsing

    # Copy original GFF to output dir
    outfile = os.path.join(out_dir, 
                           os.path.basename(gff).replace('.gff', '-with-utrs.gff'))
    shutil.copyfile(gff, outfile)

    # Open output file for appending
    fp = open(outfile, 'a')

    # Write header to output
    fp.writelines([entry + "\n" for entry in entries])

    # clean up
    fp.close()

def build_gff_utr_entry(feature, gene, chr_id):
    """Takes a single SeqFeature corresponding to either a primary SL or
    polyadenylation site, as well as its associated gene, and returns a GFF 5'
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
            start = gene.location.end
            end = feature.location.start
        else:
            # Negative strand
            start = feature.location.end
            end = gene.location.start
    else:
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

    # Description
    # ID=TcCLB.511911.98_5utr;Name=TcCLB.511911.98;description=hypothetical+protein,+conserved
    short_type = '5utr' if entry_type == 'five_prime_UTR' else '3utr'

    desc = 'ID=%s_%s;Name=%s;description=%s' % (gene.id, short_type, gene.id,
                                                gene.qualifiers['description'][0])

    # GFF parts
    strand = '+' if gene.strand == 1 else '-'
    score = feature.qualifiers['score'][0]

    return [chr_id, 'El-Sayed', entry_type, start + 1, int(end), score, strand, '.', desc]

def create_summary_csv_files(out_dir, utr5_entries, utr3_entries):
    """Creates 5' and 3'UTR summary CSV files"""
    utr5_outfile = os.path.join(out_dir, '5utr_stats.csv')
    utr3_outfile = os.path.join(out_dir, '3utr_stats.csv')

    field_names = ['gene', 'primary_site', 'length', 'num_reads', 'gc', 'ct']

    # write 5'UTR summary statistics
    with  open(utr5_outfile, 'w') as fp:
        writer = csv.writer(fp)
        writer.writerow(field_names)
        writer.writerows(utr5_entries)

    # write 3'UTR summary statistics
    with  open(utr3_outfile, 'w') as fp:
        writer = csv.writer(fp)
        writer.writerow(field_names)
        writer.writerows(utr3_entries)

