"""
I/O functions
"""
import csv
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


def write_output_gff(output):
    # Iterate over inter-CDS regions and find ORFs of at least the specified
    # length in any of the six possible reading frames and output as GFF entries
    fp = open(output, 'w')

    # Write csv header
    fp.write("##gff-version\t3\n")
    fp.write("##feature-ontology\tsofa.obo\n")
    fp.write("##attribute-ontology\tgff3_attributes.obo\n")

    # Write header to output
    writer = csv.writer(fp, delimiter='\t')

    # Iterate through inter_cds regions
    for orf in orfs:
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

    fp = open(gff)
    for entry in GFF.parse(fp):
        if len(entry.features) > 0 and entry.features[0].type == 'chromosome':
            annotations[entry.id] = entry
    fp.close()

    return annotations
