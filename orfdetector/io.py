"""
I/O functions
"""
import csv
from argparse import ArgumentParser

def parse_args():
    """Parses input and returns arguments"""
    parser = ArgumentParser(description='Detect novel ORFs using RNA-Seq data.')

    # parser.add_argument('-b', '--bam', help='Input bam file', required=True)
    parser.add_argument('-c', '--coverage', 
                            help='Single nucleotide resolution genome coverage map',
                            required=True)
    parser.add_argument('-f', '--fasta', help='Input FASTA file',
                        required=True)
    parser.add_argument('-g', '--gff', help='Input GFF file', required=True)
    parser.add_argument('-l', '--min-protein-length', required=True, 
                        type='int', default=30,
                        help=('Minimum size in amino acids allowed for '
                                'novel ORFs. (default=30)'))
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


def write_output_gff(:
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
