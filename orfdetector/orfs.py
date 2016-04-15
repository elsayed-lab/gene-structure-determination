"""
Functions for detecting ORFs in a genome.
"""
def search_genome_for_orfs(genome_sequence, inter_cds_regions, min_protein_length=30):
    # Create a dictionary to store output; for ease of use, we will store
    # each set of detected ORFs using the same keys as the inter_cds_regions
    # dictionary (chromosome, strand, and region number)
    orfs = {}

    # Iterate over chromosomes
    for i, chr_id in enumerate(inter_cds_regions, 1):
        orfs[chr_id] = {}
        print("Processing %s (%d/%d)" % (chr_id, i, len(genome_sequence)))

        # Iterate over strands
        for strand in [-1, 1]:
            orfs[chr_id][strand] = []

            # Iterate through inter-CDS regions
            for j, region in enumerate(inter_cds_regions[chr_id][strand]):
                # Get sequence record for the range
                record = genome_sequence[chr_id][region[0]:region[1]]

                # Find ORFs in each of the six possible reading frames
                # that are at least the specified length
                inter_cds_orfs = find_orfs(record.seq,
                                           min_protein_length, strand)

                # Choose the most like ORF among the possibilities, based
                # on support from coverage and predicted SL acceptor sites
                # and polyadenylation sites
                # TODO

                # Write GFF entries for each match
                orfs_detected = []

                for k, orf in enumerate(inter_cds_orfs, 1):
                    # Convert coordinates to chromosomal position
                    #if orf[2] == "+":
                    start = orf[0] + region[0] + 1
                    stop = orf[1] + region[0]
                    str_strand = orf[2]

                    orfs_detected.append({
                        'id': "%s.ORF.%d.%d" % (chr_id, j, k),
                        'chr': chr_id,
                        'start': start,
                        'stop': stop,
                        'strand': str_strand
                    })

                # add to output
                orfs[chr_id][strand].append(orfs_detected)

    return orfs


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

