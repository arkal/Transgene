from __future__ import print_function
import collections
import csv
import logging
import re
from cStringIO import StringIO

import swalign

from common import read_fasta


# Standard Genetic Code from NCBI
amino = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
base1 = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
base2 = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
base3 = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
genetic_code = dict((''.join([b1, b2, b3]), aa) for aa, b1, b2, b3 in zip(amino, base1, base2, base3))


def get_transcriptome_data(infile):
    """
    Parses Gencode transcript FASTA file and returns CDS sequences keyed by the transcript ID

    :param file infile: Gencode FASTA file object
    :return: Dictionary of transcripts keyed by transcript ID and dictionary of transcript IDs keyed by gene ID
    :rtype: tuple
    """
    regex = r"(?P<transcript_id>ENST[0-9A-Z]+.\d+)\|(?P<gene_id>ENSG[0-9A-Z]+.\d+).*CDS:(?P<start>\d+)-(?P<stop>\d+)"
    gene_transcripts = collections.defaultdict(list)
    transcript_cds = {}
    for header, comment, seq in read_fasta(infile, 'ACGT'):
        match = re.search(regex, header)
        if match:
            # Remove the version number on gene ID
            gene_id, version = match.group("gene_id").split('.')
            transcript_id = match.group("transcript_id")
            # GTF is one-based
            start = int(match.group("start")) - 1
            stop = int(match.group("stop"))
            cds = seq[start: stop]

            gene_transcripts[gene_id].append(transcript_id)
            transcript_cds[transcript_id] = cds
    return transcript_cds, gene_transcripts


def translate(seq):
    """
    Translates DNA sequence into protein sequence using globally defined genetic code

    :param str seq: DNA sequence
    :returns: Translated sequence
    :rtype: str

    >>> translate('ATGTTTCGTT')
    'MFR'
    """
    start = 0
    n = len(seq)
    codons = (seq[i: i+3] for i in range(start, n - n % 3, 3))
    protein = [genetic_code[codon] for codon in codons]
    return ''.join(protein)


def read_fusions(fusion_file):
    """
    Reads in gene fusion predictions in modified BEDPE format. In addition to the basic BEDPE features, this
    function requires the fusion junction sequences and HUGO names for the donor and acceptor genes.

    :param file fusion_file: Fusion calls in BEDPE format
    :returns: list of BEDPE namedtuples
    :rtype: list

    Modified BEDPE format
    chrom1:         Chromosome of first feature
    start1:         Zero-based starting position of the first feature or '.'
    end1:           One-based ending position of the first feature -- 5' fusion breakpoint
    chrom2:         Chromosome of second feature
    start2:         Zero-based starting position of the second feature -- 3' fusion breakpoint
    end2:           One-based ending position of thh second feature or '.'
    name:           Hyphenated Ensembl gene IDs (i.e. ENSG00000168172-ENSG00000165731)
    score:          Optional fusion score
    strand1:        Strand of first feature
    strand2:        Strand of second feature
    junctionSeq1:   Fusion junction sequence in first feature
    junctionSeq2:   Fusion junction sequence in second feature
    hugo1:          HUGO name for first feature
    hugo2:          HUGO name for second feature
    """
    BEDPE = collections.namedtuple('BEDPE',
                                   'chrom1, start1, end1, '
                                   'chrom2, start2, end2, '
                                   'name, score, '
                                   'strand1, strand2, '
                                   'junctionSeq1, junctionSeq2, '
                                   'hugo1, hugo2')

    calls = []
    for line in csv.reader(fusion_file, delimiter='\t'):
        if line[0].startswith('#'):
            continue
        try:
            calls.append(BEDPE(*line))

        except TypeError:
            raise ValueError("ERROR: fusion file is malformed.\n{}".format(read_fusions.__doc__))

    return calls

def align_filter(ref, query, mode):
    alignment_stats = collections.namedtuple('AlignStats', 'qstart, qstop, rstart, rstop, insertions, deletions')

    bounds_regex = re.compile(r'Query\s*:\s*(?P<qstart>\d*)\s*\w*\s*(?P<qstop>\d*)\s*[\|\s]*\s*Ref\s*:\s*(?P<rstart>\d*)\s*\w*\s*(?P<rstop>\d*)')
    match_regex = re.compile(r'Matches: \d+\s\((?P<percent>\d*)')

    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
    alignment = sw.align(ref, query)

    # First check that the donor sequence is in frame
    insertions = 0
    deletions = 0
    for chr, num in alignment.cigar:
        if chr == 'I':
            insertions += num

        elif chr == 'D':
            deletions += num

    # Next grab the alignment statistics
    string = StringIO()
    alignment.dump(out=string)
    dump = string.getvalue()
    string.close()

    # If it's not a near perfect match, then the quality of the assembly may not be good
    m = match_regex.search(dump)
    if m:
        percent = int(m.group('percent'))
        if percent < 99:
            # print('Percent matching %d' % percent)
            # print(dump)
            return

    # If the fusion transcript passes these filters, then grab the bounds of the alignment
    s = bounds_regex.search(dump)
    if s:
        qstart = int(s.group('qstart')) - 1    # Make zero-based
        qstop = int(s.group('qstop'))

        # If the end of the fusion transcript doesn't align, then skip this transcript
        if mode == 'donor' and qstop != len(query):
            # print('Donor doesn\'t end at assembled transcript end')
            # print('Query len', len(query))
            # print('qstop', qstop)
            # print(dump)
            return

        elif mode == 'acceptor' and qstart != 0:
            # print('Acceptor doesn\'t start at one')
            # print(dump)
            return

        rstart = int(s.group('rstart')) - 1    # Make zero-based
        rstop = int(s.group('rstart'))

        return alignment_stats(qstart, qstop, rstart, rstop, insertions, deletions)

    else:
        return

def scan_frame(reference_start):
    """
    Find the frame of a sequence relative to an in frame sequence

    :param reference_start:
    :return: Number of bases to slice sequence to make it in-frame
    :rtype: int
    """
    # Find the position of the transcript where the predicted donor sequence starts
    in_frame_adjustment = 0
    while (reference_start + in_frame_adjustment) % 3 != 0:
        in_frame_adjustment += 1

        if in_frame_adjustment > 3:
            return

    return in_frame_adjustment


def insert_fusions(transcriptome, fusion_calls, gene_transcripts, peplen, outfile):
    """
    Generates fusion peptides for MHC binding prediction.

    :param dict transcriptome: Dictionary of transcripts keyed by transcript ID
    :param list[BEDPE] fusion_calls: List of BEDPE namedtuples
    :param dict[list] gene_transcripts: Dictionary of transcript IDs keyed by gene ID
    :param int peplen: Length of n-mer peptide
    :param file outfile: File object to write n-mer peptides in FASTA format
    """
    logging.info('Generating fusion peptides of length %d' % peplen)

    found = False

    # Iterate over fusion calls
    for call in fusion_calls:

        # Check to see if the junction sequence is large enough for neoepitope prediction
        if len(call.junctionSeq1) // 3 < peplen - 1 or len(call.junctionSeq2) // 3 < 1:
            continue

        logging.info('Processing fusion: %s' % repr(call.name))

        donor_name, acceptor_name = call.name.split('-')
        for donor_transcript_id in gene_transcripts[donor_name]:
            # print(call.hugo1, call.hugo2)
            # Generate the fusion transcript using the GENCODE transcript annotation
            try:
                donor_transcript = transcriptome[donor_transcript_id]
                bounds = align_filter(donor_transcript, call.junctionSeq1, 'donor')
                if bounds is None:
                    # print("DONOR BOUNDS")
                    continue

            except KeyError:
                continue

            # Find the position of the transcript where the predicted donor sequence starts
            in_frame_adjustment = scan_frame(bounds.rstart)
            if in_frame_adjustment is None:
                continue

            # Slice the donor transcript such that the first base is in-frame
            # May have a overhang base
            donor_start = bounds.qstart + in_frame_adjustment
            in_frame_donor_seq = call.junctionSeq1[donor_start:]
            in_frame_donor_seq_len = len(in_frame_donor_seq)

            # If acceptor Ensembl ID is not in the GENCODE annotation, then label transcript as not applicable. This
            # may occur for intergenic fusions.
            if acceptor_name in gene_transcripts:
                acceptor_transcript_ids = gene_transcripts[acceptor_name]

            else:
                acceptor_transcript_ids = ['NA']

            # Now iterate over all possible acceptor transcript sequences
            for acceptor_transcript_id in acceptor_transcript_ids:
                acceptor_in_frame_adjustment = None
                try:
                    acceptor_transcript = transcriptome[acceptor_transcript_id]

                    # Find the position of the transcript where the predicted donor sequence starts
                    acceptor_bounds = align_filter(acceptor_transcript, call.junctionSeq2, 'acceptor')

                    # If the acceptor sequence could not be found, move onto the next transcript
                    if acceptor_bounds is None:
                        # print('ACCEPTOR BOUNDS')
                        continue

                    acceptor_in_frame_adjustment = scan_frame(acceptor_bounds.rstart)
                    if acceptor_in_frame_adjustment is None:
                        # print('Acceptor frame adjustment')
                        continue

                # If we can't find a transcript for the acceptor, then continue without a transcript ID
                except KeyError:
                    pass

                acceptor_frame_shift = False

                # Figure out if frame-shift occurred as a result of the fusion
                # If we don't have a reference transcript, then just continue without this information
                if acceptor_in_frame_adjustment is None:
                    pass

                # TODO need to develop this more
                # There could be an indel that corrects the frame downstream
                elif (in_frame_donor_seq_len % 3 + acceptor_in_frame_adjustment) % 3 != 0:
                    acceptor_frame_shift = True

                # Use the donor and acceptor sequences to generate the fusion transcript
                fusion_seq = in_frame_donor_seq + call.junctionSeq2

                # The fusion peptide for MHC binding prediction includes all sliding windows that would generate a
                # novel peptide sequence. Therefore, the donor sequence starts at the position of the last donor codon
                # minus the fusion peptide length plus two codons. The acceptor portion consists of the last donor codon
                # plus the peptide length.
                #
                #  Transcript Start: 0  3  6  9  12 15
                #    Transcript End: 2  5  8  11 14 17
                #       Codon Index: 0  1  2  3  4  5
                #    Fusion Protein: D  D  D  A  A  A
                #
                #                       D  D  A
                #                          D  A  A
                #              Result:  D  D  A  A
                #
                #   start = 6 - 9 + 6 = 3
                #   end   = 6 + 9     = 14

                last_donor_codon = in_frame_donor_seq_len - in_frame_donor_seq_len % 3 - 3

                start = last_donor_codon - 3 * peplen + 6

                end = last_donor_codon + 3 * peplen

                in_frame_junction_seq = fusion_seq[start: end]

                fusion_peptide = translate(in_frame_junction_seq)

                # This may occur if the fusion sequence occurs near the end of the transcript. Warn the user, but
                # still output the fusion peptide because there is at least one fusion sequence that can be generated
                # from it.
                if len(fusion_peptide) != 2 * peplen - 2:
                    logging.warning('Fusion peptide %s-%s is malformed.\n'
                                    'Expected peptide of length %d, got %d' % (donor_name,
                                                                               acceptor_name,
                                                                               2 * peplen - 2,
                                                                               len(fusion_peptide)))

                # print(call.hugo1, call.hugo2)
                # print(donor_name, acceptor_name)
                # print(donor_transcript_id, acceptor_transcript_id)
                # print(fusion_peptide)

                # Write fusion neoepitope sequence to FASTA file
                fasta = '>{donor_gene}-{acceptor_gene}_' \
                        '{donor_transcript}-{acceptor_transcript}_' \
                        '{donor_hugo}-{acceptor_hugo}_' \
                        'FUSION\n{sequence}\n'.format(donor_gene=donor_name,
                                                      acceptor_gene=acceptor_name,
                                                      donor_transcript=donor_transcript_id,
                                                      acceptor_transcript=acceptor_transcript_id,
                                                      donor_hugo=call.hugo1,
                                                      acceptor_hugo=call.hugo2,
                                                      sequence=fusion_peptide)
                outfile.write(fasta)
                found = True

    if not found:
        logging.warn("Did not generate any fusion peptides from data")
