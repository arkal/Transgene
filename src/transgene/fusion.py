from __future__ import print_function

import collections
import csv
import logging
import pickle
import re
from cStringIO import StringIO

import swalign
from transgene.common import read_fasta, trans, BEDPE, translate


def get_transcriptome_data(infile, drop_transcript_version=False):
    """
    Parses Gencode transcript FASTA file and returns CDS sequences keyed by the transcript ID

    :param file infile: Gencode FASTA file object
    :param bool drop_transcript_version: Drop the version part of transcript ids.
    :return: {transcript ID: transcripts} and {gene ID: transcript IDs}
    :rtype: tuple
    """
    regex = r"(?P<transcript_id>ENST[0-9A-Z]+.\d+)\|(?P<gene_id>ENSG[0-9A-Z]+.\d+)" \
            r".*CDS:(?P<start>\d+)-(?P<stop>\d+)"
    gene_transcripts = collections.defaultdict(list)
    transcript_cds = {}
    for header, comment, seq in read_fasta(infile, 'ACGT'):
        match = re.search(regex, header)
        if match:
            # Remove the version number on gene ID
            gene_id, version = match.group("gene_id").split('.')
            transcript_id = match.group("transcript_id")
            if drop_transcript_version:
                transcript_id = transcript_id.split('.')[0]
            # GTF is one-based. Make interval [zero, one-based)
            start = int(match.group("start")) - 1
            stop = int(match.group("stop"))
            cds = seq[start: stop]
            gene_transcripts[gene_id].append(transcript_id)
            # Save transcript to gene ID mapping
            gene_transcripts[transcript_id] = match.group("gene_id")
            transcript_cds[transcript_id] = cds
    return transcript_cds, gene_transcripts


def rna_gene_in_bedpe(record):
    """
    Determine if one of the two candidates in a BEDPE line is an rna gene.

    :param BEDPE record: A BEDPE line from the input file
    :returns: True if one of the candidates is an RNA gene and False if not
    :rtype: bool
    """
    #  We will accept fusions that have an RP11- (lncRNA) 3' partner since they can still be
    # translated. This is a heuristic.
    return 'RP11-' in record.hugo1


def readthrough_in_bedpe(record, annotation, rt_threshold):
    """
    Determine if the two genes in the record are within `rt_threshold` bp of each other on the same
    chromosome.

    :param BEDPE record: A BEDPE line from the input file
    :param dict(str, GTFRecord) annotation: see `read_fusions:gene_annotations`
    :param rt_threshold: The genomic distance on the same chromosome below which we will call a
           candidate fusion a readthrough.
    :returns: True if the pair is considered a readthrough and False if not
    :rtype: bool
    """
    return (record.chrom1 == record.chrom2 and
            ((annotation[record.hugo1].start <= annotation[record.hugo2].start <=
                annotation[record.hugo1].end + rt_threshold) or
             (annotation[record.hugo2].start <= annotation[record.hugo1].start <=
                annotation[record.hugo2].end + rt_threshold)))


def read_fusions(fusion_file, gene_annotations, filter_mt, filter_ig, filter_rg, filter_rt,
                 rt_threshold, out_bedpe):
    """
    Reads in gene fusion predictions in modified BEDPE format.
    In addition to the basic BEDPE features, this function requires the fusion
    junction sequences and HUGO names for the donor and acceptor genes.

    :param file fusion_file: Fusion calls in BEDPE format
    :param dict(str, GTFRecord) gene_annotations: The gene annotations from the gtf
    :param bool filter_mt: Filter mitochondrial events?
    :param bool filter_ig: Filter immunoglobulin pairs?
    :param bool filter_rg: Filter RNA-Gene events?
    :param bool filter_rt: Filter transcriptional read-throughs?
    :param int rt_threshold: Distance threshold to call a readthrough
    :param file out_bedpe: A file handle to an output BEDPE file
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

    calls = []

    for line in csv.reader(fusion_file, delimiter='\t'):
        if line[0].startswith('#'):
            print('\t'.join(line), file=out_bedpe)
            continue
        try:
            record = BEDPE(*line)
        except TypeError:
            raise ValueError("ERROR: fusion file is malformed.\n{}".format(read_fusions.__doc__))

        if filter_mt and ('M' in record.chrom1 or 'M' in record.chrom2):
            logging.warning("Rejecting %s-%s for containing a Mitochondrial gene.", record.hugo1,
                            record.hugo2)
            continue
        elif filter_ig and record.hugo1.startswith('IG') and record.hugo2.startswith('IG'):
            # This will drop some Insulin-like growth factor (IGF) proteins but they have a lot of
            # homology too so its ok.
            logging.warning("Rejecting %s-%s an an Immunoglobulin gene pair.", record.hugo1,
                            record.hugo2)
            continue
        elif filter_rg and rna_gene_in_bedpe(record):
            logging.warning("Rejecting %s-%s for containing a 5' RNA gene.", record.hugo1,
                            record.hugo2)
            continue
        elif filter_rt and readthrough_in_bedpe(record, gene_annotations, rt_threshold):
            logging.warning("Rejecting %s-%s as a potential readthrough.", record.hugo1,
                            record.hugo2)
            continue
        else:
            logging.info("Accepting %s-%s for further study.", record.hugo1, record.hugo2)
            print('\t'.join(line), file=out_bedpe)
            calls.append(record)

    if not calls:
        logging.warning('Input bedpe file was empty or had no actionable mutations.')
    return calls

# Namedtuple for storing alignment metrics
# Needs to be global for pickling
AlignStats = collections.namedtuple('AlignStats',
                                    'qstart, qstop, rstart, rstop, insertions, deletions')


def align_filter(ref, query, mode, mismatches_per_kb=1):
    """
    Aligns query to reference CDS sequence using the Smith-Waterman algorithm.
    Returns None if the alignment is clipped at the fusion boundary.

    :param str ref: In-frame reference transcript
    :param str query: Query transcript
    :param str mode: 'donor' or 'acceptor'
    :param int mismatches_per_kb: Allowed number of mismatches per kilobase
                                  of the alignment
    :return: Alignment features
    :rtype: namedtuple
    """

    bound_regex_str = r'Query\s*:\s*(?P<qstart>\d*)\s*[\w-]*\s*(?P<qstop>\d*)\s*[\|\s]*\s*' \
                      r'Ref\s*:\s*(?P<rstart>\d*)\s*[\w-]*\s*(?P<rstop>\d*)'
    bounds_regex = re.compile(bound_regex_str)
    mismatch_regex = re.compile(r'Mismatches: (?P<mismatches>\d+)')

    # Use default swalign parameters
    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
    alignment = sw.align(ref, query)

    # Count the number of insertions and deletions
    insertions = 0
    deletions = 0
    for chr, num in alignment.cigar:
        if chr == 'I':
            insertions += num

        elif chr == 'D':
            deletions += num

    # Next grab the alignment statistics
    sw_output = StringIO()
    alignment.dump(out=sw_output)
    dump = sw_output.getvalue()
    sw_output.close()

    num_mismatches = None
    m = mismatch_regex.search(dump)
    if m:
        num_mismatches = int(m.group('mismatches'))

    s = bounds_regex.search(dump)
    if s:
        qstart = int(s.group('qstart')) - 1    # Make zero-based
        qstop = int(s.group('qstop'))

        # Filter alignments that have more than the allowed number of mismatches per kilobase
        if num_mismatches > int( mismatches_per_kb * (qstop - qstart) ):
            logging.debug("Mismatch filter: %d > %d" % (num_mismatches,
                                                        int(mismatches_per_kb * (qstop - qstart))))
            return

        # Filter alignments that do not include the donor breakpoint
        if mode == 'donor' and qstop != len(query):
            logging.debug('Donor alignment does not reach fusion boundary')
            return

        # Filter alignments that do not include the acceptor breakpoint
        elif mode == 'acceptor' and qstart != 0:
            logging.debug('Acceptor alignment does not reach fusion boundary')
            return

        rstart = int(s.group('rstart')) - 1    # Make zero-based
        rstop = int(s.group('rstart'))

        return AlignStats(qstart, qstop, rstart, rstop, insertions, deletions)


def scan_frame(reference_start):
    """
    Find the frame of a sequencing using the alignment starting position

    :param int reference_start: Alignment start position
    :return: Number of bases to slice sequence to make it in-frame
    :rtype: int
    """
    in_frame_adjustment = 0
    while (reference_start + in_frame_adjustment) % 3 != 0:
        in_frame_adjustment += 1

        if in_frame_adjustment > 3:
            return

    return in_frame_adjustment


def get_donor_junction_exon(breakpoint, exons):
    """
    Finds the first exon before the fusion breakpoint

    :param str transcript_id: donor Ensembl transcript ID
    :param int breakpoint: Genomic position of donor breakpoint
    :param list exons: Sorted list of GTFRecord objects keyed
    :return: Donor junction sequence
    :rtype: str
    """

    # Donor
    # We want the first exon that is upstream of the breakpoint
    #            exon.start < breakpoint        exon.start < breakpoint < exon.end
    #            exon.end  <= breakpoint
    #          _________     |      ___________           ______|_____
    # -->--->--|___1____|----|-----|__________|-----------|_____|_____|--------
    #                        |                                  |
    #                                                                       Traverse <--

    if exons[0].strand == '+':
        for exon in exons[::-1]:
            if exon.start < breakpoint and exon.end <= breakpoint:
                return exon.sequence

            elif exon.start < breakpoint < exon.end:
                # Gencode and breakpoint annotations are one-based
                return exon.sequence[: breakpoint - exon.start]


    # Donor
    # Want first exon that is downstream of breakpoint
    # exon.start < breakpoint < exon.end                  exon.start >= breakpoint
    #                                                     exon.end > breakpoint
    #          ____|_____          ___________      |      ___________
    # --<---<--|___|____|---------|__________|------|-----|_____1_____|--------
    #              |                                |
    # Traverse -->

    elif exons[0].strand == '-':
        for exon in exons[::-1]:
            if exon.start >= breakpoint and exon.end > breakpoint:
                return exon.sequence.translate(trans)[::-1]

            elif exon.start < breakpoint < exon.end:
                seq = exon.sequence[breakpoint - exon.start:]
                return seq.translate(trans)[::-1]


def get_acceptor_junction(breakpoint, exons):
    """
    Finds the first exon after the fusion breakpoint

    :param str transcript_id: donor Ensembl transcript ID
    :param int breakpoint: Genomic position of donor breakpoint
    :param dict exons:
    :return:
    """

    # Acceptor
    # Want first exon downstream of breakpoint
    #               exon.start >= breakpoint      exon.start < breakpoint < breakpoint
    #               exon.end  > breakpoint
    #          _________     |      ___________           ______|_____
    # -->--->--|_______ |----|-----|_____x_____|----------|_____|_____|--------
    #                        |                                  |
    # Traverse -->

    if exons[0].strand == '+':
        for exon in exons:
            if exon.start >= breakpoint and exon.end > breakpoint:
                return exon.sequence

            elif exon.start < breakpoint < exon.end:
                seq = exon.sequence[breakpoint - exon.start:]
                return seq.translate(trans)[::-1]

    # Acceptor
    # Want first exon downstream of breakpoint
    #               exon.start >= breakpoint      exon.start < breakpoint < breakpoint
    #               exon.end  > breakpoint
    #          _________     |      ___________           ______|_____
    # --<---<--|_______ |----|-----|_____x_____|----------|_____|_____|--------
    #                        |                                  |
    #                                                                   <-- Traverse

    elif exons[0].strand == '-':
        for index, exon in enumerate(exons):
            if exon.start < breakpoint and exon.end <= breakpoint:
                return exon.sequence.translate(trans)[::-1]

            elif exon.start < breakpoint < exon.end:
                seq = exon.sequence[: breakpoint - exon.start]
                return seq.translate(trans)[::-1]


def insert_fusions(transcriptome, fusion_calls, gene_transcripts, peplen, outfile, exons=None, output_frameshift=False):
    """
    Generates fusion peptides for MHC binding prediction. If the junction sequence is not
    provided for the donor and acceptor, a theoretical sequence is inferred using a
    genome annotation file.

    :param dict transcriptome: Dictionary of transcripts keyed by transcript ID
    :param list[BEDPE] fusion_calls: List of BEDPE namedtuples
    :param dict[list] gene_transcripts: Dictionary of transcript IDs keyed by gene ID
    :param int peplen: Length of n-mer peptide
    :param file outfile: File object to write n-mer peptides in FASTA format
    :param bool output_frameshift: Output peptides that result from a frameshift in the acceptor sequence
    :param dict exons: Dictionary mapping transcript IDs to ordered list of GTFRecord objects
    """
    logging.info('Generating fusion peptides of length %d' % peplen)

    found = False

    # Load previous alignments
    try:
        with open('transgene_fusion_alignments.pkl') as f:
            previous_alignments = pickle.load(f)

    except (OSError, IOError):
        previous_alignments = {}

    # Iterate over fusion calls
    for call in fusion_calls:

        logging.info('Processing fusion: %s' % repr(call.name))

        donor_name, acceptor_name = call.name.split('-')

        junction_seq1 = call.junctionSeq1
        end1 = int(call.end1)

        # Iterate over all possible transcripts for the donor gene
        for donor_transcript_id in gene_transcripts[donor_name]:

            # If a donor sequence is not provided, find the theoretical one
            if exons and exons[donor_transcript_id] and call.junctionSeq1 == '.':
                logging.debug("Inferring donor junction sequence using breakpoint")
                junction_seq1 = get_donor_junction_exon(end1, exons[donor_transcript_id])

            if junction_seq1 is None:
                continue

            # Check to see if the junction sequence is large enough for neoepitope prediction
            if len(junction_seq1) // 3 < peplen - 1 :
                logging.warn('Donor sequence is not long enough for neoepitope prediction'
                             '\n%s' % junction_seq1)
                continue

            # Find the frame in the donor sequence using a reference transcriptome
            try:
                donor_transcript = transcriptome[donor_transcript_id]

                # Check if we have already aligned this sequence (ie for another peptide length)
                try:
                    bounds = previous_alignments[(donor_transcript, junction_seq1)]

                # If not, then align the sequence and save the alignment information
                except KeyError:
                    bounds = align_filter(donor_transcript, junction_seq1, 'donor')
                    previous_alignments[(donor_transcript, junction_seq1)] = bounds

                # Donor sequence may not align to some isoforms
                if bounds is None:
                    logging.debug('Donor sequence did not align to reference transcripts')
                    continue

            except KeyError:
                continue

            # Find the position of the transcript where the predicted donor sequence starts
            in_frame_adjustment = scan_frame(bounds.rstart)
            if in_frame_adjustment is None:
                continue

            # Slice the donor transcript such that the first base is in-frame
            # May have an overhang base
            donor_start = bounds.qstart + in_frame_adjustment
            in_frame_donor_seq = junction_seq1[donor_start:]
            in_frame_donor_seq_len = len(in_frame_donor_seq)

            if acceptor_name in gene_transcripts:
                acceptor_transcript_ids = gene_transcripts[acceptor_name]

            else:
                logging.debug('Could not find acceptor transcript ID')
                break

            # Now iterate over all possible acceptor transcript sequences
            start2 = int(call.start2)
            for acceptor_transcript_id in acceptor_transcript_ids:
                junction_seq2 = call.junctionSeq2

                # Predict acceptor sequence if none was provided
                if exons and acceptor_transcript_id in exons and junction_seq2 == '.':
                    junction_seq2 = get_acceptor_junction(start2, exons[acceptor_transcript_id])

                if junction_seq2 is None:
                    logging.debug('Could not infer theoretical acceptor junction sequence')
                    continue

                acceptor_in_frame_adjustment = None
                try:
                    # Get reference coding transcript for acceptor gene
                    acceptor_transcript = transcriptome[acceptor_transcript_id]

                    try:
                        # Find the position of the transcript where the predicted acceptor sequence starts
                        acceptor_bounds = previous_alignments[(acceptor_transcript, junction_seq2)]

                    except KeyError:
                        acceptor_bounds = align_filter(acceptor_transcript, junction_seq2, 'donor')
                        previous_alignments[(acceptor_transcript, junction_seq2)] = bounds

                    # If the acceptor sequence could not be found, move onto the next transcript
                    if acceptor_bounds is None:
                        continue

                    acceptor_in_frame_adjustment = scan_frame(acceptor_bounds.rstart)
                    if acceptor_in_frame_adjustment is None:
                        continue

                # We can still generate a fusion transcript without a reference for the acceptor gene
                except KeyError:
                    pass

                # Figure out if frame-shift occurred as a result of the fusion
                acceptor_frame_shift = False

                # If we don't have a reference transcript, then just continue without this information
                if acceptor_in_frame_adjustment is None:
                    pass

                # There could be an indel that corrects the frameshift downstream
                elif (in_frame_donor_seq_len % 3 + acceptor_in_frame_adjustment) % 3 != 0:
                    logging.warn('%s--%s has a 3\' frame shift' % (call.hugo1, call.hugo2))
                    acceptor_frame_shift = True

                # Use the donor and acceptor sequences to generate the fusion transcript
                fusion_seq = in_frame_donor_seq + junction_seq2
                logging.debug(fusion_seq)

                # The fusion peptide for MHC binding prediction includes all sliding windows that would
                # generate a novel peptide sequence. Therefore, the donor sequence starts at the position
                # of the last donor codon minus the fusion peptide length plus two codons. The acceptor
                # portion consists of the last donor codon plus the peptide length.
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
                #   peplen = 9
                #
                #   start = 6 - 9 + 6 = 3
                #   end   = 6 + 9     = 14

                last_donor_codon = in_frame_donor_seq_len - in_frame_donor_seq_len % 3 - 3

                start = last_donor_codon - 3 * peplen + 6

                end = last_donor_codon + 3 * peplen

                in_frame_junction_seq = fusion_seq[start: end]

                fusion_peptide = translate(in_frame_junction_seq)
                logging.debug("Fusion peptide:\n%s" % fusion_peptide)

                if '*' in fusion_peptide:
                    logging.warn('Skipping fusion epitope containing stop codon')
                    continue


                # This may occur if the fusion sequence occurs near the end of the transcript. Warn the user, but
                # still output the fusion peptide because there is at least one fusion sequence that can be generated
                # from it.
                if len(fusion_peptide) < peplen:
                    logging.warning('Fusion peptide %s-%s is malformed.\n'
                                    'Expected peptide of length %d, got %d' % (donor_name,
                                                                               acceptor_name,
                                                                               2 * peplen - 2,
                                                                               len(fusion_peptide)))
                    break

                gencode_donor_id = gene_transcripts[donor_transcript_id]
                gencode_acceptor_id = gene_transcripts[acceptor_transcript_id]

                # Write fusion neoepitope sequence to FASTA file
                fasta = '>{donor_gene}-{acceptor_gene}_' \
                        '{donor_transcript}-{acceptor_transcript}_' \
                        '{donor_hugo}-{acceptor_hugo}_' \
                        'FUSION_{score}\n' \
                        '{sequence}\n'.format(donor_gene=gencode_donor_id,
                                              acceptor_gene=gencode_acceptor_id,
                                              donor_transcript=donor_transcript_id,
                                              acceptor_transcript=acceptor_transcript_id,
                                              donor_hugo=call.hugo1,
                                              acceptor_hugo=call.hugo2,
                                              score=call.score,
                                              sequence=fusion_peptide)
                outfile.write(fasta)
                found = True

                # A frameshift in the acceptor sequence also generates
                # candidate neoepitope sequences
                if output_frameshift and acceptor_frame_shift:
                    for i in range(end, len(fusion_seq) - 3 * peplen, 3):
                        shift_sequence = fusion_seq[i: i + 3 * peplen]
                        shift_peptide = translate(shift_sequence)

                        if '*' in shift_peptide:
                            break

                        if len(shift_peptide) != peplen:
                            break

                        logging.debug('Created fusion frameshift peptide:\n%s' % shift_peptide)


                        # Write frameshift peptide to FASTA file
                        fasta = '>{donor_gene}-{acceptor_gene}_' \
                                '{donor_transcript}-{acceptor_transcript}_' \
                                '{donor_hugo}-{acceptor_hugo}_' \
                                'FUSION_ACCEPTOR_FRAMESHIFT\n' \
                                '{sequence}\n'.format(donor_gene=gencode_donor_id,
                                                      acceptor_gene=gencode_acceptor_id,
                                                      donor_transcript=donor_transcript_id,
                                                      acceptor_transcript=acceptor_transcript_id,
                                                      donor_hugo=call.hugo1,
                                                      acceptor_hugo=call.hugo2,
                                                      sequence=shift_peptide)
                        outfile.write(fasta)

    if not found:
        logging.warn("Did not generate any fusion peptides from data!")

    # Save alignments for different peptide lengths
    with open('transgene_fusion_alignments.pkl', 'w') as f:
        pickle.dump(previous_alignments, f)

