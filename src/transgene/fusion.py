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


def get_donor_junction_exon(breakpoint, exons, cds_start, peplen):
    """
    Finds the first exon before the fusion breakpoint

    :param int breakpoint: Genomic position of donor breakpoint
    :param list exons: Sorted list of exons for the current transcript as GTFRecord objects
    :param int cds_start: The CDS start for the transcript
    :param int peplen: The length of peptides that will be estimated from this IAR.
    :return: In-frame donor junction sequence with (n-1)*3 + overhang nucleotides (for n-1 normal
             AAs and < 3 bp that will make the first modifed AA, and a flag indicating that the
             breakpoint was in an intron.
    :rtype: str, bool
    """
    sequence = []
    intron_junction = True
    if exons[0].strand == '+':
        if breakpoint < cds_start:
            logging.warning('Predicted breakpoint %s is in in the 5\' UTR of %s. Skipping',
                            breakpoint, exons[0].transcript_id)
            return None, None
        for i, exon in enumerate(exons):
            if breakpoint < exon.start:
                # The breakpoint was in the previous intron
                if breakpoint <= exons[i - 1].end + 23:
                    logging.warning('Breakpoint within 23bp of exon end. Splicing may be affected.')
                break
            if cds_start > exon.end:
                # This handles an exon of just 5' UTR
                continue
            start = max(0, cds_start - exon.start)
            end = exon.end - exon.start if breakpoint > exon.end else breakpoint - exon.start
            sequence.append(exon.sequence[start:end + 1])
            if exon.start <= breakpoint <= exon.end:
                intron_junction = False
                break
    else:
        if breakpoint > cds_start:
            logging.debug('Predicted breakpoint in the 5\' UTR. Skipping')
            return None, None
        for i, exon in enumerate(exons):
            if breakpoint > exon.end:
                # The breakpoint was in the previous intron
                if breakpoint >= exons[i - 1].start - 23:
                    logging.warning('Breakpoint within 23bp of exon end. Splicing may be affected.')
                break
            if cds_start < exon.start:
                # This handles an exon of just 5' UTR
                continue
            start = 0 if breakpoint < exon.start else breakpoint - exon.start
            end = min(cds_start + 2, exon.end) - exon.start
            sequence.append(exon.sequence[start:end + 1][::-1].translate(trans))
            if exon.start <= breakpoint <= exon.end:
                intron_junction = False
                break
    sequence = ''.join(sequence)
    overhang = len(sequence) % 3
    if overhang == 0:
        return sequence[-(peplen - 1) * 3:], intron_junction
    else:
        return sequence[-((peplen - 1) * 3 + overhang):], intron_junction


def get_acceptor_junction(breakpoint, exons, cds_start, chrom_seq, peplen):
    """
    Finds the first exon after the fusion breakpoint

    :param int breakpoint: Genomic position of donor breakpoint
    :param list exons: Sorted list of exons for the current transcripts as GTFRecord objects
    :param int cds_start: The CDS start for the transcript
    :param list chrom_seq: The sequence for the current chromosome
    :param int peplen: The length of peptides that will be estimated from this IAR.
    :return: An in-frame acceptor sequence of len 2 * peplen * 3 with a 5' overhang,
             A genomic sequence of length 99(if breakpoint intronic else None,
             and a flag indicating if the breakpoint was in an intron
    :rtype: str, str|None, bool
    """
    sequence = []
    genomic_sequence = None
    intron_junction = False
    acceptor_exon_reached = False
    skipped_nucs = 0
    if exons[0].strand == '+':
        #                 E1                        E2
        #    -------|~~~~~=========|-----------|==========|----ETC
        #                 C
        # B1:          ^                                                  : Handled outside the fn
        # B2:                 ^                                           : B < E1e
        #                                                                   E1s < B < E1e
        # B3:                            ^                                : B > E1e, B < E2e
        #                                                                   B < E2s => INTRON
        # B4:                                       ^                     : B > E1e, B < E2e
        #                                                                   E2s < B < E2e
        for i, exon in enumerate(exons):
            if breakpoint > exon.end:
                # Breakpoint is after this exon
                if exon.start <= cds_start <= exon.end:
                    skipped_nucs += exon.end - cds_start + 1
                else:
                    skipped_nucs += exon.end - exon.start + 1
                continue
            if not acceptor_exon_reached:
                acceptor_exon_reached = True
                if breakpoint < exon.start:
                    intron_junction = True
                    logging.info('Introninc breakpoint detected.')
                    genomic_sequence = get_genomic_sequence(breakpoint, chrom_seq, 33,
                                                            exon.strand)
                    if breakpoint >= exon.start - 23:
                        logging.warning('Breakpoint within 23bp of exon start. Splicing may be '
                                        'affected.')
                else:
                    skipped_nucs += (breakpoint - exon.start)
            start = breakpoint - exon.start if breakpoint >= exon.start else 0
            sequence.append(exon.sequence[start:])
    else:
        #                       E2                     E1
        #      CTE----------|=======|-----------|=======~~~~~|-------
        #                                              C
        # B1:                                             ^              : Handled outside the fn
        # B2:                                     ^                      :  B > E1s
        #                                                                       E1s < B < E1e
        # B3:                             ^                              :  B < E1s, B > E2s
        #                                                                       B > E2e => INTRON
        # B4:                 ^                                          :  B < E1s, B > E2s
        #                                                                       E2s < B < E2e
        for i, exon in enumerate(exons):
            if breakpoint < exon.start:
                # Breakpoint is after this exon
                if exon.start <= cds_start + 2 <= exon.end:
                    skipped_nucs += cds_start + 2 - exon.start + 1
                else:
                    skipped_nucs += exon.end - exon.start + 1
                continue
            if not acceptor_exon_reached:
                acceptor_exon_reached = True
                if breakpoint > exon.end:
                    intron_junction = True
                    logging.info('Introninc breakpoint detected.')
                    genomic_sequence = get_genomic_sequence(breakpoint, chrom_seq, 33,
                                                            exon.strand)
                    if breakpoint <= exon.end + 23:
                        logging.warning('Breakpoint within 23bp of exon start. Splicing may be '
                                        'affected.')
                else:
                    skipped_nucs += (exon.end - breakpoint)
            end = breakpoint - exon.start + 1 if breakpoint <= exon.end else None
            sequence.append(exon.sequence[:end][::-1].translate(trans))
    sequence = ''.join(sequence)

    if skipped_nucs % 3 == 0:
        # The first base in the junction is in-frame on the acceptor. Return 2x peplen
        sequence = sequence[:min(len(sequence), peplen * 3 * 2)]
    else:
        # The first base was not in-frame. Return 2x peplen + the overhang
        sequence = sequence[:min(len(sequence),  3 - (skipped_nucs % 3) + peplen * 3 * 2)]
    return sequence, genomic_sequence, intron_junction


def insert_fusions_from_junctions(transcriptome, call, gene_transcripts, peplen, outfile,
                                  previous_alignments):
    """
    Generates fusion peptides for MHC binding prediction. If the junction sequence is not
    provided for the donor and acceptor, a theoretical sequence is inferred using a
    genome annotation file.

    :param dict transcriptome: Dictionary of transcripts keyed by transcript ID
    :param call: A single BEDPE fusion call
    :param dict[list] gene_transcripts: Dictionary of transcript IDs keyed by gene ID
    :param int peplen: Length of n-mer peptide
    :param file outfile: File object to write n-mer peptides in FASTA format
    :param dict previous_alignments: A dict of previous alignments from handling another read length
    :returns: True if an epitope was fond else False
    :rtype: bool
    """
    found = False
    donor_name, acceptor_name = call.name.split('-')
    logging.info('Processing using 5\' and 3\' junction sequences')
    junction_seq1 = call.junctionSeq1

    # Iterate over all possible transcripts for the donor gene
    for donor_transcript_id in gene_transcripts[donor_name]:
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
        # Now retain only 8 codons + overhang
        overhang = len(in_frame_donor_seq) % 3
        if overhang == 0:
            in_frame_donor_seq = junction_seq1[-(peplen - 1) * 3:]
        else:
            in_frame_donor_seq = junction_seq1[-((peplen - 1) * 3 + overhang):]

        donor_overhang = len(in_frame_donor_seq) % 3

        if donor_overhang == 0:
            acceptor_seq = call.junctionSeq2[:(peplen - 1) * 3]
        else:
            acceptor_seq = call.junctionSeq2[:3 - donor_overhang + (peplen - 1) * 3]

        # Use the donor and acceptor sequences to generate the fusion transcript
        junction_seq = [in_frame_donor_seq, acceptor_seq]

        logging.debug('Inferred junction seq: %s %s', junction_seq[0], junction_seq[1])

        fusion_peptide = translate(''.join(junction_seq))
        logging.debug("Inferred fusion peptide: %s" % fusion_peptide)

        if '*' in fusion_peptide:
            fusion_peptide = fusion_peptide[:fusion_peptide.find('*')]

        if len(fusion_peptide) < peplen:
            logging.warning('Fusion peptide %s from %s-%s is is less than %s '
                            'residues.', fusion_peptide, donor_name, acceptor_name, peplen)
            continue

        # Write fusion neoepitope sequence to FASTA file
        fasta = '>{donor_gene}-{acceptor_gene}_' \
                '{donor_transcript}-{acceptor_transcript}_' \
                '{donor_hugo}-{acceptor_hugo}_' \
                '{donor_breakpoint}-{acceptor_breakpoint}_' \
                'FUSION_{score}\n' \
                '{sequence}\n'.format(donor_gene=donor_name,
                                      acceptor_gene=acceptor_name,
                                      donor_transcript=donor_name,
                                      acceptor_transcript=acceptor_name,
                                      donor_hugo=call.hugo1,
                                      acceptor_hugo=call.hugo2,
                                      donor_breakpoint=call.end1,
                                      acceptor_breakpoint=call.start2,
                                      score=call.score,
                                      sequence=fusion_peptide)
        outfile.write(fasta)
        found = True

    return found


def insert_fusions(transcriptome, fusion_calls, gene_transcripts, peplen, outfile, chroms=None,
                   exons=None, cds_starts=None, extend_length=10):
    """
    Generates fusion peptides for MHC binding prediction. If the junction sequence is not
    provided for the donor and acceptor, a theoretical sequence is inferred using a
    genome annotation file.

    :param dict transcriptome: Dictionary of transcripts keyed by transcript ID
    :param list[BEDPE] fusion_calls: List of BEDPE namedtuples
    :param dict[list] gene_transcripts: Dictionary of transcript IDs keyed by gene ID
    :param int peplen: Length of n-mer peptide
    :param file outfile: File object to write n-mer peptides in FASTA format
    :param dict chroms: Contains the chromosomal dna in the form of a dict where keys hold the
           chromosome name and the values holds the sequence.
    :param dict exons: See return value of `get_exons`
    :param dict cds_starts: See return value of `transgene::get_exons`
    :param int extend_length: The number of codons downstream of a fusion to process
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

        if call.junctionSeq1 != '.' and call.junctionSeq2 != '.':
            # If the user provides BOTH junctions, just use their sequences.
            f = insert_fusions_from_junctions(transcriptome, call, gene_transcripts, peplen,
                                              outfile, previous_alignments)
            if f:
                found = True
            continue

        # Iterate over all possible transcripts for the donor gene
        for donor_transcript_id in gene_transcripts[donor_name]:
            logging.debug('Processing donor transcript %s', donor_transcript_id)
            # A Flag indicating if the donor junction was in an intron
            donor_in_intron = False
            # If a donor sequence is not provided, find the theoretical one
            if exons and exons[donor_transcript_id]:
                logging.debug("Inferring donor junction sequence using breakpoint")
                if donor_transcript_id not in cds_starts:
                    logging.warning('Cowardly refusing to process %s as the 5\' donor sequence '
                                    'since it does not have a defined CDS start.',
                                    donor_transcript_id)
                    continue
                in_frame_donor_seq, donor_in_intron = get_donor_junction_exon(
                    int(call.end1), exons[donor_transcript_id], cds_starts[donor_transcript_id],
                    peplen)
                logging.debug('Inferred 5\' Junction sequence with breakpoint %s: %s', call.end1,
                              in_frame_donor_seq)
            else:
                if call.junctionSeq1 != '.':
                    logging.warning('Cannnot process %s as the 5\' donor sequence since it does '
                                    'not have an entry in the input GTF. A junction sequence was '
                                    'provided but not used since the 3\' partner has no junction '
                                    'sequence .', donor_transcript_id)
                else:
                    logging.warning('Cannnot process %s as the 5\' donor sequence since it does '
                                    'not have an entry in the input GTF and a junction sequence '
                                    'was not provided.', donor_transcript_id)
                continue

            if in_frame_donor_seq is None:
                # Skip 5' UTR
                continue
            # Check to see if the junction sequence is large enough for neoepitope prediction
            if len(in_frame_donor_seq) < 1:
                logging.warn('Donor sequence is not long enough for neoepitope prediction'
                             '\n%s' % in_frame_donor_seq)
                continue

            # Now process the acceptor
            denovo_acceptor = False
            if acceptor_name in gene_transcripts:
                acceptor_transcript_ids = gene_transcripts[acceptor_name]
                genomic_seq = None
            else:
                logging.debug('Could not find acceptor transcript ID for %s. Attempting to infer '
                              'an acceptor sequence from the genome.', acceptor_name)
                # If there is no downstream transcript, it could be a mapping to an RNA Gene. We
                # will have to use sequence from the genome.
                strand = call.strand2
                chrom = call.chrom2
                denovo_acceptor = True
                genomic_seq = get_genomic_sequence(int(call.start2), chroms[chrom],
                                                   extend_length, strand)
                logging.debug('Inferred acceptor sequence with breakpoint %s: %s', call.start2,
                              genomic_seq)
                acceptor_transcript_ids = [acceptor_name]

            # Now iterate over all possible acceptor transcript sequences
            start2 = int(call.start2)
            donor_overhang = len(in_frame_donor_seq) % 3
            for acceptor_transcript_id in acceptor_transcript_ids:
                logging.debug('Processing acceptor transcript %s', acceptor_transcript_id)
                acceptor_seq = ''
                acceptor_in_intron = False
                acceptor_in_utr = False
                if not denovo_acceptor:
                    # Predict acceptor sequence if none was provided
                    if exons and acceptor_transcript_id in exons:
                        if (acceptor_transcript_id not in cds_starts or
                                (call.strand2 == '+' and
                                 start2 < cds_starts[acceptor_transcript_id]) or
                                (call.strand2 == '-' and
                                 start2 > cds_starts[acceptor_transcript_id])):
                            # If the 3' acceptor doesn't have a defnied start codon or if the
                            # breakpoint is in the gene but before the CDS start (in the UTR)
                            genomic_seq = get_genomic_sequence(start2, chroms[call.chrom2],
                                                               extend_length, call.strand2)
                            acceptor_in_utr = True
                        else:
                            acceptor_seq, genomic_seq, acceptor_in_intron = get_acceptor_junction(
                                start2, exons[acceptor_transcript_id],
                                cds_starts[acceptor_transcript_id], chroms[call.chrom2], peplen)
                    else:
                        if call.junctionSeq2 != '.':
                            logging.warning('Cannnot process %s as the 3\' donor sequence since it '
                                            'does not have an entry in the input GTF. A junction '
                                            'sequence was provided but not used since the 5\' '
                                            'partner has no junction sequence .',
                                            acceptor_transcript_id)
                        else:
                            logging.warning('Cannnot process %s as the 3\' donor sequence since '
                                            'it does not have an entry in the input GTF and a '
                                            'junction sequence was not provided.',
                                            acceptor_transcript_id)
                        continue

                in_frame_acceptor_seq = ''
                if donor_in_intron and not acceptor_in_intron:
                    if acceptor_in_utr:
                        # If the donor breakpoint was in an intron, the acceptor must be too else we
                        # cannot properly figure out splicing
                        logging.warning('%s-%s was predicted to have a 5\' breakpoint in an intron '
                                        'but the 3\' beakpoint was in the 5\' UTR of the acceptor '
                                        'gene. Cannot denovo assess splicing. Skipping.',
                                        donor_transcript_id, acceptor_transcript_id)
                    elif denovo_acceptor:
                        # If the donor breakpoint was in an intron, the acceptor must be too else we
                        # cannot properly figure out splicing
                        logging.warning('%s-%s was predicted to have a 5\' breakpoint in an intron '
                                        'but the 3\' beakpoint was in the 5\' UTR of the acceptor '
                                        'gene. Cannot denovo assess splicing. Skipping.',
                                        donor_transcript_id, acceptor_transcript_id)
                    else:
                        # If the donor breakpoint was in an intron, the acceptor must be too else we
                        # cannot properly figure out splicing
                        logging.warning('%s-%s was predicted to have a 5\' breakpoint in an intron '
                                        'but the 3\' beakpoint was not. Cannot denovo assess '
                                        'splicing. Skipping.', donor_transcript_id,
                                        acceptor_transcript_id)
                    continue
                elif (not donor_in_intron and
                          (acceptor_in_intron or acceptor_in_utr or denovo_acceptor)):
                    if donor_overhang == 0:
                        in_frame_acceptor_seq = genomic_seq[:(extend_length - 1) * 3]
                    else:
                        in_frame_acceptor_seq = genomic_seq[:(3 - donor_overhang +
                                                              (extend_length - 1) * 3)]
                else:
                    # donor in intron and acceptor in intron
                    # donor not in intron and acceptor not in intron)
                    if donor_overhang == 0:
                        in_frame_acceptor_seq = acceptor_seq[:(peplen - 1) * 3]
                    else:
                        in_frame_acceptor_seq = acceptor_seq[:(3 - donor_overhang +
                                                               (peplen - 1) * 3)]

                logging.debug('Inferred 3\' junction seq with breakpoint %s: %s ', call.start2,
                              in_frame_acceptor_seq)
                # Use the donor and acceptor sequences to generate the fusion transcript
                junction_seq = [in_frame_donor_seq, in_frame_acceptor_seq]

                logging.debug('Inferred junction seq: %s %s', junction_seq[0], junction_seq[1])

                fusion_peptide = translate(''.join(junction_seq))
                logging.debug("Inferred fusion peptide: %s" % fusion_peptide)

                if '*' in fusion_peptide:
                    fusion_peptide = fusion_peptide[:fusion_peptide.find('*')]

                if len(fusion_peptide) < peplen:
                    logging.warning('Fusion peptide %s from %s-%s (%s-%s) is is less than %s '
                                    'residues.', fusion_peptide, donor_name, acceptor_name,
                                    donor_transcript_id, acceptor_transcript_id, peplen)
                    continue

                gencode_donor_id = gene_transcripts[donor_transcript_id]
                if acceptor_name in gene_transcripts:
                    gencode_acceptor_id = gene_transcripts[acceptor_transcript_id]
                else:
                    gencode_acceptor_id = acceptor_transcript_id

                # Write fusion neoepitope sequence to FASTA file
                fasta = '>{donor_gene}-{acceptor_gene}_' \
                        '{donor_transcript}-{acceptor_transcript}_' \
                        '{donor_hugo}-{acceptor_hugo}_' \
                        '{donor_breakpoint}-{acceptor_breakpoint}_' \
                        'FUSION_{score}\n' \
                        '{sequence}\n'.format(donor_gene=gencode_donor_id,
                                              acceptor_gene=gencode_acceptor_id,
                                              donor_transcript=donor_transcript_id,
                                              acceptor_transcript=acceptor_transcript_id,
                                              donor_hugo=call.hugo1,
                                              acceptor_hugo=call.hugo2,
                                              donor_breakpoint=call.end1,
                                              acceptor_breakpoint=call.start2,
                                              score=call.score,
                                              sequence=fusion_peptide)
                outfile.write(fasta)
                found = True

    if not found:
        logging.warn("Did not generate any fusion peptides from data!")

    # Save alignments for different peptide lengths
    with open('transgene_fusion_alignments.pkl', 'w') as f:
        pickle.dump(previous_alignments, f)


def get_genomic_sequence(breakpoint, chrom_seq, extend_length, strand):
    breakpoint -= 1
    if strand == '+':
        genomic_sequence = chrom_seq[breakpoint:breakpoint + extend_length * 3]
    else:
        genomic_sequence = chrom_seq[breakpoint - extend_length * 3:breakpoint + 1]
        genomic_sequence = genomic_sequence[::-1].translate(trans)
    return genomic_sequence

