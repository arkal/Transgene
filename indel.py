import logging

import pysam
from common import reject_decision


def reject_indel(indel, rna_bam=None, reject_threshold=None, rna_min_alt_freq=None):
    """
    Decide whether the mutation should be rejected based on the expression filter.

    :param dict indel: A single vcf INDEL record split by the tab character and keyed by the fields
    :param str rna_bam: See reject_mutations:`rna_bam`
    :param int reject_threshold: See reject_mutations:`reject_threshold`
    :param float rna_min_alt_freq: See reject_mutations:`rna_min_alt_freq`
    :return: A named tuple of the True/False if the snv should be rejected or not, and a reason
    for rejection
    :rtype: reject_decision
    """
    if rna_bam is None:
        logging.warning('REJECTING INDEL %s:%s>%s.', indel['POS'] + 1, indel['REF'], indel['ALT'])
        return reject_decision(reject=True, reason='NoGenomeFile', coverage='.,./.,.', vaf=0.0)
    else:
        output_counts, reads = get_indel_alignment_info(rna_bam, indel)
        reject = None
        coverage = '/'.join([','.join([str(reads[x]), '.']) for x in ('covering', 'spanning')])
        low_rna_coverage = False
        if output_counts['ALT'] == 0:
            if reads['covering'] < reject_threshold and reads['spanning'] < reject_threshold:
                logging.debug('Mutation at position %s:%s has no evidence of existence in the '
                              'RNA-seq bam. However the coverage at the region (%s/%s '
                              'reads::Covering/spanning) is below the threshold for rejecting '
                              '(%s). Accepting.', indel['CHROM'], indel['POS'] + 1,
                              reads['covering'], reads['spanning'], reject_threshold)
                reject = reject_decision(reject=False, reason='LowRNACoverage', coverage=coverage,
                                         vaf=0.0)
            else:
                logging.warning('Mutation at position %s:%s has no evidence of existence in the '
                                'RNA-seq bam. Coverage = %s/%s reads (Covering/spanning). '
                                'Rejecting.', indel['CHROM'], indel['POS'] + 1, reads['covering'],
                                reads['spanning'])
                reject = reject_decision(reject=True, reason='NoAltDetected', coverage=coverage,
                                         vaf=0.0)
        else:
            # There was some evidence of the indel
            alt_freq = 1.0 * output_counts['ALT'] / sum([output_counts[x] for x in output_counts])
            if alt_freq < rna_min_alt_freq:
                logging.warning('Mutation at position %s:%s has less than %s ALT allele '
                                'frequency (%s) in the RNA-seq bam. Rejecting.',
                                indel['CHROM'], indel['POS'] + 1, rna_min_alt_freq,
                                round(alt_freq, 2))
                reject = reject_decision(reject=True, reason='LowAltFrequency', coverage=coverage,
                                         vaf=0.0)
            else:
                logging.debug('Accepted mutation at position %s:%s with %s read coverage and %s '
                              'VAF.',
                              indel['CHROM'], indel['POS'] + 1, coverage, round(alt_freq, 2))
                reject = reject_decision(reject=False, reason='.', coverage=coverage, vaf=alt_freq)
        return reject


def get_indel_alignment_info(rna_bam, vcf_record):
    """
    Get information about the spanning and covering reads at a given genomic locus for an indel.

    :param str rna_bam: Path to the RNA-Seq bam file
    :param dict vcf_record: A single vcf record split by the tab character and keyed by the fields
    :return: see get_alignment_info:`return`
    :rtype: tuple(collections.Counter|collections.defaultdict(Collections.Counter), dict(dict))
    """

    indel_length = len(vcf_record['REF']) - len(vcf_record['ALT'])
    deletion = indel_length > 0  # If deletion is false, it has to be an insertion
    indel_length = abs(indel_length)

    samfile = pysam.Samfile(rna_bam, 'rb')

    reads = {'spanning': 0, 'covering': 0}
    output_counts = {'REF': 0, 'ALT': 0}

    pileups = samfile.pileup(vcf_record['CHROM'], vcf_record['POS'],
                             vcf_record['POS'] + 1, truncate=True)
    for pileup in pileups:
        for read in pileup.pileups:
            reads['spanning'] += 1
            if read.is_del and read.is_refskip:
                # This is a spliced alignment. We will not consider it
                continue
            reads['covering'] += 1
            if read.alignment.seq[read.query_position] != vcf_record['REF'][0]:
                logging.warning('read %s has %s at pos %s. expected %s', read.alignment.query_name,
                                read.alignment.seq[read.query_position], vcf_record['POS'] + 1,
                                vcf_record['REF'][0])
                continue
            if read.indel == 0:
                # No indel on this read
                output_counts['REF'] += 1
                continue
            if deletion:
                if read.indel == -indel_length:
                    if (ord(read.alignment.qual[read.query_position]) >= 63 and
                                ord(read.alignment.qual[read.query_position+1]) >= 63):
                        output_counts['ALT'] += 1
                else:
                    # TODO figure out left indel stuff
                    # Not a deletion of the required length (we will handle non-left-aligned
                    # events in the future.
                    continue
            else:  # Insertion
                if read.indel == indel_length:
                    quals = [ord(read.alignment.qual[read.query_position + x])
                             for x in range(0, indel_length+1)]
                    if sum(quals)/(indel_length+1) >= 63:
                        # On average, the bases must be Q30 or higher
                        output_counts['ALT'] += 1
                else:
                    # TODO figure out left indel stuff
                    # Not an insertion of the required length (we will handle non-left-aligned
                    # events in the future.
                    continue
        break  # Use only the first pileup position
    samfile.close()
    return output_counts, reads


def get_exon_start_pos(transcript_name, peptide_pos, exons, cds_starts):
    """
    Get the start position (1-based) of the exon containing `peptide_pos` in `peptide_name`.

    :param str transcript_name: The peptide name
    :param int peptide_pos: A protein position in a peptide
    :param dict exons: See return value of `transgene::get_exons`
    :param dict cds_starts: See return value of `transgene::get_exons`
    :return: The start position and a value identifying if the gene is on the positive strand
    :rtype: int, bool
    """
    positive_strand = exons[transcript_name][0].strand == '+'
    if positive_strand:
        transcript_exons = [(exon.start, exon.end) for exon in exons[transcript_name]
                            if exon.end >= cds_starts[transcript_name]]
        transcript_exons[0] = (cds_starts[transcript_name], transcript_exons[0][1])
    else:
        transcript_exons = [(exon.start, exon.end) for exon in exons[transcript_name]
                            if exon.start <= cds_starts[transcript_name]]
        transcript_exons[0] = (transcript_exons[0][0], cds_starts[transcript_name]+2)

    frame = 0  # identifies the starting frame for the current exon
    total = 0  # tracks the total number of AAs accounted for at the end of the prev exon
    for idx, exon in enumerate(transcript_exons):
        exon_len = frame + exon[1] - (exon[0] - 1)
        if total + exon_len/3 >= peptide_pos:
            # The codon lies in this exon
            diff = peptide_pos - total
            if diff == 1 and frame != 0:
                # This means the start is in the previous exon
                if positive_strand:
                    return transcript_exons[idx-1][1] - frame + 1, True
                else:
                    return transcript_exons[idx-1][0] + frame - 1, False
            else:
                if frame != 0:
                    diff -= 1  # codon starting previous exon
                if positive_strand:
                    return exon[0] + (diff - 1) * 3 + ((3 - frame) if frame else 0), True
                else:
                    return exon[1] - (diff - 1) * 3 - ((3 - frame) if frame else 0), False
        else:
            total += exon_len/3
            frame = exon_len % 3
