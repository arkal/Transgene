import collections
import logging

import pysam
from common import reject_decision


def reject_snv(snv, rna_bam=None, reject_threshold=None, rna_min_alt_freq=None, dna_bam=None,
               oxog_min_alt_freq=None):
    """
    Decide whether the mutation should be rejected based on the expression filter.

    :param dict snv: A single vcf SNV record split by the tab character and keyed by the fields
    :param str rna_bam: See reject_mutations:`rna_bam`
    :param int reject_threshold: See reject_mutations:`reject_threshold`
    :param float rna_min_alt_freq: See reject_mutations:`rna_min_alt_freq`
    :param str dna_bam: See reject_mutations:`dna_bam`
    :param float oxog_min_alt_freq: See reject_mutation:`oxog_min_alt_freq`

    :return: A named tuple of the True/False if the snv should be rejected or not, and a reason
    for rejection
    :rtype: reject_decision
    """
    bams = {'rna': rna_bam,
            'dna': dna_bam}
    output_counts, reads = get_snv_alignment_info(rna_bam, dna_bam, snv)
    reject = None
    possible_oxog_artefact = False
    low_rna_coverage = None
    rna_alt_freq = None
    coverage = '/'.join([','.join([str(reads[x]['rna']),
                                   str(reads[x]['dna'])]) for x in ('covering', 'spanning')])
    # Very importantly, we first look at dna and then rna
    for bam in 'dna', 'rna':
        # index
        if bams[bam] is None:
            continue
        if reads['spanning'][bam] == 0:
            logging.warning('Mutation at position %s:%s has no coverage in the %s-seq bam. '
                            'Rejecting.', snv['CHROM'], snv['POS'] + 1, bam.upper())
            reject = reject_decision(reject=True, reason='NoReadCoverage', coverage=coverage,
                                     vaf=0.0)
            break
        elif reads['covering'][bam] == 0 and bam == 'rna':
            # This concept of covering vs spanning reads is only in RNA-Seq.
            # I.e. you can have spanning reads but none of them covered the mutation (split).
            logging.warning('Mutation at position %s:%s appears to be in a region that has no '
                            'coverage in the RNA-seq bam, but has %s spanning reads. Possible '
                            'splicing event. Rejecting.', snv['CHROM'], snv['POS'] + 1,
                            reads['spanning'][bam])
            reject = reject_decision(reject=True, reason='SpliceDetected', coverage=coverage,
                                     vaf=0.0)
            break
        elif output_counts[bam][snv['ALT']]['counts'] == 0:
            # This branch will only be hit in RNA-seq cases since 0 Alt will probably never be
            # called by a mutations caller.  If we don't have too many reads over the location, it
            # could either be an undersampled location, or a splice.
            assert bam == 'rna'  # Sanity
            if (reads['covering'][bam] < reject_threshold and
                    reads['spanning'][bam] < reject_threshold):
                if possible_oxog_artefact:
                    # Flagged for RNA rescue but could not be rescued.
                    logging.debug('Mutation at position %s:%s has no evidence of existence in the '
                                  'RNA-seq bam. Possible OxoG variant. rejecting',
                                  snv['CHROM'], snv['POS'] + 1)
                    reject = reject_decision(reject=True, reason='OxoGArtefactLowAlt',
                                             coverage=coverage, vaf=0.0)
                    break
                else:
                    logging.debug('Mutation at position %s:%s has no evidence of existence in the '
                                  'RNA-seq bam. However the coverage at the region (%s/%s '
                                  'reads::Covering/spanning) is below the threshold for rejecting '
                                  '(%s). Accepting.', snv['CHROM'], snv['POS'] + 1,
                                  reads['covering']['rna'], reads['spanning']['rna'],
                                  reject_threshold)
                    low_rna_coverage = True
            else:
                logging.warning('Mutation at position %s:%s has no evidence of existence in the '
                                '%s-seq bam. Coverage = %s/%s reads (Covering/spanning). '
                                'Rejecting.', snv['CHROM'], snv['POS'] + 1, bam.upper(),
                                reads['covering'][bam], reads['spanning'][bam])
                reject = reject_decision(reject=True, reason='NoAltDetected', coverage=coverage,
                                         vaf=0.0)
                break
        else:
            alt_freq = (1.0 * output_counts[bam][snv['ALT']]['counts'] /
                        sum([output_counts[bam][x]['counts'] for x in output_counts[bam]]))
            if bam == 'dna':
                # This can only mean we want OxoG filtering
                if (snv['REF'], snv['ALT']) in (('C', 'A'), ('G', 'T')):
                    if alt_freq <= oxog_min_alt_freq:
                        logging.warning('Mutation at position %s:%s has less than %s ALT allele '
                                        'frequency (%s) in the %s-seq bam. Possible oxoG artefact. '
                                        'Flagging for RNA rescue.', snv['CHROM'], snv['POS'] + 1,
                                        oxog_min_alt_freq, round(alt_freq, 2), bam.upper())
                        possible_oxog_artefact = True
                    elif 0 in [output_counts[bam][snv['ALT']][x] for x in 'fr']:
                        # If either number of 'f'orward or 'r'everse reads is a 0 it means all reads
                        # came from only one strand and hence is possibly an oxoG artefact.
                        if output_counts[bam][snv['ALT']]['f']:
                            dominant_strand = 'forward'
                        else:
                            dominant_strand = 'reverse'
                        logging.warning('Mutation at position %s:%s is called from reads '
                                        'originating only from %s reads in the %s-seq bam. '
                                        'Possible oxoG artefact. Rejecting.', snv['CHROM'],
                                        snv['POS'] + 1, dominant_strand, bam.upper())
                        if dominant_strand == 'forward':
                            reject = reject_decision(reject=True, reason='OxoGArtefactAllFwd',
                                                     coverage=coverage, vaf=0.0)
                        else:
                            reject = reject_decision(reject=True, reason='OxoGArtefactAllRev',
                                                     coverage=coverage, vaf=0.0)
                        break
            else:  # rna
                if alt_freq < rna_min_alt_freq:
                    logging.warning('Mutation at position %s:%s has less than %s ALT allele '
                                    'frequency (%s) in the %s-seq bam. Rejecting.',
                                    snv['CHROM'], snv['POS'] + 1, rna_min_alt_freq,
                                    round(alt_freq, 4), bam.upper())
                    reject = reject_decision(reject=True, reason='LowAltFrequency',
                                             coverage=coverage, vaf=0.0)
                else:
                    rna_alt_freq = alt_freq
    if reject is None:
        # This means the call has not been rejected for any reason. It may have still have been
        # flagged so we need to handle that.
        if possible_oxog_artefact:
            if rna_bam:
                assert rna_alt_freq is not None
                logging.debug('Mutation at position %s:%s with has evidence in the RNA-Seq '
                              '(ALT frequency = %s) despite being flagged as a possible OxoG '
                              'artefect. Accepting', snv['CHROM'], snv['POS'] + 1, rna_alt_freq)
                reject = reject_decision(reject=False, reason='PossibleOxoGArtefactLowAltRNARescue',
                                         coverage=coverage, vaf=rna_alt_freq)
            else:
                logging.warning('Mutation at position %s:%s was flagged as possible oxoG artefact '
                                'and could not be rescued without matching RNA-Seq. Rejecting',
                                snv['CHROM'], snv['POS'] + 1)
                reject = reject_decision(reject=True, reason='OxoGArtefactLowAlt',
                                         coverage=coverage, vaf=0.0)
        else:
            reason = 'LowRNACoverage' if low_rna_coverage else '.'
            logging.debug('Accepted mutation at position %s:%s with %s read coverage and %s VAF.',
                          snv['CHROM'], snv['POS'] + 1, coverage,
                          round(rna_alt_freq, 2) if rna_alt_freq is not None else 'NA')
            reject = reject_decision(reject=False, reason=reason, coverage=coverage,
                                     vaf=0.0 if rna_alt_freq is None else round(rna_alt_freq, 2))
    return reject


def get_snv_alignment_info(rna_bam, dna_bam, vcf_record):
    """
    Get information about the spanning and covering reads at a given genomic locus for an snv.

    :param str rna_bam: Path to the RNA-Seq bam file
    :param str dna_bam: Path to the DNA-Seq bam file
    :param dict vcf_record: A single vcf record split by the tab character and keyed by the fields
    :return: information about the alignment at the position in the vcf record
             Output_counts- The counts of all nucleotides at the position in 'rna' and 'dna', and
             the counts of reads mapping to the 'f'orward and 'r'everse strands for each nucleotide.
             reads - The counts of reads mapping at ('covering') and across ('spanning') the
             position in 'rna' and 'dna'
    :rtype: tuple(collections.Counter|collections.defaultdict(Collections.Counter), dict(dict))
    """
    bams = {'rna': rna_bam, 'dna': dna_bam}
    output_counts = {'rna': collections.Counter(), 'dna': collections.Counter()}
    reads = {'spanning': {'rna': '.', 'dna': '.'},
             'covering': {'rna': '.', 'dna': '.'}}
    for bam in bams:
        if not bams[bam]:
            continue
        samfile = pysam.Samfile(bams[bam], 'rb')
        output_counts_ = collections.defaultdict(collections.Counter)
        covering_reads = 0
        spanning_reads = 0
        pileups = samfile.pileup(vcf_record['CHROM'], vcf_record['POS'], vcf_record['POS'] + 1,
                                 truncate=True)
        for pileup in pileups:
            for read in pileup.pileups:
                if not (read.is_refskip or read.is_del):
                    base = read.alignment.seq[read.query_position]
                    if ord(read.alignment.qual[read.query_position]) >= 63:  # (33PHRED + 30 QUAL)
                        output_counts_[base]['counts'] += 1
                        if read.alignment.is_read2:
                            output_counts_[base]['r'] += 1
                        else:
                            output_counts_[base]['f'] += 1
                    covering_reads += 1
                spanning_reads += 1
        samfile.close()
        reads['spanning'][bam] = spanning_reads
        reads['covering'][bam] = covering_reads
        output_counts[bam] = output_counts_
    return output_counts, reads
