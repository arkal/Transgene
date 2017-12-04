#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : Transgene/transgene


Program info can be found in the docstring of the main function
"""
from __future__ import division, print_function

import argparse
import collections
import json
import logging
import os
import pysam
import random
import re
import shutil
import sys

from functools import partial
from multiprocessing import Manager, Pool
from operator import itemgetter
from tempfile import mkstemp

from common import read_fasta, chrom_sort, reject_decision, trans
from fusion import get_transcriptome_data, read_fusions, insert_fusions, get_exons
from snv import reject_snv


def reject_mutation(call, rna_bam=None, reject_threshold=None, rna_min_alt_freq=None, dna_bam=None,
                    oxog_min_alt_freq=None):
    """
    Decide whether the mutation should be rejected based on the expression filter.

    :param dict call: A single vcf record split by the tab character and keyed by the fields
    :param str rna_bam: The path to the RNA-seq file. Can be None suggesting that expression
           filtering should not be carried out.
    :param int reject_threshold: The read coverage above which we reject a no-ALT-coverage event
    :param float rna_min_alt_freq: The ALT allele frequency below which calls will be rejected as
           being having insufficient RNA evidence.
    :param str dna_bam:  The path to the tumor DNA-seq file. Can be None suggesting that oxoG
           filtering should not be carried out.
    :param float oxog_min_alt_freq: The ALT allele frequency below which calls will be rejected as
           being oxoG artifacts.

    :return: A named tuple of the True/False if the mutation should be rejected or not, and a reason
    for rejection
    :rtype: tuple(reject_decison)
    """
    # We need either a rna bam or a dna bam to operate
    assert rna_bam or dna_bam, "Cannot filter without at least one bam file."
    logging.info('Processing %s>%s mutation at position %s:%s', call['REF'], call['ALT'],
                 call['CHROM'], call['POS'] + 1)
    if len(call['REF']) == len(call['ALT']):
        return reject_snv(call, rna_bam, reject_threshold, rna_min_alt_freq, dna_bam,
                          oxog_min_alt_freq),
    else:
        if ',' in call['ALT']:
            # MAV.
            alts = call['ALT'].split(',')
            decisions = []
            for i, alt in enumerate(alts):
                # Modify this for each call to reject_snv
                call['ALT'] = alt
                decisions.append(reject_snv(call, rna_bam, reject_threshold, rna_min_alt_freq,
                                            dna_bam, oxog_min_alt_freq))
            # Return it to the original state
            call['ALT'] = ','.join(alts)
            return tuple(decisions)
        else:
            # Indel
            return reject_decision(reject=True, reason='INDEL', coverage='.,./.,.', vaf=0.0),


def write_to_vcf(data, out_vcf, header=False):
    """
    Write results to a vcf file
    :param dict data: The data line to write to the vcf
    :param file out_vcf: An open file handle to the output vcf
    :param header: Is this the field header
    :return: None
    """
    assert out_vcf.mode == 'w'
    # Hard coding the output order
    fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    if header:
        print('#', end='', sep='', file=out_vcf)
        print('\t'.join(fields), file=out_vcf)
    else:
        print('\t'.join([str(data[field]) if field in data else '.' for field in fields]),
              file=out_vcf)


def read_snpeffed_vcf(snpeff_file, rna_file=None, out_vcf=None, reject_threshold=None,
                      rna_min_alt_freq=None, dna_file=None, oxog_min_alt_freq=None, processes=1):
    """
    This module reads in the SNVs from the SnpEffed vcf file. It assumes that the SnpEff results
    have been writen to the INFO column. If the header contains a line that starts with  `#CHROM`,
    we will attempt to dynamically find the index of the `CHROM`, `POS`, `REF`, `ALT` and `INFO`
    columns.  If not, indexes 0, 1, 3, 4 and 7 (columns 1, 2, 4, 5 and 8) will be used.
    Example vcf column header:
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

    :param file snpeff_file: An open file handle to the snpeffed vcf file.
    :param str rna_file: The path to the RNA-seq file. Can be None suggesting that expression
           filtering should not be carried out.
    :param file out_vcf: An open file handle to the output vcf file
    :param int reject_threshold: The read coverage above which we reject a no-ALT-coverage event
    :param float rna_min_alt_freq: The ALT allele frequency below which calls will be rejected as
           being having insufficient RNA evidence.
    :param str dna_file:  The path to the tumor DNA-seq file. Can be None suggesting that oxoG
           filtering should not be carried out.
    :param float oxog_min_alt_freq: The ALT allele threshold below which calls will be rejected as
           being oxoG artifacts.
    :param int processes: How many processes should handle the parsing?
    :returns: A parsed dictionary of snvs to be processed.
    :rtype: tuple(collections.Counter, collections.Counter)
    """
    if out_vcf:
        assert out_vcf.mode == 'w', 'Output vcf file is not open for writing'
    logging.info('Reading in Mutations from %s', snpeff_file.name)
    indexes = {'CHROM': 0,
               'POS': 1,
               'REF': 3,
               'ALT': 4,
               'INFO': 7}
    # Set up the Manager, and the shared dict and queue
    manager = Manager()
    vcf_queue = manager.Queue()
    mutations = manager.dict()
    # Read in the vcf header and identify the columns in the file
    for line in snpeff_file:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                # This is the last line of the header.
                line = line.strip().lstrip('#').split('\t')
                temp_indexes = {val: pos for pos, val in enumerate(line)}
                for column in indexes.keys():
                    if column in temp_indexes:
                        logging.debug('Dynamically identified %s column from vcf header to be '
                                      'column %s.', column, temp_indexes[column] + 1)
                        indexes[column] = temp_indexes[column]
                    else:
                        logging.warning('Could not dynamically identify %s column from vcf header. '
                                        'Using column %s.', column, indexes[column] + 1)
                if rna_file or dna_file:
                    print('##INFO=<ID=reason,Number=.,Type=String,Description=Reason for rejecting/'
                          'accepting a mutation.', file=out_vcf)
                    print('##INFO=<ID=coverage,Number=1,Type=String,Description=Reads supporting '
                          'the decision (Covering_RNA,Covering_DNA/Spanning_RNA,Spanning_DNA).',
                          file=out_vcf)
                    print('##INFO=<ID=vaf,Number=1,Type=String,Description=The >Q30 variant '
                          'allele frequency for the mutation in the RNA-Seq bam.', file=out_vcf)
                    write_to_vcf(indexes, out_vcf, header=True)
            elif rna_file or dna_file:
                print(line, end='', file=out_vcf)
            continue
        elif len(line.strip()) == 0:
            continue
        else:
            # Populate the queue.
            vcf_queue.put(line)
    for _ in range(0, processes):
        # Add the graceful worker shutdown signal to the queue for every worker.
        vcf_queue.put(None)
    # Start the workers
    logging.info('Starting %s workers for vcf parsing.' % processes)
    pool = Pool(processes=processes)
    parse_vcf_line_partial = partial(parse_vcf_line,
                                     queue=vcf_queue,
                                     indexes=indexes,
                                     out_calls=mutations,
                                     rna_file=rna_file,
                                     reject_threshold=reject_threshold,
                                     rna_min_alt_freq=rna_min_alt_freq,
                                     dna_file=dna_file,
                                     oxog_min_alt_freq=oxog_min_alt_freq)
    pool.map(parse_vcf_line_partial, range(0, processes))
    pool.close()
    pool.join()
    # Now that all the processes have completed, we need to process the worker specific vcfs
    worker_vcfs = {i: '.worker_%s.vcf' % i for i in range(0, processes)}
    if rna_file or dna_file:
        if processes == 1:
            # The queue is ordered and the lone process would have run the lines in order so we can
            # just copy the contents over as-is without spending compute on sorting, etc.
            with open(worker_vcfs[0]) as v_file:
                for line in v_file:
                    print(line, end='', file=out_vcf)
        else:
            chrom_lines = collections.defaultdict(dict)
            for vcf in worker_vcfs:
                with open(worker_vcfs[vcf]) as v_file:
                    for line in v_file:
                        line_ = line.strip().split()
                        chrom_lines[line_[indexes['CHROM']]][int(line_[indexes['POS']])] = line
            for chrom in chrom_sort(chrom_lines.keys()):
                for pos in sorted(chrom_lines[chrom].keys()):
                    print(chrom_lines[chrom][pos], end='', file=out_vcf)
    for vcf in worker_vcfs:
        # Delete the temp vcfs now that we're done with them.
        os.remove(worker_vcfs[vcf])

    # Merge the dict structure into a dict of dicts
    out_muts = collections.defaultdict(dict)
    non_syn_seen = False

    for transcript, mut, tlen in mutations.keys():
        out_muts[transcript][mut] = mutations[(transcript, mut, tlen)]
        if (not non_syn_seen and (mut[0] != mut[-1] and mut[-1] != '*')):
            # If we haven't seen a non_synonymous mutation, and this is non synonymous, update the
            # flag
            non_syn_seen = True
        # probably overkill
        if 'len' in out_muts[transcript]:
            assert out_muts[transcript]['len'] == int(tlen)
        else:
            out_muts[transcript]['len'] = int(tlen)

    if not (out_muts and non_syn_seen):
        logging.error('Input snpeffed mutations file was empty or had no actionable mutations.')
        raise RuntimeError('Input snpeffed mutations file was empty or had no actionable '
                           'mutations.')

    return dict(out_muts)


def get_ref_pos_alt_aa(aa_change):
    """
    Determine the reference AA(s), the position in the protein, and the alternate AA(s) for a given
    mutation in the form A123B

    :param aa_change: A String describing a mutation in protein space
    :return: reference AA(s), position, Alternamte AA(s)
    :rtype: tuple(str, int, str)
    """
    # We expect the mutation to be 3-parts, non-numeric, numeric, non-numeric or NA
    expected_isdigit = False
    parts = ['']
    for char in aa_change:
        if char.isdigit() is expected_isdigit:
            parts[-1] += char
        else:
            parts.append(char)
            expected_isdigit = not expected_isdigit
    parts[1] = int(parts[1])
    if len(parts) == 2:
        parts.append('')
    return parts


def get_codon_aa_change(codon_change):
    """
    For a given codon change E.g. aAt/aGt, return the changed base ([A,G] in teh example).
    :param str codon_change: Thecodon change
    :return: The modified ref and alt
    :rtype: tuple
    """
    ref, alt = codon_change.split('/')
    for i in range(3):
        # range(3) is ok because this will only happen with MAV
        if ref[i] == alt[i]:
            continue
        else:
            return ref[i], alt[i]
    assert False


def parse_vcf_line(worker_id, queue, indexes, out_calls, rna_file=None, reject_threshold=None,
                   rna_min_alt_freq=None, dna_file=None, oxog_min_alt_freq=None):
    """
    Parse one mutation-containing line in the input vcf.

    :param int worker_id: The id for this worker
    :param multiprocessing.Manager.Queue queue: A queue that will supply a line of the vcf as every
           item
    :param dict indexes: A dict indicating the fields in the vcf line.
    :param multiprocessing.Manager.dict out_calls: The shared output vcf that will contain the
           output calls.
    :param str rna_file: The path to the RNA-seq file. Can be None suggesting that expression
           filtering should not be carried out.
    :param int reject_threshold: The read coverage above which we reject a no-ALT-coverage event
    :param float rna_min_alt_freq: The ALT allele frequency below which calls will be rejected as
           being having insufficient RNA evidence.
    :param str dna_file:  The path to the tumor DNA-seq file. Can be None suggesting that oxoG
           filtering should not be carried out.
    :param float oxog_min_alt_freq: The ALT allele threshold below which calls will be rejected as
           being oxoG artifacts.
    """
    logging.info('VCF parsing worker %s is up and running.' % worker_id)
    with open('.worker_%s.vcf' % worker_id, 'w') as out_vcf:
        while True:
            line = queue.get(timeout=30)
            if line is None:
                break
            line = line.strip().split('\t')
            line = {x: line[y] for x, y in indexes.items()}

            alt_alleles = line['ALT'].split(',')
            changes = line['INFO'].split(';')[-1]
            if rna_file or dna_file:
                # NOTE: vcf snvs are 1-based and not 0-based!
                line['POS'] = int(line['POS']) - 1
                if changes.startswith('EFF'):
                    decisions = reject_mutation(line, rna_file, reject_threshold, rna_min_alt_freq,
                                                dna_file, oxog_min_alt_freq)
                else:
                    logging.warning('Mutation at position %s:%s corresponds to a non-exonic '
                                    'region. Rejecting.', line['CHROM'], line['POS'] + 1)
                    decisions = (reject_decision(reject=True, reason='NonExonic', coverage='NA',
                                                 vaf='0.0'),)
                # Print the line to the output vcf
                line['INFO'] = ';'.join([line['INFO'],
                                         'reason=' + ','.join(d.reason for d in decisions),
                                         'coverage=' + str(decisions[0].coverage),
                                         'VAF=' + ','.join(str(d.vaf) for d in decisions)])
                line['FILTER'] = 'REJECT' if all(d.reject for d in decisions) else 'PASS'
                # Temporarily change back to 1-based for printing the vcf record
                line['POS'] += 1
                write_to_vcf(line, out_vcf)
                line['POS'] -= 1
                # If all ALT alleles were rejected, we can skip to the next call
                if all(d.reject for d in decisions):
                    continue
            else:
                decisions = tuple([reject_decision(False, None, None, None) for _ in alt_alleles])
            changes = re.sub('EFF=', '', changes)
            changes = [x for x in changes.split(',')
                       if x.startswith((
                            'NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING', 'STOP_GAINED',  # SNVs
                            )) and 'protein_coding' in x]

            for di, decision in enumerate(decisions):
                if decision.reject:
                    # happens if one of the ALT alleles was rejected but the other wasn't
                    continue
                for i in range(0, len(changes)):
                    temp = changes[i].split('|')
                    # MAVs will have a mix of both ALTs in the EFF line so we need to ensure
                    # this change corresponds to the current ALT
                    codon_aa_change = get_codon_aa_change(temp[2])
                    if codon_aa_change not in [(line['REF'], alt_alleles[di]),
                                               (line['REF'].translate(trans),
                                                alt_alleles[di].translate(trans))]:
                        continue
                    ref_aa, pos, alt_aa = get_ref_pos_alt_aa(temp[3])
                    if not alt_aa:
                        # Synonymous change
                        alt_aa = ref_aa
                    if alt_aa.endswith('?'):
                        assert ref_aa == '?' and temp[-1] == 'WARNING_TRANSCRIPT_INCOMPLETE)', \
                            '? seen in an SNV for the ALT AA'

                    out_calls[(temp[8], temp[3], temp[4])] = {
                        'AA': {'REF': ref_aa,
                               'ALT': alt_aa,
                               'change': temp[2],
                               'POS': pos},
                        'NUC': {'REF': line['REF'],
                                'ALT': alt_alleles[di],
                                'POS': line['POS'],
                                'CHROM': line['CHROM']},
                    }
    logging.info('VCF parsing worker %s received signal to go down.' % worker_id)


def build_groups(mut_list, base_mut=None, out_list=None):
    """
    Recursively propagates a binary tree of muts

    :param list base_mut: The base for the mutation
    :param list mut_list: The remaining mutations to add to the tree
    :param list out_list: The list to append values to
    :return: a list of lists containing all the possible mutations

    >>> build_groups(['a', 'b'])
    [('a', 'b'), ('a',), ('b',)]
    >>> build_groups(['a'], ['b'])
    [('b', 'a'), ('b',)]
    """
    # Since we can't make mutables the default parameters
    if base_mut is None:
        base_mut = []
    if out_list is None:
        out_list = []
    if not mut_list:
        # base case
        if base_mut:
            out_list.append(tuple(base_mut))
        return None
    else:
        new_base = [mut_list[0]]
        build_groups(mut_list[1:], base_mut + new_base, out_list)
        build_groups(mut_list[1:], base_mut, out_list)
        return out_list


def fix_group(grp, phased_groups):
    """
    This is used to fix a group that has a 0 in it, using the other groups
    :param tuple grp: The group to be fixed
    :param set phased_groups: a set of tuples of other groups
    :return:
    """
    index = grp.index(0)
    left_index = index - 1 if index > 0 else None
    right_index = index + 1 if index < len(grp)-1 else None
    possible_values = set([x[index]
                           for x in phased_groups
                           if ((left_index is not None and
                                x[left_index] == grp[left_index]) or
                               left_index is None) and
                           ((right_index is not None and
                             x[right_index] == grp[right_index]) or
                            right_index is None)])
    if not possible_values:
        logging.debug('Could not fix a group.', str(grp))
        return None
    else:
        grp = list(grp)
        for val in possible_values:
            if tuple(grp[:index] + [val] + grp[index+1:]) in phased_groups:
                return None
        # If we reach here, we haven't found a perfect match in phased_groups so we create one using
        # one of the possible values for the mutation.
        grp = tuple(grp[:index] + [random.choice(list(possible_values))] + grp[index+1:])
    if 0 in grp:
        return fix_group(grp, phased_groups)
    else:
        return grp


def get_mutation_groups(p_snvs, peplen, pfasta_name, rna_bam=None):
    """
    Produces a list of mutation groups for the given peptide

    :param dict p_snvs: parsed snvs from the input file
    :param int peplen: The length of the peptide
    :param str pfasta_name: The full name of the protein fasta (for debugging purposes)
    :param str rna_bam: path to the bam file
    :return: mutation groups
    :rtype: dict
    """
    mutations = sorted(p_snvs.keys(), key=lambda p_mut: int(p_mut[1:-1]))
    groups = []
    while True:
        group_start = mutations[0]
        groups.append([group_start] + [y for y in mutations[1:]
                                       if int(group_start[1:-1]) + peplen >
                                       int(y[1:-1]) >=
                                       int(group_start[1:-1])])
        for i in groups[-1][1:]:
            mutations.remove(i)
        mutations.remove(group_start)
        if not mutations:
            break
    mutation_groups = [build_groups(x) for x in groups]
    if rna_bam is None:
        logging.debug('Peptide %s detected to have the following mutation groups %s', pfasta_name,
                      mutation_groups)
        return mutation_groups
    else:
        out_groups = []
        for group in mutation_groups:
            # We will parse for co-expressed mutations
            if len(group) == 1:
                # This is guaranteed to have coverage because we have already done a coverage filter
                out_groups.append([list(group[0])])  # List to make it consistent with bigger groups
                continue
            samfile = pysam.AlignmentFile(rna_bam, 'rb')
            mutations = group[0]
            if p_snvs[mutations[0]]['NUC']['POS'] > p_snvs[mutations[-1]]['NUC']['POS']:
                reverse_strand = True
            else:
                reverse_strand = False
            # This will become a dict with each mutation being a key and value being an empty list.
            # The list items will correspond to the event seen at the mutational position in each
            # read (0 = No coverage; 1 = REF; 2 = ALT)
            mut_relations = {x: [] for x in mutations}
            # Get the region spanning the group
            if reverse_strand:
                alignment = samfile.fetch(p_snvs[mutations[0]]['NUC']['CHROM'],
                                          p_snvs[mutations[-1]]['NUC']['POS'],
                                          p_snvs[mutations[0]]['NUC']['POS'])
            else:
                alignment = samfile.fetch(p_snvs[mutations[0]]['NUC']['CHROM'],
                                          p_snvs[mutations[0]]['NUC']['POS'],
                                          p_snvs[mutations[-1]]['NUC']['POS'])
            num_spanning_reads = 0
            for read in alignment:
                spanned_muts = [x for x in mutations if
                                read.positions[0] <= p_snvs[x]['NUC']['POS'] <= read.positions[-1]]
                for mut in spanned_muts:
                    pos = [x for x, y in read.get_aligned_pairs() if
                           y == p_snvs[mut]['NUC']['POS'] and x is not None]
                    # Make this 0 by default and then modify if we see otherwise. The else clause
                    # for the following 3 if's are handled by this line.
                    mut_relations[mut].append(0)
                    if pos:
                        pos = pos[0]
                        if read.seq[pos] is not None and read.qual[pos] > 30:
                            if read.seq[pos] == p_snvs[mut]['NUC']['ALT']:
                                mut_relations[mut][-1] = 2
                            elif read.seq[pos] == p_snvs[mut]['NUC']['REF']:
                                mut_relations[mut][-1] = 1
                for mut in mutations:
                    # Add 0s to the unspanned positions
                    if mut in spanned_muts:
                        continue
                    else:
                        mut_relations[mut].append(0)
                num_spanning_reads += 1
            assert len(mut_relations[mutations[0]]) == num_spanning_reads
            if num_spanning_reads == 0:
                logging.warning('Could not find spanning reads for the mutations (%s). Using all '
                                'combinations', ','.join(mutations))
                out_groups.append(group)
            else:
                # Get a list of how the mutations were phased per read and then process.
                phased_groups = set([tuple([mut_relations[mut][x] for mut in mutations])
                                     for x in xrange(0, num_spanning_reads)])
                # Throw away the useless entries. Useless entries are those that are all 0's (no
                # coverage in the region) and all 1's (No ALT seen in the read)
                useless_groups = [grp for grp in phased_groups
                                  if sum(grp) == 0 or (0 not in grp and 2 not in grp)]
                for grp in useless_groups:
                    phased_groups.remove(grp)
                # Now cleanup the ones with 0's in them.
                # We can fix 0's by phasing with the nearest neighbor from other groups
                final_groups = []
                for grp in phased_groups:
                    if 0 in grp:
                        grp = fix_group(grp, phased_groups)
                        if grp is None:
                            continue
                    final_groups.append([mutations[x] for x, y in enumerate(grp) if y != 1])
                out_groups.append(final_groups)
        logging.debug('Peptide %s detected to have the following mutation groups %s', pfasta_name,
                      out_groups)
        return out_groups


def merge_adjacent_snvs(snvs, muts):
    """
    Given 2-3 mutations that modify the same codon, find the actual modified AA.

    :param dict snvs: The collection of snvs in the protein
    :param list muts: The group of mutations found to be coexpressed
    :return: The corrected tuple of mutations in the group
    :rtype: tuple

    >>> snvs = {'R148R': {'AA': {'change': 'Cgg/Agg'}}, 'R148Q': {'AA': {'change': 'cGg/cAg'}}}
    >>> merge_adjacent_snvs(snvs, ['R148R', 'R148Q'])
    'R148K
    >>> snvs = {'A123B': {'NUC': {'ALT': 'C'}}, \
                'C123D': {'NUC': {'ALT': 'A'}}, \
                'E123F': {'NUC': {'ALT': 'T'}}}
    >>> merge_adjacent_snvs(snvs, ['A123B', 'C123D','E123F'])
    'A123H'
    """
    from fusion import genetic_code
    if len(muts) == 3:
        return muts[0][:-1] + genetic_code[''.join([snvs[x]['NUC']['ALT'] for x in muts])]
    else:
        # Has to be 2 mutations
        changes = zip(snvs[muts[0]]['AA']['change'][4:], snvs[muts[1]]['AA']['change'][4:])
        out_mut = ''
        for change in changes:
            if change[0] == change[1]:
                out_mut += change[0]
            else:
                out_mut += change[0] if change[0].isupper() else change[1]
        return muts[0][:-1] + genetic_code[out_mut.upper()]


def correct_for_same_codon(snvs, mut_group):
    """
    If two mutations in the group modify the same codon, snpeff calls them separately and the naive
    rna filter in Transgene also handles them independently. This method corrects the calls if they
    are found to be on the same allele.

    :param dict(str, dict(str, str|int)) snvs: The collection of snvs in the protein
    :param tuple(stri) mut_group: The group of mutations found to be coexpressed
    :return: The corrected tuple of mutations in the group
    :rtype: tuple(str)

    >>> snvs = {'A123B': {'AA': {'change': 'Xxx/Yxx', 'POS': 123}}, \
                'C124D': {'AA': {'change': 'Zzz/Azz', 'POS': 124}}}
    >>> correct_for_same_codon(snvs, ('A123B',))
    ('A123B',)
    >>> correct_for_same_codon(snvs, ('A123B', 'C124D'))
    ('A123B', 'C124D')
    >>> snvs = {'R148R': {'AA': {'change': 'Cgg/Agg', 'POS': 148}, 'NUC': {'POS': 12345}}, \
                'R148Q': {'AA': {'change': 'cGg/cAg', 'POS': 148}, 'NUC': {'POS': 12346}}, \
                'R149X': {'AA': {'change': 'xXx/xYx', 'POS': 149}, 'NUC': {'POS': 12349}}}
    >>> correct_for_same_codon(snvs, ('R148R', 'R148Q', 'R149X'))
    ('R148K', 'R149X')
    >>> snvs = {'R148R': {'AA': {'change': 'Cgg/Agg', 'POS': 148}, 'NUC': {'POS': 12345}}, \
                'R148W': {'AA': {'change': 'Cgg/Tgg', 'POS': 148}, 'NUC': {'POS': 12345}}, \
                'R149X': {'AA': {'change': 'xXx/xYx', 'POS': 149}, 'NUC': {'POS': 12349}}}
    >>> correct_for_same_codon(snvs, ('R148R', 'R148W', 'R149X'))
    ()
    """
    if len(mut_group) == 1 or len(mut_group) == len({snvs[x]['AA']['POS'] for x in mut_group}):
        # If there is only 1 mutation or if the positions are unique
        return mut_group
    else:
        out_list = []
        # A  set of positions in this group
        positions = {snvs[x]['AA']['POS'] for x in mut_group}
        # These groups are going to be small so we aren't losing much by brute-forcing this.
        for pos in sorted(positions):
            muts = [mut for mut in mut_group if str(pos) in mut]
            if len(muts) == 1:
                out_list.extend(muts)
            elif len({snvs[mut]['NUC']['POS'] for mut in muts}) == 1:
                # MAVs affect the same codon AND the same nucleotide. They CANNOT be in the same
                # group
                logging.info('Skipping a group (%s) containing a Multi-Allelic Variant.', muts)
                return ()
            else:
                logging.info('Identified co-expressed mutations (%s)that affect the same '
                             'codon/AA.  Merging.', muts)
                merged_group = merge_adjacent_snvs(snvs, muts)
                out_list.append(merged_group)
                logging.info('Merged (%s) into %s.', muts, merged_group)
        return tuple(out_list)


def insert_mutations(protein_fa, mutations, tumfile, normfile, peplen, rna_bam=None, chroms=None):
    """
    This module uses the mutation data contained in `mutations` and inserts them into the
    genome contained in `chroms`.

    :param dict protein_fa: Contains the peptides in the form of a dict where keys hold the protein
           name and the values holds the sequence.
    :param dict mutations: Contains the snvs parsed from teh input SnpEFF file.
    :param file tumfile: The file to write tumor output to.
    :param file normfile: The file to write normal output to.
    :param int peplen: Length of peptides which will be generated from the output file.
    :param str rna_bam: The path to the RNA-Seq bam file
    :param dict chroms: Contains the chromosomal dna int he form of a dict where keys hold the
           chromosome name and the values holds the sequence.
    """
    logging.info('Inserting Mutations into IARs')
    for pept in mutations.keys():
        # First, grab the positions of all mutations in the peptide.  If there are multiple
        # mutations, then for each mutation, list the other mutations that are within peplen
        # positions of it IN THE POSITIVE DIRECTION.  The mutations are potentially co-expressed on
        # the same peptide. The reason can be visually shown for a 9-mer peptide with three
        # mutations, X, Y and Z as below (O = non-mutated AA).  The region depicts the extended IAR
        # for the 3. Possible mutational combinations are shown. The WT (O,O,O) isn't shown since it
        # isn't useful.
        #
        # POSITIONS:    a b c d e f g h i j k l m n o p q r s t u v w x y
        #
        # X+Y+Z        O O O O O O O O X O O O Y O O O Z O O O O O O O O
        # X + Z        O O O O O O O O X O O O O O O O Z O O O O O O O O
        # X + Y        O O O O O O O O X O O O Y O O O O O O O O O O O O
        # X only       O O O O O O O O X O O O O O O O O O O O O O O O O
        #
        # Y + Z        O O O O O O O O O O O O Y O O O Z O O O O O O O O
        # Y only       O O O O O O O O O O O O Y O O O O O O O O O O O O
        #
        # Z only       O O O O O O O O O O O O O O O O Z O O O O O O O O
        #
        # The data for what mutations are grouped together, and in what combinations is stored in
        # the mutation_groups list object. Each group in in the list is a list of all mutations in
        # the group.  If there are no neighboring mutations, the group contains only the mutation
        # itself.
        pfasta_name = [x for x in protein_fa.keys() if ''.join([pept, '_']) in x][0]
        # We have all the mutations and the protein sequence here. Do a sanity check on the results.
        protein = protein_fa[pfasta_name]
        if 'len' in mutations[pept]:
            # If 'len' is not in this dict, it means a previous process has already cleaned up the
            # dict. This will need to change if this module is multithreaded.
            expected_tlen = mutations[pept].pop('len')
            for mutation in mutations[pept].keys():
                mut_pos = int(mutation[1:-1])
                if protein[mut_pos - 1] != mutations[pept][mutation]['AA']['REF']:
                    if mutations[pept][mutation]['AA']['REF'] == mutations[pept][mutation]['AA']['ALT']:
                        # Let synonymous mutations pass through
                        continue
                    logging.warning('%s seen at position %s in %s.... %s expected.',
                                    protein_fa[pfasta_name][mut_pos - 1], mut_pos, pept,
                                    mutations[pept][mutation]['AA']['REF'])
                    if len(protein_fa[pfasta_name]) - expected_tlen in (1, 2):
                        # This is a case seen often with SNPEff where the n+1 position will have
                        # the mutation.  Check the n+1 residue
                        if protein[mut_pos] == mutations[pept][mutation]['AA']['REF']:
                            logging.info('Mutation at position %s in %s detected to exist at '
                                         'position %s.', mut_pos, pept, mut_pos + 1)
                            new_mutation = ''.join([mutations[pept][mutation]['AA']['REF'],
                                                    str(mut_pos + 1),
                                                    mutations[pept][mutation]['AA']['ALT']])
                            mutations[pept][new_mutation] = mutations[pept][mutation]
                            mutations[pept].pop(mutation)
                    else:
                        logging.warning('Cannot handle mutation at position %s in %s. '
                                        'Disregarding.', mut_pos, pept)
                        mutations[pept].pop(mutation)
        if not mutations[pept]:
            continue
        all_mutation_groups = get_mutation_groups(mutations[pept], peplen, pfasta_name, rna_bam)

        for mut_group in all_mutation_groups:
            for group_muts in mut_group:
                group_muts = correct_for_same_codon(mutations[pept], group_muts)
                # After correcting, if there is only ONE mutation and it is a stop gain, disregard
                # this call
                if len(group_muts) == 1 and group_muts[0][-1] == '*':
                    continue
                group_muts = tuple([x for x in group_muts if x[0] != x[-1]])
                if len(group_muts) == 0:
                    # If the group only consisted of one synonymous mutation, or was rejected for
                    # having an MAV.
                    continue
                out_pept = {'pept_name': [pfasta_name],
                            'tum_seq': [],
                            'norm_seq': []
                            }
                mut_pos = int(group_muts[0][1:-1])
                prev = max(mut_pos - peplen, 0)  # First possible AA in the IAR
                for group_mut in group_muts:
                    mut_pos = int(group_mut[1:-1])
                    out_pept['pept_name'].append(group_mut)
                    out_pept['tum_seq'].extend(protein[prev:mut_pos - 1])
                    out_pept['norm_seq'].extend(protein[prev:mut_pos - 1])
                    if group_mut[-1] == '*':
                        prev = None
                        break
                    out_pept['tum_seq'].append(group_mut[-1])
                    out_pept['norm_seq'].append(group_mut[0])
                    prev = mut_pos
                if prev is not None:
                    out_pept['tum_seq'].extend(protein[prev:min(mut_pos + peplen - 1,
                                                                len(protein))])
                    out_pept['norm_seq'].extend(protein[prev:min(mut_pos + peplen - 1,
                                                                 len(protein))])
                write_pepts_to_file(out_pept, tumfile, normfile, peplen)
    return None


def write_pepts_to_file(peptide_info, tumfile, normfile, peplen):
    """
    This module will accept the info for a given peptide including the name, the tumor peptide
    sequence, and the corresponding normal peptide sequence, and will print then to tumfile and
    normfile
    """
    blacklist = set(list('BJOUXZ'))
    peptide_name = '_'.join(peptide_info['pept_name'])
    tum_peptide_seq = ''.join(peptide_info['tum_seq'])
    norm_peptide_seq = ''.join(peptide_info['norm_seq'])
    truncated = False
    if not blacklist.isdisjoint(tum_peptide_seq):
        # If the final peptide contains a blacklisted amino acid,
        # handle it.
        logging.warning('Blacklisted amino acid(s) %s seen in IAR %s... Truncating.',
                        list(set(list(tum_peptide_seq)).intersection(blacklist)), peptide_name)
        matches = [i for i, x in enumerate(tum_peptide_seq) if x in blacklist]
        truncated = True
        if len(matches) > 1:
            # Out of 94285 transcripts in gencode 25, Only 14 have multiple blacklisted AAs. Of
            # these 14, only 6 transcripts from 4 genes (3 + 1 + 1 + 1) have the possibility of >2
            # blacklisted AAs being in the same single-mutation 9-, 10- and 15-mer, and double 9-
            # and 10-mers. Another transcript from another gene (total 7) has a possibility of >2
            # blacklisted AAs in a  double-mutation 15-mer.  It seems logically ok to just disregard
            # this edge case for now.
            logging.warning('Cannot handle multiple blacklisted AAs per IAR.')
            return None
        else:
            # If there is one truncation point, you get 2 piece of either the same size (necessarily
            # both are less than peplen in size), or dissimilar size where if the larger piece is
            # greater than or equal to peplen in size, it must necessarily contain the mutant.
            # Pieces smaller than peplen are neglected.
            pos = matches[0]
            if pos > len(tum_peptide_seq)/2:
                tum_peptide_seq = tum_peptide_seq[:pos]
                norm_peptide_seq = norm_peptide_seq[:pos]
            else:
                tum_peptide_seq = tum_peptide_seq[pos+1:]
                norm_peptide_seq = norm_peptide_seq[pos+1:]

    if len(tum_peptide_seq) < peplen:
        prefix = 'Truncated ' if truncated else ''
        logging.warning((prefix + 'IAR %s is less than %s residues (%s).'), peptide_name, peplen,
                        tum_peptide_seq)
    else:
        print('>', peptide_name, sep='', file=tumfile)
        print('>', peptide_name, sep='', file=normfile)
        print(tum_peptide_seq, file=tumfile)
        print(norm_peptide_seq, file=normfile)
    return None


def parse_peptides(tumfile, normfile, prefix, suffix):
    """
    This module takes in a peptides file and squashes it into the minimum number of peptides
    required to describe the potential neo-immunopeptidome.  It's main function is to take
    transcript-level mutation calls and merge them by gene if the sequences they contain are
    identical (to correct for splicing information).
    """
    with open(normfile, 'r') as i_f:
        normal_peptides = {}
        for pep_name, _, pep_seq in read_fasta(i_f, 'ARNDCQEGHILKMFPSTWYVBZJUOX'):
            normal_peptides[pep_name] = pep_seq
    with open(tumfile, 'r') as i_f:
        peptides = collections.Counter()
        for full_pep_name, _, pep_seq in read_fasta(i_f, 'ARNDCQEGHILKMFPSTWYVBZJUOX'):
            pep_name = full_pep_name.split('_')
            gene_name = pep_name[0]
            hugo_gene = pep_name[2]
            transcript_mutation = '_'.join([pep_name[1]]+pep_name[3:])
            if peptides[(gene_name, hugo_gene)] == 0:
                try:
                    peptides[(gene_name, hugo_gene)] = {transcript_mutation:
                                                            (pep_seq, normal_peptides[full_pep_name])}
                except KeyError:
                    peptides[(gene_name, hugo_gene)] = {transcript_mutation: (pep_seq, None)}
            else:
                try:
                    peptides[(gene_name, hugo_gene)][transcript_mutation] = (pep_seq, normal_peptides[full_pep_name])
                except KeyError:
                    peptides[(gene_name, hugo_gene)][transcript_mutation] = (pep_seq, None)
    outmap = {}
    out_tum_file = '_'.join([prefix, 'tumor', suffix])
    out_norm_file = '_'.join([prefix, 'normal', suffix])
    with open(out_tum_file, 'w') as t_f, open(out_norm_file, 'w') as n_f:
        # Mhc predictors can't handle the giant peptide names created by merging the names of
        # peptide emerging from transcripts bearing similar mutations hence we use an incremental
        # system to name the peptides and store the peptide number to actual peptide name
        # information in a map file.
        peptide_number = 1
        for gene in peptides.keys():
            unique_seqs = set(peptides[gene].values())
            for group in [[x for x, y in peptides[gene].items() if y[0] == z[0]]
                          for z in unique_seqs]:
                pepname = ''.join(['neoepitope_', str(peptide_number)])
                group_info = '\t'.join(list(gene) + [','.join(group)])
                group_tum_seq = peptides[gene][group[0]][0]
                group_norm_seq = peptides[gene][group[0]][1]
                # Print to faa files
                print('>', pepname, sep='', file=t_f)
                print(group_tum_seq, file=t_f)

                if group_norm_seq is not None:
                    print('>', pepname, sep='', file=n_f)
                    print(group_norm_seq, file=n_f)
                # Save to map dict
                outmap[pepname] = group_info
                peptide_number += 1
    #  Json dump the results to file
    with open(''.join([out_tum_file, '.map']), 'w') as mapfile:
        json.dump(outmap, mapfile)
    return None


def get_proteome_data(infile):
    """
    Loads GENCODE translated transcript sequences into a dictionary

    :param file infile: GENCODE file object
    :return: Mapping of concatenated Ensembl gene IDs, Ensembl transcript IDs, and HUGO names to protein sequences
    :rtype: dict
    """
    # Read the proteomic fasta
    logging.info('Reading GENCODE proteome fasta file')
    peptides = collections.Counter()
    for fa_seq in read_fasta(infile, 'ARNDCQEGHILKMFPSTWYVBZJUOX'):
        #  Fastq headers are ALWAYS 7 or 8 |-separated fields long. The fields are
        #  0a. Ensembl Peptide (optional) -- e.g ENSP00000334393.3
        #  0b. Ensembl Transcript         -- e.g ENST00000335137.3
        #  1. Ensembl Gene                -- e.g ENSG00000186092.4
        #  2. Havana Gene                 -- e.g OTTHUMG00000001094.2
        #  3. Havana Transcript           -- e.g OTTHUMT00000003223.2
        #  4. HGNC gene splice variant    -- e.g OR4F5-001
        #  5. HUGO name / HGNC symbol     -- e.g OR4F5
        #  6. Length in AA residues       -- e.g 305
        #  We need columns 0b, 1,  and 5
        try:
            record_name = '_'.join(itemgetter(-6, -7, -2)(fa_seq[0].split('|')))
        except IndexError:
            raise RuntimeError('Was the input peptides file obtained from Gencode?')
        peptides[record_name] = list(fa_seq[2])
    return peptides


def main(params):
    """
    Transgene accepts a SnpEffed vcf of translated mutations -- i.e. a file that has information of
    the mutation in genomic space (e.g. chr1:47403668) and in proteomic space (e.g.
    ENST00000371904.4:D113N) -- and outputs files that contain n-mer peptides (n>0) that can be
    passed to MHC:peptide binding prediction tools. Additionally, this tool accepts gene fusion
    annotations in BEDPE format and outputs n-mer peptides across the gene fusion boundary.

    Transgene can be run in an RNA-aware mode where it will only accept small mutations that have
    evidence in an RNA-seq bam (Mapped to the genome reference). It does not look for evidence of
    fusion events in the RNA-seq bam.
    """
    if params.logfile is not None:
        params.logfile = os.path.abspath(os.path.expanduser(params.logfile))
        assert os.path.exists(os.path.dirname(params.logfile)), 'Invalid Directory to write logfile'
        print('Logging to file: %s' % params.logfile)
    logging.basicConfig(filename=params.logfile,
                        level=getattr(logging, params.log_level),
                        format='%(levelname)s: %(message)s')

    # Process the arguments
    if params.rna_file is not None:
        params.rna_file = os.path.abspath(os.path.expanduser(params.rna_file))
        if not os.path.exists(params.rna_file):
            raise RuntimeError('Could not file the RNA-seq bam file at the provided location.')
        if not os.path.exists(params.rna_file + '.bai'):
            raise RuntimeError('Could not file the RNA-seq bai file at the provided location.')

    if params.dna_file is not None:
        params.dna_file = os.path.abspath(os.path.expanduser(params.dna_file))
        if not os.path.exists(params.dna_file):
            raise RuntimeError('Could not file the input tumor DNA-seq bam file at the provided '
                               'location.')
        if not os.path.exists(params.dna_file + '.bai'):
            raise RuntimeError('Could not file the input tumor DNA-seq bai file at the provided '
                               'location.')

    # Load proteomic data
    proteome_data = None
    if params.peptide_file:
        proteome_data = get_proteome_data(params.peptide_file)
        params.peptide_file.close()

    # Read in snpeffed vcf files
    mutations = None
    if params.snpeff_file:
        if params.rna_file or params.dna_file:
            out_vcf = open('_'.join([params.prefix, 'transgened.vcf']), 'w')
        else:
            out_vcf = None
        try:
            mutations = read_snpeffed_vcf(snpeff_file=params.snpeff_file,
                                          rna_file=params.rna_file,
                                          out_vcf=out_vcf,
                                          reject_threshold=params.reject_threshold,
                                          rna_min_alt_freq=params.rna_min_alt_freq,
                                          dna_file=params.dna_file,
                                          oxog_min_alt_freq=params.oxog_min_alt_freq,
                                          processes=params.cores)
        finally:
            if params.rna_file or params.dna_file:
                out_vcf.close()
            params.snpeff_file.close()

    # Load data for generating fusion peptides
    if params.transcript_file:
        transcriptome, gene_transcript_ids = get_transcriptome_data(params.transcript_file)
    else:
        transcriptome, gene_transcript_ids = None, None

    # Load data from fusion file
    if params.fusion_file:
        fusions = read_fusions(params.fusion_file)
    else:
        fusions = None

    if params.genome_file and params.annotation_file:
        exons = get_exons(params.genome_file, params.annotation_file)
    else:
        exons = None

    for peplen in params.pep_lens.split(','):
        logging.info('Processing %s-mers', peplen)
        outfile1, tumfile_path = mkstemp()
        outfile2, normfile_path = mkstemp()
        os.close(outfile1)
        os.close(outfile2)

        with open(tumfile_path, 'w') as tumfile, open(normfile_path, 'w') as normfile:
            if proteome_data and mutations:
                insert_mutations(proteome_data, mutations, tumfile, normfile, int(peplen),
                                 params.rna_file)
            if transcriptome and gene_transcript_ids and fusions:
                insert_fusions(transcriptome, fusions, gene_transcript_ids, int(peplen), tumfile,
                               exons=exons)

        if params.no_json_dumps:
            shutil.move(tumfile_path, '_'.join([params.prefix, 'tumor', peplen,
                                                'mer_snpeffed.faa']))
            shutil.move(normfile_path, '_'.join([params.prefix, 'normal', peplen,
                                                 'mer_snpeffed.faa']))
        else:
            parse_peptides(tumfile_path, normfile_path, params.prefix,
                           '_'.join([peplen, 'mer_snpeffed.faa']))
            os.remove(tumfile_path)
            os.remove(normfile_path)

    # Temporary file used during fusion alignment steps
    if os.path.exists('transgene_fusion_alignments.pkl'):
        os.remove('transgene_fusion_alignments.pkl')


def run_transgene():
    """
    This will try to run transgene from system arguments
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--peptides', dest='peptide_file',
                        type=argparse.FileType('r'),
                        help='Path to GENCODE translation FASTA file')
    parser.add_argument('--transcripts', dest='transcript_file',
                        type=argparse.FileType('r'),
                        help='Path to GENCODE transcript FASTA file')
    parser.add_argument('--snpeff', dest='snpeff_file',
                        type=argparse.FileType('r'),
                        help='Path to snpeff file')
    parser.add_argument('--cores', dest='cores', type=int,
                        help='Number of cores to use for the filtering step.', required=False,
                        default=1)
    parser.add_argument('--fusions', dest='fusion_file',
                        help='Path to gene fusion file',
                        type=argparse.FileType('r'))
    parser.add_argument('--genome', dest='genome_file',
                        help='Path to reference genome file',
                        type=argparse.FileType('r'))
    parser.add_argument('--annotation', dest='annotation_file',
                        help='Path to gencode annotation file',
                        type=argparse.FileType('r'))
    parser.add_argument('--prefix', dest='prefix', type=str,
                        help='Prefix for output file names', required=True)
    parser.add_argument('--pep_lens', dest='pep_lens', type=str,
                        help='Desired peptide lengths to process. '
                             'The argument should be in the form of comma separated values.  '
                             'E.g. 9,15', required=False, default='9,10,15')
    parser.add_argument('--no_json_dumps', action='store_true',
                        help='Do not educe peptide fasta record names in the output by dumping the '
                             'mapping info into a .map json file.', required=False, default=False)
    # RNA-Aware options
    parser.add_argument('--rna_file', dest='rna_file', help='The path to an RNA-seq bam file. If '
                        'provided, the vcf will be filtered for coding mutations only. The file '
                        'must be indexed with samtools index.',
                        required=False, default=None)
    parser.add_argument('--reject_threshold', dest='reject_threshold', help='The minimum number of '
                        'reads containing the REF allele required to reject an event where no '
                        'reads contain the ALT allele.  If the value is lower than this, we will '
                        'accept the mutation.', type=int, required=False, default=5)
    parser.add_argument('--min_rna_alt_freq', dest='rna_min_alt_freq', help='The ALT allele '
                        'frequency (as a fraction) in the RNA-Seq below which we will reject the '
                        'mutation.', type=float, required=False, default=0.1)
    # OxoG filtering options
    parser.add_argument('--filterOxoG', dest='filter_oxog', action='store_true', help='Filter the '
                        'calls for OxoG artifacts. This feature requires a tumor dna bam as input.',
                        required=False, default=False)
    parser.add_argument('--dna_file', dest='dna_file', help='The path to an tumor DNA-seq bam '
                        'file. This is required for OxoG artifact filtering. The file must be '
                        'indexed with samtools index.', required=False, default=None)
    parser.add_argument('--min_OxoG_variant_freq', dest='oxog_min_alt_freq', help='The ALT '
                        'allele frequency (as a fraction) in the DNA-Seq below which we will flag'
                        'the mutation as being an OxoG variant.',
                        type=float, required=False, default=0.1)
    parser.add_argument('--log_level', dest='log_level', help='The level of logging above which '
                        'messages should be printed.', required=False, choices={'DEBUG', 'INFO',
                                                                                'WARNING', 'ERROR'},
                        default='INFO')
    parser.add_argument('--log_file', dest='logfile', help='A path to a logfile.', type=str,
                        required=False, default=None)
    params = parser.parse_args()

    if params.snpeff_file and not params.peptide_file:
        raise ValueError('VCF file requires GENCODE translation FASTA file')

    if params.fusion_file and not params.transcript_file:
        raise ValueError('Fusion file requires GENCODE transcripts FASTA file')

    if params.filter_oxog:
        if not params.dna_file:
            raise ValueError('OxoG artifact filtering required a tumor DNA-Seq bam file.')
        if params.oxog_min_alt_freq > 1.0 or params.oxog_min_alt_freq < 0.0:
            raise ValueError('--min_OxoG_variant_threshold must be a fraction between 0 and 1.')

    if params.rna_file:
        if params.rna_min_alt_freq > 1.0 or params.rna_min_alt_freq < 0.0:
            raise ValueError('--min_rna_alt_threshold must be a fraction between 0 and 1.')

    return main(params)


if __name__ == '__main__':
    sys.exit(run_transgene())
