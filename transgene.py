#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : Transgene/transgene


Program info can be found in the docstring of the main function
"""
from __future__ import division, print_function

import random
from operator import itemgetter
from tempfile import mkstemp

import argparse
import collections
import json
import logging
import os
import pysam
import re
import shutil
import sys


def read_fasta(input_file, alphabet):
    """
    This module will accept an input fasta and yield individual sequences
    one at a time.
    """
    regexp_string = ''.join(['[^', alphabet, alphabet.lower(), ']+'])
    nucs_regexp = re.compile(regexp_string)  # Regexp to identify sequence lines
    first_seq = True  # Used to bypass first seq
    seq, seq_id, seq_comments = None, None, None
    for line in input_file:
        line = line.lstrip()
        if len(line) == 0:  # blank line
            continue
        if line[0] == '>':  # id line
            if first_seq:
                first_seq = False
                #  >seq_id comments becomes ['', 'seq_id','comments']
                temp_id = re.split(r'[>\s,]+', line, maxsplit=2)
                seq_id = temp_id[1]
                #  right strip to remove trailing newline character
                if len(temp_id) == 3:
                    seq_comments = temp_id[2].rstrip()
                seq = ''
                continue
            if nucs_regexp.findall(seq):
                seq = re.sub(nucs_regexp, '', seq)
            yield [seq_id, seq_comments, seq.upper()]
            #  >seq_id comments becomes ['', 'seq_id','comments']
            temp_id = re.split(r'[>\s,]+', line, maxsplit=2)
            seq_id = temp_id[1]
            if len(temp_id) == 3:
                seq_comments = temp_id[2].rstrip()
            seq = ''
            continue
        #  Remove whitespaces in sequence string
        line = ''.join(line.split())
        seq += line
    # If seq hasn't been initialized, the input file was empty.
    try:
        seq
    except UnboundLocalError:
        raise RuntimeError('Input peptides file was empty.')
    if nucs_regexp.findall(seq):
        seq = re.sub(nucs_regexp, '', seq)
    yield [seq_id, seq_comments, seq.upper()]


# A named tuple describing the result of running reject_mutation on a mutation.
reject_decision = collections.namedtuple('reject_decision', (
    # The decision on whether to reject the mutation or not
    'reject',
    # The reason for rejection, if decision == True
    'reason',
    # The number of reads at the position
    'coverage'))


def reject_mutation(snv, rna_bam, reject_threshold):
    """
    Decide whether the mutation should be rejected based on the expression filter.

    :param dict snv: A single vcf record split by the tab character and keyed by the fields
    :param str rna_bam: Path to the RNA-seq bam file
    :param int reject_threshold: The read coverage above which we reject a no-ALT-coverage event
    :return: A named tuple of the True/False if the mutation should be rejected or not, and a reason
    for rejection
    :rtype: reject_decision
    """
    # NOTE: vcf snvs are 1-based and not 0-based!
    snv['POS'] = int(snv['POS']) - 1
    logging.info('Processing mutation at position %s:%s', snv['CHROM'], snv['POS'] + 1)
    samfile = pysam.Samfile(rna_bam, 'rb')
    alignment = samfile.fetch(snv['CHROM'], snv['POS'], snv['POS'] + 1)
    output_counts = collections.Counter()
    total_reads = 0
    spanning_reads = 0
    for read in alignment:
        read_pos = [x for x, y in read.get_aligned_pairs() if y == snv['POS'] and x is not None]
        if read_pos:
            read_pos = read_pos[0]
            # This means read_pos is not an empty dict and it should have EXACTLY 1 entry
            if read.seq[read_pos] is not None and read.qual[read_pos] > 30:
                output_counts[read.seq[read_pos]] += 1
            total_reads += 1
        spanning_reads += 1
    if spanning_reads == 0:
        logging.warning('Mutation at position %s:%s has no coverage in the RNA-seq file. '
                        'Rejecting.', snv['CHROM'], snv['POS'] + 1)
        reject = reject_decision(reject=True, reason='NoReadCoverage', coverage='0/0')
    elif not output_counts:
        logging.warning('Mutation at position %s:%s appears to be in a region that has no coverage,'
                        ' but has %s spanning reads. Possible Splicing event. Rejecting.',
                        snv['CHROM'], snv['POS'] + 1, spanning_reads)
        reject = reject_decision(reject=True, reason='SpliceDetected',
                                 coverage='0/' + str(spanning_reads))
    elif output_counts[snv['ALT']] == 0:
        if output_counts[snv['REF']] < 5:
            logging.debug('Mutation at position %s:%s has no evidence of existence in the RNA seq. '
                          'However the coverage at the region (%s/%s reads::Covering/spanning) '
                          'is below the threshold for rejecting (%s). Accepting.',
                          snv['CHROM'], snv['POS'] + 1, total_reads, spanning_reads,
                          reject_threshold)
            reject = reject_decision(reject=False, reason='LowReadCoverage',
                                     coverage='/'.join([str(total_reads), str(spanning_reads)]))
        else:
            logging.warning('Mutation at position %s:%s has no evidence of existence in the RNA '
                            'seq. Coverage = %s/%s reads (Covering/spanning). Rejecting.',
                            snv['CHROM'], snv['POS'] + 1, total_reads, spanning_reads)
            reject = reject_decision(reject=True, reason='NoAltDetected',
                                     coverage='/'.join([str(total_reads), str(spanning_reads)]))
    else:
        logging.debug('Accepted mutation at position %s:%s with %s read coverage.', snv['CHROM'],
                      snv['POS'] + 1, total_reads)
        reject = reject_decision(reject=False, reason='',
                                 coverage='/'.join([str(total_reads), str(spanning_reads)]))
    samfile.close()
    return reject


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


def read_snvs(snpeff_file, rna_file=None, out_vcf=None, reject_threshold=None):
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
    :returns: A parsed dictionary of snvs to be processed.
    :rtype: collections.Counter
    """
    if out_vcf:
        assert out_vcf.mode == 'w', 'Output vcf file is not open for writing'
    logging.info('Reading in SNVs')
    snvs = collections.Counter()
    indexes = {'CHROM': 0,
               'POS': 1,
               'REF': 3,
               'ALT': 4,
               'INFO': 7}
    for line in snpeff_file:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
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
                if rna_file:
                    print('##INFO=<ID=reason,Number=.,Type=String,Description=Reason for rejecting/'
                          'accepting a mutation.', file=out_vcf)
                    print('##INFO=<ID=coverage,Number=1,Type=String,Description=Reads supporting '
                          'the decision (Covering/Spanning).', file=out_vcf)
                    write_to_vcf(indexes, out_vcf, header=True)
            elif rna_file:
                print(line, end='', file=out_vcf)
            continue
        line = line.strip().split('\t')
        line = {x: line[y] for x, y in indexes.items()}
        changes = line['INFO'].split(';')[-1]
        if changes.startswith('EFF'):
            if rna_file is not None:
                decision = reject_mutation(line, rna_file, reject_threshold)
                # Print the line to the output vcf
                line['INFO'] = ';'.join([line['INFO'], 'reason=' + decision.reason,
                                         'coverage=' + str(decision.coverage)])
                line['FILTER'] = 'REJECT' if decision.reject else 'PASS'
                write_to_vcf(line, out_vcf)
                # Act on the decision
                if decision.reject:
                    continue
            changes = re.sub('EFF=', '', changes)
            changes = [x for x in changes.split(',') if
                       (x.startswith('NON_SYNONYMOUS_CODING') or
                        x.startswith('STOP_GAINED')) and 'protein_coding' in x]
            for i in range(0, len(changes)):
                temp = changes[i].split('|')
                if snvs[temp[8]] == 0:
                    snvs[temp[8]] = collections.Counter({
                        temp[3]: {
                            'AA': {'REF': temp[3][0], 'ALT': temp[3][-1]},
                            'NUC': {'REF': line['REF'], 'ALT': line['ALT'], 'POS': line['POS'],
                                    'CHROM': line['CHROM']}
                        }
                    })
                else:
                    snvs[temp[8]][temp[3]] = {
                            'AA': {'REF': temp[3][0], 'ALT': temp[3][-1]},
                            'NUC': {'REF': line['REF'], 'ALT': line['ALT'], 'POS': line['POS'],
                                    'CHROM': line['CHROM']}
                    }
        else:
            pass
    if len(snvs) == 0:
        raise RuntimeError('Input snpeffed mutations file was empty or had no actionable '
                           'mutations.')
    return snvs


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
                                       int(y[1:-1]) >
                                       int(group_start[1:-1])])
        for i in groups[-1][1:]:
            mutations.remove(i)
        mutations.remove(group_start)
        if not mutations:
            break
    mutation_groups = [build_groups(x) for x in groups]
    if rna_bam is None:
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


def insert_snvs(chroms, snvs, tumfile, normfile, peplen, rna_bam=None):
    """
    This module uses the snv data contained in snvs and inserts them into the genome contained in
    chroms.
    :param dict chroms: Contains the peptides in the form of dicts where keys hold the protein name
    and the values holds the sequence.
    :param dict snvs: Contains the snvs parsed from teh input SnpEFF file.
    :param file outfile: The file to write output to.
    :param int peplen: Length of peptides which will be generated from the output file.
    :param str rna_bam: The path to the RNA-Seq bam file
    """
    logging.info('Inserting SNVs into IARs')
    for pept in snvs.keys():
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
        pfasta_name = [x for x in chroms.keys() if ''.join([pept, '_']) in x][0]
        all_mutation_groups = get_mutation_groups(snvs[pept], peplen, pfasta_name, rna_bam)
        protein = chroms[pfasta_name]
        for mut_group in all_mutation_groups:
            for group_muts in mut_group:
                # If the first mutation in the group is a stop gain, it does not yield any
                # neo-epitopes.  This also handles the possibility of the ONLY mutation being a
                # stop gain
                if group_muts[0][-1] == '*':
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
                    if protein[mut_pos-1] != snvs[pept][group_mut]['AA']['REF']:
                        logging.warning('%s seen at position %s in %s.... %s expected.',
                                        chroms[pfasta_name][mut_pos-1], mut_pos, pept,
                                        snvs[pept][group_mut]['AA']['REF'])
                    if snvs[pept][group_mut]['AA']['ALT'] == '*':
                        prev = None
                        break
                    out_pept['tum_seq'].append(snvs[pept][group_mut]['AA']['ALT'])
                    out_pept['norm_seq'].append(snvs[pept][group_mut]['AA']['REF'])
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
    # If the final peptide contains a blacklisted amino acid,
    # discard it.
    peptide_name = '_'.join(peptide_info['pept_name'])
    tum_peptide_seq = ''.join(peptide_info['tum_seq'])
    norm_peptide_seq = ''.join(peptide_info['norm_seq'])
    if not blacklist.isdisjoint(tum_peptide_seq):
        logging.warning('Blacklisted amino acids %s seen in %s.',
                        list(set(list(tum_peptide_seq)).intersection(blacklist)), peptide_name)
    elif len(tum_peptide_seq) < peplen:
        logging.warning('IAR %s is less than %s resides (%s).', peplen, peptide_name,
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
                peptides[(gene_name, hugo_gene)] = {transcript_mutation:
                                                        (pep_seq, normal_peptides[full_pep_name])}
            else:
                peptides[(gene_name, hugo_gene)][transcript_mutation] = \
                    (pep_seq, normal_peptides[full_pep_name])
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
                print('>', pepname, sep='', file=n_f)
                print(group_tum_seq, file=t_f)
                print(group_norm_seq, file=n_f)
                # Save to map dict
                outmap[pepname] = group_info
                peptide_number += 1
    #  Json dump the results to file
    with open(''.join([out_tum_file, '.map']), 'w') as mapfile:
        json.dump(outmap, mapfile)
    return None


def main(params):
    """
    This tool accepts a vcf of translated mutations -- i.e. a file that has information of the
    mutation in genomic space (e.g. chr1:47403668) and in proteomic space (e.g.
    ENST00000371904.4:D113N) -- and outputs files that contain n-mer peptides (n>0) that can be
    passed to MHC:peptide binding prediction tools.

    Transgene can be run in an RNA-aware mode where it will only accept mutations that have evidence
    in an RNA-seq bam (Mapped to the genome reference).
    """
    logging.basicConfig(level=getattr(logging, params.log_level), format='%(levelname)s: '
                                                                         '%(message)s')

    # Process the arguments
    if params.rna_file is not None:
        params.rna_file = os.path.abspath(os.path.expanduser(params.rna_file))
        if not os.path.exists(params.rna_file):
            raise RuntimeError('Could not file the RNA-seq bam file at the provided location.')
        if not os.path.exists(params.rna_file + '.bai'):
            raise RuntimeError('Could not file the RNA-seq bai file at the provided location.')

    # Read the proteomic fasta
    logging.info('Reading the Input fasta file')
    chroms = collections.Counter()
    for fa_seq in read_fasta(params.input_file, 'ARNDCQEGHILKMFPSTWYVBZJUOX'):
        #  Fastq headers are ALWAYS 7 |-separated fields long. The fields are
        #  0. Ensembl Transcript         -- e.g ENST00000511116.1
        #  1. Ensembl Gene               -- e.g ENSG00000113658.12
        #  2. Havana Gene                -- e.g OTTHUMG00000163212.3
        #  3. Havana Transcript          -- e.g OTTHUMT00000372098.1
        #  4. HGNC gene splice variant   -- e.g SMAD5-003
        #  5. HUGO name / HGNC symbol    -- e.g SMAD5
        #  6. Length in AA residues      -- e.g 134
        #  We need columns 1, 2,  and 6
        try:
            record_name = '_'.join(itemgetter(1, 0, 5)(fa_seq[0].split('|')))
        except IndexError:
            raise RuntimeError('Was the input peptides file obtained from Gencode?')
        chroms[record_name] = list(fa_seq[2])

    # Read in snpeff file
    # The naming convention of the variables and functions comes from a previous functionality
    out_vcf = None
    try:
        if params.rna_file:
            out_vcf = open('_'.join([params.prefix, 'transgened.vcf']), 'w')
        snvs = read_snvs(params.snpeff_file, params.rna_file, out_vcf, params.reject_threshold)
    finally:
        if out_vcf is not None:
            out_vcf.close()
    params.input_file.close()
    params.snpeff_file.close()
    for peplen in params.pep_lens.split(','):
        logging.info('Processing %s-mers', peplen)
        outfile1, tumfile_path = mkstemp()
        outfile2, normfile_path = mkstemp()
        os.close(outfile1)
        os.close(outfile2)

        with open(tumfile_path, 'w') as tumfile, open(normfile_path, 'w') as normfile:
            insert_snvs(chroms, snvs, tumfile, normfile, int(peplen), params.rna_file)
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


def run_transgene():
    """
    This will try to run transgene from system arguments
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--peptides', dest='input_file',
                        type=argparse.FileType('r'), help='Input peptide' +
                                                          ' FASTA file', required=True)
    parser.add_argument('--snpeff', dest='snpeff_file',
                        type=argparse.FileType('r'),
                        help='Input snpeff file name', required=True)
    parser.add_argument('--prefix', dest='prefix', type=str,
                        help='Output FASTQ file prefix', required=True)
    parser.add_argument('--pep_lens', dest='pep_lens', type=str,
                        help='Desired peptide lengths to process.  The ' +
                             'argument should be in the form of comma separated ' +
                             'values.  E.g. 9,15', required=False, default='9,10,15')
    parser.add_argument('--no_json_dumps', action='store_true',
                        help='Do not educe peptide fasta record names in the ' +
                             'output by dumping the mapping info into a .map json ' +
                             'file.', required=False, default=False)
    parser.add_argument('--rna_file', dest='rna_file', help='The path to an RNA-seq bam file. If '
                        'provided, the vcf will be filtered for coding mutations only. The file '
                        'must be indexed with samtools index.',
                        required=False, default=None)
    parser.add_argument('--reject_threshold', dest='reject_threshold', help='The minimum number of '
                        'reads containing the REF allele required to reject an event where no '
                        'reads contain the ALT allele.  If the value is lower than this, we will '
                        'accept the mutation.', type=int, required=False, default=5)
    parser.add_argument('--log_level', dest='log_level', help='The level of logging above which '
                        'messages should be printed.', required=False, choices={'DEBUG', 'INFO',
                                                                                'WARNING', 'ERROR'},
                        default='INFO')
    params = parser.parse_args()
    return main(params)

if __name__ == '__main__':
    sys.exit(run_transgene())
