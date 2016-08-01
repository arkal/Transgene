#!/usr/bin/env python2.7
'''
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : Transgene/transgene


Program info can be found in the docstring of the main function
'''
from __future__ import division, print_function
from operator import itemgetter
from tempfile import mkstemp

import argparse
import collections
import json
import os
import re
import shutil
import sys



def read_fasta(input_file, alphabet):
    '''
    This module will accept an input fasta and yield individual sequences
    one at a time.
    '''
    regexp_string = ''.join(['[^', alphabet, (alphabet).lower(), ']+'])
    nucs_regexp = re.compile(regexp_string)  # regexp to identify sequence lines
    first_seq = True  # used to bypass first seq
    temp_id = None  # variable to hold the id string before processing
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
                #  rstrip to remove trailing newline character
                if len(temp_id) == 3:
                    seq_comments = temp_id[2].rstrip()
                seq = ''
                continue
            if nucs_regexp.findall(seq) != []:
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
        seq = seq + line
    # If seq hasn't been initialized, the input file was empty.
    try:
        seq
    except UnboundLocalError:
        raise RuntimeError('Input peptides file was empty.')
    if nucs_regexp.findall(seq) != []:
        seq = re.sub(nucs_regexp, '', seq)
    yield [seq_id, seq_comments, seq.upper()]


def read_snvs(snpeff_file, chromnames):
    '''
    This module reads in the SNVs from the snpeff_file provided to the program.
    It assumes that the file has a non #-beginning header containing the words
    contig/chrom, pos*, ref*, alt*. In the absence of a header, it will take
    the following as default
    [1] - chrom/contig
    [2] - position (1-based)
    [4] - ref allele
    [5] - alt allele
    '''
    snvs = collections.Counter()
    for line in snpeff_file:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        changes = line[7].split(';')[-1]
        if changes.startswith('EFF'):
            changes = re.sub('EFF=', '', changes)
            changes = [x for x in changes.split(',') if
                       (x.startswith('NON_SYNONYMOUS_CODING') or
                        x.startswith('STOP_GAINED')) and 'protein_coding' in x]
            for i in range(0, len(changes)):
                temp = changes[i].split('|')
                if snvs[temp[8]] == 0:
                    snvs[temp[8]] = collections.Counter({temp[3]:(temp[3][0],
                                                                  temp[3][-1])})
                else:
                    snvs[temp[8]][temp[3]] = (temp[3][0], temp[3][-1])
        else:
            pass
    if len(snvs) == 0:
        raise RuntimeError('Input snpeffed mutations file was empty or had no' +
                           ' actionable mutations.')
    return snvs


def insert_snvs(chroms, snvs, outfile, peplen):
    '''
    This module uses the snv data contained in snvs and inserts them into the
    genome contained in chroms. The arguments are
    chroms  (DICT) - contains the peptides in the form of dicts where keys
                     hold the protein name and the values holds the sequence.
    snvs    (DICT) - contains the snvs parsed from teh input SnpEFF file.
    outfile (FILE) - the file to write output to.
    peplen  (INT)  - length of peptides which will be generated from the ouput
                     file.
    '''
    for pept in snvs.keys():
        #  First, grab the positions of all mutations in the peptide.  If there
        #  are multiple mutations, then for each mutation, list the other
        #  mutations that are within peplen positions of it IN THE POSITIVE
        #  DIRECTION.  The mutations are potentially coexpressed on the same
        #  peptide. The reason can be visually shown for a 9-mer peptive with
        #  three mutations, X, Y and Z as below (O = non-mutated AA).  The
        #  region depicts the extended IAR for the 3. Possible mutational
        #  combinations are shown. The WT (O,O,O) isn't shown since it isn't
        #  useful.
        #
        #  POSITONS:    a b c d e f g h i j k l m n o p q r s t u v w x y
        #
        #  X+Y+Z        O O O O O O O O X O O O Y O O O Z O O O O O O O O
        #  X + Z        O O O O O O O O X O O O O O O O Z O O O O O O O O
        #  X + Y        O O O O O O O O X O O O Y O O O O O O O O O O O O
        #  X only       O O O O O O O O X O O O O O O O O O O O O O O O O
        #
        #  Y + Z        O O O O O O O O O O O O Y O O O Z O O O O O O O O
        #  Y only       O O O O O O O O O O O O Y O O O O O O O O O O O O
        #
        #  Z only       O O O O O O O O O O O O O O O O Z O O O O O O O O
        #
        #  It is clear that if we only look forward, all combinations of mutated
        #  and non mutated peptides can be seen.  For example, the peptide from
        #  positions g through o. The three combinations are X and Y, X and O,
        #  and O and Y.  The first 2 will be met through X and forward, and the
        #  3rd will be met with Y and forward.
        #
        #  The data is stored in the mutation_groups list object.  Each group in
        #  in the list is a list of all mutations in the group.  If there are no
        #  neighboring mutations, the group contains only the mutation itself.
        mutation_groups = {x: [y for y in snvs[pept].keys() if
                               int(y[1:-1]) < int(x[1:-1]) + peplen and
                               int(y[1:-1]) > int(x[1:-1])] for x in
                           snvs[pept].keys()}
        for mutation, mut_group in mutation_groups.items():
            mut_pos = int(mutation[1:-1])
            pfasta_name = [x for x in chroms.keys() if
                           ''.join([pept, '_']) in x][0]
            #  The keys for snvs[pept] are mutations in the form of
            #  <REF><POS><ALT> (Eg A521C, V98F).  If the reference
            #  protein stored in chroms has the same ref allele as
            #  described in the mutation, then continue - ELSE print
            #  an error and skip.
            if chroms[pfasta_name][mut_pos-1] == snvs[pept][mutation][0]:
                #  If the firs mutation is a stop gain, it does not yield any
                #  neo epitopes
                if snvs[pept][mutation][1] == '*':
                    continue
                #  Sort grouped mutations by position.
                mut_group = sorted(mut_group, key=lambda x: int(x[1:-1]))
                num_iars = 2**len(mut_group)
                protein = chroms[pfasta_name]
                #  n total mutations in a peplen sequence will yield 2^(n-1)
                #  combinations.  These will be stored in a dict of dicts object
                #  as
                #  out_pepts
                #       |- XXXXXXXXXX
                #       |  |-pept_name: []
                #       |  +-pept_seq: []
                #       |  ..
                #       |  ..
                #       |- n
                #       |  |-pept_name: []
                #       |  +-pept_seq: []
                #  The pept_name for each seq will be updated as new mutations
                #  are added, as will the sequences.  The dict will grow as in
                # a binary tree-like fashion
                #  If a stop is seen, the peptide is printed and the value is
                #  popped from the dict.
                #  The keys of the outer dict serve no real purpose
                first = max(mut_pos-peplen, 0)  # First possible AA in the IAR
                #  The dict of dicts is initialized with the sequence before and
                #  upto the first mutation in the IAR.  This is ok because all
                #  sequence in all outputs are identical.
                out_pepts = {
                    os.urandom(10):{'pept_name': [pfasta_name, mutation],
                                    'pept_seq': protein[first:mut_pos-1] +
                                        [snvs[pept][mutation][1]]}}
                last = mut_pos  #  The last AA that was procesed
                for mut in mut_group:
                    for iar in out_pepts.keys():
                        child_iars = grow_out_pept(out_pepts[iar], last,
                                                   mut, protein, outfile,
                                                   peplen)
                        out_pepts.pop(iar)
                        out_pepts.update(child_iars)
                        last = int(mut[1:-1])
                for iar in out_pepts.values():
                    iar['pept_seq'] += protein[last:min(mut_pos + peplen - 1,
                                                       len(protein))]
                write_pepts_to_file(out_pepts, outfile, peplen)
            else:
                print('ERROR: ', chroms[pfasta_name][mut_pos-1], ' seen at ',
                      'position ', mut_pos, ' in ', pept, '.  ',
                      snvs[pept][mutation][0], ' expected.', sep='',
                      file=sys.stderr)
    return None


def grow_out_pept(curr_pept, last, mutation, protein, outfile, peplen):
    '''
    Accept a growing IAR and add in infor for the region spanning until the next
    mutation.
    curr_pept is a dict as
        |- pept_name: list of identifiers for the name
        +- pept_seq: list of characters
    last is an int
    mutation is a string like A980V
    protein is a list of characters for the protein seq
    '''
    mut_pos = int(mutation[1:-1])
    mut_pept = curr_pept['pept_seq'] + protein[last:mut_pos-1]
    unmut_pept = curr_pept['pept_seq'] + protein[last:mut_pos-1]
    unmut_pept.append(mutation[0])
    if mutation[-1] == '*':
        #  Write the mutated STOP containing peptide to the output
        write_pepts_to_file({os.urandom(10):{
                                'pept_name':curr_pept['pept_name']+[mutation],
                                'pept_seq': mut_pept}}, outfile, peplen)
        #  Return only the non mutated one
        return {os.urandom(10):{'pept_name':curr_pept['pept_name'],
                                'pept_seq': unmut_pept}}
    elif protein[mut_pos-1] != mutation[0]:
        print('ERROR :', chroms[pept][int(pos)-1], 'seen at position',
                      pos, 'in', pept, '.', snvs[pept][mutation][0],
                      'expected.', sep=' ', file=sys.stderr)
        #  Return only the non mutated one
        return {os.urandom(10):{'pept_name':curr_pept['pept_name'],
                                'pept_seq': unmut_pept}}
    else:
        mut_pept.append(mutation[-1])
        return {os.urandom(10):{'pept_name':curr_pept['pept_name']+[mutation],
                                'pept_seq': mut_pept},
                os.urandom(10):{'pept_name':curr_pept['pept_name'],
                                'pept_seq': unmut_pept}}


def write_pepts_to_file(pept_dict, outfile, peplen):
    '''
    This module will accept a dict of dicts object as
    pept_dict
        out_pepts
            |- XXXXXXXXXX
            |  |-pept_name: []
            |   +-pept_seq: []
            |  ..
            |  ..
            |- n
            |  |-pept_name: []
            |  +-pept_seq: []
    and will print then to outfile
    '''
    blacklist = set(list('BJOUXZ'))
    for pept_info in pept_dict.values():
        peptide_name = '_'.join(pept_info['pept_name'])
        peptide_sequence = ''.join(pept_info['pept_seq'])
        # If the final peptide contains a blacklisted amino acid,
        # discard it.
        if not blacklist.isdisjoint(peptide_sequence):
            print('INFO: Blacklisted amino acids seen in ', peptide_name, '.',
                  sep='', file=sys.stderr)
            continue
        if len(peptide_sequence) < peplen:
            print('INFO: IAR less than %s' % peplen + ' resides seen at\n',
                  peptide_name, '\n', peptide_sequence, sep='', file=sys.stderr)
            continue
        print('>', peptide_name, sep='', file=outfile)
        print(peptide_sequence, file=outfile)
    return None

def parse_peptides(infile, outfile):
    '''
    This module takes in a peptides file and squashes it into the miminum number
    of peptides required to describe the potential neo-immunopeptidome.  It's
    main function is to take transcript-level mutation calls and merge them by
    gene if the sequences they contain are identical (to correct for splicing
    information).
    '''
    with open(infile, 'r') as i_f:
        peptides = collections.Counter()
        for pep_name, _, pep_seq in read_fasta(i_f,
                                               'ARNDCQEGHILKMFPSTWYVBZJUOX'):
            pep_name =  pep_name.split('_')
            gene_name = pep_name[0]
            hugo_gene = pep_name[2]
            transcript_mutation = '_'.join([pep_name[1]]+pep_name[3:])
            if peptides[(gene_name, hugo_gene)] == 0:
                peptides[(gene_name, hugo_gene)] = {transcript_mutation:
                                                    pep_seq}
            else:
                peptides[(gene_name, hugo_gene)][transcript_mutation] = pep_seq
    outmap = {}
    with open(outfile, 'w') as o_f:
        #  mhc predictors can't handle the giant peptide names created by
        #  merging the names of peptide emerging from transcripts bearing
        #  similar mutations hence we use an incremental system to naem the
        #  peptides and store the peptide number to actual peptide name
        #  information in a map file.
        peptide_number=1
        for gene in peptides.keys():
            unique_seqs = set(peptides[gene].values())
            for group in [[x for x, y in peptides[gene].items() if y == z] for
                          z in unique_seqs]:
                pepname = ''.join(['neoepitope_', str(peptide_number)])
                group_info = '\t'.join(list(gene) +[','.join(group)])
                groupseq = peptides[gene][group[0]]
                #  Print to faa files
                print('>', pepname, sep='', file=o_f)
                print(groupseq, file=o_f)
                #  Save to map dict
                outmap[pepname] = group_info
                peptide_number+=1
    #  Json dump the results to file
    with open(''.join([outfile, '.map']), 'w') as mapfile:
        json.dump(outmap, mapfile)
    return None


def main(params):
    '''
    This tool is used to take a vcf of translated mutations -- i.e. a file that
    has information of the mutation in genomic space (e.g. chr1:47403668) and in
    proteomic space (e.g. ENST00000371904.4:D113N) -- and obtain files that
    contain 9-, 10- and 15-mer peptides that can be passed to MHC:peptide
    binding prediction tools.
    '''

    # Read the proteomic fasta
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
            raise RuntimeError('Was the input peptides file obtained from ' +
                               'Gencode?')
        chroms[record_name] = list(fa_seq[2])

    # Read in snpeff file
    # The naming convention of the variables and functions comes from a previous
    # functionality
    snvs = read_snvs(params.snpeff_file, chroms.keys())
    params.input_file.close()
    params.snpeff_file.close()
    for peplen in params.pep_lens.split(','):
        outfile, outfile_path = mkstemp()
        os.close(outfile)
        with open(outfile_path, 'w') as outfile:
            insert_snvs(chroms, snvs, outfile, int(peplen))
        if params.no_json_dumps:
            shutil.move(outfile_path, '_'.join([params.prefix, 'tumor',
                                                peplen, 'mer_snpeffed.faa']))
        else:
            parse_peptides(outfile_path, '_'.join([params.prefix, 'tumor',
                                                   peplen, 'mer_snpeffed.faa']))
            os.remove(outfile_path)

if __name__ == '__main__':
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
    params = parser.parse_args()
    sys.exit(main(params))
