import collections
import gzip
import os
import re


# A translation table use to complement string sequences
import string

forward = 'ACGTN'
reverse = 'TGCAN'
trans = string.maketrans(forward, reverse)


def read_fasta(input_file, alphabet):
    """
    This module will accept an input fasta and yield individual sequences
    one at a time.

    :param file input_file: FASTA file object
    :param str alphabet: Allowed characters in FASTA sequence
    :return: FASTA ID, Comment, and Sequence
    :rtype: generator
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
    if seq:
        if nucs_regexp.findall(seq):
            seq = re.sub(nucs_regexp, '', seq)
        yield [seq_id, seq_comments, seq.upper()]


def chrom_sort(in_chroms):
    """
    Sort a list of chromosomes in the order 1..22, X, Y, M.

    :param list in_chroms: Input chromsomes
    :return: Sorted chromosomes
    :rtype: list[str]
    """
    if not in_chroms:
        return in_chroms
    chr_prefix = False
    if in_chroms[0].startswith('chr'):
        in_chroms = [x.lstrip('chr') for x in in_chroms]
        chr_prefix = True
    assert in_chroms[0] in [str(x) for x in range(1, 23)] + ['X', 'Y', 'M']
    in_chroms = sorted(in_chroms, key=lambda c: int(c) if c not in ('X', 'Y', 'M') else c)
    try:
        m_index = in_chroms.index('M')
    except ValueError:
        pass
    else:
        in_chroms.pop(m_index)
        in_chroms.append('M')
    # At this point it should be nicely sorted
    if chr_prefix:
        in_chroms = [''.join(['chr', x]) for x in in_chroms]
    return in_chroms


class GTFRecord(object):
    def __init__(self, line):
        """
        Converts a GTF record into a python object

        :param line str|list line: GTF record
        """
        if isinstance(line, str):
            line = line.strip().split('\t')
        if len(line) != 9:
            msg = '\n'.join(self.attributes)
            raise ValueError('Malformed GTF line.'
                             ' Must have the following attributes: {}'.format(msg))
        for (attr, func), value in zip(self.attributes, line):
            setattr(self, attr, func(value))
        for feature in self.attribute.split(';'):
            try:
                attr, value = feature.split()
                if attr in self.__dict__:
                    if isinstance(getattr(self, attr), list):
                        setattr(self, attr, getattr(self, attr) + [value.replace('"', '')])
                    else:
                        setattr(self, attr, [getattr(self, attr), value.replace('"', '')])
                else:
                    setattr(self, attr, value.replace('"', ''))
            except ValueError:
                pass

    attributes = [('seqname', str), ('source', str), ('feature', str), ('start', int),
                  ('end', int), ('score', str), ('strand', str), ('frame', str), ('attribute', str)]

    def __repr__(self):
        if self.feature == 'gene':
            return 'GTF({}, gene:{}, start:{}, length:{})'.format(self.seqname,
                                                                  self.gene_name,
                                                                  self.start,
                                                                  self.end - self.start + 1)
        elif self.feature == 'start_codon':
            return 'GTF({}, start codon:{}, start:{}, length:{})'.format(self.seqname,
                                                                         self.gene_name,
                                                                         self.start,
                                                                         self.end - self.start + 1)
        elif self.feature == 'transcript':
            return 'GTF({}, transcript:{}, start:{}, length:{})'.format(self.seqname,
                                                                        self.transcript_name,
                                                                        self.start,
                                                                        self.end - self.start + 1)
        elif self.feature == 'exon':
            return 'GTF({}, exon:{}, start:{}, length:{}, id:{})'.format(self.seqname,
                                                                         self.exon_number,
                                                                         self.start,
                                                                         self.end - self.start + 1,
                                                                         self.exon_id)
        elif self.feature == 'CDS':
            return 'GTF({}, CDS:{}, start:{}, length:{}, id:{})'.format(self.seqname,
                                                                         self.exon_number,
                                                                         self.start,
                                                                         self.end - self.start + 1,
                                                                         self.exon_id)

    def __hash__(self):
        if self.feature == 'gene':
            return hash(self.gene_id)
        elif self.feature == 'transcript':
            return hash(self.transcript_id)
        elif self.feature == 'CDS':
            return hash('CDS' + self.exon_id)
        elif self.feature == 'exon':
            return hash(self.exon_id)

    def __eq__(self, other):
        try:
            if self.feature == 'gene':
                return self.gene_id == other.gene_id
            elif self.feature == 'transcript':
                return self.transcript_id == other.transcript_id
            elif self.feature == 'CDS':
                if other.feature == 'CDS':
                    return self.exon_id == other.exon_id
                else:
                    return False
            elif self.feature == 'exon':
                return self.exon_id == other.exon_id
            else:
                return False
        except AttributeError:
            return False

# A named tuple describing the result of running reject_mutation on a mutation.
reject_decision = collections.namedtuple('reject_decision', (
    # The decision on whether to reject the mutation or not
    'reject',
    # The reason for rejection, if decision == True
    'reason',
    # The number of reads at the position
    'coverage',
    # The variant allele frequency
    'vaf'))

# Standard Genetic Code from NCBI
amino = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
base1 = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
base2 = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
base3 = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
genetic_code = {''.join([b1, b2, b3]): aa
                for aa, b1, b2, b3 in zip(amino, base1, base2, base3)}

EFF = collections.namedtuple('EFF',
                             'annotation, '  # The annotation of the event
                             'codon_change, '  # A tuple of the changed nucleotide in the codon
                             'aa_change, '
                             'gene_name, '
                             'transcript_id, '
                             'transcript_len, '
                             'warnings'
                             )

BEDPE = collections.namedtuple('BEDPE',
                               'chrom1, start1, end1, '
                               'chrom2, start2, end2, '
                               'name, '
                               'score, '
                               'strand1, strand2, '
                               'junctionSeq1, junctionSeq2, '
                               'hugo1, hugo2')

# A dictionary to translate used legacy snpeff annotations to sequence ontology
snpeff_to_so = {'NON_SYNONYMOUS_CODING': 'missense_variant',
                'STOP_GAINED': 'stop_gained',
                'SYNONYMOUS_CODING': 'synonymous_variant',
                'FRAME_SHIFT': 'frameshift_variant',
                'CODON_CHANGE_PLUS_CODON_DELETION': 'disruptive_inframe_deletion',
                'CODON_CHANGE_PLUS_CODON_INSERTION': 'disruptive_inframe_insertion',
                'CODON_DELETION': 'inframe_deletion',
                'CODON_INSERTION': 'inframe_insertion',
                }


amino_letter = {'Ala': 'A',
                'Arg': 'R',
                'Asn': 'N',
                'Asp': 'D',
                'Cys': 'C',
                'Glu': 'E',
                'Gln': 'Q',
                'Gly': 'G',
                'His': 'H',
                'Ile': 'I',
                'Leu': 'L',
                'Lys': 'K',
                'Met': 'M',
                'Phe': 'F',
                'Pro': 'P',
                'Ser': 'S',
                'Thr': 'T',
                'Trp': 'W',
                'Tyr': 'Y',
                'Val': 'V',
                'Ter': '*',
                '*': '*'
                }


def get_exons(genome_file, annotation_file, genes_of_interest, drop_transcript_version=False):
    """
    Generates list of GTFRecord objects for each transcript and the position of the start codon
    in a separate dict

    :param file genome_file: Reference genome FASTA file
    :param file annotation_file: Genome annotation file (GTF)
    :param set(str) genes_of_interest: The genes in this sample that might need translation.
    :param bool drop_transcript_version: Drop the version part of transcript ids.
    :return: list GTFRecord of exons for each transcript, metadata for each transcript
    :rtype: dict(str, list(GTFRecord)), dict(str, int|bool)
    """
    annotation_file.seek(0)
    chroms = {}
    exons = collections.defaultdict(list)
    cds_starts = {}
    for header, comment, seq in read_fasta(genome_file, 'ACGTN'):
        chroms[header] = seq

    for line in annotation_file:
        if line.startswith('#'):
            continue
        else:
            gtf = GTFRecord(line)
            if gtf.feature == 'exon' and gtf.gene_name in genes_of_interest:
                transcript_id = gtf.transcript_id
                if drop_transcript_version:
                    transcript_id = transcript_id.split('.')[0]
                gtf.sequence = chroms[gtf.seqname][gtf.start - 1: gtf.end]
                exons[transcript_id].append(gtf)
            elif gtf.feature == 'start_codon' and gtf.gene_name in genes_of_interest:
                transcript_id = gtf.transcript_id
                if drop_transcript_version:
                    transcript_id = transcript_id.split('.')[0]
                cds_starts[transcript_id] = gtf.start
    return exons, cds_starts


def read_genes_from_gtf(gtf_file):
    """
    Read the gene annotations into a dict

    :param file gtf_file: A file handle ot the annotation file.
    :returns:  A dict of a gtf record for each gene
    :rtype: dict(string, GTFRecord)
    """
    gene_annotations = {}
    for line in gtf_file:
        if line.startswith('#'):
            continue
        else:
            gtf = GTFRecord(line)
            if gtf.feature == 'gene':
                gene_annotations[gtf.gene_name] = gtf
    return gene_annotations


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


def first_mismatch(str1, str2):
    """
    Returns the position of the first mismatch between 2 strings.

    :param str1: The first string
    :param str2: The second string
    :return: Position of the first mismatch else -1 if no mismatches

    >>> first_mismatch('abc', 'def')
    0
    >>> first_mismatch('abc', 'abe')
    2
    >>> first_mismatch('abc', 'abc')
    -1
    >>> first_mismatch('abc', 'abcde')
    -1
    """
    for i, j in enumerate(zip(str1, str2)):
        if j[0] != j[1]:
            return i
    else:
        return -1


def is_gzipfile(filename):
    """
    Attempt to ascertain the gzip status of a file based on the "magic signatures" of the file.

    This was taken from the stack overflow post
    http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type\
        -and-uncompress

    :param str filename: A path to a file
    :return: True if the file appears to be gzipped else false
    :rtype: bool
    """
    assert os.path.exists(filename), 'Input {} does not point to a file.'.format(filename)
    with open(filename, 'rb') as in_f:
        start_of_file = in_f.read(3)
        if start_of_file == '\x1f\x8b\x08':
            return True
        else:
            return False


def file_type(filename):
    """
    This module is used to open an input file in the appropriate method dependiing on whether it is a
    regular file, a gzipped
    file or a url.

    :param str filename: The path to the file
    :return: an open file handle to the file
    :rtype: file
    """
    if is_gzipfile(filename):
        return gzip.GzipFile(filename, 'r')
    else:
        return open(filename, 'r')


def get_snpeff_3_6_changes(eff):
    """
    Get the changes for the call from a snpEff 3.6 annotated vcf.

    :param str eff: The EFF field in the INFO column
    :return: The Effects of the mutation on the given protein
    :rtype: list(EFF)
    """
    out_effs = []
    eff = re.sub('EFF=', '', eff)
    for e in eff.split(','):
        if 'protein_coding' in e:
            if e.startswith(('NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING', 'STOP_GAINED')):
                e = e.strip().split('|')
                out_effs.append(EFF(annotation=snpeff_to_so[e[0].split('(')[0]],
                                    codon_change=e[2],
                                    aa_change=e[3],
                                    gene_name=e[5],
                                    transcript_id=e[8],
                                    transcript_len=e[4],
                                    warnings=e[-1].rstrip(')') if len(e) == 12 else ''))
            elif e.startswith(('CODON', 'FRAME')):
                e = e.strip().split('|')
                out_effs.append(EFF(annotation=snpeff_to_so[e[0].split('(')[0]],
                                    codon_change='',
                                    aa_change=e[3],
                                    gene_name=e[5],
                                    transcript_id=e[8],
                                    transcript_len=e[4],
                                    warnings=e[-1].rstrip(')') if len(e) == 12 else ''))
    return out_effs


def get_snpeff_4_0_changes(eff):
    """
    Get the changes for the call from a snpEff 4.0 annotated vcf.

    :param str eff: The EFF field in the INFO column
    :return: The Effects of the mutation on the given protein
    :rtype: list(EFF)
    """
    out_effs = []
    eff = re.sub('EFF=', '', eff)
    for e in eff.split(','):
        if 'protein_coding' in e:
            if e.startswith(('missense_variant', 'synonymous_variant', 'stop_gained')):
                e = e.strip().split('|')
                out_effs.append(EFF(annotation=e[0].split('(')[0],
                                    codon_change=e[2],
                                    aa_change=translate_aa_change_4_0(e[3]),
                                    gene_name=e[5],
                                    transcript_id=e[8],
                                    transcript_len=e[4],
                                    warnings=e[-1].rstrip(')') if len(e) == 12 else ''))
            elif e.startswith(('frame', 'inframe', 'disruptive')):
                e = e.strip().split('|')
                # Sinve SNPEff 4.0 does things differently, we will handle all indels similarly to
                # frameshifts (otherwise we will predict the wrong things based on the annotation.
                # The only downside to this is that we will predict excessively long peptides if the
                # user specified --indel_extend_length > n
                aa_change = list(translate_aa_change_4_0(e[3]))
                while True:
                    temp = aa_change.pop()
                    if temp.isdigit():
                        aa_change.append(temp + '?')
                        aa_change = ''.join(aa_change)
                        break
                annotation = 'frameshift_variant'
                out_effs.append(EFF(annotation=annotation,
                                    codon_change='',
                                    aa_change=aa_change,
                                    gene_name=e[5],
                                    transcript_id=e[8],
                                    transcript_len=e[4],
                                    warnings=e[-1].rstrip(')') if len(e) == 12 else ''))
    return out_effs


def translate_aa_change_4_0(aa_change_field):
    """
    Translate an amino acid change in snpeff 4.0 to 3.6 form.

    :param str aa_change_field: The field containing the aa change
    :return: The changed amino acid(s) in the for <REF><POS><ALT>
    :rtype: str

    >>> translate_aa_change_4_0('p.Gly144_Cys145del/c.431_433delGCT')
    'GC144'
    >>> translate_aa_change_4_0('p.Asp479_Ser480fs')
    'DS479?'
    >>> translate_aa_change_4_0('p.His87_Ala88insGlnLeu')
    'HA87QL'
    >>> translate_aa_change_4_0('p.Gly100Arg')
    'G100R'
    >>> translate_aa_change_4_0('p.Trp102fs')
    'W102?'
    >>> translate_aa_change_4_0('p.Tyr17*')
    'Y17*'
    """
    if '/' in aa_change_field:
        aa_change_field = aa_change_field.split('/')[0]

    aa_change_field = aa_change_field[2:]  # remove the leading "p."

    if any(x in aa_change_field for x in ['ins', 'del', 'fs']):
        # Indel
        aa_change_field = aa_change_field.split('_')

        # handle the last one first, then move towards the first
        prefix, aa_change_field[-1] = amino_letter[aa_change_field[-1][:3]], aa_change_field[-1][3:]
        pos = aa_change_field[-1].translate(None, string.letters)
        aa_change_field[-1] = aa_change_field[-1].translate(None, string.digits)

        if aa_change_field[-1] == 'fs':
            suffix = '?'
        else:
            suffix = ''.join([amino_letter[aa_change_field[-1][x:x + 3]]
                              for x in range(0, len(aa_change_field[-1]), 3)
                              if aa_change_field[-1][x:x + 3] not in ['del', 'ins']])

        aa_change_field.pop()
        for change in aa_change_field:
            prefix = amino_letter[change[:3]] + prefix
            pos = change[3:]
    else:
        prefix, aa_change_field = amino_letter[aa_change_field[:3]], aa_change_field[3:]
        if aa_change_field.endswith('*'):
            pos, suffix = aa_change_field[:-1], '*'
        else:
            pos, suffix = aa_change_field[:-3], amino_letter[aa_change_field[-3:]]

    return prefix + pos + suffix


def get_snpeff_4_x_changes(eff):
    """
    Get the changes for the call from a snpEff 4.1+ annotated vcf.

    :param str eff: The ANN field in the INFO column
    :return: The Effects of the mutation on the given protein
    :rtype: list(EFF)
    """
    out_effs = []
    eff = re.sub('ANN=', '', eff)
    for e in eff.split(','):
        if 'protein_coding' in e:
            e = e.strip().split('|')
            if e[1] in ('missense_variant', 'synonymous_variant', 'stop_gained'):
                out_effs.append(EFF(annotation=e[1],
                                    codon_change=e[9],
                                    aa_change=translate_aa_change_4_x(e[10]),
                                    gene_name=e[3],
                                    transcript_id=e[6],
                                    transcript_len=e[13].split('/')[1],
                                    warnings=e[-1].rstrip(')') if len(e) == 16 else ''))
            elif any(x in e[1] for x in ('frame', 'inframe', 'disruptive')):
                out_effs.append(EFF(annotation=e[1],
                                    codon_change=e[9],
                                    aa_change=translate_aa_change_4_x(e[10]),
                                    gene_name=e[3],
                                    transcript_id=e[6],
                                    transcript_len=e[13].split('/')[1],
                                    warnings=e[-1].rstrip(')') if len(e) == 16 else ''))
    return out_effs


def translate_aa_change_4_x(aa_change_field):
    """
    Translate an amino acid change in snpeff 4.1+ to 3.6 form.

    :param str aa_change_field: The field containing the aa change
    :return: The changed amino acid(s) in the for <REF><POS><ALT>
    :rtype: str

    >>> translate_aa_change_4_x('p.Gly144_Cys145del/c.431_433delGCT')
    'GC144'
    >>> translate_aa_change_4_x('p.Asp479_Ser480fs')
    'DS479?'
    >>> translate_aa_change_4_x('p.His87_Ala88insGlnLeu')
    'A88QLA'
    >>> translate_aa_change_4_x('p.Gly100Arg')
    'G100R'
    >>> translate_aa_change_4_x('p.Trp102fs')
    'W102?'
    >>> translate_aa_change_4_x('p.Tyr17*')
    'Y17*'
    >>> translate_aa_change_4_x('p.Tyr174_Arg175delinsTrp')
    'YR174W'
    >>> translate_aa_change_4_x('p.Arg175delinsTrpLeu')
    'R175WL'
    """
    if '/' in aa_change_field:
        aa_change_field = aa_change_field.split('/')[0]

    aa_change_field = aa_change_field[2:]  # remove the leading "p."

    if any(x in aa_change_field for x in ['ins', 'del', 'fs']):
        # Indel
        aa_change_field = aa_change_field.split('_')

        # handle the last one first, then move towards the first
        prefix, aa_change_field[-1] = amino_letter[aa_change_field[-1][:3]], aa_change_field[-1][3:]
        pos = aa_change_field[-1].translate(None, string.letters)
        aa_change_field[-1] = aa_change_field[-1].translate(None, string.digits)

        if aa_change_field[-1] == 'fs':
            suffix = '?'
        else:
            suffix = ''
            if aa_change_field[-1].startswith('delins'):
                aa_change_field[-1] = aa_change_field[-1][6:]  # Strip "delins"
            elif aa_change_field[-1] == 'del':
                aa_change_field[-1] = ''  # Strip "del"
            else:
                assert aa_change_field[-1].startswith('ins')
                aa_change_field = [aa_change_field[-1][3:]]  # Strip "ins"
                suffix = prefix  # Make it so the results will be C>ABC
            suffix = ''.join([amino_letter[aa_change_field[-1][x:x + 3]]
                              for x in range(0, len(aa_change_field[-1]), 3)]) + suffix
        aa_change_field.pop()
        for change in aa_change_field:
            prefix = amino_letter[change[:3]] + prefix
            pos = change[3:]
    else:
        prefix, aa_change_field = amino_letter[aa_change_field[:3]], aa_change_field[3:]
        if aa_change_field.endswith('*'):
            pos, suffix = aa_change_field[:-1], '*'
        else:
            pos, suffix = aa_change_field[:-3], amino_letter[aa_change_field[-3:]]

    return prefix + pos + suffix


def get_vep_changes(eff):
    """
    Get the changes for the call from a VEP annotated vcf.

    :param str eff: The CSQ field in the INFO column
    :return: The Effects of the mutation on the given protein
    :rtype: list(EFF)
    """
    out_effs = []
    eff = re.sub('ANN=', '', eff)
    for e in eff.split(','):
        if 'protein_coding' in e:
            e = e.strip().split('|')
            if e[1] in ('missense_variant', 'synonymous_variant', 'stop_gained'):
                out_effs.append(EFF(annotation=e[1],
                                    codon_change=e[16],
                                    aa_change=translate_aa_change_vep(e[14], e[15]),
                                    gene_name=e[3],
                                    transcript_id=e[6],
                                    transcript_len=-1,
                                    warnings=''))
            elif any(x in e[1] for x in ('frame', 'inframe', 'protein_altering')):
                out_effs.append(EFF(annotation=e[1],
                                    codon_change=e[16],
                                    aa_change=translate_aa_change_vep(e[14], e[15]),
                                    gene_name=e[3],
                                    transcript_id=e[6],
                                    transcript_len=-1,
                                    warnings=''))
    return out_effs


def translate_aa_change_vep(pos_field, aa_change_field):
    """
    Translate an amino acid change in snpeff 4.1+ to 3.6 form.

    :param str aa_change_field: The field containing the aa change
    :return: The changed amino acid(s) in the for <REF><POS><ALT>
    :rtype: str

    >>> translate_aa_change_vep('144-145', 'GC/-')
    'GC144'
    >>> translate_aa_change_vep('479-480', 'DS/X')
    'DS479?'
    >>> translate_aa_change_vep('87-88', '-/QL')
    '-87QL'
    >>> translate_aa_change_vep('100', 'G/R')
    'G100R'
    >>> translate_aa_change_vep('102-103', 'W/X')
    'W102?'
    >>> translate_aa_change_vep('17', 'Y/*')
    'Y17*'
    >>> translate_aa_change_vep('174', 'YR/W')
    'YR174W'
    >>> translate_aa_change_vep('175', 'R/WL')
    'R175WL'
    >>> translate_aa_change_vep('175', 'R')
    'R175R'
    """
    aa_change_field = aa_change_field.split('/')
    if len(aa_change_field) == 1:
        # Synonymous mutation
        aa_change_field = aa_change_field + aa_change_field

    pos_field = pos_field.split('/')[0].split('-')[0]
    if aa_change_field[1].endswith('X'):
        aa_change_field[1] = '?'
    elif aa_change_field[1] == '-':
        aa_change_field[1] = ''
    prefix, pos, suffix = aa_change_field[0], pos_field, aa_change_field[1]
    return prefix + pos + suffix
