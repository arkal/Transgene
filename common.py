import re


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

