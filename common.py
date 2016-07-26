import re


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
    if seq:
        if nucs_regexp.findall(seq):
            seq = re.sub(nucs_regexp, '', seq)
        yield [seq_id, seq_comments, seq.upper()]
