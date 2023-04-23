import re

gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}
dna_reverse_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def is_valid_dna(seq):
    """The function accepts a string and returns True if all the letters in the string are valid nucleotides (ie A or
    C or G or T). Otherwise, the function returns False. """
    for Nucleotide in seq.upper():
        if Nucleotide != 'A' and Nucleotide != 'T' and Nucleotide != 'G' and Nucleotide != 'C':
            return False
    return True


def get_gc_content(seq):
    """The function receives a DNA sequence and returns %GC, i.e. the percentage of nucleotides that are G or C"""
    if not is_valid_dna(seq):
        return -1
    count_G_C = 0
    for Nucleotide in seq.upper():
        if Nucleotide == 'G' or Nucleotide == 'C':
            count_G_C += 1
    return (count_G_C / len(seq)) * 100


def reverse_complement(seq):
    """The function receives a string (a DNA sequence) and returns the sequence that is the reverse complement for
    it. For example, for the parameter 'ACGTTG' the function will return 'CAACGT'. """
    dna_reverse_complement_str = ''
    if not is_valid_dna(seq):
        return dna_reverse_complement_str
    for Nucleotide in seq.upper():
        dna_reverse_complement_str += dna_reverse_complement[Nucleotide]
    return dna_reverse_complement_str[::-1]


def get_transcription(seq, strand):
    """The function accepts two parameters: A. A string (DNA sequence) on the positive strand that we want to
    reproduce. B. strand (of type int: the positive strand is represented by 1 and the negative by -1) The function
    returns the RNA sequence that is copied from the seq sequence. If the strand is -1 this means that the coded
    sequence is on the complementary strand """
    if not is_valid_dna(seq):
        return ''

    def dna_strand_positive_to_mRNA(_seq):
        mRNA = ''
        for Nucleotide in _seq.upper():
            if Nucleotide == 'T':
                mRNA += 'U'
            else:
                mRNA += Nucleotide
        return mRNA

    if strand == 1:
        return dna_strand_positive_to_mRNA(seq)
    elif strand == -1:
        reverse_complement_srt = reverse_complement(seq)
        return dna_strand_positive_to_mRNA(reverse_complement_srt)
    else:
        return ''


def stats_amino_acids():
    """The genetic code is a table that maps codons to amino acids. There are 64 codons and 20 different amino acids.
    We would like to count for each amino acid, how many codons code for it. A function acids_amino_stats which
    returns a dictionary with these counts. The key in the dictionary will be the name of an amino acid (one letter)
    and the value will be the corresponding count. It also counts how many termination codons there are, indicating
    this by the key END. """
    dic_count = {}
    for value in gencode.values():
        _value = ''
        if value == '_':
            _value = 'END'
        else:
            _value = value
        if _value in dic_count.keys():
            dic_count[_value] += 1
        else:
            dic_count[_value] = 1
    return dic_count


def translate_seq(seq, reading_frame=None):
    """the function  translates a DNA sequence into a protein. The function receives a string that is a DNA sequence
    (checks that the sequence is a valid DNA sequence, if not, return an empty list) and translates the sequence into
    a protein. If frame_reading is None , the function returns a list of the three possible protein sequences that
    the sequence may encode. Otherwise, if frame_reading is a number between 1 and 3, the function will calculate the
    corresponding protein when reading starts at position frame_reading. In this case, the function will return a
    list of length one, containing the protein sequence. """
    if not is_valid_dna(seq):
        return []

    def codons_to_proteins(_seq):
        list_codons = re.findall('...', _seq)
        protein = ''
        for codon in list_codons:
            protein += gencode[codon.upper()]
        return protein

    if reading_frame == 1:
        return [codons_to_proteins(seq[:])]
    elif reading_frame == 2:
        return [codons_to_proteins(seq[1:])]
    elif reading_frame == 3:
        return [codons_to_proteins(seq[2:])]
    elif reading_frame is None:
        return [codons_to_proteins(seq[:]), codons_to_proteins(seq[1:]), codons_to_proteins(seq[2:])]
    else:
        return []


def count_codons(seq):
    """The function receives a string that is a DNA sequence. The function counts for each codon in the sequence,
    how many times it appears in the sequence, and returns a dictionary with the corresponding counts (the key is a
    codon, the value is the number of occurrences). For example, for the sequence 'ACGTGTTGT', the dictionary is
    returned: {' ACG': 1, 'TGT': 2 } ) That is, codons must be counted without overlap. """
    dic = {}
    list_codons = re.findall('...', seq.upper())
    for codon in list_codons:
        if codon in dic.keys():
            dic[codon] += 1
        else:
            dic[codon] = 1
    return dic
