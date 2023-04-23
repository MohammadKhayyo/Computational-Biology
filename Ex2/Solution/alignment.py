import numpy as np
from io import StringIO
from Bio import Entrez, SeqIO, pairwise2, Align
from Bio.pairwise2 import format_alignment
import os


def get_record_name(record, xml=True):
    """get record name"""
    if xml:
        return record["TSeq_orgname"]
    else:
        return record.name


def get_record_seq(record, xml=True):
    """get record seq"""
    if xml:
        return record["TSeq_sequence"].upper()
    else:
        return record.seq.upper()


def percentage_Calculator(records, xml=False, only_one=True, m=None, s=None):
    """this function Calculate Matching between list of records and return dictionary that contain the result"""
    dic = {}
    first_record = None
    last_record = None
    for index, record in enumerate(records):
        if index == 0:
            first_record = record
        if last_record is not None:
            dic[f'{get_record_name(last_record, xml)} --> {get_record_name(record, xml)}'] = \
                max_best_matching(get_record_seq(last_record, xml), get_record_seq(record, xml), only_one, m, s)
            dic[f'{get_record_name(record, xml)} --> {get_record_name(last_record, xml)}'] = \
                max_best_matching(get_record_seq(record, xml), get_record_seq(last_record, xml), only_one, m, s)
        last_record = record
    if first_record is not None and last_record is not None:
        dic[f'{get_record_name(last_record, xml)} --> {get_record_name(first_record, xml)}'] = \
            max_best_matching(get_record_seq(last_record, xml), get_record_seq(first_record, xml), only_one, m, s)
        dic[f'{get_record_name(first_record, xml)} --> {get_record_name(last_record, xml)}'] = \
            max_best_matching(get_record_seq(first_record, xml), get_record_seq(last_record, xml), only_one, m, s)
    return dic


def load_data_from_genbank(email="jcecomputationalbiology@gmail.com",
                           mitochondrial_s=None):
    """load data from genbank"""
    if mitochondrial_s is None:
        mitochondrial_s = ["AF451972", "AF176731", "X90314"]
    Entrez.email = email  # Provide an email address
    with Entrez.efetch(db="nucleotide", id=mitochondrial_s,
                       rettype="fasta", retmode="xml") as handle:
        features = Entrez.read(handle)
    return features


def global_xx_mx_ms(x, y, m=None, s=None):
    """call global functions"""
    if not m:
        return pairwise2.align.globalxx(x, y)
    if m and not s:
        return pairwise2.align.globalmx(x, y, m[0], m[1])
    if not m and s:
        return pairwise2.align.globalxs(x, y, s[0], s[1])
    if m and s:
        return pairwise2.align.globalms(x, y, m[0], m[1], s[0], s[1])


def max_best_matching(x, y, only_one=True, m=None, s=None):
    """characterizes the various mutations that have occurred based on the resulting page"""
    max_matching = -1
    best_matching = ''
    alignments = global_xx_mx_ms(x, y, m, s)
    if not only_one:
        for i, current_align in enumerate(alignments):
            align1, align2, score, begin, end = current_align
            str_format_alignment = format_alignment(align1=align1, align2=align2, score=score, begin=begin, end=end)
            count_matching = str_format_alignment.count('|')
            if count_matching > max_matching:
                best_matching = str_format_alignment
                max_matching = count_matching
    if only_one:
        align1, align2, score, begin, end = alignments[0]
        best_matching = format_alignment(align1=align1, align2=align2, score=score, begin=begin, end=end)
        max_matching = best_matching.count('|')
    np_best_matching = np.loadtxt(StringIO(best_matching), dtype=str, delimiter=',')
    gaps_1 = np_best_matching[0].count('-')
    gaps_2 = np_best_matching[2].count('-')
    Transition = 0
    Transversion = 0
    dic_Transition = {'A': 'G', 'G': 'A', 'T': 'C', 'C': 'T'}
    for first, second, third in zip(np_best_matching[0], np_best_matching[1], np_best_matching[2]):
        if second == '.':
            if dic_Transition[first] == third:
                Transition += 1
            else:
                Transversion += 1
    dic_info = {'matches': max_matching, 'transitions': Transition, 'transversions': Transversion, 'gaps_1': gaps_1,
                'gaps_2': gaps_2}
    _best_matching = np_best_matching[0] + '\n' + np_best_matching[1] + '\n' + np_best_matching[2] + '\n'
    return {'Seq with Score': best_matching, 'Seq without Score': _best_matching,
            'matching': (max_matching / len(np_best_matching[0])) * 100, 'dic_info': dic_info}


def info(_dic):
    """print the results"""
    for key, value in _dic.items():
        print(f'{key}:')
        print(f'matching: {value["matching"]}%')
        print(f'{value["dic_info"]}')
        print(f'{value["Seq with Score"]}')


def main(_m=None, _s=None):
    """Calling functions"""
    fasta_path = 'ex2_sequences_a.fasta'
    assert (os.path.isfile(fasta_path))
    _dic_genbank = percentage_Calculator(load_data_from_genbank(), xml=True, only_one=True, m=_m, s=_s)
    info(_dic_genbank)
    _dic_from_fasta = percentage_Calculator(SeqIO.parse(fasta_path, "fasta"), only_one=True, m=_m, s=_s)
    info(_dic_from_fasta)


if __name__ == "__main__":
    print("count matches for globalxx:")
    main()
    _m = [2, -1]
    _s = [-2.5, -2.1]
    print("#" * 500)
    print(f"count matches for globalms: with match: {_m[0]}, mismatch: {_m[1]}, open: {_s[0]}, gap: {_s[1]}")
    main(_m=_m, _s=_s)
