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


def count_mutation_by_type(position, type):
    if position <= 0 or position > 3:
        print("position must be 1 or 2 or 3")
        return 0
    count = 0
    if type == 'nonsense':
        gencode_stop = {key: value for key, value in gencode.items() if value == '_'}
        for key_all, value_all in gencode.items():
            if value_all == '_':
                continue
            for key_stop, value_stop in gencode_stop.items():
                if key_all[:position - 1] + key_all[position:] == key_stop[:position - 1] + key_stop[position:]:
                    count += 1
    if type == 'synonymous':
        Nucleotides = ['A', 'G', 'T', 'C']
        dic_rest = {}
        for Nucleotide in Nucleotides:
            dic_rest[Nucleotide] = [item for item in Nucleotides if item != Nucleotide]
        for key_all, value_all in gencode.items():
            for key_with_other_Nucleotide in dic_rest[key_all[position - 1]]:
                key_without_position = key_all[:position - 1] + key_with_other_Nucleotide + key_all[position:]
                if gencode[key_without_position] == gencode[key_all]:
                    count += 1
    return count


print(count_mutation_by_type(3, 'synonymous'))
