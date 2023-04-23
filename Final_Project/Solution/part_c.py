from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from Bio.Data import CodonTable
import warnings

global_gencode = {
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


# This function counts the synonymous substitutions for each codon in the given gencode dictionary
def count_synonymous_substitutions(gencode):
    nuc_options = ["A", "C", "G", "T"]
    result_dict = dict()

    # looping through the genbank codons
    for key, value in gencode.items():

        not_stop_codons_counter = 0
        s_substitutions_counter = 0

        # looping through the three possible positions of the nucleotides
        for position in range(0, 3):

            rest_nucs = list(filter(lambda nuc: nuc != str(key[position]), nuc_options))

            # looping through the rest nucleotides options
            # and modifying the counters according to the new_codon and mutation
            for curr_nuc in rest_nucs:
                new_codon = key[:position] + curr_nuc + key[position + 1:]

                s_substitutions_counter += 1 if gencode[new_codon] == value else 0

                not_stop_codons_counter += 1 if gencode[new_codon] != '_' else 0

        # calculating the synonymous substitutions according to the low from the lectures
        result_dict[key] = 3 * s_substitutions_counter / not_stop_codons_counter

    return result_dict


# This function returns genes coding statistics of the given genbank_record
# if as_sets=True it returns data structures of: 1- genes names, 2- genes that could be converted to protine,
# 3- genes that could not be converted to protein, 4- genes that can be converted to multi proteins [as dictionary gene_name: number_of_additional_proteins_option]
# if as_sets=False it returns the lengths of the mentioned above data_structures [numeric statistics]
def get_gb_coding_statistics(genbank_record, as_sets=False):
    genes_set = set()
    to_protein_genes_set = set()
    multi_protein_genes = dict()

    for feature in genbank_record.features:

        if feature.type.lower() == "gene":
            genes_set.add(feature.qualifiers["gene"][0])

        elif feature.type.lower() == "cds":

            if feature.qualifiers["gene"][0] in to_protein_genes_set:
                if feature.qualifiers["gene"][0] in multi_protein_genes.keys():
                    multi_protein_genes[feature.qualifiers["gene"][0]] += 1
                else:
                    multi_protein_genes[feature.qualifiers["gene"][0]] = 1

            else:
                to_protein_genes_set.add(feature.qualifiers["gene"][0])

    not_to_protein_genes_set = genes_set.difference(to_protein_genes_set)

    if as_sets:
        return genes_set, to_protein_genes_set, not_to_protein_genes_set, multi_protein_genes

    return len(genes_set), len(to_protein_genes_set), len(not_to_protein_genes_set), len(multi_protein_genes)


# This function prints the genes coding statistics that the function get_gb_coding_statistics() returns
def print_gb_coding_statistics(file_path, as_sets=False):
    with open(file_path, "r") as input_handle:
        gen = SeqIO.parse(input_handle, "genbank")
        genbank_record = next(gen)

        all_genes, to_protein_genes, not_to_protein_genes, multi_protein_genes = get_gb_coding_statistics(
            genbank_record, as_sets)

        print(f"File Path: {file_path}")
        print(f"All Genes: {all_genes}")
        print(f"To Protine Genes: {to_protein_genes}")
        print(f"Not To Protine Genes: {not_to_protein_genes}")
        print(f"Multi Protine Genes: {multi_protein_genes}")


# This function returns genbank records comparison statistics of the given two genbank records
# if as_sets=True it returns data structures of: 1- shared genes, 2- genes that appear just in the first record
# 3- genes that appear just in the second record
# if as_sets=False it returns the lengths of the mentioned above data_structures [numeric statistics]
def get_gbs_comparison_statistics(first_genbank_record, second_genbank_record, as_sets=False):
    first_all_genes, first_to_protein_genes, first_not_to_protein_genes, first_multi_protein_genes = get_gb_coding_statistics(
        first_genbank_record, as_sets=True)

    second_all_genes, second_to_protein_genes, second_not_to_protein_genes, second_multi_protein_genes = get_gb_coding_statistics(
        second_genbank_record, as_sets=True)

    shared_genes = first_all_genes.intersection(second_all_genes)
    just_in_first_genes = first_all_genes.difference(second_all_genes)
    just_in_second_genes = second_all_genes.difference(first_all_genes)

    if as_sets:
        return shared_genes, just_in_first_genes, just_in_second_genes

    return len(shared_genes), len(just_in_first_genes), len(just_in_second_genes)


# This function prints the genes comparison statistics that the function get_gbs_comparison_statistics() returns
def print_gb_comparison_statistics(first_genbank_file_path, second_genbank_file_path, as_sets=False):
    with open(first_genbank_file_path, "r") as first_file, open(second_genbank_file_path, "r") as second_file:
        first_gen = SeqIO.parse(first_file, "genbank")
        first_genbank_record = next(first_gen)

        second_gen = SeqIO.parse(second_file, "genbank")
        second_genbank_record = next(second_gen)

        shared_genes, just_in_first_genes, just_in_second_genes = get_gbs_comparison_statistics(first_genbank_record,
                                                                                                second_genbank_record,
                                                                                                as_sets)

        print(f"Shared Genes: {shared_genes}")
        print(f"Just In First Genes: {just_in_first_genes}")
        print(f"Just In Second Genes: {just_in_second_genes}")


# This function returns the translation for the given genes names to_protine_genes_names
# it gets these genes from the given genbank file path then it translates them
# after that it compares its translation with the translations of the gen bank file for the same genes
# if there is a different it will keep its translation, but it will return the names of the genes that has different translations
# at the end it returns the translations of the genes as dictionary [ gene_name: translation ]
# a set of the genes names for the genes that have different translation in the genbank file compared to our translation
def get_genes_translation(genbank_file_path, to_protine_genes_names):
    with open(genbank_file_path, "r") as genbank_file:
        genbank_parser = SeqIO.parse(genbank_file, "genbank")
        genbank_record = next(genbank_parser)
        full_genbank_seq = genbank_record.seq
        translations_result_dict = dict()
        different_translation_genes = set()

        for feature in genbank_record.features:

            if feature.type.lower() == "gene":
                if feature.qualifiers["gene"][0] in to_protine_genes_names:
                    curr_gene_seq = full_genbank_seq[feature.location.start: feature.location.end]

                    if feature.location.strand == 1:
                        curr_gene_translation = curr_gene_seq.translate(table=11, cds=True)
                    else:
                        curr_gene_seq = curr_gene_seq.reverse_complement()
                        curr_gene_translation = curr_gene_seq.translate(table=11, cds=True)

                    translations_result_dict[feature.qualifiers["gene"][0]] = curr_gene_translation

        for feature in genbank_record.features:

            if feature.type.lower() == "cds":
                if feature.qualifiers["gene"][0] in translations_result_dict.keys():
                    if translations_result_dict[feature.qualifiers["gene"][0]] != feature.qualifiers["translation"][0]:
                        different_translation_genes.add(feature.qualifiers["gene"][0])

        return translations_result_dict, different_translation_genes


# This function returns the aligned protine sequence to nucleotides sequence based on the given original nucleotides sequence
def back_to_nucleotides(original_nucleotides_seq, aligned_codons_seq):
    nucleotides_str = ""
    original_pointer = 0

    for codon in aligned_codons_seq:
        if codon == "-":
            nucleotides_str += "---"
        else:
            nucleotides_str += original_nucleotides_seq[original_pointer: original_pointer + 3]
            original_pointer += 3

    return nucleotides_str


# This function returns the genes sequences from the genbank file of the given genes names [genes_list]
def get_genes_sequences(genbank_file_path, genes_list):
    genes_sequences_dict = dict()

    with open(genbank_file_path, "r") as genbank_file:
        genbank_parser = SeqIO.parse(genbank_file, "genbank")
        genbank_record = next(genbank_parser)
        full_genbank_seq = genbank_record.seq

        for feature in genbank_record.features:
            if feature.type.lower() == "gene":
                if feature.qualifiers["gene"][0] in genes_list:
                    genes_sequences_dict[feature.qualifiers["gene"][0]] = full_genbank_seq[
                                                                          feature.location.start: feature.location.end]

        return genes_sequences_dict


# This function makes comparison between the genes that given in the list shared_genes_to_compare
# from the two given files so for each pair of genes (the same gene but from two different files)
# it creates a tuple (dN, dS, dN_dS_ratio) it returns dictionary like [gene_name: (dN, dS, dN_dS_ratio)]
# ofcourse it calls the other functions that we wrote to accomplish this task
def shared_genes_comparison(first_genbank_file_path, second_genbank_file_path, shared_genes_to_compare):
    comparison_results_dict = dict()

    first_translations_dict, first_different_translation_genes = get_genes_translation(first_genbank_file_path,
                                                                                       shared_genes_to_compare)
    second_translations_dict, second_different_translation_genes = get_genes_translation(second_genbank_file_path,
                                                                                         shared_genes_to_compare)

    first_genes_sequences_dict = get_genes_sequences(first_genbank_file_path, shared_genes_to_compare)
    second_genes_sequences_dict = get_genes_sequences(second_genbank_file_path, shared_genes_to_compare)

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    for gene in shared_genes_to_compare:
        alignment = aligner.align(first_translations_dict[gene], second_translations_dict[gene])[0]
        alignment_strs_list = alignment.__str__().split('\n')

        first_nucleotides_alignment = str(back_to_nucleotides(first_genes_sequences_dict[gene], alignment_strs_list[0]))
        second_nucleotides_alignment = str(
            back_to_nucleotides(second_genes_sequences_dict[gene], alignment_strs_list[2]))

        first_seq = CodonSeq(first_nucleotides_alignment)
        second_seq = CodonSeq(second_nucleotides_alignment)

        dN, dS = cal_dn_ds(first_seq, second_seq, codon_table=CodonTable.generic_by_id[11])

        if dS != 0:
            dN_dS_ratio = float(dN / dS)
        else:
            dN_dS_ratio = 0

        comparison_results_dict[gene] = (dN, dS, dN_dS_ratio)

    return comparison_results_dict


# This function prints the shared genes comparison that the function shared_genes_comparison() returns
def print_shared_genes_comparison(first_genbank_file_path, second_genbank_file_path, shared_genes_to_compare):
    comparison_results_dict = shared_genes_comparison(first_genbank_file_path, second_genbank_file_path,
                                                      shared_genes_to_compare)

    for key, value in comparison_results_dict.items():
        print(f"Gene_Name: {key}")
        print(f"dN: {value[0]}")
        print(f"dS: {value[1]}")
        print(f"dN_dS_ratio: {value[2]}")
        print(f"-----------------------")


# All what we do in the main is just calling the function that we wrote above
# to get the answers of part 3 of the Final_Job [Homework]
if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    first_genbank_file_path = "corona_2021.gb"
    second_genbank_file_path = "corona_2022.gb"

    print("--------------------- Synonymous Substitutions ---------------------\n")
    print(count_synonymous_substitutions(global_gencode))
    print()

    print("--------------------- Corona 2021 Coding Statistics ---------------------\n")
    print_gb_coding_statistics(first_genbank_file_path, as_sets=False)
    print()

    print("--------------------- Corona 2022 Coding Statistics ---------------------\n")
    print_gb_coding_statistics(second_genbank_file_path, as_sets=False)
    print()

    print("--------------------- Corona 2021/2022 Genes [Shared / Not Shared] Statistics ---------------------\n")
    print_gb_comparison_statistics(first_genbank_file_path, second_genbank_file_path, as_sets=False)
    print()

    print(
        "--------------------- Corona 2021/2022 Shared Genes Comparison [dN, dS, dN_dS_ratio] ---------------------\n")
    shared_to_protein_genes_to_compare = ['E', 'M', 'ORF8', 'N', 'ORF7b']
    print_shared_genes_comparison(first_genbank_file_path, second_genbank_file_path, shared_to_protein_genes_to_compare)
