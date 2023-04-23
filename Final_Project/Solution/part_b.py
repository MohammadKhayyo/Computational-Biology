from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings


def print_comparison(result_gb_diff, result_uni_diff, result_gb_uni):
    """ This function pritn all comparison for Section A"""

    print("comparison: ")
    print("common proteins number between GenBank and UniProt is:", len(result_gb_uni))
    print()

    print(len(result_gb_diff), "proteins are in GenBank file but not in UniProt file")
    print()

    print(len(result_uni_diff), 'proteins are in UniProt file but not in GenBank file')


def get_transmembrane_seq(xl_content):
    """This function get transmembrane sequences from UniProt file """

    # get all proteins sequences
    seq_list = xl_content["Sequence"].dropna().values.tolist()

    # get all proteins Transmembrane data if they have
    transmembrane_column = xl_content["Transmembrane"].dropna()

    # create index list
    index_list = transmembrane_column.index.values.tolist()

    # convert data to list
    transmembranes_list = transmembrane_column.values.tolist()

    # Initialise list
    transmembrane_seq_list = []

    # this loop get the range of transmembrane sequence in protein sequence and get the sequence
    for i in range(len(transmembranes_list)):
        str = transmembranes_list[i]
        for j in range(transmembranes_list[i].count('TRANSMEM')):
            start = str.find('TRANSMEM') + len('TRANSMEM')
            end = str.find(';')
            star_end = str[start:end].split('..')
            protien_seq = seq_list[int(index_list[i])]

            transmembrane_seq_list.append(protien_seq[int(star_end[0]) - 1: int(star_end[1])])
            str = str[end:]
            str = str[str.find('TRANSMEM'):]

    return transmembrane_seq_list


def get_data(xl_file_path):
    """This function get data from GeneBank and uniProt files"""
    # read from uniprot file
    xl_content = pd.read_excel(xl_file_path)

    # put the data in list
    copy_gene_name_list_uni = xl_content["Gene Names (primary)"].dropna().values.tolist()

    # Initialise list
    gene_name_list_uni = []

    # if there are more than one name for gene
    for gene_name in copy_gene_name_list_uni:
        if (';' in gene_name):
            gene_name_list_uni.extend(gene_name.split('; '))

        else:
            gene_name_list_uni.append(gene_name)

    # Initialise list
    gene_name_list_gb = []

    # read from GenBank file
    with open("BS168.gb", "r") as input_handle:
        gen = SeqIO.parse(input_handle, "genbank")
        record_gb = next(gen)  # content of 1st record
        features_list = record_gb.features
    for feature in features_list:
        if feature.type == 'gene' and 'gene' in feature.qualifiers.keys():
            gene_name_list_gb.extend(feature.qualifiers['gene'])

    return gene_name_list_uni, gene_name_list_gb, xl_content, features_list, record_gb


def Average(lst):
    """This function  caculate the Avarge"""

    return sum(lst) / len(lst)


def show_plot_transmembrane(transmembrane_seq_len):
    """This function show figure histogram of distribution of lengths for transmembrane sequences """

    # print statistic
    print("Max Transmembrane Sequence length:", max(transmembrane_seq_len))

    print("Min Transmembrane Sequence length:", min(transmembrane_seq_len))

    average = Average(transmembrane_seq_len)
    print("Average Transmembrane Sequence:", average)

    # show the figure
    plt.figure(1)
    plt.title("The histogram of distribution of lengths for transmembrane sequences")
    plt.ylabel("number of sequences")
    plt.xlabel("lengths")
    plt.hist(transmembrane_seq_len)


def calculate_at_percentage(sequence):
    """Calculates and returns the AT percentage of a DNA sequence."""
    num_at = sequence.count("A") + sequence.count("T")
    at_percentage = num_at / len(sequence) * 100
    return at_percentage


def create_B_Group(gene_name_list_uni, gene_name_list_gb, xl_content):
    """This function create group B,  The group of gene sequences found in the intersection between UniProt and GenBank,
    locate the genes that contain at least one transmembrane region"""

    # do intersection between UniProt and GenBank genes
    result_uni_and_gb = set(gene_name_list_uni) & set(gene_name_list_gb)

    # Initialise list
    B_list_gene_name = []

    # get genes that contain at least one transmembrane region and gene sequences found in the intersection between UniProt and GenBank files
    for gene_name in result_uni_and_gb:
        if len(xl_content[(xl_content["Gene Names (primary)"] == gene_name) & (
        xl_content["Transmembrane"].notnull())].values) > 0:
            B_list_gene_name.append(gene_name)

    return B_list_gene_name


def create_A_Group(features_list, record_gb):
    """This function create group A, The group of genes that are CDS"""

    # Initialise two list
    A_list = []
    A_list_name = []

    # This loop get name and sequences for genes that are CDS
    for feature in features_list:
        if feature.type == 'CDS':

            A_list.append(feature.location.extract(record_gb.seq))
            if 'gene' in feature.qualifiers.keys():
                A_list_name.extend(feature.qualifiers["gene"])

    return A_list, A_list_name


def get_group_seq(group_gene_name, features_list, record_gb):
    """This function receive list of gene name then return list of proteins sequences for this genes names from
    GeneBank file """
    # Initialise string
    seq_group = ''

    # Initialise list
    list_translation = []

    # tThis loop get the proteins sequences from GeneBank file
    for gene_name in group_gene_name:
        for feature in features_list:
            if feature.type == 'CDS' and 'gene' in feature.qualifiers.keys() and feature.qualifiers['gene'][
                0] == gene_name:
                # location = feature.location[]

                list_translation.append(feature.location.extract(record_gb.seq))

                # list_translation.append(feature.qualifiers['translation'][0])
                # seq_group +=feature.qualifiers['translation'][0]

    return list_translation


def seq_list_to_len_list(seq_list):
    """This function calculate length of sequences that receive them from list and return them in list"""

    # Initialise list
    len_list = []
    for i in range(len(seq_list)):
        len_list.append(len(seq_list[i]))

    return len_list


def seq_list_to_string(seq_list):
    """ this function convert list of sequences to string by marge them """

    # Initialise string
    str = ''

    # marge the sequences from list to string
    for i in range(len(seq_list)):
        str += seq_list[i]
    return str


def calculate_statistic(group_seq):
    """ this function calculate max, min, average, Standard deviation"""

    # calculate length of sequences and put them in list
    group_seq_len = seq_list_to_len_list(group_seq)

    # calculate maximum
    max_group = np.max(group_seq_len)

    # calculate minimum
    min_group = np.min(group_seq_len)

    # calculate average
    ava_group = np.mean(group_seq_len)

    # calculate Standard deviation
    std_group = np.std(group_seq_len)

    return max_group, min_group, ava_group, std_group


def convert_statistic_to_dataframe(A_group_seq, B_group_seq):
    """This function convert statistic for group A and group B to data frame"""

    # calculate max, min, avarage, Standard deviation
    max_A_group, min_A_group, ava_A_group, std_A_group = calculate_statistic(A_group_seq)
    max_B_group, min_B_group, ava_B_group, std_B_group = calculate_statistic(B_group_seq)

    # create data frame
    data = {'GroupA': [max_A_group, min_A_group, ava_A_group, std_A_group],
            'GroupB': [max_B_group, min_B_group, ava_B_group, std_B_group]}

    # Creates pandas DataFrame.
    df = pd.DataFrame(data, index=['Maximum', 'Minimum', 'Average', 'Standard deviation'])

    print(df)


def calculate_percentages_per_group(group_seq):
    """This function calculate percentages of AT in sequence"""

    # Initialise list
    list_group_percentages = []

    # calculate percentages of AT and put in list
    for i in range(len(group_seq)):
        list_group_percentages.append(calculate_at_percentage(group_seq[i]))

    return list_group_percentages


def plot_percentages(list_of_percentage_group_A, list_of_percentage_group_B, list_of_percentage_group_A_notIn_B):
    """This function draw a figures of histogram AT percentage for group A, group B , For groupA not in groupB ,
    group B and group A not in groupB """

    # Initialise the subplot function using number of rows and columns
    figure, axis = plt.subplots(2, 2)
    figure.suptitle('The histogram for AT%')
    # For group A
    axis[0, 0].hist(list_of_percentage_group_A)
    axis[0, 0].set_title("group A")

    # For group B
    axis[0, 1].hist(list_of_percentage_group_B)
    axis[0, 1].set_title("group B")

    # For group A not in groupB
    axis[1, 0].hist(list_of_percentage_group_A_notIn_B)
    axis[1, 0].set_title("group A not B ")

    # For group B and groupA not in group B
    axis[1, 1].hist(list_of_percentage_group_A_notIn_B)
    axis[1, 1].hist(list_of_percentage_group_B)
    axis[1, 1].set_title("group B and group A not B ")

    # do a labels for x and y axis
    for ax in axis.flat:
        ax.set(xlabel='AT%', ylabel='number of sequences')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axis.flat:
        ax.label_outer()
    plt.legend(['A not B', 'B'])


def calculate_hydrophobic_distribution(sequence):
    """This function calculate the hydrophobic distribution in sequence that receive the return the percentage """

    # Define a list of hydrophobic amino acids
    hydrophobic_amino_acids = ['A', 'I', 'L', 'M', 'F', 'W', 'P', 'V']

    # Iterate through the sequence and count the number of each amino acid
    total_count = 0
    for amino_acids in hydrophobic_amino_acids:
        total_count += sequence.count(amino_acids)

    # Calculate the distribution as the percentage of each amino acid in the sequence
    percentage = total_count / len(sequence) * 100

    return percentage


def Transmembrane_sequence_hydrophobic_statistic(list_Transmembrane_sequence):
    """This function calculate distribution of prcentge hydrophobic amino acid in Transmembrane sequences and
    caculate the avarage """

    # Initialise list
    percentage_of_hydrophobic_list = []

    # calculate hydrophobic distribution and put it in list
    for seq in list_Transmembrane_sequence:
        percentage_of_hydrophobic_list.append(calculate_hydrophobic_distribution(seq))

    str_all_Transmembrane_seq = seq_list_to_string(list_Transmembrane_sequence)
    print()
    print("The percentage of hydrophobic amino acid in all Transmembrane sequences",
          calculate_hydrophobic_distribution(str_all_Transmembrane_seq))
    print("The average value of percentage hydrophobicover amino acid in Transmembrane sequences",
          np.mean(percentage_of_hydrophobic_list))


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    xl_file_path = 'UniProt_data.xlsx'

    # Section A: the Cross between the proteins from the GeneBank file and the proteins from the UniP file

    print("============================================= Section A ==============================================")
    gene_name_list_uni, gene_name_list_gb, xl_content, features_list, record_gb = get_data(xl_file_path)
    result_uni_diff = set(gene_name_list_uni) - set(gene_name_list_gb)
    result_gb_diff = set(gene_name_list_gb) - set(gene_name_list_uni)
    result_gb_uni = set(gene_name_list_gb) & set(gene_name_list_uni)

    print_comparison(result_gb_diff, result_uni_diff, result_gb_uni)

    # Section B: pulled out the transmembrane sequences: 1- Caculate distribution of their lengths and show it in
    # histogram, and calculate statistic maximum , minumum , avarage ,Standard deviation. 2- Caculate the
    # distribution of the percentage of hydrophobic amino acids in transmembrane sequences, and the average value
    # across all these sequences.

    print("\n============================================= Section B ==============================================\n")
    transmembrane_seq_list = get_transmembrane_seq(xl_content)
    transmembrane_seq_len = seq_list_to_len_list(transmembrane_seq_list)
    show_plot_transmembrane(transmembrane_seq_len)
    Transmembrane_sequence_hydrophobic_statistic(transmembrane_seq_list)

    # Section C: divied the proteins fot two groups A and B, calculate statistic and show distribution of the percentage AT in histogram

    print("\n============================================= Section C ==============================================\n")
    B_group_gene_name = create_B_Group(gene_name_list_uni, gene_name_list_gb, xl_content)
    A_group_seq, A_group_gene_name = create_A_Group(features_list, record_gb)
    B_group_seq = get_group_seq(B_group_gene_name, features_list, record_gb)

    A_notIn_B_group_gene_name = list(set(A_group_gene_name) - set(B_group_gene_name))
    A_notIn_B_group_seq = get_group_seq(A_notIn_B_group_gene_name, features_list, record_gb)

    convert_statistic_to_dataframe(A_group_seq, B_group_seq)
    list_of_percentage_group_A = calculate_percentages_per_group(A_group_seq)
    list_of_percentage_group_B = calculate_percentages_per_group(B_group_seq)
    list_of_percentage_group_A_notIn_B = calculate_percentages_per_group(A_notIn_B_group_seq)
    group_B_str = seq_list_to_string(B_group_seq)
    print("AT% in the gene sequences in group B", calculate_at_percentage(group_B_str))
    plot_percentages(list_of_percentage_group_A, list_of_percentage_group_B, list_of_percentage_group_A_notIn_B)

    plt.show()
