import os
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import warnings


def parse_genbank(gene_bank_file):
    """The parse_genbank() function takes in a file path gene_bank_file and returns a single record of a genome in
    the GenBank format. It first checks that the file path is valid, then opens and parses the file using the
    SeqIO.parse() method and extracts the first record of the genome using the next() function. This record is then
    returned as the output of the function. This function uses the python's built-in os and Bio library's SeqIO
    modules. """
    assert (os.path.exists(gene_bank_file))  # check that the path is valid
    with open(gene_bank_file, "r") as input_handle:
        genome = SeqIO.parse(input_handle, "genbank")
        record_gb = next(genome)
    return record_gb


def create_dataframe(record_gb):
    """The create_dataframe() function takes in a single GenBank record record_gb and returns a Pandas dataframe
    containing information about the features in the record. """
    dictionary_before_df = {'type': [], 'location': [], 'start': [], 'end': [], 'strand': []}
    name_of_genes = {}
    for index, feature in enumerate(record_gb.features):
        dictionary_before_df['type'].append(feature.type)
        dictionary_before_df['start'].append(feature.location.start)
        dictionary_before_df['end'].append(feature.location.end)
        dictionary_before_df['strand'].append(feature.strand)
        dictionary_before_df['location'].append(feature.location)
        if feature.type == 'source':
            name_of_genes[index] = ''
            continue
        if 'gene' in feature.qualifiers:
            name_of_genes[index] = feature.qualifiers['gene'][0]
        else:
            name_of_genes[index] = feature.qualifiers['locus_tag'][0]
    df = pd.DataFrame.from_dict(dictionary_before_df)
    df['length'] = df['end'] - df['start']
    df['name'] = pd.Series(name_of_genes)
    df = df.sort_values(by="start")
    return df


def count_elements(df):
    """The count_elements() function takes in a Pandas dataframe df and returns a dictionary containing the count of
    occurrences of each unique element in the 'type' column of the dataframe. """
    dic = df['type'].value_counts().to_dict()
    dic.pop('source')
    return dic


def draw_Multiple_histogram(all_genes_coding, protein_coding, other_coding, title, xlabel, ylabel):
    """The draw_multiple_histogram() function takes in four lists of numbers all_genes_coding, protein_coding,
    other_coding and three strings title, xlabel, ylabel and creates a 2x2 subplot of histograms using the Matplotlib
    library. """
    figure, axis = plt.subplots(2, 2)
    figure.suptitle(title)
    axis[0, 0].hist(all_genes_coding, bins=250)
    axis[0, 0].set_title("All genes coding")
    axis[0, 1].hist(protein_coding, bins=250)
    axis[0, 1].set_title("Protein coding")
    axis[1, 0].hist(other_coding, bins=25)
    axis[1, 0].set_title("Other coding")
    axis[1, 1].hist(protein_coding, bins=250)
    axis[1, 1].hist(other_coding, bins=25)
    axis[1, 1].set_title("Protein coding and Other coding")
    for ax in axis.flat:
        ax.set_ylim([0, 350])
        ax.set_xlim([0, 3000])
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.label_outer()
    plt.legend(['protein_coding', 'other_coding'])
    plt.show()


def calculate_AT_composition_df(df, sequence):
    """The calculate_AT_composition_df() function takes in a Pandas dataframe df and a DNA sequence ,
    and calculates the AT composition of each feature in the dataframe. """
    at_percentages = {}
    for index in range(df.shape[0]):
        seq = df.iloc[index]['location'].extract(sequence)
        at_percentages[index] = ((seq.count('A') + seq.count('T')) / len(seq)) * 100
    df['AT_composition'] = pd.Series(at_percentages)
    return df


def divide_genes(df):
    """The divide_genes() function takes in a Pandas dataframe df and returns 3 separate dataframes: one containing
    only the rows with 'gene' in the 'type' column, another containing only the rows with 'CDS' in the 'type' column,
    and the last one containing all the remaining rows that are not 'gene' or 'CDS' in the 'type' column. """
    return df[df['type'] == 'gene'], df[df['type'] == 'CDS'], df[
        (df['type'] != 'gene') & (df['type'] != 'CDS') & (df['type'] != 'source')]


def calculate_statistics(series):
    """Calculates the minimum, maximum, average, and standard deviation of a series."""
    return {
        "min": series.min(),
        "max": series.max(),
        "mean": series.mean(),
        "std": series.std()
    }


def draw_histogram(series, title, xlabel, ylabel):
    """Draws a histogram of a series."""
    plt.hist(series)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()


def find_richest_in_AT(df, n):
    """Finds the n genes with the highest AT composition in a dataframe."""
    df = df.sort_values(by="AT_composition", ascending=False)
    return df.head(n)


def find_poorest_in_AT(df, n):
    """Finds the n genes with the lowest AT composition in a dataframe."""
    df = df.sort_values(by="AT_composition")
    return df.head(n)


def find_cell_wall_genes(df, record_gb):
    """The find_cell_wall_genes() function takes in a Pandas dataframe df and a GenBank record record_gb, and adds a
    new column 'cell wall' to the dataframe which indicates if a gene is a cell wall gene or not. The function
    returns the modified dataframe. """
    description = {}
    for index, feature in enumerate(record_gb.features):
        if 'product' in feature.qualifiers and 'cell wall' in feature.qualifiers['product'][0]:
            description[index] = feature.qualifiers['product'][0]
        else:
            description[index] = ''
    df['cell wall'] = pd.Series(description)
    return df


def check_for_exceptions(df_CDS, record_gb):
    """The check_for_exceptions() function takes in a dataframe df_CDS and a GenBank record record_gb, and checks for
    exceptions such as if the length of the CDS is not a multiple of 3 or if the translation of the CDS is not the
    same as the translation in the genbank file. It returns a dictionary containing the name of the gene as key and
    the error message as value. """
    error_dic = {}

    def find_exceptions_by_strand(df_strand, strand=True):
        for index in range(df_strand.shape[0]):
            if df_strand.iloc[index]['length'] % 3 != 0:
                error_dic[df_strand.iloc[index]['name']] = 'Cant be divided into 3'
                continue
            start = df_strand.iloc[index]['start']
            end = df_strand.iloc[index]['end']
            seq = record_gb.seq[start:end]
            if not strand:
                seq = seq.reverse_complement()
            try:
                translation_seq = str(
                    seq.translate(table=int(record_gb.features[df_strand.index[index]].qualifiers['transl_table'][0]),
                                  cds=True))
                translation_form_file = record_gb.features[df_strand.index[index]].qualifiers['translation'][0]
                if translation_seq == translation_form_file:
                    continue
                else:
                    error_dic[df_strand.iloc[index]['name']] = 'dna seq is not compatible with the protein seq'
            except Exception as e:
                error_dic[df_strand.iloc[index]['name']] = e

    plus = df_CDS[df_CDS['strand'] == 1]
    negative = df_CDS[df_CDS['strand'] == -1]
    find_exceptions_by_strand(plus, strand=True)
    find_exceptions_by_strand(negative, strand=False)
    error_df = pd.DataFrame(error_dic.items())
    error_df.columns = ['Name Gene', 'Error message']
    error_df.to_csv('gene_exceptions.csv')
    return error_dic


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    _record_gb = parse_genbank("BS168.gb")
    df_all = create_dataframe(_record_gb)
    print("\n============================================= Section 1 ==============================================\n")
    print("Number of each type of region:")
    regions = count_elements(df_all)
    print(regions)
    print("\n============================================= Section 2 ==============================================\n")
    _df_gene, _df_CDS, _df_other = divide_genes(df_all)
    all_genes_coding_length_stats = calculate_statistics(_df_gene["length"])
    print("All genes:")
    print("Minimum length:", all_genes_coding_length_stats['min'])
    print("Maximum length:", all_genes_coding_length_stats['max'])
    print("Average length:", all_genes_coding_length_stats['mean'])
    print("Standard deviation:", all_genes_coding_length_stats['std'])
    print()
    protein_coding_length_stats = calculate_statistics(_df_CDS["length"])
    print("Protein-coding genes:")
    print("Minimum length:", protein_coding_length_stats['min'])
    print("Maximum length:", protein_coding_length_stats['max'])
    print("Average length:", protein_coding_length_stats['mean'])
    print("Standard deviation:", protein_coding_length_stats['std'])
    print()
    other_length_stats = calculate_statistics(_df_other["length"])
    print("Non-protein-coding genes:")
    print("Minimum length:", other_length_stats['min'])
    print("Maximum length:", other_length_stats['max'])
    print("Average length:", other_length_stats['mean'])
    print("Standard deviation:", other_length_stats['std'])
    draw_Multiple_histogram(_df_gene["length"], _df_CDS["length"], _df_other["length"], "histograms for length",
                            "Length (bp)",
                            "Number of Genes")
    print("\n============================================= Section 3 ==============================================\n")
    print('calculate AT composition')
    df_all = calculate_AT_composition_df(df_all, _record_gb.seq)
    _df_gene, _df_CDS, _df_other = divide_genes(df_all)
    genome_at_percentage = df_all[df_all['type'] == 'source']['AT_composition'].values[0]
    protein_coding_gene_at_percentages = calculate_statistics(_df_CDS["AT_composition"])
    non_protein_coding_gene_at_percentages = calculate_statistics(_df_other["AT_composition"])
    print("Genome AT percentage:", genome_at_percentage)
    print("Protein-coding gene AT percentage:", protein_coding_gene_at_percentages['mean'])
    print("Non-protein-coding gene AT percentage:", non_protein_coding_gene_at_percentages['mean'])
    if genome_at_percentage > protein_coding_gene_at_percentages["mean"]:
        print("The average AT composition for the protein-coding genes is lower than the average for the entire genome")
    else:
        print(
            "The average AT composition for the protein-coding genes is higher than the average for the entire genome")
    print()
    # Draw the histogram for the AT composition of the protein-coding genes
    draw_histogram(_df_CDS["AT_composition"], "Protein-Coding Genes (AT Composition)", "% AT",
                   "Number of Genes")
    print("The 5 top AT-rich genes are:")
    at_rich_genes = find_richest_in_AT(_df_gene, 5)
    at_rich_genes = at_rich_genes.drop(columns=['type', 'location'])
    print(at_rich_genes)
    print()
    print("The 5 top AT-poor genes are:")
    at_poor_genes = find_poorest_in_AT(_df_gene, 5)
    at_poor_genes = at_poor_genes.drop(columns=['type', 'location'])
    print(at_poor_genes)
    print("\n============================================= Section 4 ==============================================\n")
    df_all = find_cell_wall_genes(df_all, _record_gb)
    _df_wall_cell = df_all[df_all['cell wall'].str.contains('cell wall') == True]
    print(f"Number of genes containing 'cell wall' in the description: {_df_wall_cell.shape[0]}")
    print()
    wall_cell_length_stats = calculate_statistics(_df_wall_cell["length"])
    print("Gene lengths:")
    print("Minimum length:", wall_cell_length_stats['min'])
    print("Maximum length:", wall_cell_length_stats['max'])
    print("Average length:", wall_cell_length_stats['mean'])
    print("Standard deviation:", wall_cell_length_stats['std'])
    print()
    draw_histogram(_df_wall_cell["length"], "Genes with 'cell wall' in Description (Length)", "Length (bp)",
                   "Number of Genes")
    wall_cell_AT_stats = calculate_statistics(_df_wall_cell["AT_composition"])
    print("AT percentages:")
    print("Minimum:", wall_cell_AT_stats['min'])
    print("Maximum:", wall_cell_AT_stats['max'])
    print("Average:", wall_cell_AT_stats['mean'])
    print("Standard deviation:", wall_cell_AT_stats['std'])
    print()
    draw_histogram(_df_wall_cell["AT_composition"], "Genes with 'cell wall' in Description (AT Composition)", "% AT",
                   "Number of Genes")
    print("\n============================================= Section 5 ==============================================\n")
    _error_dic = check_for_exceptions(_df_CDS, _record_gb)
    if _error_dic:
        print()
        print(_error_dic)
    else:
        print()
        print("No exceptions found.")
    df_all = df_all.drop([0])
    df_all = df_all.drop(columns=['location'])
    df_all = df_all.sort_values(by=['start'], ascending=True)
    df_all.to_csv('./part_a.csv')
