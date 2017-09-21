from Bio import Phylo
from ete3 import Tree
from io import StringIO
from numpy import floor
from collections import Counter
#todo: write as class object with methods to get sites, consensus sequence, etc.
###############################################################################################
#Find number of homoplasies that exist in the sequences. Run as follows:
#python homoplasies.py <sequences.fasta> <clonal_frame.nwk>
# written by Tom Brown retrieved 01/20/2017; added lines to midpoint root FastTree generated
# Newick tree; steps: 1. create binary versions of the sequences based on differences with con-
# sensus sequences; 2. create all possible non-homoplastic binary patterns of snp variation based
# on rooted tree; 3 verify if observed patterns are found in non-homoplastic pattern collection; if
# yes, score as non-homoplasy, no, score as homoplasy
###############################################################################################

#Read variable sites of sequence and create a binary object for each nucleotide if the current leaf
#has a variable site at this location

# sequence_file =open(argv[1],'r')
# sequences = []
# labels = []
# for line in sequence_file:
#     if line[0] != '>':
#         sequences.append(line)
#         sequences[-1] = sequences[-1][:-1]
#     else:
#         labels.append(line)
#         labels[-1] = labels[-1][1:-1]
#
#
#
# #Create a consensus sequence with the most present nucleotide at each site
# dummy_list = [None]*len(sequences)
# consensus_seq = []
# for nuc in range(len(sequences[0])):
#     for seq in range(len(sequences)):
#         dummy_list[seq] = sequences[seq][nuc]
#     consensus_seq.append(Counter(dummy_list).most_common(1)[0][0])
#
# print(consensus_seq)
#
# seq_binary = []
# for nuc in range(len(sequences[0])):
#     local_binary = []
#     count=0
#     for seq in range(len(sequences)):
#         if consensus_seq[nuc] == sequences[seq][nuc]:
#             local_binary.append(0)
#         else:
#             local_binary.append(1)
#     if sum(local_binary) > (len(sequences)/2):
#         for i in range(len(local_binary)):
#             local_binary[i] = 1 - local_binary[i]
#     seq_binary.append(local_binary)
#
# print(len(seq_binary))
#
# #Construct a list of binary objects where each node has a 1 or 0 representing whether a leaf is
# #a descendent of the current node
#
# #midpoint root tree with ETE 3
# treeETE = Tree(self.output_prefix)
# R = treeETE.get_midpoint_outgroup()
# treeETE.set_outgroup(R)
# print(treeETE)
# t = StringIO(treeETE.write())
# #original script
# tree = Phylo.read(t,"newick")
# internal_nodes = tree.get_nonterminals()
# binary_list = []
# for node in tree.find_clades():
#     print(node)
#     local_binary = []
#     for leaf in labels:
#         print(leaf)
#         if node in tree.get_path(next(tree.find_clades(leaf))):
#             local_binary.append(1)
#         else:
#             local_binary.append(0)
#
#     binary_list.append(local_binary)
# print(binary_list[3])
# print(len(binary_list[3]))
#
# #Include the opposite of each number to account for the other half of the tree at each node
# reverse_list = []
# for i in range(len(binary_list)):
#     local_binary = []
#     for j in range(len(binary_list[i])):
#         if (binary_list[i][j]):
#             local_binary.append(0)
#         else:
#             local_binary.append(1)
#     reverse_list.append(local_binary)
# for i in range(len(reverse_list)):
#     binary_list.append(reverse_list[i])
#
# print(reverse_list[3])
# del reverse_list
#
# non_homoplasies = []
# homoplasies = []
# for i in range(int(floor(len(sequences)/2))+1):
#     non_homoplasies.append(0)
#     homoplasies.append(0)
#
# for nuc in seq_binary:
#     if nuc in binary_list:
#         non_homoplasies[sum(nuc)] += 1
#     else:
#         homoplasies[sum(nuc)] += 1
# homoplasy_file = open("homoplasy.txt","w")
# non_homoplasy_file = open("non-homoplasy.txt","w")
# for i in range(len(homoplasies)):
#     homoplasy_file.write(str(homoplasies[i]))
#     homoplasy_file.write('\n')
#     non_homoplasy_file.write(str(non_homoplasies[i]))
#     non_homoplasy_file.write('\n')
# homoplasy_file.close()
# non_homoplasy_file.close()

def find_homoplasious_sites(alignment, newick_tree):
    # todo adjust to deal with missing data; done for now by excluding missing data
    # either delete missing states or impute most likely state....
    #implement deletion first...
    labels = []
    sequences = []
    for label, sequence in alignment:
        labels.append(label)
        sequences.append(sequence)
    # Create a consensus sequence with the most present nucleotide at each site
    dummy_list = [None] * len(sequences)
    consensus_seq = []
    for nuc in range(len(sequences[0])):
        for seq in range(len(sequences)):
            dummy_list[seq] = sequences[seq][nuc]
        cnt = Counter(dummy_list)
        del cnt['-']
        consensus_seq.append(cnt.most_common(1)[0][0])
    #Create binary sequences
    seq_binary = []
    for nuc in range(len(sequences[0])):
        local_binary = []
        count = 0
        for seq in range(len(sequences)):
            if consensus_seq[nuc].upper() == sequences[seq][nuc].upper():
                local_binary.append(0)
            elif sequences[seq][nuc] == '-':
                local_binary.append('-')
            else:
                local_binary.append(1)
        # if sum(local_binary) > (len(sequences) / 2):
        #     for i in range(len(local_binary)):
        #         local_binary[i] = 1 - local_binary[i]
        seq_binary.append(local_binary)

    # midpoint root tree with ETE 3
    treeETE = newick_tree
    R = treeETE.get_midpoint_outgroup()
    treeETE.set_outgroup(R)
    t = StringIO(treeETE.write())
    # original script
    tree = Phylo.read(t, "newick")
    internal_nodes = tree.get_nonterminals()
    binary_list = []
    for node in tree.find_clades():
        local_binary = []
        for leaf in labels:
            if node in tree.get_path(next(tree.find_clades(leaf))):
                local_binary.append(1)
            else:
                local_binary.append(0)

        binary_list.append(local_binary)

    # Include the opposite of each number to account for the other half of the tree at each node
    reverse_list = []
    for i in range(len(binary_list)):
        local_binary = []
        for j in range(len(binary_list[i])):
            if (binary_list[i][j]):
                local_binary.append(0)
            else:
                local_binary.append(1)
        reverse_list.append(local_binary)
    for i in range(len(reverse_list)):
        binary_list.append(reverse_list[i])
    del reverse_list
    #here the homoplasy site detection starts
    non_homoplasies =[]
    homoplasies = []
    for i  in range(len(seq_binary)):
        if '-' in seq_binary[i]:
            super_count = 0
            for j in binary_list:
                count = 0
                for index in range(len(seq_binary[i])):
                    if (seq_binary[i][index] == j[index]) or seq_binary[i][index] == '-':
                        continue
                    else: count += 1
                if count == 0:
                    non_homoplasies.append(i)
                    break
                else: super_count += 1
            if super_count == len(binary_list):
                homoplasies.append(i)
        else:
            if seq_binary[i] in binary_list:
                non_homoplasies.append(i)
            else:
                homoplasies.append(i)
    return non_homoplasies, homoplasies


