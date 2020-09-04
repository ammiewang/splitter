#!/usr/bin/env pypy
from ete3 import Tree
from Bio import Phylo, SeqIO
import multiprocessing as mp
from consts import THRESHHOLD
from nodes import SplitterNode
from tree_funcs import *

def ecotype_test(nodes):
    print("Ecotype Test")
    nodes_for_splitting = []
    ecotypes = []
    for node in nodes:
        #any node with 2 or less species is automatically binned as an ecotype
        if len(node.seqs) <= 2:
            ecotypes.append(node)
        else:
            node.find_in_fasta(sequences_list)
            true_freq = node.frequencies_fasta()
            true_ent = node.entropy(true_freq)
            #node.true_ent = true_ent
            print('Actual Entropy', true_ent)
            if true_ent == 0.0:
                ecotypes.append(node)
            else:
                sim_results = node.check()
                #node.sim_results = sim_results
                #if the actual entropy exceeds the bottom 5% of the simulated entropies, the node will be binned as an ecotype
                #otherwise, the entropy will be tagged for further splitting
                if true_ent < sim_results[THRESHHOLD] or true_ent > sim_results[-1]:
                    nodes_for_splitting.append(node)
                else:
                    ecotypes.append(node)
    return (ecotypes, nodes_for_splitting)


#splits each designated node into 2 sub-nodes for further testing
def splitter(nodes):
    print("Splitting Ecotypes")
    ecotypes = nodes[0]
    splitting = nodes[1]
    while splitting != []:
        i = 0
        k = len(splitting)
        cleared = []
        add_to_split = []
        while i < k:
            t = splitting[i].tree
            if len(t.children) > 1:
                a = t.children
                #if t only has 2 immediate desendants, we split the tree at the root node into two new nodes
                #these nodes will then be sent back to the ecotype_test function
                if len(a) == 2:
                    left = a[0]
                    right = a[1]
                    left2 = []
                    right2 = []
                    if not left.is_leaf() or not right.is_leaf():
                        left2 = get_names(left)
                        right2 = get_names(right)
                        cleared.append(splitting[i])
                        add_to_split.append(SplitterNode(left2, left, read_length))
                        add_to_split.append(SplitterNode(right2, right, read_length))
                else:
                    branch_lengths = []
                    #measures and sorts the descendants by their amount of mutation
                    for child in t.children:
                        child_size = mini_tree_size(child)/len(child)
                        branch_lengths.append((child, child_size),)
                    #sorted_child_sizes = sorted(branch_lengths, key=lambda x: x[1], reverse=True)
                    new_out = min(branch_lengths, key=lambda x: x[1])
                    #the most mutated branch is chosen to be a new node
                    #the rest of the branches are grouped into a second node
                    #both nodes are then sent back to the ecotype_test function and retested
                    if not new_out.is_leaf():
                        new_node = []
                        new_node = get_names(new_out)
                        add_to_split.append(SplitterNode(new_node, new_out, read_length))
                        for out in new_node:
                            splitting[i].seqs.remove(out)
                        new_tree = prune_tree(splitting[i].seqs, nwk_file)
                        add_to_split.append(SplitterNode(splitting[i].seqs, new_tree, read_length))
                        cleared.append(splitting[i])
                    else:
                        add_to_split.append(SplitterNode([new_out.name], new_out, read_length))
                        splitting[i].seqs.remove(new_out.name)
                        new_tree = prune_tree(splitting[i].seqs)
                        add_to_split.append(SplitterNode(splitting[i].seqs, new_tree, read_length))
                        cleared.append(splitting[i])
            i += 1
        for read in cleared:
            splitting.remove(read)
        for read in add_to_split:
            splitting.append(read)
        split_further = ecotype_test(splitting)
        for read in split_further[0]:
            ecotypes.append(read)
            splitting.remove(read)
    return ecotypes

    #parse new ecotypes into a text file
def to_log(ecotypes):
    new_ecotypes = open(out_file, "w")
    count = 1

    for ecotype in ecotypes:
        new_ecotypes.write("Ecotype " + str(count) + ": " + "[" + ", ".join(x for x in ecotype.seqs) + "]" + "\n")
        count += 1


if __name__ == '__main__':
    nwk_file = input("Input your newick tree file: ")
    algn_file = input("Input your alignment file: ")
    out_file = input("Input the desired name of your output file: ")
    sequences = SeqIO.parse(algn_file, "fasta")
    sequences_list = list(sequences)
    read_length = len(sequences_list[0])
    t = Tree(nwk_file)
    node = [seq.id for seq in sequences_list]
    nodes = [SplitterNode(node, t, read_length)]

    s = ecotype_test(nodes)
    t = splitter(s)
    for ecotype in t:
        print(ecotype.seqs)
    to_log(t)
