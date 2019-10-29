#!/usr/bin/env pypy
from ete3 import Tree
from Bio import Phylo
from Bio import SeqIO
import random
import math
import operator
import numpy

nwk_file = input("Input your newick tree file: ")
algn_file = input("Input your alignment file: ")
sequences = SeqIO.parse(algn_file, "fasta")
sequences_list = list(sequences)
read_length = len(sequences_list[0])
t = Tree(nwk_file)
node = []
for seq in sequences_list:
    node.append(seq.id)
nodes = [(node,t)]

#prunes newick tree for a subtree with the given species IDs
def prune_tree(node):
    t = Tree(nwk_file)
    t.prune(node)
    return t

#measures the size of a tree by estimating the number of mutations
def tree_size(t):
    t.write(format=1, outfile = "small_tree.nwk")
    newT = Phylo.read("small_tree.nwk", "newick")
    size = newT.total_branch_length()
    return int(round(size*read_length))

#measures the size of smaller subtrees by their individual branch lengths
def mini_tree_size(t):
    t.write(format=1, outfile = "small_tree.nwk")
    newT = Phylo.read("small_tree.nwk", "newick")
    size = newT.total_branch_length()
    return size

def ecotype_test(nodes):
    nodes_for_splitting = []
    ecotypes = []
    for node in nodes:
        #any node with 2 or less species is automatically binned as an ecotype
        if len(node[0]) <= 2:
            ecotypes.append(node)
        else:
            num_reads = len(node[0])
            true_freq = frequencies_fasta(find_in_fasta(node[0]))
            true_ent = entropy(true_freq, num_reads)
            sim_results = run_sim(node, num_reads)
            #if the actual entropy exceeds the bottom 5% of the simulated entropies, the node will be binned as an ecotype
            #otherwise, the entropy will be tagged for further splitting
            if true_ent < sim_results[0] or true_ent > sim_results[-1]:
                nodes_for_splitting.append(node)
            else:
                count = 0
                threshhold = 1
                for result in sim_results:
                    if result <= true_ent:
                        count += 1
                    if count > threshhold:
                        ecotypes.append(node)
                        break
                if count < threshhold:
                    nodes_for_splitting.append(node)
    return (ecotypes, nodes_for_splitting)

def find_in_fasta(node):
    seq_for_decomp = []
    for sequence in sequences_list:
        if sequence.id in node:
            seq_for_decomp.append(sequence)
    return seq_for_decomp

#produces the results of 20 simulations on one node
def run_sim(node, num_reads):
    sim_results = []
    num_mutations = tree_size(node[1])
    for i in range(0,20):
        print("Running Simulation " + str(i+1) + "/20")
        start_matrix = fill_sim(num_reads)
        switched_nodes = switcher(start_matrix, node, num_reads, num_mutations)
        freqs = frequencies(switched_nodes)
        ent = entropy(freqs, num_reads)
        sim_results.append(ent)
    sorted_sim_results = sorted(sim_results)
    return sorted_sim_results

#calculates how many times the same sequence occurs in one simulated node

def frequencies(node):
    print("Gathering Frequencies")
    frequencies = []
    existing_reads = []
    new_reads = []
    for read in node:
        x = "".join(read)
        new_reads.append(x)
    for read in new_reads:
        if read not in existing_reads:
            freq = new_reads.count(read)
            frequencies.append((read, freq),)
            existing_reads.append(read)
    return frequencies

def frequencies_fasta(node):
    print("Gathering Frequencies")
    existing_reads = []
    frequencies = []
    new_reads = []
    for read in node:
        new_reads.append(read.seq)
    for read in new_reads:
        if read not in existing_reads:
            freq = new_reads.count(read)
            frequencies.append((read, freq),)
            existing_reads.append(read)
    sorted_freqs = sorted(frequencies, key=operator.itemgetter(1), reverse=True)
    return sorted_freqs

#if a node has n reads of length m, fill_sim creates n homogenous reads of length m
def fill_sim(num_reads):
    simulation_nodes = []
    for i in range(0, num_reads):
        inner_nodes = []
        for j in range(0, read_length):
            inner_nodes.append("0")
        simulation_nodes.append(inner_nodes)
    return simulation_nodes

#uses 2 random number generators to select the sequences and nucleotide sites which will be mutated in a simulation
def switcher(sim_nodes, node, num_reads, num_mutations):
    print("Mutating Sequences")
    print(num_mutations)
    x = random.sample(range(0, num_reads*num_mutations), num_mutations)
    y = random.sample(range(0, read_length*num_mutations), num_mutations)
    for i in range(0, num_mutations):
        sim_nodes[x[i]%num_reads][y[i]%read_length] = "1" if sim_nodes[x[i]%num_reads][y[i]%read_length] == "0" else "0"
        if i%99==0 or i == num_mutations-1:
            print("Current Number of Mutations: " + str(i+1))

    return sim_nodes

#calculates the Shannon entropy of one simulated node
def entropy(frequencies, num_reads):
    print("Calculating Entropy")
    if len(frequencies) == 0:
        return 0.0
    else:
        total = []
        for freq in frequencies:
            ind = (freq[1] * 1.0 / num_reads) + 0.0000000000000000001
            total.append(ind * math.log2(ind))
        return -(sum(total))

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
            t = splitting[i][1]
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
                        for clade in left.traverse():
                            if clade.name:
                                left2.append(clade.name)
                        for clade in right.traverse():
                            if clade.name:
                                right2.append(clade.name)
                        cleared.append(splitting[i])
                        add_to_split.append((left2,left),)
                        add_to_split.append((right2,right),)
                else:
                    children = []
                    #finds all the immediate descendants of t
                    for child in t.children:
                        children.append(child)
                    branch_lengths = []
                    #measures and sorts the descendants by their amount of mutation
                    for child in children:
                        child_size = mini_tree_size(child)/len(child)
                        branch_lengths.append((child, child_size),)
                    sorted_child_sizes = sorted(branch_lengths, key=operator.itemgetter(1), reverse=True)
                    new_out = sorted_child_sizes[0][0]
                    #the most mutated branch is chosen to be a new node
                    #the rest of the branches are grouped into a second node
                    #both nodes are then sent back to the ecotype_test function and retested
                    if not new_out.is_leaf():
                        new_node = []
                        for clade in new_out.traverse():
                            if clade.name:
                                new_node.append(clade.name)
                        add_to_split.append((new_node, new_out),)
                        for out in new_node:
                            splitting[i][0].remove(out)
                        new_tree = prune_tree(splitting[i][0])
                        add_to_split.append((splitting[i][0], new_tree),)
                        cleared.append(splitting[i])
                    else:
                        add_to_split.append(([new_out.name], new_out),)
                        splitting[i][0].remove(new_out.name)
                        new_tree = prune_tree(splitting[i][0])
                        add_to_split.append((splitting[i][0], new_tree),)
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

s = ecotype_test(nodes)
t = splitter(s)
print(t)

#parse new ecotypes into a text file
def to_log(ecotypes):
    new_ecotypes = open("new_ecotypes.txt", "w")
    count = 1
    ecotype_list = []

    for ecotype in ecotypes:
        ecotype_list.append(ecotype[0])

    for ecotype in ecotype_list:
        new_ecotypes.write("Ecotype " + str(count) + ": " + "[" + ", ".join(x for x in ecotype) + "]" + "\n")
        count += 1

to_log(t)
