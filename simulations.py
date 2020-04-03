#!/usr/bin/env pypy
from ete3 import Tree
from Bio import Phylo
from Bio import SeqIO
import random
import math
import operator
import numpy
import multiprocessing as mp


def run_sim(node, tree, num_reads, read_length, num_mutations):
    #print("Running Simulation")
    start_matrix = fill_sim(num_reads, read_length)
    switched_nodes = switcher(start_matrix, node, num_reads, num_mutations, read_length)
    freqs = frequencies(switched_nodes)
    ent = entropy(freqs, num_reads)
    return ent

#calculates how many times the same sequence occurs in one simulated node

def frequencies(node):
    #print("Gathering Frequencies")
    new_reads = []
    for read in node:
        x = "".join(read)
        new_reads.append(x)
    freq_dict = {}
    for read in new_reads:
        if read not in freq_dict:
            freq_dict[read] = 1
        else:
            freq_dict[read] += 1
    return freq_dict

#if a node has n reads of length m, fill_sim creates n homogenous reads of length m
def fill_sim(num_reads, read_length):
    simulation_nodes = []
    for i in range(0, num_reads):
        inner_nodes = []
        for j in range(0, read_length):
            inner_nodes.append("0")
        simulation_nodes.append(inner_nodes)
    return simulation_nodes

#uses 2 random number generators to select the sequences and nucleotide sites which will be mutated in a simulation
def switcher(sim_nodes, node, num_reads, num_mutations, read_length):
    #print("Mutating Sequences")
    x = random.sample(range(0, num_reads*num_mutations), num_mutations)
    y = random.sample(range(0, read_length*num_mutations), num_mutations)
    for i in range(0, num_mutations):
        sim_nodes[x[i]%num_reads][y[i]%read_length] = "1" if sim_nodes[x[i]%num_reads][y[i]%read_length] == "0" else "0"

    return sim_nodes

#calculates the Shannon entropy of one simulated node
def entropy(frequencies, num_reads):
    #print("Calculating Entropy")
    if len(frequencies) == 1:
        return 0.0
    else:
        total = []
        for freq in frequencies:
            ind = (frequencies[freq] * 1.0 / num_reads)
            total.append(ind * math.log2(ind))
        return -(sum(total))

def check(node, tree, num_reads, read_length, num_mutations):
    pool = mp.Pool(mp.cpu_count())
    result_objects = [pool.apply_async(run_sim, args=(node, tree, num_reads, read_length, num_mutations)) for i in range(20)]
    results = [r.get() for r in result_objects]
    sorted_results = sorted(results)
    pool.close()
    print("Simulation Results: ", sorted_results)
    return sorted_results
