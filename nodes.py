from ete3 import Tree
from Bio import Phylo, SeqIO
import random
import math
import multiprocessing as mp
from consts import SIMULATIONS
from tree_funcs import tree_size

class SplitterNode():
    def __init__(self, seqs, tree, read_length):
        self.seqs = seqs
        self.tree = tree
        self.read_length = read_length
        self.num_reads = len(self.seqs)
        self.num_mutations = tree_size(self.tree, self.read_length)

    def find_in_fasta(self, sequences_list):
        self.seq_for_decomp = []
        seq_set = set(self.seqs)
        for sequence in sequences_list:
            if sequence.id in seq_set:
                self.seq_for_decomp.append(sequence)
        #return self.seq_for_decomp

    def frequencies_fasta(self):
        print("Gathering Frequencies")
        new_reads = []
        for read in self.seq_for_decomp:
            x = str(read.seq)
            new_reads.append(x)
        freq_dict = {}
        for read in new_reads:
            if read not in freq_dict:
                freq_dict[read] = 1
            else:
                freq_dict[read] += 1
        return freq_dict

    def run_sim(self):
        print("Running Simulation")
        start_matrix = self.fill_sim()
        mutated = self.switcher(start_matrix)
        freqs = self.frequencies(mutated)
        ent = self.entropy(freqs)
        return ent

    #calculates how many times the same sequence occurs in one simulated node

    def frequencies(self, sim_nodes):
        print("Gathering Frequencies")
        new_reads = []
        for read in sim_nodes:
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
    def fill_sim(self):
        return [['0' for _ in range(self.read_length)] for _ in range(self.num_reads)]

    #uses 2 random number generators to select the sequences and nucleotide sites which will be mutated in a simulation
    def switcher(self, sim_nodes):
        print("Mutating Sequences")
        #print(num_mutations)
        x = random.sample(range(self.num_reads*self.num_mutations), self.num_mutations)
        y = random.sample(range(self.read_length*self.num_mutations), self.num_mutations)
        for i in range(self.num_mutations):
            sim_nodes[x[i]%self.num_reads][y[i]%self.read_length] = "1" if sim_nodes[x[i]%self.num_reads][y[i]%self.read_length] == "0" else "0"
            if i == self.num_mutations-1:
                print("Current Number of Mutations: " + str(i))

        return sim_nodes

    #calculates the Shannon entropy of one simulated node
    def entropy(self, frequencies):
        print("Calculating Entropy")
        if len(frequencies) == 1:
            return 0.0
        else:
            total = []
            for freq in frequencies:
                ind = (frequencies[freq] * 1.0 / self.num_reads)
                total.append(ind * math.log2(ind))
            return -(sum(total))

    def check(self):
        pool = mp.Pool(mp.cpu_count())
        result_objects = [pool.apply_async(self.run_sim) for i in range(SIMULATIONS)]
        results = [r.get() for r in result_objects]
        self.sorted_results = sorted(results)
        pool.close()
        print("Results: ", self.sorted_results)
        return self.sorted_results
