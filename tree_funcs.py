from ete3 import Tree
from Bio import Phylo

def tree_size(t, read_length):
      t.write(format=1, outfile = "small_tree.nwk")
      newT = Phylo.read("small_tree.nwk", "newick")
      size = newT.total_branch_length()
      return round(size*read_length)

#prunes newick tree for a subtree with the given species IDs
def prune_tree(node, nwk_file):
  t = Tree(nwk_file)
  t.prune(node)
  return t

#measures the size of smaller subtrees by their individual branch lengths
def mini_tree_size(t):
  t.write(format=1, outfile = "small_tree.nwk")
  newT = Phylo.read("small_tree.nwk", "newick")
  size = newT.total_branch_length()
  return size

def get_names(t):
  names = [] #set({})
  for clade in t.traverse():
      if clade.name:
          names.append(clade.name)#names.add(clade.name)
  return names
