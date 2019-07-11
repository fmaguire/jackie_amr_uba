from ete3 import Tree
from ete3 import TreeStyle

import sys

def get_arguments():
    return open(sys.argv[1])

def read_arguments_as_text(arguments):
    return arguments.readlines()

def car(array):
    return array[0]

#def cdr(array):
#    return array[0]

def ete_treeify(newick_string):
    return Tree(newick_string, format=0)

tree = (ete_treeify(car(read_arguments_as_text(get_arguments()))))

print(car(read_arguments_as_text(get_arguments())))

ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
ts.complete_branch_lines_when_necessary = True 
tree.render("mytree2.png", w=2048, units="px", tree_style=ts )
