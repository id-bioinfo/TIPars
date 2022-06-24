#!/usr/bin/env python

#####Copyright Owned by Marcus Shum, 01/01/2019######

import os
import argparse
import sys
import ete3
from ete3 import Tree
import datetime

def main():
    i=0
    print("Tree input for Innode Name Addition is: "+str(sys.argv[1]))
    t = Tree(sys.argv[1], format=1)
    ##Adding node name to internal node
    print("Innode Name Addition In Progress......")
    for node in t.traverse("preorder"):
        if node.is_root() != 1:
            if not node.is_leaf():
            #if node.name == "":
                name = "INNODE"+str(i)
                i=i+1
                node.name=name
        else:
            node.name = "ROOT"
    file_name = sys.argv[1]+"_InnodeNameAdded"
    tree_file = open(file_name, "w")			
    tree_file.write(t.write(format=1)[:-1]+"ROOT:0;")
    tree_file.close()
    print("Innode Name Added Tree created as "+"\""+str(file_name)+"\"")
    ##print ancestor list 
    ancestor_file = sys.argv[1]+"_ancestor"
    ancestor_outfile=open(ancestor_file, 'w')
    today = datetime.date.today()
    ancestor_outfile.writelines("# created on " + today.ctime() + "\n")
    ancestor_outfile.writelines("name\tparent\n")
    for node in t.traverse("postorder"):
        if node.is_root() != 1:
                parent = node.up
                ancestor_outfile.writelines(node.name + "\t" + parent.name + "\t" + str(node.dist) + "\n")
    ancestor_outfile.close()
    print("Ancestor list created as "+"\""+str(ancestor_file)+"\"")
    sys.stderr.write("\nTree Innode Name Addtion COMPLETED\n")

main()


