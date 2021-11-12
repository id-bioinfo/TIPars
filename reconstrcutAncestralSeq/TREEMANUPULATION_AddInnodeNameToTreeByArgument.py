#!/usr/bin/env python

#####Copyright Owned by Marcus Shum, 01/01/2019######

import os
import argparse
import sys
import ete3
from ete3 import Tree

def main():
	i=0
	print "Tree input for Innode Name Addition is: "+str(sys.argv[1])
	t = Tree(sys.argv[1], format=1)
	##Adding node name to internal node
	print "Innode Name Addition In Progress......"
	for node in t.traverse("preorder"):
		if node.is_root() != 1:
			if node.name == "":
				name = "INNODE"+str(i)
				i=i+1
				node.name=name
	file_name = sys.argv[1]+"_InnodeNameAdded"
	tree_file = open(file_name, "w")			
	tree_file.write(t.write(format=1)[:-1]+"ROOT:0;")
	print "Innode Name Added Tree created as "+"\""+str(file_name)+"\""
	sys.stderr.write("\nTree Innode Name Addtion COMPLETED\n")

main()


