# TIPars

Taxa Insertion by Parsimony (TIPars) is Java program for inserting taxa into a reference phylogeny based on parsimony criterion. 
Ancestral sequences of all internal nodes are required as an input. Reference tree is maintained unchanged (except a new branch added) in the process. 
Length of the new branch(taxon) is calculated simply as P-distance (Can use JC69, K81, ... given the known substitution model, pending to be implemented). 
This program uses a specific and modified set of BEAST library (included in the repository; the official BEAST version cannot be used by ours).
This is still undergoing development.

# Quick Usage
[Compilation]
javac -classpath .:beast.jar TIPars.java

[Run]\n
java -classpath .:beast.jar TIPars\n
... answer the following question interactively ...
Enter your input nexus tree file: RefTree.nex
Enter your input taxa seq file [fasta name is taxaname]: RefSeq.fas
Enter your input ancestral seq file [fasta name is nid]: RefAncSeq.fas
Enter your input query seq file [fasta name is taxaname]: SeqToBeInserted.fas
Enter your output nexus tree file: inserted_RefTree.nex
Want to output ABQdis, Bnid and genotype info? (0=no|1=yes): 0

Note: classpath separator varies between OS.

# Acknowledgements
- This project is supported in part by Hong Kong Research Grants Council General Research Fund (17150816).
