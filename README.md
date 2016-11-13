# TIPars

Taxa Insertion by PARSimony (TIPars) is Java program for inserting taxa into a reference phylogeny based on parsimony criterion. 
Ancestral sequences of all internal nodes are required as an input. Reference tree is maintained unchanged (except a new branch added) in the process. 
Length of the new branch(taxon) is calculated simply as P-distance (Can use JC69, K81, ... given the known substitution model, pending to be implemented). 
This program uses a specific and modified version of BEAST library (included in the repository; the official BEAST version cannot be used by ours).<br />
This is still undergoing development.

# Quick Usage
[Compilation]<br />
javac -classpath .:beast.jar TIPars.java <br />
<br />
[Run] <br />
java -classpath .:beast.jar TIPars <br />
... answer the following question interactively ... <br />
Enter your input nexus tree file: RefTree.nex <br />
Enter your input taxa seq file [fasta name is taxaname]: RefSeq.fas <br />
Enter your input ancestral seq file [fasta name is nid]: RefAncSeq.fas <br />
Enter your input query seq file [fasta name is taxaname]: SeqToBeInserted.fas <br />
Enter your output nexus tree file: inserted_RefTree.nex <br />
Want to output ABQdis, Bnid and genotype info? (0=no|1=yes): 0 <br />
<br />
Note: classpath separator varies between OS. <br />

# Acknowledgements
- This project is supported in part by Hong Kong Research Grants Council General Research Fund (17150816).
