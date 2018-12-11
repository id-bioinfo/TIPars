# TIPars

Taxa Insertion by Parsimony (`TIPars`) is Java program for inserting taxa into a reference phylogeny based on parsimony criterion.
Ancestral sequences of all internal nodes are required as an input. Reference tree is maintained unchanged (except a new branch added) in the process.
Length of the new branch(taxon) is calculated simply as P-distance (Can use JC69, K81, ... given the known substitution model, pending to be implemented).
This program uses a specific and modified version of BEAST library
(included in the repository; the official BEAST version cannot be used
by ours).


This is still undergoing development.

# Authors


Tommy Tsan-Yuk LAM and Guangchuang YU


# Quick Usage

```bash
./tipars -t tree \
	 -s aligned_taxa_sequence \
         -a aligned_ancestral_sequence \
	 -q query_sequence \
	 -m model \
	 -g gap \
	 -p type \	 
	 -o output_file
```

# Option Details

(input files)
-t: Provide the name of reference tree file for insertion
-s: Provide the name of alignment file of the taxa sequences; Their names should match the taxa names in the reference tree)
-a: Provide the name of alignment file of the ancestral sequences; Their names should match the node labels in the reference tree)
-q: Provide the name of alignment file of query sequence(s). Their nucleotides should be in aligned positions with the other input alignments.

(model option)
-m: Choose the distance model for estimating the branch length (pendant length) of the inserted query sequence. 
-g: Choose gap treatment method: (1) as an extra character or just 


# Acknowledgements

- This project is supported by Hong Kong Research Grants Council General Research Fund (17150816).
