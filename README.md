# TIPars

Taxa Insertion by Parsimony (`TIPars`) is Java program for inserting taxa into a reference phylogeny based on parsimony criterion.
Ancestral sequences of all internal nodes are required as an input. Reference tree is maintained unchanged (except a new branch added) in the process.
Length of the new branch(taxon) is calculated simply as P-distance (Can use JC69, K81, ... given the known substitution model, pending to be implemented).
This program uses a modified version of BEAST library (included in the repository).


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
	 -o output_file \
```

# Option Details

## input files

+ `-t`: Provide the file name of reference tree for insertion, in newick format.
+ `-s`: Provide the file name of alignment of the taxa sequences; The names of sequences should match the taxa names in the reference tree. 
+ `-a`: Provide the file name of alignment of the ancestral sequences; The names of sequences should match the node labels in the reference tree. 
+ `-q`: Provide the file name of alignment of query sequence(s). The nucleotides should be in aligned positions with the other input alignments. 

## model option

+ `-m`: Choose the distance model for estimating the branch length (pendant length) of the inserted query sequence
  + JC69 (Juke-Cantor 1969)
  + K2P (Kimura two-parameter)
  + LE (Local estimation using branches of the triplet tree, default)
+ `-g`: Choose gap treatment method:
  + `ignore`: ignore all gaps (default)
  + `inner`: only count inner gaps
  + `all`: count all gaps as sequence characters

## output file

+ `-o`: output file name (`TIPars_output.tree` by default)
+ `-p`: algorithm type
  + `insertion` (default) for query sequence(s) insertion
    - output `newick` tree file with query sequence(s) inserted
    - mainly for updating tree
    - progressively insertions
  + `placement` for query sequence(s) placement
    - output `jplace` tree file that incorporates original tree with placement information
    - mainly for query sequence(s) classification
    - independent placement

# Acknowledgements

This project is supported by Hong Kong Research Grants Council General Research Fund (17150816).
