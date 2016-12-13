# TIPars

Taxa Insertion by Parsimony (`TIPars`) is Java program for inserting taxa into a reference phylogeny based on parsimony criterion.
Ancestral sequences of all internal nodes are required as an input. Reference tree is maintained unchanged (except a new branch added) in the process.
Length of the new branch(taxon) is calculated simply as P-distance (Can use JC69, K81, ... given the known substitution model, pending to be implemented).
This program uses a specific and modified version of BEAST library
(included in the repository; the official BEAST version cannot be used
by ours).


This is still undergoing development.

# Authors


Tommy Tsan-Yuk LAM and Guangchuang Yu


# Quick Usage

```bash
./tipars -t tree \
	     -s aligned_taxa_sequence \
         -a aligned_ancestral_sequence \
		 -q query_sequence \
		 -m model \
		 -o output_file
```

# Acknowledgements

- This project is supported in part by Hong Kong Research Grants Council General Research Fund (17150816).
