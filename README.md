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
	 -q aligned_query_sequence \
	 -o output_file \
	 -f sequence_fileFormat (optional) \
	 -m bestplacement or not (optional) \
	 -p insertion or placement (optional)\
	 -d print to screen or not (optional) \
	 -x java Xmx setting (optional) \
	 
```
## demo
./tipars -t test/testdata/ref.tree -s test/testdata/taxa.fasta -a test/testdata/ancseq.fasta -q test/testdata/query.fasta -o test/testdata/tipars.tree


# Option Details

## input 

+ `-t`: tree file, in Newick format
+ `-s`: fasta/vcf file contains aligned taxa sequences
+ `-a`: fasta/vcf file contains aligned ancestral sequences
+ `-q`: fasta/vcf file contains one or multiple query seqence(s)
+ `-f`: sequences file format, one of 'fasta' and 'vcf', default(fasta)

## output

+ `-o`: output tree/jplace file name, default('TIPars_output.tree')
+ `-m`: choose multiplacement('true')(default) or multiplacement('false') for user notices
+ `-p`: algorithm type
  + `insertion` (default) for query sequence(s) insertion
    - output `newick` tree file with query sequence(s) inserted
    - mainly for updating tree
    - sequentially insertion
  + `placement` for query sequence(s) placement
    - output `jplace` tree file that incorporates original tree with placement information
    - mainly for query sequence(s) classification
    - independent placement
+ `-d`: print placement result to screen or not (default true)

# Acknowledgements

This project is supported by Hong Kong Research Grants Council General Research Fund (17150816).
