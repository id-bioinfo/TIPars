# Reconstrcut Ancestral Sequences using PastML

reconstructAncestralSeq.pl is a perl script to reconstrcut ancestral sequences using PastML
that includes two additional in-house scripts for adding names (INNODE1...INNODEN) for internal nodes of the tree 
and generating the annotation files as input for PastML.


# Authors

Marcus Shum and Yongtao Ye

# Dependency

To run reconstructAncestralSeq.pl, please first install [PastML](https://github.com/evolbioinfo/pastml) [^1] and [ETE3](http://etetoolkit.org/new_download/) [^2]
using Python3 as well as OpenMP for our in-house c++ script.

# How It Works 

reconstructAncestralSeq.pl includes three steps.

1) Generate an annotation table specifying tip states (nucleotide or amino acid) by extracting each column in the MSA of taxa, 
using the in-house c++ program `splitEachColumn` that uses OpenMP for parallelization.

2) Add names as INNODE1 to INNODEN (N is the number of internal nodes) for the input tree (should be rooted)
by the in-house python3 script `TREEMANUPULATION_AddInnodeNameToTreeByArgument.py` using ETE3.

3) Reconstrcut the ancestral sequences using PastML parallelly.

# Quick Usage

`perl reconstructAncestralSeq.pl <tree.nwk> <taxaMSA.fas> <outdir> <numberthreads>`

## input 
1) tree.nwk: Newick format tree file
  
2) taxaMSA.fas: fasta file contains aligned taxa sequences
  
3) outdir: output directory for reconstructed ancestral sequences and working space for PastML
  
4) numberthreads: an integer for number of threads to use

## output 
All output files are in the input path of <outdir>.
  
1) ancestral_sequence.fasta: ancestral sequences constructed by PastML in <outdir>
  
2) <tree.nwk>\_InnodeNameAdded: the tree file with added internal node names in the same directory of <tree.nwk>
  
3) <tree.nwk>\_ancestor: a tsv file indicating two ends of each branch and its branch length in the same directory of <tree.nwk>

## toy test

```bash
perl reconstructAncestralSeq.pl trial.tree trial.fasta outdir 4
```

Messages printed out to screen should be as follow.
 
PLESAE CONFIRMED THAT YOU HAVE INSTALLED:
1. PastML (https://pastml.pasteur.fr/)
2. ETE3 (http://etetoolkit.org/download/)
3. Python3, Perl5 (or above) and OpenMP
4. The in-house script for adding internal node names (TREEMANUPULATION)
5. The in-house script for generating statefiles for PastML (splitEachColumn)
Both above in-house scripts should be under the SAME directory.

Step1: GENERATING TABLE FILE FOR PASTML...

Processing used 0.00259685516357422 seconds

File Generation for PastML Completed!!!

Step2: ADDITION OF INNODE NAME TO TREE FILE...

Processing used 0.487257957458496 seconds

Addition of Name Completed!!!

Step3: START RUNNING PASTML......

Processing used 2.05268001556396 seconds

Finished Running PASTML!!!

Step4: START TO COMBINE PASTML RESULT...

Processing used 0.0245828628540039 seconds

Printing Out the Ancestral Sequence...

Finished All Process!!! Thank you for using!!!

# Reference
[^1]: Ishikawa, S.A., et al., A Fast Likelihood Method to Reconstruct and Visualize Ancestral Scenarios. Molecular Biology and Evolution, 2019. 36(9): p. 2069-2085.
[^2]: Huerta-Cepas J, Serra F, Bork P. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. Mol Biol Evol. 2016;33(6):1635-1638. 
