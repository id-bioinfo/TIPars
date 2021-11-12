# Reconstrcut Ancestral Sequences using PastML

reconstrcutAncestralSeq.pl is an in-house perl script to reconstrcut ancestral sequences using PastML.

# Authors

Marcus Shum

# Dependency

To run reconstrcutAncestralSeq.pl, please first install [PastML](https://github.com/evolbioinfo/pastml) [^1] and [ETE3](http://etetoolkit.org/new_download/) [^2]
using Python2.

# Quick Usage

perl reconstrcutAncestralSeq.pl

## toy test

```bash
perl reconstrcutAncestralSeq.pl
```
PLESAE CONFIRMED THAT YOU HAVE INSTALLED:
1. PastML (https://pastml.pasteur.fr/)
2. ETE3 (http://etetoolkit.org/download/)
3. Python2 (NOT!!3!!!)
4. The In-house python script of the TREEMANUPULATION is under the SAME directory

Beforing using this script
If you have confirmed that you have downloaded the above programs, please click "Enter" to proceed...
```bash
press 'Enter' button
```
Please type in the name of/full path to the tree file(in newick format):
```bash
trial.tree
```
Please type in the name of/full path to the sequence file(in single-lined fasta format):
```bash
trial.fasta
```
Step1: READING IN FASTA FILE...

FINISH READING IN FASTA FILE!!!

Step2: NOW GENERATING TABLE FILE FOR PASTML...

position: 1/8
position: 2/8
position: 3/8
position: 4/8
position: 5/8
position: 6/8
position: 7/8
position: 8/8
FILE GENERATION FOR PASTML COMPLETED!!!

Step3: ADDITION OF INNODE NAME TO TREE FILE...


Tree Innode Name Addtion COMPLETED
Addition of Name Completed!!!

Step4: START RUNNING PASTML......

Finished Running PASTML!!!

Step5: START TO COMBINE PASTML RESULT...

1
2
3
4
5
6
7
8
Printing Out the Ancestral Sequence...

Finished All Process!!! Thank you for using!!!

# Reference
[^1]:
Ishikawa, S.A., et al., A Fast Likelihood Method to Reconstruct and Visualize Ancestral Scenarios. Molecular Biology and Evolution, 2019. 36(9): p. 2069-2085.
[^2]:
Huerta-Cepas J, Serra F, Bork P. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. Mol Biol Evol. 2016;33(6):1635-1638. 