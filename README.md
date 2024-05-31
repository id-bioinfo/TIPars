# WebApp

TIPars is available online at [www.tipars.hku.hk](https://tipars.hku.hk/), with Influenza A/H5 reference tree ready for insertion [here](https://tipars.hku.hk/reference/Influenza-A-H5)

# Robust expansion of phylogeny for fast-growing genome sequence data

TIPars is JAVA program to do fast and accurate insertion or placement of new samples onto 
a reference phylogenetic tree based on parsimony criterion and utilized the pre-computed ancestral sequences. 
It is reliable for phylogenies comprise both densely sampled sequences (very short tree branches), e.g. SARS-CoV-2 genomes and Influenza viruses,
and divergent sequences (long tree branches), e.g. bacterial 16S ribosomal RNA sequences and Newcastle disease virus. 
It uses BEAST library [^1] and requires taxa sequences, the reference tree, 
and ancestral sequences of all its internal nodes as input data. 
Reference tree is maintained unchanged (except that a new branch added) in the insertion process. 
TIPars can do insertion of single or multiple new sequences with Newick format tree file as output, 
but can also obtain phylogenetic placement of the query sequences by generating Jplace format file [^2] for other downstream analyses. 

If there is any question about using TIPars, please send email to tipars@d24h.hk.

# Authors

Yongtao Ye, Marcus Shum, Joseph Tsui, Guangchuang Yu, Tommy Lam

# How It Works 

Given the multiple sequence alignments of taxa and ancestral sequences for an existing reference phylogenetic tree, 
TIPars computes the substitution scores of the query sequence against all branches in the tree using a specific substitution scoring table based on the IUPAC nucleotide ambiguity codes 
and searches for the minimal branch as the best insertion position. Length of new branches will be recalculated based on 
a local estimation model. An example of single query insertion is illustrated in Figure 1. 
Query Q(ACG**T**) differs from both two ends of branch G(ACC**G**)-D(ACG**C**) with one mutation at the fourth site,
resulting in the minimal insertion branch.
In case where multiple branches will result into the same minimal substitution score, 
TIPars applies simple yet practical rules to filter them. For details of the algorithm, 
please refer to our preprint of this work ([link](https://www.biorxiv.org/content/10.1101/2021.12.30.474610v1)).
For the required input data ancestral sequence reconstruction would be done using an in-house 
script with PastML [^3] ([link](https://github.com/id-bioinfo/TIPars/tree/master/reconstructAncestralSeq)).
Other methods such as ML joint or marginal methods (such as those available in HYPHY [^4]) 
are also acceptable. TIPars accepts both Fasta and Vcf file formats for input sequences. 
To convert a Fasta file to a Vcf file, it is suggested to use a Python package PoMo/FastaToVCF.py [^5] ([link](https://github.com/pomo-dev/PoMo/blob/master/scripts/FastaToVCF.py)).

<img src="/img/illustration.png" width="600">

# Installation

A precompiled executable program is available as TIPars.jar (required Java 11 or above). 

For users to compile TIPars source code from GitHub, 
```bash
git clone https://github.com/id-bioinfo/TIPars.git
cd TIPars
make
```

# Quick Usage

Given a reference tree (-t) and alignments of taxa and ancestral sequences (-s and -a),
TIPars would placement a set of **aligned** query samples (-q) to `jplace` placement file
or insert them to `newick` tree file according to user setting model (-p).

TIPars expects nucleotides by default, please use `-aa` for protein sequences that uses Blosum62 scoring matrix instead.

```bash
./tipars -aa (optional) \
         -t tree \
	 -s aligned_taxa_sequence \
         -a aligned_ancestral_sequence \
	 -q aligned_query_sequence \
	 -o output_file \
	 -f sequence_fileFormat (optional) \
	 -m multiplacement or not (optional) \
	 -p insertion or placement (optional)\
	 -d print to screen or not (optional) \
	 -x java Xmx setting (optional) \ 	 
```

## toy test

A small SARS-CoV-2 dataset with 1340 sequences is provided for a toy test. 
Due to GISAID's data sharing policy [^6], only Accession Numbers are provided for the sequences downloaded from GISAID (https://www.gisaid.org/).

If you just want to have a try on TIPars, regardness of SARS2, we recommend you to test on our NDV completed benchmark dataset at the folder `Benchmark_datasets/NDV` ([link](https://github.com/id-bioinfo/TIPars/tree/master/Benchmark_datasets/NDV)). 

Any problems about the usage of TIPars, please send email to tipars@d24h.hk.

To run on `test/sars2_1k`, please make sure you have downloaded the sequences from GISAID corresponding to the provided Accession Numbers!
```bash
./tipars -t test/sars2_1k/ref.tree -s test/sars2_1k/taxa.fasta -a test/sars2_1k/ancseq.fasta -q test/sars2_1k/query.fasta -o test/sars2_1k/tipars.tree
```

To run on `Benchmark_datasets/NDV`
```bash
cd Benchmark_datasets/NDV
../../tipars -t NDV_tree.nwk -s NDV_taxa.fas -a NDV_anc.fas -q NDV_query.fas -o tipars.tree
```

# Option Details

## input 

+ `-aa`: only set it when analyzing protein sequences
+ `-t`: tree file, in Newick format
+ `-s`: fasta/vcf file contains aligned taxa sequences
+ `-a`: fasta/vcf file contains aligned ancestral sequences
+ `-q`: fasta/vcf file contains one or multiple query seqence(s)
+ `-f`: sequences file format, one of 'fasta' and 'vcf', default (fasta)
+ `-x`: java Xmx setting, e.g.,1G,4G,8G, default (8G)

## output

+ `-o`: output tree/jplace file name, default ('TIPars_output.tree')
+ `-m`: choose bestplacement ('true') (default) or single best placement ('false') for user notices
+ `-p`: algorithm type
  + `insertion` (default) for query sequence(s) insertion
    - output `newick` tree file with query sequence(s) inserted
    - mainly for updating tree
    - sequentially insertion
  + `placement` for query sequence(s) placement
    - output `jplace` placement file that incorporates original tree with placement information
    - mainly for query sequence(s) classification
    - independent placement

# Benchmark datasets

The benchmark datasets used for this study can be referred on the folder `Benchmark_datasets`([link](https://github.com/id-bioinfo/TIPars/tree/master/Benchmark_datasets)),
including 16S, H3N2, NDV, SARS2-100k and SARS2-660k. Both the tree file and alignment files of taxa and ancestral sequences are available except SARS-CoV-2 datasets.
Due to GISAID's data sharing policy, only Accession Numbers are provided for the sequences downloaded from GISAID.
For the reference tree of SARS2-660k, please refer to the phylogeny (dated on 6 September 2021) under Audacity from GISAID.

# How to reconstruct ancestral sequences

We provided a perl script `reconstructAncestralSeq.pl` to reconstruct ancestral sequences using PastML[^3] parallelly. 
Input with a rooted tree and corresponding multiple sequence alignment of taxa,
the script ouputs the reconstructed ancestral sequences to fasta file and the tree with all internal node named as "INNODEXXX" to newick file. 
More details can be check in ([link](https://github.com/id-bioinfo/TIPars/tree/master/reconstructAncestralSeq)).

# Docker setup
We provided a Dockerfile for building Docker image, based on Ubuntu 22.04. The Dockerfile installed all nessesary software and libraries needed to run TIPars and ancestral sequence reconstruction using PastML (reconstructAncestralSeq). Here is how to use it:
1. Make sure you have Docker installed and running.
2. Set a shared directory in your host computer (shared with the docker container) and put all input files required to run TIPars or reconstructAncestralSeq to it. 
3. To run TIPars or reconstructAncestralSeq in any directory of your host computer.

+ Tipars

`sudo docker run --rm -v ${MY_PATH}:/home ghcr.io/id-bioinfo/tipars:1.1.1 /tipars/tipars -t /home/<tree file name> -s /home/<taxa file name> -a /home/<anc file name> -q /home/<query file name> -o /home/<output file name>`

Example (A toy test of NDV dataset in the Benchmark_datasets):
```bash
MY_PATH=/home/ytye/TIPars/Benchmark_datasets/NDV 
cd $MY_PATH
sudo docker run --rm -v $MY_PATH:/home ghcr.io/id-bioinfo/tipars:1.1.1 /tipars/tipars -t /home/NDV_tree.nwk -s /home/NDV_taxa.fas -a /home/NDV_anc.fas -q /home/NDV_query.fas -o /home/tipars.tree
```
+ reconstructAncestralSeq

create a folder ${outdir} at the shared directory to store the ancestral sequecnes, and then run reconstructAncestralSeq `sudo docker run --rm -v ${MY_PATH}:/home -w /tipars/reconstructAncestralSeq ghcr.io/id-bioinfo/tipars:1.1.1 perl reconstructAncestralSeq.pl /home/<tree file name> /home/<taxa file name> /home/${outdir} <number of parallel processes>`

Example (a small trial data in the reconstructAncestralSeq directory):
```bash
MY_PATH=/home/ytye/TIPars/reconstructAncestralSeq/
cd $MY_PATH && mkdir outdir
sudo docker run --rm -v $MY_PATH:/home -w /tipars/reconstructAncestralSeq ghcr.io/id-bioinfo/tipars:1.1.1 perl reconstructAncestralSeq.pl /home/trial.tree /home/trial.fasta /home/outdir 4
``` 

${MY_PATH} is the absolute path of shared directory created in step 4.

4. The output will be in the shared directory ${MY_PATH}.

# How to Cite

Ye Y, Shum MH, Tsui JL, Yu G, Smith DK, et al. (2024) Robust expansion of phylogeny for fast-growing genome sequence data. PLOS Computational Biology 20(2): e1011871. https://doi.org/10.1371/journal.pcbi.1011871

# Acknowledgements

This project is supported by the Hong Kong Research Grants Council General Research Fund (17150816), the NSFC Excellent Young Scientists Fund (Hong Kong and Macau) (31922087),
the Health and Medical Research Fund (COVID1903011-549 WP1) and the Innovation and Technology Commission’s InnoHK funding (D<sup>2</sup>4H).

# Reference
[^1]: Suchard, M.A., et al. Bayesian phylogenetic and phylodynamic data integration using BEAST 1.10. Virus evolution, 2018.
[^2]: Matsen, F.A., R.B. Kodner, and E.V. Armbrust, pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics, 2010. 11(1): p. 538.
[^3]: Ishikawa, S.A., et al., A Fast Likelihood Method to Reconstruct and Visualize Ancestral Scenarios. Molecular Biology and Evolution, 2019. 36(9): p. 2069-2085.
[^4]: Kosakovsky Pond, S.L., et al., HyPhy 2.5-A Customizable Platform for Evolutionary Hypothesis Testing Using Phylogenies. Molecular Biology and Evolution, 2020. 37(1): p. 295-299.
[^5]: Schrempf, D., et al., Reversible polymorphism-aware phylogenetic models and their application to tree inference. J Theor Biol, 2016. 407: p. 362-370.
[^6]: Shu Y, McCauley J. GISAID: Global initiative on sharing all influenza data - from vision to reality. Euro Surveill. 2017. 22(13): p. 30494.
