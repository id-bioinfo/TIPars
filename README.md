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
a local estimation model. An example of a single query sequence is illustrated in Figure 1. 
In case where multiple branches will result into the same minimal substitution score, 
TIPars applies simple yet practical rules to filter them. For details of the algorithm, 
please refer to our preprint of this work ([link](https://www.biorxiv.org/content/10.1101/2021.12.30.474610v1)).
For the required input data ancestral sequence reconstruction would be done using an in-house 
script with PastML [^3] ([link](https://github.com/id-bioinfo/TIPars/tree/master/reconstrcutAncestralSeq)).
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

```bash
./tipars -t tree \
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

A H3N2 dataset with 800 HA genes is provided for a toy test. 
Any problems about the usage of TIPars, please send email to tipars@d24h.hk.

```bash
./tipars -t Benchmarkdatasets/H3N2/H3N2_tree.nwk -s Benchmarkdatasets/H3N2/H3N2_taxa.fas -a Benchmarkdatasets/H3N2/H3N2_anc.fas -q Benchmarkdatasets/H3N2/H3N2_query.fas -o Benchmarkdatasets/H3N2/tipars.tree
```

# Option Details

## input 

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
    - output `jplace` tree file that incorporates original tree with placement information
    - mainly for query sequence(s) classification
    - independent placement

# Benchmark datasets

The benchmark datasets used for this study can be referred on the folder `Benchmark datasets`([link](https://github.com/id-bioinfo/TIPars/tree/master/Benchmark%20datasets)),
including 16S, H3N2, NDV, SARS2-100k and SARS2-660k. Both the tree file and alignment files of taxa and ancestral sequences are available except SARS-CoV-2 datasets.
Due to GISAID's data sharing policy, only Accession Numbers are provided for the sequences downloaded from GISAID.
For the reference tree of SARS2-660k, please refer to the phylogeny (dated on 6 September 2021) under Audacity from GISAID.

# How to Cite

Yongtao Ye, Marcus Shum, Joseph Tsui, Guangchuang Yu, David Smith, Huachen Zhu, Joseph Wu, Yi Guan, Tommy Tsan-Yuk Lam. Robust expansion of phylogeny for fast-growing genome sequence data.
bioRxiv 2021.12.30.474610; doi: https://doi.org/10.1101/2021.12.30.474610

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
