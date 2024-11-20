# codes to generate multationSequenceMap, seqIdxMap, and ref_sequences

This branch is used for generating three serialized input, multationSequenceMap, seqIdxMap, and ref_sequences, so as to preprocess the taxa and ancestral sequences.

## Step to reproduce

1. Build the dev.Dockerfile as development environment container (in my case it is convert as singularity .sif file)
2. Prepare the taxa vcf and anc vcf, using faToVcf and mafft(if needed)
3. Build the .jar by running `make` in the development environment container created in step 1
4. Use the container in step 1 to run `./tipars`, which is a python program to start the TIPars.jar. Enter the right parameters (taxa, anc, output_dir)
5. three files are generated in the output_dir
