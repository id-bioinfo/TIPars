#!/bin/sh

echo "Job Start  Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

taxonNameFile=~/yytao/16S/16S.names.txt
cd ~/yytao/16S/remove1
less $taxonNameFile | while read line
do
    query=$line
	queryfile=query.$query
	cat $queryfile taxa.$line > taxafile
	~/yytao/package/usher/build/faToVcf taxafile taxaVCF
    treefile=remove1.nwk.$query
    outputfile=usher.tree.$query
	~/yytao/package/usher/build/usher -t $treefile -v taxaVCF --write-uncondensed-final-tree -d ~/yytao/16S/remove1/usher
	mv ~/yytao/16S/remove1/usher/uncondensed-final-tree.nh $outputfile

    echo $query" finished"
done

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

