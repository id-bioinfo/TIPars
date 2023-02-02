#!/bin/sh

echo "Job Start  Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

taxonNameFile=~/yytao/16S/16S.names.txt
cd ~/yytao/16S/remove1
less $taxonNameFile | while read line
do
    query=$line
	queryfile=query.$query
	taxafile=taxa.$line
    treefile=remove1_binary.nwk.$query
	~/yytao/package/pagan2-msa/bin/pagan2 --ref-treefile $treefile --ref-seqfile $taxafile --queryfile $queryfile --guidetree --one-placement-only
    echo $query" finished"
done

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

