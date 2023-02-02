#!/bin/sh

echo "Job Start  Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

taxonNameFile=~/yytao/16S/16S.names.txt
cd ~/yytao/16S/remove1
less $taxonNameFile | while read line
do
    query=$line
	queryfile=query.$query
	cat $queryfile taxa.$line > taxafile
    treefile=remove1.nwk.$query
	~/yytao/package/iqtree-2.1.3/bin/iqtree2 -g $treefile -s taxafile -n 0 -m GTR -fixbr -pre $query
    echo $query" finished"
done

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

