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
	~/yytao/package/epa-ng/bin/epa-ng --ref-msa $taxafile --tree $treefile --query $queryfile --model GTR -T 0
    echo $query" finished"
done

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

