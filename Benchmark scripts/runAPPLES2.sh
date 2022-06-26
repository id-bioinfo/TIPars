#!/bin/sh

echo "Job Start  Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

taxonNameFile=~/yytao/16S/16S.names.txt
cd ~/yytao/16S/remove1
less $taxonNameFile | while read line
do
    query=$line
    queryfile=query.$query
	taxafile=taxa.$line
    treefile=remove1.nwk.$query
    outputfile=apples_jplace.$query
    python ~/yytao/package/apples/run_apples.py -D -s $taxafile -q $queryfile -t $treefile > $outputfile
    echo $query" finished"
done

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

