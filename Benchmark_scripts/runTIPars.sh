#!/bin/sh

echo "Job Start  Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

taxonNameFile=~/yytao/16S/16S.names.txt
cd ~/yytao/16S/remove1
less $taxonNameFile | while read line
do
    query=$line
    queryfile=query.$query
	taxafile=taxa.$line
	ancFile=~/yytao/16S/16S_anc.fas
    treefile=remove1.nwk.$query
    outputfile=tipars.tree.$query
    ~/yytao/package/TIPars/tipars -t $treefile -s taxafile -a $ancFile -q $queryfile -o $outputfile
    echo $query" finished"
done

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

