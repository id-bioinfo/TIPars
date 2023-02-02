#!/bin/sh

echo "Job Start  Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

taxonNameFile=~/yytao/16S/16S.names.txt
refTreefile=~/yytao/16S/16S_besttree500-labelled-T2.tree
cd ~/yytao/package/TreeCmp
less $taxonNameFile | while read line
do
    treefile=~/yytao/16S/remove1/tipars.tree.$query
	outputfile=~/yytao/16S/remove1/rf_tipars.tree.$query
	java -jar bin/TreeCmp.jar -d rf -i $treefile -r $refTreefile -o $outputfile
    echo $query" finished"
done

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

