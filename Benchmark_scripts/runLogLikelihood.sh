#!/bin/sh

echo "Job Start  Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

taxonNameFile=~/yytao/16S/16S.names.txt
cd ~/yytao/16S/remove1
less $taxonNameFile | while read line
do
    query=$line
	queryfile=query.$query
	cat $queryfile taxa.$line > taxafile
    treefile=tipars.binary_tree.$query
	treefile_opt=tipars.binary_tree.opt.$query
	logfile=tipars.log.$query
	~/yytao/package/fasttree/FastTreeDbl -gamma -nt -gtr -nome -mllen -intree $treefile -log $logfile taxafile > $treefile_opt
    echo $query" finished"
done

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

