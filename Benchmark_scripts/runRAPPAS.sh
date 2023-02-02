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
	workplace=rappas/$line
	databasefile=$workplace"/DB_session_k8_o1.5.union"
	
	java -Xmx64G -Xss2m -jar ~/yytao/package/RAPPAS/dist/RAPPAS.jar -p b -s nucl -b ~/yytao/package/raxml-ng/bin/raxml-ng -w $workplace -r $taxafile -t $treefile --use_unrooted
	java -Xmx64G -Xss2m -jar ~/yytao/package/RAPPAS/dist/RAPPAS.jar -p p -d $databasefile -q $queryfile
	
    echo $query" finished"
done

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

