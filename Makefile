TIPars: TIPars.java
	javac -classpath .:beast.jar ./tippack/TIPars.java;\
	jar cvfm TIPars.jar MANIFEST.MF ./tippack/*class;\
	rm -rf ./tippack/*class

test: TIPars
	java -classpath .:beast.jar TIPars FourExamples/Try1_KC542905/NDV_fullg.tree FourExamples/Try1_KC542905/NDV_fullg.fas FourExamples/Try1_KC542905/anc_NDV_fullg.fas FourExamples/Try1_KC542905/KC542905.fas test.tree 1 label GenName

testRscript: TIPars
	./tipars.rsh -t FourExamples/Try1_KC542905/NDV_fullg.tree -s FourExamples/Try1_KC542905/NDV_fullg.fas -a FourExamples/Try1_KC542905/anc_NDV_fullg.fas -q FourExamples/Try1_KC542905/KC542905.fas -o test.nex


testscript: TIPars
	./tipars -t FourExamples/rm2.tree -s FourExamples/taxa_seq.fas -a FourExamples/anc_seq.fas -q FourExamples/KC542905_KM885167.fas -o test.nwk


testReadFASTA:
	cd test;\
	javac -classpath ..:../beast.jar testReadFASTA.java;\
	echo "";\
	echo "--> sample fasta file:";\
	echo "";\
	cat testData/sampleFASTA.fas;\
	echo "";\
	echo "--> file contains parsed by readFastaAlignmentFile:";\
	echo "";\
	java -classpath ..:../beast.jar:../TIPars.jar testReadFASTA

testTreeReader:
	cd test;\
	javac -classpath ..:../beast.jar:../TIPars.jar testTreeReader.java;\
	java -classpath ..:../beast.jar:../TIPars.jar testTreeReader

