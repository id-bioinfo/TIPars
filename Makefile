TIPars: TIPars.java
	javac -classpath .:beast.jar TIPars.java;\
	jar cvfm TIPars.jar MANIFEST.MF *class;\
	rm -rf *class

test: TIPars
	java -classpath .:beast.jar TIPars FourExamples/Try1_KC542905/NDV_fullg.tree FourExamples/Try1_KC542905/NDV_fullg.fas FourExamples/Try1_KC542905/anc_NDV_fullg.fas FourExamples/Try1_KC542905/KC542905.fas test.nex 1 label GenName

testRscript: TIPars
	./tipars -t FourExamples/Try1_KC542905/NDV_fullg.tree -s FourExamples/Try1_KC542905/NDV_fullg.fas -a FourExamples/Try1_KC542905/anc_NDV_fullg.fas -q FourExamples/Try1_KC542905/KC542905.fas -o test.nex


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

