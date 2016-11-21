TIPars:
	javac -classpath .:beast.jar TIPars.java;\
	java -classpath .:beast.jar TIPars FourExamples/Try1_KC542905/NDV_fullg.tree FourExamples/Try1_KC542905/NDV_fullg.fas FourExamples/Try1_KC542905/anc_NDV_fullg.fas FourExamples/Try1_KC542905/KC542905.fas test.nex 1 label GenName

testReadFASTA:
	javac -classpath .:beast.jar testReadFASTA.java;\
	echo "";\
	echo "--> sample fasta file:";\
	echo "";\
	cat testData/sampleFASTA.fas;\
	echo "";\
	echo "--> file contains parsed by readFastaAlignmentFile:";\
	echo "";\
	java -classpath .:beast.jar testReadFASTA

testTreeReader:
	javac -classpath .:beast.jar testTreeReader.java;\
	java -classpath .:beast.jar testTreeReader

