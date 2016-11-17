TIPars:
	javac -classpath .:beast.jar TIPars.java

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
