import java.io.*;
import java.util.*;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;


public class testReadFASTA {
    public static void main(String args[]) {
	String fn = "testData/sampleFASTA.fas";
	TIPars  x = new TIPars();
	//SimpleAlignment sa = x.readSingleLineFastaAlignmentFile(fn);
	SimpleAlignment sa = x.readFastaAlignmentFile(fn);

	for (int i=0; i<sa.getSequenceCount(); i++) {
	    Sequence seqseq = sa.getSequence(i);
	    String desc = seqseq.getTaxon().getId();
	    String seq = seqseq.getSequenceString();
	    System.out.println(">" + desc + "\n" + seq + "\n");
	}
    }

}
