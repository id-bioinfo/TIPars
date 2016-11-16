import java.io.*;
import java.util.*;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;


public class testReadFASTA {
    public static void main(String args[]) {
	String fn = "testData/sampleFASTA.fas";
	TIPars  x = new TIPars();
	SimpleAlignment sa = x.readFastaAlignmentFile(fn);
    }

}
