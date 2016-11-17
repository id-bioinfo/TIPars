import java.util.*;
import java.io.*;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;

public class testTreeReader {
    public static void main(String args[]) {
	String intfn = "FourExamples/Try1_KC542905/NDV_fullg.nex";
	String insfn = "FourExamples/Try1_KC542905/NDV_fullg.fas";
	String nwk = "FourExamples/Try1_KC542905/NDV_fullg.tree";

	TIPars  x = new TIPars();
	SimpleAlignment taxa_align = x.readFastaAlignmentFile(insfn);
	try{
	    NexusImporter tni = new NexusImporter(new FileReader(intfn));
	    Tree tree = tni.importTree(taxa_align);
	    System.out.println(tree.getNodeCount());


	    NewickImporter tnwk = new NewickImporter(new FileReader(nwk));
	    Tree tree2 = tnwk.importTree(taxa_align);
	    System.out.println(tree2.getNodeCount());

	} catch(Exception e) {
	    e.printStackTrace();
	}
    }
}
