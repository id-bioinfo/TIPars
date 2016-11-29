import java.util.*;
import java.io.*;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evolution.tree.FlexibleNode;

public class testTreeReader {
    public static void main(String args[]) {
	String intfn = "../FourExamples/Try1_KC542905/NDV_fullg.nex";
	String insfn = "../FourExamples/Try1_KC542905/NDV_fullg.fas";
	String ancfn = "../FourExamples/Try1_KC542905/anc_NDV_fullg.fas";
	String nwk = "../FourExamples/Try1_KC542905/NDV_fullg.tree";

	TIPars  x = new TIPars();
	SimpleAlignment taxa_align = x.readFastaAlignmentFile(insfn);
	SimpleAlignment anc_align = x.readFastaAlignmentFile(ancfn);

	try{
	    NexusImporter tni = new NexusImporter(new FileReader(intfn));
	    Tree tree = tni.importTree(taxa_align);
	    tree = x.cleanStringAttributeInTree(tree);
	    System.out.println(tree.getInternalNodeCount());


	    NewickImporter tnwk = new NewickImporter(new FileReader(nwk));
	    Tree tree2 = tnwk.importTree(taxa_align);
	    tree2 = x.cleanStringAttributeInTree(tree2);
	    System.out.println(tree2.getInternalNodeCount());

	    //tree2 = tree;
	    for (int i=0; i<tree2.getInternalNodeCount(); i++) {
		FlexibleNode n = (FlexibleNode) tree2.getInternalNode(i);
		Iterator attnames = n.getAttributeNames();
		while(attnames != null && attnames.hasNext()) {
		    String ahname = (String) attnames.next();
		    //System.out.println(ahname);
		}
		//		int t = anc_align.getTaxonIndex((String)(n.getAttribute("nid")));
		System.out.println("node number:" + i + "\t" + "nid:" + n.getAttribute("label"));
		//System.out.println("node-index$" + t + "    " + "inode-attribute="+(String)(n.getAttribute("nid")));
	    }



	} catch(Exception e) {
	    e.printStackTrace();
	}
    }
}
