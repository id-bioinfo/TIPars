import dr.evolution.sequence.Sequence;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.parsimony.FitchParsimony;
import dr.evolution.tree.Tree;
import dr.evolution.tree.FlexibleTree;
import dr.app.tools.NexusExporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.FastaImporter;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.stats.DiscreteStatistics;
import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import dr.evolution.alignment.SitePatterns;

// Author: Tommy Tsan-Yuk LAM, Guangchuang YU
// Description: Insert taxa into a reference phylogeny by parsimony.

public class TIPars{
    private boolean DEBUG = false;
    private SimpleAlignment taxaseq = null;
    private SimpleAlignment ancseq = null;
    private Tree mytree = null;
    private String internalnode_nidname = "nid";
    private HashMap <FlexibleNode, String> node2Sseq = null;
    private HashMap <String, String> name2seq = null;
    private FlexibleNodeBranch<FlexibleNode, Double> myBandBPbranch = new FlexibleNodeBranch<FlexibleNode, Double>(); // for evaluation
    public TIPars(){

    }

    public TIPars(SimpleAlignment taxaseq, SimpleAlignment ancseq, Tree mytree, String internalnode_nidname){
	this.mytree = mytree;
	this.taxaseq = taxaseq;
	this.ancseq = ancseq;
	this.internalnode_nidname = internalnode_nidname;
	setupHashtableOfNode2Seq();
    }

    private void setupHashtableOfNode2Seq(){
	node2Sseq = new HashMap <FlexibleNode, String>();
	for(int i=0; i<mytree.getInternalNodeCount(); i++){
	    FlexibleNode n = (FlexibleNode)mytree.getInternalNode(i);
	    int t = ancseq.getTaxonIndex((String)(n.getAttribute(this.internalnode_nidname)));
	    if(t >= 0){
		node2Sseq.put(n, ancseq.getAlignedSequenceString(t));
	    }
	    else if(DEBUG){ System.out.println("inode-attribute="+(String)(n.getAttribute(this.internalnode_nidname))+" not found in alignment."); }
	}
	for(int i=0; i<mytree.getExternalNodeCount(); i++){
	    FlexibleNode n = (FlexibleNode)mytree.getExternalNode(i);
	    int t = taxaseq.getTaxonIndex(n.getTaxon().getId());
	    if(t >= 0){
		node2Sseq.put(n, taxaseq.getAlignedSequenceString(t));
	    }
	    else if(DEBUG){ System.out.println("enode-taxa="+n.getTaxon().getId()+" not found in alignment."); }
	}
    }

    public void setDebug(boolean d){
	this.DEBUG = d;
    }

    public Tree addQuerySequence(String qname, String nodeQseq, String qid, String pid, boolean printDisInfoOnScreen, String nidname, String attname, double[] ABQ_brlen){
	// Travel all internal nodes for all possible nodeA-nodeB pairs
	//HashMap <FlexibleNode, int[]> nodeB2score = new HashMap<FlexibleNode, int[]>();   // map the *nodeB* to the score array (nodeA score, nodeB score, nodeQ score)
	//HashMap <FlexibleNode, String> nodeB2nodePseq = new HashMap<FlexibleNode, String>();  // map the *nodeB* to the nodeP sequence
	String nodePseq = "";
	int[] selectedScores = new int[3];
	int minQScore = 9999999;
	int selectedNodeBIndex = -1;// the index of selected nodeB (with min score) in the tree's internal .
	String selectedBnid = "";
	String selectedBatt = "";
	for(int i=0; i<mytree.getNodeCount(); i++){
	    // We use this annotation:
	    //    nodeA
	    //       |
	    //     nodeP
	    //       |  \
	    //       |   \
	    //    nodeB nodeQ

	    FlexibleNode nodeB = (FlexibleNode)mytree.getNode(i);
	    if(!nodeB.isRoot()){
		// nodeA-nodeB pairs
		FlexibleNode nodeA = (FlexibleNode)nodeB.getParent();
		/////////// IMPORTANT NOTE: Ignore the case when nodeA is root first.
		///////////                 Reason: 1) BEAST cannot intake root attribute; 2) HYPHY will not reconstruct root anc seq.
		///////////                 My be try to solve this problem later.
		////////////////////////////////////////////////////
		if(nodeA.isRoot()){
		    if(DEBUG){ System.out.println("Root attribute - "+ nodeA.getAttribute("nid")); }
		}
		else{
		    String nodeAseq = getSequenceByNode(nodeA);
		    String nodeBseq = getSequenceByNode(nodeB);
		    int[] scores = new int[3];
		    if(nodeAseq == null){
			if(DEBUG){ System.out.println("nodeAseq cannot find seq - "+nodeA.getAttribute("nid")); }
		    }
		    else if(nodeBseq == null){
			if(DEBUG){ System.out.println("nodeBseq cannot find seq - "+nodeB.getAttribute("nid")); }
		    }
		    else{
			Integer[] scores1 = new Integer[3];
			String tmp_nodePseq = getStringAndScoreFromNodeABQSeq(nodeAseq, nodeBseq, nodeQseq, scores1);
			scores[0] = scores1[0].intValue();
			scores[1] = scores1[1].intValue();
			scores[2] = scores1[2].intValue();
			if(scores[2] < minQScore){
			    if(DEBUG) System.out.println("Score: "+scores[2]);
			    String tmp_nodePseq3 = getStringAndScoreFromNodeABQSeq(nodeAseq, nodeBseq, nodeQseq, scores1);
			    minQScore = scores[2];
			    selectedScores = scores; // the mutations(scores) at A-P, B-P and P-Q branch after taxa insertion.
			    nodePseq = tmp_nodePseq;
			    selectedNodeBIndex = i;
			    selectedBnid = (String)nodeB.getAttribute(nidname);
			    selectedBatt = (String)nodeB.getAttribute(attname);
			}
		    }
		}
	    }
	}


	// add the query node to the Tree copy.
	double[] afterscores = {0.0, 0.0, 0.0};
	MyFlexibleTree mynewTree = new MyFlexibleTree(mytree, true);
	this.copyAttributeFromOneNodeToAnother((FlexibleNode)mytree.getRoot(), (FlexibleNode)mynewTree.getRoot());
	mynewTree.beginTreeEdit();
	FlexibleNode selected_nodeB = (FlexibleNode)mynewTree.getNode(selectedNodeBIndex);
	double original_B = selected_nodeB.getLength();
	FlexibleNode selected_nodeA = selected_nodeB.getParent();
	Taxon qtaxon = new Taxon(qname);
	FlexibleNode selected_nodeQ = new FlexibleNode(qtaxon);
	selected_nodeQ.setLength(selectedScores[2]/((double)getAlignmentLength()));
	FlexibleNode selected_nodeP = new FlexibleNode();
	// set the attributes of newly added node.
	selected_nodeP.setAttribute("nid", pid);
	selected_nodeQ.setAttribute("nid", qid);

	if(selectedScores[0] == 0){ // A-P is zero branch length, meaning that Q is inserted into A directly.
	    selected_nodeP = selected_nodeA;
	    selected_nodeA.addChild(selected_nodeQ);  //NOTE: should node A nid change to node P nid? AWAITING TO SOLVE!!
	    afterscores[0] = 0.0;
	    afterscores[1] = selected_nodeB.getLength();
	    afterscores[2] = selected_nodeQ.getLength();
	}
	else if(selectedScores[1] == 0){  // P-B is zero branch length, meaning that Q is inserted into B directly.
	    if(selected_nodeB.isExternal()){ // If B is leaf, cannot add the Q directly there, must add P node.
		double newNodeBLength = 0.0;
		double newNodePLength = selected_nodeB.getLength();
		selected_nodeP.setLength(newNodePLength);
		selected_nodeB.setLength(newNodeBLength);
		selected_nodeP.addChild(selected_nodeQ);
		selected_nodeA.removeChild(selected_nodeB);
		selected_nodeA.addChild(selected_nodeP);
		selected_nodeP.addChild(selected_nodeB);
	    }
	    else{
		selected_nodeB.addChild(selected_nodeQ);
	    }
	    afterscores[0] = selected_nodeB.getLength();
	    afterscores[1] = 0.0;
	    afterscores[2] = selected_nodeQ.getLength();
	}
	else{
	    double Pratio = selectedScores[0]/((double)(selectedScores[0]+selectedScores[1]));
	    double newNodeBLength = selected_nodeB.getLength()*(1.0-Pratio);
	    double newNodePLength = selected_nodeB.getLength()*Pratio;
	    selected_nodeP.setLength(newNodePLength);
	    selected_nodeB.setLength(newNodeBLength);
	    selected_nodeP.addChild(selected_nodeQ);
	    selected_nodeA.removeChild(selected_nodeB);
	    selected_nodeA.addChild(selected_nodeP);
	    selected_nodeP.addChild(selected_nodeB);
	    afterscores[0] = selected_nodeP.getLength();
	    afterscores[1] = selected_nodeB.getLength();
	    afterscores[2] = selected_nodeQ.getLength();
	}
	ABQ_brlen[0] = afterscores[0];
	ABQ_brlen[1] = afterscores[1];
	ABQ_brlen[2] = afterscores[2];
	myBandBPbranch = new FlexibleNodeBranch<FlexibleNode, Double>((FlexibleNode)(mytree.getNode(selectedNodeBIndex)), new Double(afterscores[1])); // For EVALUATION
	// myBandBPbranch.setNode(myTree.getNode(selectedNodeBIndex)); // For EVALUATION - must take the original tree node.
	// myBandBPbranch.setBrlen(new Double(afterscores[1]));  // For EVALUATION - take update branch length(sub/site, not absolute mutations).

	if(printDisInfoOnScreen) System.out.println(""+selectedScores[0]+"\t"+selectedScores[1]+"\t"+original_B+"\t"+selectedScores[2]+"\t"+afterscores[0]+"\t"+afterscores[1]+"\t"+afterscores[2]+"\t"+selectedBnid+"\t"+selectedBatt+"\n");
	mynewTree.endTreeEdit();
	if(DEBUG) System.out.println("TaxonCount: "+mynewTree.getTaxonCount());
	mynewTree.toAdoptNodes((FlexibleNode)mynewTree.getRoot());
	if(DEBUG) System.out.println("TaxonCount after taxalist refresh: "+mynewTree.getTaxonCount());
	return mynewTree;
    }

    class MyFlexibleTree extends FlexibleTree
    {
	MyFlexibleTree(Tree t){
	    super(t);
	}
	MyFlexibleTree(Tree t, boolean keepAttribute){
	    super(t, keepAttribute);
	}
	public void toAdoptNodes(FlexibleNode n){
	    super.adoptNodes(n);
	}

    }


    // Note: the n1 and n2 nodes are on the same tree. Not yet tested.
    public static double[] calcDisBetweenTwoNodes(FlexibleTree t, FlexibleNode n1, FlexibleNode n2, double n1_dis, double n2_dis){
	double node_dis = 0.0; // Number of node difference
	double dis = 0.0;
	if(n1 == n2){ dis = Math.abs(n1_dis-n2_dis); } // if n1 is same as n2
	else{
	    double[] dis0 = nodeAisAncestorOfNodeB(n1, n2);
	    if(dis0[0] >= 0){ dis = dis0[0]-n2_dis+n1_dis; node_dis = dis0[1];} // if n1 is above n2;
	    else{
		dis0 = nodeAisAncestorOfNodeB(n2, n1);
		if(dis0[0] >= 0){ dis = dis0[0]-n1_dis+n2_dis; node_dis = dis0[1];} // if n2 is above n1
		else{  // if n1 and n2 not on the same pathway - use the path thr their common ancestor
		    HashMap<FlexibleNode, double[]> listOf_n1_parents = new HashMap<FlexibleNode, double[]>(20);
		    dis0[0] = n1.getLength();
		    dis0[1] = 1;
		    FlexibleNode k = n1.getParent();
		    listOf_n1_parents.put(k, dis0);
		    while(!k.isRoot()){
			dis0[0] += k.getLength();
			dis0[1] += 1;
			k = k.getParent();
			listOf_n1_parents.put(k, dis0);
		    }
		    double dis1 = n2.getLength();
		    double dis1b = 1;
		    k = n2.getParent();
		    while(!listOf_n1_parents.containsKey(k)){
			dis1 += k.getLength();
			dis1b += 1;
			k = k.getParent();
		    }
		    dis1 += ((double[])(listOf_n1_parents.get(k)))[0]; // branch-len difference
		    dis = dis1-n1_dis-n2_dis; // branch-len difference that consider the n1_br-len and n2_br-len
		    node_dis = dis1b + ((double[])(listOf_n1_parents.get(k)))[1]; // node# difference
		}
	    }
	}
	return new double[]{dis, node_dis};
    }



    // If not found return -1, -99, else return branch length & number of node distance
    public static double[] nodeAisAncestorOfNodeB(FlexibleNode a, FlexibleNode b){
	if(a == b){
	    return new double[]{0, 0};
	}
	else{
	    double outdis = b.getLength();
	    FlexibleNode k = b.getParent();
	    boolean found = false;
	    double num_node = 1;
	    while(!k.isRoot()){
		if(k == a){
		    found = true;
		    break;
		}
		else{
		    outdis += k.getLength();
		    k = k.getParent();
		    num_node += 1;
		}
	    }
	    if(k.isRoot() && k == a){ outdis += k.getLength(); }
	    if(found){ return new double[]{outdis, num_node}; }
	    else{ return new double[]{-1, -99}; } // -99 is dummy value
	}
    }



    // Remove a taxon Q, and return "FlexibleNode B" AND "B-P br.length" as the below:
    //        A
    //        |
    //        P
    //       /|
    //      / |
    //     Q  B
    // Note the removal will be effective to the input tree object
    public static FlexibleNodeBranch<FlexibleNode, Double> removeTaxon(FlexibleTree t, FlexibleNode n){
	try{
	    t.beginTreeEdit();
	    FlexibleNode p = n.getParent();
	    p.removeChild(n);
	    n.setParent(null);
	    if(p.getChildCount() == 1){ // Remove this p node if if has less than 2 child nodes
		FlexibleNode b = (FlexibleNode)(p.getChild(0));
		FlexibleNode a = (FlexibleNode)(p.getParent());
		a.addChild(b);
		a.removeChild(p);
		b.setParent(a);
		double oldb2p = b.getLength();
		b.setLength(b.getLength()+p.getLength());
		p.setParent(null);
		t.endTreeEdit();
		return new TIPars().new FlexibleNodeBranch<FlexibleNode, Double>(b, new Double(oldb2p));
	    }
	    else{
		t.endTreeEdit();
		return new TIPars().new FlexibleNodeBranch<FlexibleNode, Double>(p, new Double(0.0));
	    }
	}
	catch(Exception e){
	    e.printStackTrace();
	    return null;
	}
    }

    // copy attribute from n to m
    private void copyAttributeFromOneNodeToAnother(FlexibleNode n, FlexibleNode m){
	Iterator attnames = n.getAttributeNames();
	while(attnames != null && attnames.hasNext()){
	    String ahname = (String)attnames.next();
	    m.setAttribute(ahname,  n.getAttribute(ahname));
	}
    }

    private String getSequenceByNode(FlexibleNode a){
	String seq = node2Sseq.get(a);
	return seq;
    }

    // Clean the "i102" problem in attribute whose value is String type. I have mentioned this in email.
    public static Tree cleanStringAttributeInTree(Tree t){
	for(int i=0; i<t.getNodeCount(); i++){
	    FlexibleNode nr = (FlexibleNode)t.getNode(i);
	    Iterator attnames = nr.getAttributeNames();
	    while(attnames != null && attnames.hasNext()){
		String ahname = (String)attnames.next();
		if(nr.getAttribute(ahname) instanceof String){
		    String v = (String)nr.getAttribute(ahname);
		    v = v.replaceAll("\"", "");
		    nr.setAttribute(ahname,  v);
		}
		else{
		    nr.setAttribute(ahname,  nr.getAttribute(ahname));
		}
	    }
	}
	return t;
    }

    // Suppose to contain the "FlexibleNode B" AND "B-P br.length" as the below:
    //        A
    //        |
    //        P
    //       /|
    //      / |
    //     Q  B
    public class FlexibleNodeBranch<FlexibleNode,Double> {
	public FlexibleNode a;
	public Double b;
	public FlexibleNodeBranch(FlexibleNode a, Double b) {
	    this.a = a;
	    this.b = b;
	}
	public FlexibleNodeBranch() {
	    this.a = null;
	    this.b = null;
	}
	public void setNode(FlexibleNode a2){
	    this.a = a2;
	}
	public void setBrlen(Double b2){
	    this.b = b2;
	}
	public FlexibleNode getNode(){
	    return this.a;
	}
	public Double getBrlen(){
	    return this.b;
	}
    }

    public FlexibleNodeBranch<FlexibleNode, Double> getNodeBandBPbrlen(){
	return myBandBPbranch;
    }

    public String getStringAndScoreFromNodeABQSeq(String a, String b, String c, Integer[] scores1){
	int[] scores = new int[3]; scores[0] = 0; scores[1] = 0; scores[2] = 0;
	a = a.toUpperCase();
	b = b.toUpperCase();
	c = c.toUpperCase(); // c is Q
	String p = "";
	for(int i=0; i<a.length(); i++){
	    char ai = a.charAt(i);
	    char bi = b.charAt(i);
	    char ci = c.charAt(i);
	    if(ai == ci && ai == bi){
		// do nothing
		if(DEBUG) System.out.print("AAA");
	    }
	    else if(ai != ci && ci != bi && ai != bi){   // ATC
		p += ai;     // NOTE: biased to parent node char; should try alternating ai and bi to balance the node p position
		scores[2]++;
		if(DEBUG) System.out.print("ATC");
	    }
	    else if(ai == bi && ci != bi){   // AAT
		p += ai;
		scores[2]++;
		if(DEBUG) System.out.print("AAT");
	    }
	    else if(ai != bi && ci == bi){   // ATT
		p += bi;
		scores[0]++;
		if(DEBUG) System.out.print("ATT");
	    }
	    else if(ai == ci && ci != bi){   // TAT
		p += ai;
		scores[1]++;
		if(DEBUG) System.out.print("TAT");
	    }
	    else if(DEBUG){ System.out.println("Unmatched ABQ type"); }
	}
	if(DEBUG) System.out.println("score: "+scores[0]+"|"+scores[1]+"|"+scores[2]);
	scores1[0] = new Integer(scores[0]);
	scores1[1] = new Integer(scores[1]);
	scores1[2] = new Integer(scores[2]);
	return p;
    }


    public int getAlignmentLength(){
	if(taxaseq != null){
	    return (taxaseq.getAlignedSequenceString(0)).length();
	}
	else{
	    return -1;
	}
    }


    public static SimpleAlignment readFastaAlignmentFile(String fn){
	/////////// Read the sequence/state file (expect fasta format)
	SimpleAlignment sa = null;
	try{ // This part is a parser for general data type. Also, it parser *manually* into SimpleAlignment. There is a FastaImporter in new BEAST package (below), but not make it work yet.
	    sa = new SimpleAlignment();
	    char[] agg = {'H', 'A', 'L', 'I', 'K', 'M', 'Y', 'C', 'E', 'P', 'G', 'S', 'D', 'T', 'F', 'V', 'W', 'N', 'X', '?', '-'};
	    GeneralDataType dt1 = new GeneralDataType(agg);
	    sa.setDataType(dt1);
	    BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
	    String fasline = br2.readLine();
	    String desc;
	    while(fasline != null){
		fasline = fasline.replaceAll("\n", "");
		String seqseq = "";
		if(fasline.matches("^>.+")){
		    desc = fasline;
		    fasline = br2.readLine();
		    while(fasline != null && !fasline.matches("^>.+")) {
			seqseq += fasline;
			fasline = br2.readLine();
		    }
		    desc = desc.replaceAll(">", "");
		    seqseq = seqseq.replaceAll("\n", "");
		    // manually make the alignment. Taxon is supposed to be name of the sequence
		    sa.addSequence(new dr.evolution.sequence.Sequence(new Taxon(desc), seqseq));
		}
		// fasline now store desc for next sequence
	    }
	}
	catch(Exception e){
	    e.printStackTrace();
	}
	return sa;
    }

    public static SimpleAlignment readSingleLineFastaAlignmentFile(String fn){
	/////////// Read the sequence/state file (expect fasta format)
	SimpleAlignment sa = null;
	try{ // This part is a parser for general data type. Also, it parser *manually* into SimpleAlignment. There is a FastaImporter in new BEAST package (below), but not make it work yet.
	    sa = new SimpleAlignment();
	    char[] agg = {'H', 'A', 'L', 'I', 'K', 'M', 'Y', 'C', 'E', 'P', 'G', 'S', 'D', 'T', 'F', 'V', 'W', 'N', 'X', '?', '-'};
	    GeneralDataType dt1 = new GeneralDataType(agg);
	    sa.setDataType(dt1);
	    BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
	    String line2 = br2.readLine();
	    while(line2 != null){
		line2 = line2.replaceAll("\n", "");
		if(line2.matches("^>.+")){
		    String seqseq = br2.readLine();
		    seqseq = seqseq.replaceAll("\n", "");
		    line2 = line2.replaceAll(">", "");
		    sa.addSequence(new dr.evolution.sequence.Sequence(new Taxon(line2), seqseq)); // manually make the alignment. Taxon is supposed to be name of the sequence
		}
		line2 = br2.readLine();
	    }
	}
	// 		try{ // Not yet functional.
	//  			BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
	//             FastaImporter importer = new FastaImporter(br2, Nucleotides.INSTANCE);  // Assume nucleotide sequences
	//             sa = (SimpleAlignment)(importer.importAlignment());
	// 		}
	catch(Exception e){
	    e.printStackTrace();
	}
	return sa;
    }

    // Read query fasta, only read one entry.
    public static String[] readSingleLineQueryFastaAlignmentFile(String fn){
	/////////// Read the sequence/state file (expect fasta format)
	String[] sa = new String[2];
	try{
	    BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
	    String line2 = br2.readLine();
	    while(line2 != null){
		line2 = line2.replaceAll("\n", "");
		if(line2.matches("^>.+")){
		    String seqseq = br2.readLine();
		    seqseq = seqseq.replaceAll("\n", "");
		    line2 = line2.replaceAll(">", "");
		    sa[0] = line2; sa[1] = seqseq;
		}
		line2 = br2.readLine();
		break;
	    }
	}
	catch(Exception e){
	    e.printStackTrace();
	}
	return sa;
    }


    public static void main(String[] args)
    {
	String insfn = "";
	String intfn = "";
	String inafn = "";
	String inqfn = "";
	String outfn = "";
	//String outdisfn = "";
	String nidname = "nid";
	String attname = "GenName";
	boolean outdis = true;
	try{
	    if(args.length != 8){
		intfn = getShellInput("Enter your input nexus tree file: ");
		insfn = getShellInput("Enter your input taxa seq file [fasta name is taxaname]: ");
		inafn = getShellInput("Enter your input ancestral seq file [fasta name is nid]: ");
		inqfn = getShellInput("Enter your input query seq file [fasta name is taxaname]: ");
		outfn = getShellInput("Enter your output nexus tree file: ");
		String tempstr = getShellInput("Want to output ABQdis, Bnid and genotype info? (0=no|1=yes): ");
		outdis = (tempstr.equals("0"))?false:true;
		if(outdis){
		    nidname = getShellInput("Enter your nid attribute name(e.g. nid): ");
		    attname = getShellInput("Enter your genotype attribute name (e.g. GenName): ");
		}
	    }
	    else{
		intfn = args[0];
		insfn = args[1];
		inafn = args[2];
		inqfn = args[3];
		outfn = args[4];
		String tempstr = args[5];
		outdis = (tempstr.equals("0"))?false:true;
		if(outdis){
		    nidname = args[6];
		    attname = args[7];
		}

	    }
	    /////////// Read the sequence/state file (expect fasta format)
	    SimpleAlignment taxa_align = readSingleLineFastaAlignmentFile(insfn);
	    SimpleAlignment anc_align = readSingleLineFastaAlignmentFile(inafn);
	    String[] query_align = readSingleLineQueryFastaAlignmentFile(inqfn);


	    // Read the tree from the nexus tree file. Note that nexus can also hold alignment, but not intended in this example.
	    NexusImporter tni = new NexusImporter(new FileReader(intfn));
	    Tree[] trees = (Tree[])tni.importTrees(taxa_align);
	    trees[0] = cleanStringAttributeInTree(trees[0]);
	    TIPars myAdd = new TIPars(taxa_align, anc_align, trees[0], "nid");
	    Tree outtree = myAdd.addQuerySequence(query_align[0], query_align[1], "q1", "p1", outdis, nidname, attname, new double[3]); // q1 and p1 are the attributes of nodeQ and nodeP.
	    PrintStream fw = new PrintStream(new FileOutputStream(new File(outfn)));
	    NexusExporter kne = new NexusExporter(fw); // export the tree *with attributes* to the output nexus file.
	    kne.exportTree(outtree);
	    fw.close();
	}
	catch(Exception e){
	    e.printStackTrace();
	}
    }




    public static void writeNexusTree(Tree t, String fn){
	try{
	    PrintStream fw = new PrintStream(new FileOutputStream(new File(fn)));
	    NexusExporter kne = new NexusExporter(fw); // export the tree *with attributes* to the output nexus file.
	    Tree[] ts = new Tree[1];
	    ts[0] = t;
	    kne.exportTrees(ts, true);
	    fw.close();
	}
	catch(Exception e){
	    e.printStackTrace();
	}
    }

    public static String getShellInput(String question){
	String myfile = "";
	try
	    {
		BufferedReader stdin =
		    new BufferedReader(
				       new InputStreamReader(System.in));
		System.out.print(question);
		myfile = stdin.readLine();
	    }
	catch (Exception e)
	    {
		System.out.println("Error in screen reading");
	    }
	return myfile;
    }

}
