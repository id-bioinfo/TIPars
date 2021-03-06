import dr.app.tools.NexusExporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.NodeRef;
import dr.evolution.parsimony.FitchParsimony;
import dr.evolution.tree.Tree;
import dr.evolution.tree.FlexibleTree;
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


// Author: Tommy Tsan-Yuk LAM, Guangchuang YU
// Description: Insert taxa into a reference phylogeny by parsimony.
// last update 2019-01-04

public class TIPars{
    private boolean DEBUG = false;
    private boolean OUTPUT_PSEQ = false;
    private SimpleAlignment taxaseq = null;
    private SimpleAlignment ancseq = null;
    private Tree mytree = null;
    private String internalnode_nidname = "nid";
    private HashMap <FlexibleNode, String> node2Sseq = null;
    private HashMap <String, String> name2seq = null;
    private static HashMap <FlexibleNode, Integer> node2edge = new HashMap <FlexibleNode, Integer>();
    private FlexibleNodeBranch<FlexibleNode, Double> myBandBPbranch = new FlexibleNodeBranch<FlexibleNode, Double>(); // for evaluation
    private String model;
    private String gap;
    private static String[] placements = null;
    private static int edge_number = 1;
    private String OUTPUT_FOLDER;
    public TIPars(){

    }

    public TIPars(SimpleAlignment taxaseq, SimpleAlignment ancseq, Tree mytree, String internalnode_nidname, String model, String gap, String otype, String OUTPUT_FOLDER){
        this.mytree = mytree;
        this.taxaseq = taxaseq;
        this.ancseq = ancseq;
        this.model = model;
        this.gap = gap;
        this.internalnode_nidname = internalnode_nidname;
        this.OUTPUT_FOLDER = OUTPUT_FOLDER;
        setupHashtableOfNode2Seq();
        setupHashtableOfNode2edge(otype);
    }

    private void setupHashtableOfNode2edge(String otype) {
        if (otype.equals("placement")) {

            for (int i=0; i<mytree.getExternalNodeCount(); i++) {
                Integer edge = new Integer(edge_number);
                FlexibleNode node = (FlexibleNode) mytree.getExternalNode(i);
                node2edge.put(node, edge);
                edge_number++;
            }

            for (int i=0; i<mytree.getInternalNodeCount(); i++) {
                Integer edge = new Integer(edge_number);
                FlexibleNode node = (FlexibleNode) mytree.getInternalNode(i);
                node2edge.put(node, edge);
                edge_number++;
            }
        }
    }

    private void setupHashtableOfNode2Seq(){
        node2Sseq = new HashMap <FlexibleNode, String>();
        for(int i=0; i<mytree.getInternalNodeCount(); i++){
            FlexibleNode n = (FlexibleNode)mytree.getInternalNode(i);
            int t = ancseq.getTaxonIndex((String)(n.getAttribute(this.internalnode_nidname)));
            if(t >= 0){
                node2Sseq.put(n, ancseq.getAlignedSequenceString(t));
            } else if(DEBUG){
                System.out.println("inode-attribute="+(String)(n.getAttribute(this.internalnode_nidname))+" not found in alignment.");
            }
        }
        for(int i=0; i<mytree.getExternalNodeCount(); i++){
            FlexibleNode n = (FlexibleNode)mytree.getExternalNode(i);
            int t = taxaseq.getTaxonIndex(n.getTaxon().getId());
            if(t >= 0){
                node2Sseq.put(n, taxaseq.getAlignedSequenceString(t));
            } else if(DEBUG){
                System.out.println("enode-taxa="+n.getTaxon().getId()+" not found in alignment.");
            }
        }
    }

    public void setDebug(boolean d){
        this.DEBUG = d;
    }

    public Tree addQuerySequence(String qname, String nodeQseq, String qid, String pid, boolean printDisInfoOnScreen,
                                 String nidname, String attname, double[] ABQ_brlen, String otype, int ii){
        // Travel all internal nodes for all possible nodeA-nodeB pairs
        // HashMap <FlexibleNode, int[]> nodeB2score = new HashMap<FlexibleNode, int[]>();
        // map the *nodeB* to the score array (nodeA score, nodeB score, nodeQ score)
        // HashMap <FlexibleNode, String> nodeB2nodePseq = new HashMap<FlexibleNode, String>();
        // map the *nodeB* to the nodeP sequence
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
            //       | \
            //       |  \
            //    nodeB nodeQ

            FlexibleNode nodeB = (FlexibleNode)mytree.getNode(i);

            if(!nodeB.isRoot()){
                // nodeA-nodeB pairs
                FlexibleNode nodeA = (FlexibleNode)nodeB.getParent();

                String nodeAseq = getSequenceByNode(nodeA);
                String nodeBseq = getSequenceByNode(nodeB);
                int[] scores = new int[3];
                if(nodeAseq == null){
                    if(DEBUG){ System.out.println("nodeAseq cannot find seq - "+nodeA.getAttribute(this.internalnode_nidname)); }
                } else if(nodeBseq == null){
                    if(DEBUG){ System.out.println("nodeBseq cannot find seq - "+nodeB.getAttribute(this.internalnode_nidname)); }
                } else{
                    Integer[] scores1 = new Integer[3];
                    String tmp_nodePseq = getStringAndScoreFromNodeABQSeq(nodeAseq, nodeBseq, nodeQseq, scores1);
                    scores[0] = scores1[0].intValue();
                    scores[1] = scores1[1].intValue();
                    scores[2] = scores1[2].intValue();
                    if(scores[2] < minQScore){
                        if(DEBUG) System.out.println("Score: "+scores[2]);
                        minQScore = scores[2];
                        selectedScores = scores; // the mutations(scores) at A-P, B-P and P-Q branch after taxa insertion.
                        nodePseq = tmp_nodePseq;
                        selectedNodeBIndex = i;
                        selectedBnid = (String)nodeB.getAttribute(nidname);
                        if (selectedBnid == null) {
                            selectedBnid = nodeB.getTaxon().getId();
                        }
                        selectedBatt = (String)nodeB.getAttribute(attname);
                    }
                }
                // }
            }
        }

        if (DEBUG) System.out.println("selectedScores: " + selectedScores[0] + "\t" + selectedScores[1] + "\t" + selectedScores[2]);
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
        // Q-P pendent length

        double pqlen;
        if (model.equals("JC69")) {
            // JC69
            double p = selectedScores[2]/((double)getAlignmentLength());
            pqlen = JC69(p);
        } else if (model.equals("K2P")) {
            // K2P
            pqlen = K2P(nodePseq, nodeQseq);
        } else {
            // local estimation
            pqlen = localEstimation(selectedScores, original_B);
        }
        selected_nodeQ.setLength(pqlen);

        FlexibleNode selected_nodeP = new FlexibleNode();
        // set the attributes of newly added node.
        selected_nodeP.setAttribute(this.internalnode_nidname, pid);
        selected_nodeQ.setAttribute(this.internalnode_nidname, qname);

        if(selectedScores[0] == 0){ // A-P is zero branch length, meaning that Q is inserted into A directly.
            selected_nodeP = selected_nodeA;
            selected_nodeA.addChild(selected_nodeQ);  
            afterscores[0] = 0.0;
            afterscores[1] = selected_nodeB.getLength();
            afterscores[2] = selected_nodeQ.getLength();
        } else if(selectedScores[1] == 0){  // P-B is zero branch length, meaning that Q is inserted into B directly.
            if(selected_nodeB.isExternal()) { // If B is leaf, cannot add the Q directly there, must add P node.
                double newNodeBLength = 0.0;
                double newNodePLength = selected_nodeB.getLength();
                selected_nodeP.setLength(newNodePLength);
                selected_nodeB.setLength(newNodeBLength);
                selected_nodeP.addChild(selected_nodeQ);
                selected_nodeA.removeChild(selected_nodeB);
                selected_nodeA.addChild(selected_nodeP);
                selected_nodeP.addChild(selected_nodeB);
            } else {
                selected_nodeB.addChild(selected_nodeQ);
            }
            afterscores[0] = selected_nodeB.getLength();
            afterscores[1] = 0.0;
            afterscores[2] = selected_nodeQ.getLength();
        } else {
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


        if (DEBUG)
            System.out.println("ABQ_brlen: " + ABQ_brlen[0] + "\t" + ABQ_brlen[1] + "\t" + ABQ_brlen[2] + "\t" + selectedScores[2]*(ABQ_brlen[0]+ABQ_brlen[1])/(selectedScores[0]+selectedScores[1]));

        myBandBPbranch = new FlexibleNodeBranch<FlexibleNode, Double>((FlexibleNode)(mytree.getNode(selectedNodeBIndex)), new Double(afterscores[1])); // For EVALUATION

        if(printDisInfoOnScreen)
            System.out.println(""+selectedScores[0]+"\t"+selectedScores[1]+"\t"+original_B+"\t"+selectedScores[2]+"\t"+afterscores[0]+"\t"+afterscores[1]+"\t"+afterscores[2]+"\t"+selectedBnid+"\t"+selectedBatt+"\n");

        System.out.println("\nquery sequence: " + qname);
        System.out.println("insert to edge: " + (String)selected_nodeA.getAttribute(nidname) + "-" + selectedBnid);
        System.out.println("distal_length: " + afterscores[0]);
        System.out.println("pendant_length: " + afterscores[2] + "\n");


        if (otype.equals("placement")) {
            String placeInfo = "\t{\"p\":[" + node2edge.get((FlexibleNode) mytree.getNode(selectedNodeBIndex)) + ", " + afterscores[0]  + ", " + afterscores[2] +  "]," + "\"n\":[\"" + qname + "\"]}";
            placements[ii] = placeInfo;
            if (otype.equals("placement")) {
                return mytree;
            }
        }

        mynewTree.endTreeEdit();
        if(DEBUG) System.out.println("TaxonCount: "+mynewTree.getTaxonCount());
        mynewTree.toAdoptNodes((FlexibleNode)mynewTree.getRoot());
        if(DEBUG) System.out.println("TaxonCount after taxalist refresh: "+mynewTree.getTaxonCount());

        // mynewTree (MyFlexibleTree object) is a new copy of mytree (Tree object), and need to refresh hash table of node to seq
        // add the sequence to ancseq and taxaseq for running setupHashtableOfNode2Seq to refresh the hash table
        ancseq.addSequence(new  dr.evolution.sequence.Sequence(new Taxon(pid), nodePseq));
        taxaseq.addSequence(new dr.evolution.sequence.Sequence(new Taxon(qname), nodeQseq));

        if (OUTPUT_PSEQ) writeFASTA(pid, nodePseq, OUTPUT_FOLDER);
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

    public static double JC69(double p) {
        double d = -3.0 * Math.log(1-4.0*p/3.0) / 4.0;
        return d;
    }

    public static double K2P(String nodePseq, String nodeQseq) {
        int S = 0;
        int V = 0;
        int n = nodePseq.length();

        boolean precedingGapQ = true;

        for (int i=0; i<n; i++) {
            if (nodeQseq.charAt(i) != '-') {
                precedingGapQ = false;
            }

            if (precedingGapQ && nodeQseq.charAt(i) == '-') {
                continue;
            }

            if (nodePseq.charAt(i) == nodeQseq.charAt(i)) {
                continue;
            } else if (nodePseq.charAt(i) == 'A' && nodeQseq.charAt(i) == 'G') {
                S++;
            } else if (nodePseq.charAt(i) == 'G' && nodeQseq.charAt(i) == 'A') {
                S++;
            } else if (nodePseq.charAt(i) == 'C' && nodeQseq.charAt(i) == 'T') {
                S++;
            } else if (nodePseq.charAt(i) == 'T' && nodeQseq.charAt(i) == 'C') {
                S++;
            } else {
                V++;
            }
        }

        double s = 1.0 * S / (1.0*n);
        double v = 1.0 * V / (1.0*n);
        double d = -1.0 * Math.log(1.0-2.0*s-v)/2.0 - 1.0 * Math.log(1.0-2.0*v)/4.0;
        //System.out.println(s + "\t" + v + "\t" + d);
        return d;
    }

    public static double localEstimation(int[] scores, double ABbranch) {
        if (scores[0] + scores[1] == 0) {
            return 0;
        }
        double d = ABbranch * ((double) scores[2])/((double) (scores[0]+scores[1]));

        return d;
    }


    // If not found return -1, -99, else return branch length & number of node distance
    public static double[] nodeAisAncestorOfNodeB(FlexibleNode a, FlexibleNode b){
        if(a == b){
            return new double[]{0, 0};
        } else{
            double outdis = b.getLength();
            FlexibleNode k = b.getParent();
            boolean found = false;
            double num_node = 1;
            while(!k.isRoot()){
                if(k == a){
                    found = true;
                    break;
                } else{
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
            } else{
                t.endTreeEdit();
                return new TIPars().new FlexibleNodeBranch<FlexibleNode, Double>(p, new Double(0.0)); //? p, p.getLentth ?
            }
        } catch(Exception e){
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
                } else{
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


    private void getSequenceComparisonPosition(int[] position, String a, String b, String c) {
        position[0] = 0;
        position[1] = a.length()-1;

        if (gap.equals("distinctive") || gap.equals("inner")) {
            endTrailingGap3(position, a, b, c);
        }
    }

    private static void endTrailingGap3(int[] position, String a, String b, String c) {
        endTrailingGap(position, a);
        endTrailingGap(position, b);
        endTrailingGap(position, c);
    }

    private static void endTrailingGap(int[] pos, String a) {
        int i=-1;
        while (i < a.length()) {
            ++i;
            if (a.charAt(i) != '-') {
                break;
            }
        }
        int j=a.length();
        while (j >= 0) {
            --j;
            if (a.charAt(j) != '-') {
                break;
            }
        }
        if (pos[0] < i)
            pos[0] = i;
        if(pos[1] > j)
            pos[1] = j;
    }

    public String getStringAndScoreFromNodeABQSeq(String a, String b, String c, Integer[] scores1){
        int[] scores = new int[3]; scores[0] = 0; scores[1] = 0; scores[2] = 0;
        a = a.toUpperCase();
        b = b.toUpperCase();
        c = c.toUpperCase(); // c is Q
        // String p = "";
        boolean precedingGapC = true;
        boolean precedingGapB = true;

        StringBuilder p = new StringBuilder(c);

        int[] position = new int[2];
        getSequenceComparisonPosition(position, a, b, c);
        // System.out.println(a.length() + "\t" + position[0] + "\t" + position[1] + "\n");
        boolean ingoreGap = false;
        if (gap.equals("ignore")) {
            ingoreGap = true;
        }

        for(int i=position[0]; i<=position[1]; i++){
            char ai = a.charAt(i);
            char bi = b.charAt(i);
            char ci = c.charAt(i);
            if(ai == ci && ai == bi){
                continue;
                // do nothing
            } else if(ai != ci && ci != bi && ai != bi){   // ATC
                // prefer to assign the character using the most closely seqence.
                if (scores[2] <= scores[0] && scores[2] <= scores[1]) {
                    p.setCharAt(i, ci);
                    scores[0]++;
                    scores[1]++;
                } else if (scores[1] <= scores[0] && scores[1] <= scores[2]) {
                    p.setCharAt(i, bi);
                    scores[0]++;
                    scores[2]++;
                } else {
                    p.setCharAt(i, ai);
                    scores[1]++;
                    scores[2]++;
                }


                if (ingoreGap && (ai == '-' || bi == '-' || ci == '-')) {
                    continue;
                }
            } else if(ai == bi && ci != bi){   // AAT
                //p += ai;
                p.setCharAt(i, ai);

                if (ingoreGap && (ai == '-' || bi == '-' || ci == '-')) {
                    continue;
                }
                scores[2]++;
            }
            else if(ai != bi && ci == bi){   // ATT
                //p += bi;
                p.setCharAt(i, bi);
                if (ingoreGap && (ai == '-' || bi == '-' || ci == '-')) {
                    continue;
                }
                scores[0]++;
            }
            else if(ai == ci && ci != bi){   // TAT
                p.setCharAt(i, ai);

                if (ingoreGap && (ai == '-' || bi == '-' || ci == '-')) {
                    continue;
                }

                scores[1]++;
            }
            else if(DEBUG){ System.out.println("Unmatched ABQ type"); }
        }
        if(DEBUG) System.out.println("score: "+scores[0]+"|"+scores[1]+"|"+scores[2]);
        scores1[0] = new Integer(scores[0]);
        scores1[1] = new Integer(scores[1]);
        scores1[2] = new Integer(scores[2]);
        String pseq = p.toString();
        //	System.out.println(pseq.length() + "\t" + c.length());
        return pseq;
    }


    public int getAlignmentLength(){
        if(taxaseq != null){
            return (taxaseq.getAlignedSequenceString(0)).length();
        } else{
            return -1;
        }
    }


    public static SimpleAlignment readFastaAlignmentFile(String fn){
        /////////// Read the sequence/state file (expect fasta format)
        SimpleAlignment sa = null;
        try{
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

    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();

        String insfn = "";
        String intfn = "";
        String inafn = "";
        String inqfn = "";
        String inm   = "";
        String ing   = "";
        String outfn = "";
        String otype = "";
        //String outdisfn = "";
        String nidname = "label";
        String attname = "GenName";
        boolean outdis = true;
        try{
            if(!(args.length == 11 || args.length == 9)){
                intfn = getShellInput("Enter your input tree file: ");
                insfn = getShellInput("Enter your input taxa seq file [fasta name is taxaname]: ");
                inafn = getShellInput("Enter your input ancestral seq file [fasta name is nid]: ");
                inqfn = getShellInput("Enter your input query seq file [fasta name is taxaname]: ");
                inm   = getShellInput("Enter substitution model for estimating branch ['LE', 'JC69' or 'K2P']: ");
                ing   = getShellInput("Enter gap option ['ignore', 'inner' or 'all']: ");
                outfn = getShellInput("Enter your output tree file: ");
                otype = getShellInput("Enter algorithm type [insertion or placement]: ");
                String tempstr = getShellInput("Want to output ABQdis, Bnid and genotype info? (0=no|1=yes): ");
                outdis = (tempstr.equals("0"))?false:true;
                if(outdis){
                    nidname = getShellInput("Enter your nid attribute name(e.g. nid): ");
                    attname = getShellInput("Enter your genotype attribute name (e.g. GenName): ");
                }
            } else{
                intfn = args[0];
                insfn = args[1];
                inafn = args[2];
                inqfn = args[3];
                inm   = args[4];
                ing   = args[5];
                outfn = args[6];
                otype = args[7];
                String tempstr = args[8];
                outdis = (tempstr.equals("0"))?false:true;
                if(outdis){
                    nidname = args[9];
                    attname = args[10];
                }

            }

            String OUTPUT_FOLDER = getFolder(outfn);

            /////////// Read the sequence/state file (expect fasta format)
            SimpleAlignment taxa_align = readFastaAlignmentFile(insfn);
            SimpleAlignment anc_align  = readFastaAlignmentFile(inafn);

            //System.out.println(query_align[0] + "\n\n\n" + query_align[1]);


            NewickImporter tni = new NewickImporter(new FileReader(intfn));
            Tree tree = tni.importTree(taxa_align);
            TIPars myAdd = new TIPars(taxa_align, anc_align, tree, nidname, inm, ing, otype, OUTPUT_FOLDER);

            Tree outtree = null;
            long startTime2 = System.currentTimeMillis();
            SimpleAlignment query_align = readFastaAlignmentFile(inqfn);
            if (otype.equals("insertion")) {
                int nq = query_align.getTaxonCount();
                for (int i=0; i<nq; i++) {
                    String[] tmp_query = new String[2];
                    tmp_query[0] = query_align.getTaxonId(i);
                    tmp_query[1] = query_align.getAlignedSequenceString(i);
                    String qid = "q" + (i+1);
                    String pid = "p" + (i+1);
                    outtree = myAdd.addQuerySequence(tmp_query[0], tmp_query[1], qid, pid, outdis, nidname, attname, new double[3], otype, 0);
                    // q1 and p1 are the attributes of nodeQ and nodeP.
                    myAdd.mytree = outtree;
                    myAdd.setupHashtableOfNode2Seq();
                }
            } else {
                int nq = query_align.getTaxonCount();
                placements = new String[nq];
                for (int i=0; i<nq; i++) {
                    String[] tmp_query = new String[2];
                    tmp_query[0] = query_align.getTaxonId(i);
                    tmp_query[1] = query_align.getAlignedSequenceString(i);
                    String qid = "q" + (i+1);
                    String pid = "p" + (i+1);
                    outtree = myAdd.addQuerySequence(tmp_query[0], tmp_query[1], qid, pid, outdis, nidname, attname, new double[3], otype, i);
                }
            }
            long endTime2 = System.currentTimeMillis();
            long totalTime2 = endTime2 - startTime2;

             writeToTree(outtree, outfn, otype);
            //writeToTree(myAdd.mytree, outfn, otype);

            long endTime = System.currentTimeMillis();
            long totalTime = endTime - startTime;
            System.out.println("Insertion time: " + (double) totalTime2/1000);
            System.out.println("Overall time: " + (double) totalTime/1000);
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    public static String getFolder(String filename) {
        String dir = filename.substring(0, filename.lastIndexOf("/") + 1);
        return dir;
    }

    public static void writeToTree(Tree tree, String fn, String otype) {
        try {
            StringBuilder buffer = new StringBuilder();
            if (otype.equals("placement")) {
                buffer.append("{\n\t\"tree\": \"");
            }
            toNewick(tree, (FlexibleNode) tree.getRoot(), buffer, otype);
            buffer.append(';');

            if (otype.equals("placement")) {
                buffer.append("\",\n\t\"placements\": [\n");
                for (int i=0; i<placements.length; i++) {
                    buffer.append(placements[i]);
                    if (i == placements.length -1) {
                        buffer.append("\n\t],\n");
                    } else {
                        buffer.append(",\n");
                    }
                }
                buffer.append("\t\"metadata\": {\"info\": \"placement using TIPars\"},\n");
                buffer.append("\t\"version\": 201703,\n");
                buffer.append("\t\"fields\": [\"edge_num\", \"distal_length\", \"pendant_length\"\n\t]\n}\n");
            }


            PrintStream out = new PrintStream(new FileOutputStream(new File(fn)));
            out.println(buffer.toString());
            out.close();
        } catch(Exception e){
            e.printStackTrace();
        }
    }

    public static void toNewick(Tree tree, FlexibleNode node, StringBuilder buffer, String otype) {
        if (tree.isExternal(node)) {
            String label = tree.getTaxonId(node.getNumber());
            buffer.append(label);
            appendLength(tree, node, buffer);
            appendEdgeNumber(buffer, node, otype);
        } else {
            buffer.append('(');
            int n = tree.getChildCount(node);
            for (int i=0; i<n; i++) {
                toNewick(tree, (FlexibleNode) tree.getChild(node, i), buffer, otype);
                if (i == (n-1)) {
                    buffer.append(')');
                    String label = (String) node.getAttribute("label");
                    if (label != null)
                        buffer.append(label);
                } else {
                    buffer.append(',');
                }
            }
            FlexibleNode parentNode = (FlexibleNode) tree.getParent(node);
            if (parentNode != null) {
                appendLength(tree, node, buffer);
                appendEdgeNumber(buffer, node, otype);
            }
        }
    }

    private static void appendEdgeNumber(StringBuilder buffer, FlexibleNode node, String otype) {
        if (otype.equals("placement")) {
            buffer.append('{');
            int edge = node2edge.get(node).intValue();
            buffer.append(edge);
            buffer.append('}');
        };
    }

    private static void appendLength(Tree tree, FlexibleNode node, StringBuilder buffer) {
        if (tree.hasBranchLengths()) {
            buffer.append(':');
            buffer.append(tree.getBranchLength(node));
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
        } catch(Exception e){
            e.printStackTrace();
        }
    }

    public static void writeFASTA(String desc, String seq, String OUTPUT_FOLDER) {
        String output = OUTPUT_FOLDER + desc + ".fas";
        StringBuilder buffer = new StringBuilder();
        buffer.append("> ");
        buffer.append(desc);
        buffer.append("\n");
        buffer.append(seq);
        try{
            PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
            out.println(buffer.toString());
            out.close();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    public static String getShellInput(String question){
        String myfile = "";
        try {
            BufferedReader stdin =
                new BufferedReader(new InputStreamReader(System.in)
                                   );
            System.out.print(question);
            myfile = stdin.readLine();
        } catch (Exception e) {
            System.out.println("Error in screen reading");
        }
        return myfile;
    }

}
