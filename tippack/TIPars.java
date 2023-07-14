package tippack;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.locks.ReentrantLock;
import dr.app.tools.NexusExporter;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

public class TIPars {

    private String OUTPUT_FOLDER;
    private static boolean isMultiplePlacements = false; /// is output multiple placements
    private static boolean isFastaFile = true; /// is fasta file

    private boolean DEBUG = false; /// for debug
    private boolean OUTPUT_PNODE = false; /// is output P node sequence

    private Tree mytree = null;
    private static String internalnode_nidname = "label";
    private static HashMap<FlexibleNode, Integer> node2edge = new HashMap<FlexibleNode, Integer>();
    private static String[] placements = null;
    private static int edge_number = 1;

    private static ArrayList<String> stringSequencesList = new ArrayList<String>(); /// fasta: store taxaseq and ancseq,
                                                                                    /// ordered follow the seqIdxMap
    private static HashMap<String, Integer> seqIdxMap = new HashMap<String, Integer>(); /// sequence names (taxaseq and
                                                                                        /// ancseq) map to their reading
                                                                                        /// orders from file
    private static int sequence_character_length = -1; /// alignment length
    private static HashMap<FlexibleNode, String> node2seqName = new HashMap<FlexibleNode, String>(); /// tree_nodes map
                                                                                                     /// to their
                                                                                                     /// sequences name

    private static ArrayList<ConcurrentHashMap<Integer, Byte>> multationSequencesMap = new ArrayList<ConcurrentHashMap<Integer, Byte>>(); /// vcf:
                                                                                                                                          /// store
                                                                                                                                          /// taxaseq
                                                                                                                                          /// and
                                                                                                                                          /// ancseq,
                                                                                                                                          /// ordered
                                                                                                                                          /// follow
                                                                                                                                          /// the
                                                                                                                                          /// seqIdxMap
    private static byte[] ref_sequence = new byte[100000]; /// reference sequence string

    private double minGlobalMemoryBranchScore = Double.MAX_VALUE;
    private ArrayList<FlexibleNode> minGlobalMemoryBranchScoreNodeList = new ArrayList<FlexibleNode>();

    public static char[] alphabet_nt = { 'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N',
            '-' }; /// IUPAC nucleotide codes
    public static char[] alphabet_aa = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
            'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '-' }; /// IUPAC nucleotide codes
    public static HashMap<Byte, HashSet<Byte>> _nucleotide_nomenclature = null;
    private static double[][] _nucleotide_nomenclature_scoreTable = null;
    private static HashMap<BitSet, Byte> _nucleotide_nomenclature_map2char = null;
    private static double[][] _aminoacid_scoreTable = null;
    private static double[][] _used_scoreTable = null;

    private ReentrantLock lock = new ReentrantLock();

    public static double MinDoubleNumLimit = 10 * Double.MIN_VALUE;

    public TIPars() {
    }

    // initialization
    public TIPars(Tree mytree, String otype, String output_folder) {
        this.mytree = mytree;
        this.OUTPUT_FOLDER = output_folder;
        setupHashtableOfnode2seqName();
        setupHashtableOfNode2edge(otype);
    }

    private void setupHashtableOfNode2edge(String otype) {
        node2edge.clear();
        edge_number = 1;
        if (otype.equals("placement")) {
            for (int i = 0; i < mytree.getExternalNodeCount(); i++) {
                Integer edge = new Integer(edge_number);
                FlexibleNode node = (FlexibleNode) mytree.getExternalNode(i);
                node2edge.put(node, edge);
                edge_number++;
            }

            for (int i = 0; i < mytree.getInternalNodeCount(); i++) {
                Integer edge = new Integer(edge_number);
                FlexibleNode node = (FlexibleNode) mytree.getInternalNode(i);
                node2edge.put(node, edge);
                edge_number++;
            }
        }
    }

    private void setupHashtableOfnode2seqName() {
        node2seqName.clear();
        for (int i = 0; i < mytree.getInternalNodeCount(); i++) {
            FlexibleNode n = (FlexibleNode) mytree.getInternalNode(i);
            String sequenceName = (String) (n.getAttribute(this.internalnode_nidname));
            if (seqIdxMap.containsKey(sequenceName)) {
                node2seqName.put(n, sequenceName);
            } else {
                System.out.println("internalnode=" + sequenceName + " not found in alignment.");
            }
        }
        for (int i = 0; i < mytree.getExternalNodeCount(); i++) {
            FlexibleNode n = (FlexibleNode) mytree.getExternalNode(i);
            String sequenceName = n.getTaxon().getId();
            if (seqIdxMap.containsKey(sequenceName)) {
                node2seqName.put(n, sequenceName);
            } else {
                System.out.println("externalnode=" + sequenceName + " not found in alignment.");
            }
        }
    }

    public HashMap<String, FlexibleNode> setupHashtableOfseqName2node(Tree tree) {
        HashMap<String, FlexibleNode> mySeqName2node = new HashMap<String, FlexibleNode>();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            FlexibleNode n = (FlexibleNode) tree.getInternalNode(i);
            String sequenceName = (String) (n.getAttribute(this.internalnode_nidname));
            mySeqName2node.put(sequenceName, n);
        }
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            FlexibleNode n = (FlexibleNode) tree.getExternalNode(i);
            String sequenceName = n.getTaxon().getId();
            mySeqName2node.put(sequenceName, n);
        }
        return mySeqName2node;
    }

    // Interface API
    public String addQuerySequence(String qname, ConcurrentHashMap<Integer, Byte> nodeQseq, String qid, String pid,
            boolean printDisInfoOnScreen,
            double[] ABQ_brlen, String otype, int ii) {
        return addQuerySequence_global(qname, nodeQseq, qid, pid, printDisInfoOnScreen, ABQ_brlen, otype, ii);
    }

    // for vcf file
    public String addQuerySequence_global(String qname, ConcurrentHashMap<Integer, Byte> nodeQseq, String qid,
            String pid, boolean printDisInfoOnScreen,
            double[] ABQ_brlen, String otype, int ii) {

        minGlobalMemoryBranchScore = Double.MAX_VALUE;
        minGlobalMemoryBranchScoreNodeList.clear();
        /// if parallel computes branchscore is bigger than minGlobalMemoryBranchScore,
        /// the calculation will stop.
        ArrayList<Integer> nodeIdxAndScoreList = new ArrayList<Integer>(mytree.getNodeCount());
        //// initialize the nodeIdxAndScoreList to store nodeIdx
        for (int i = 0; i < mytree.getNodeCount(); ++i) {
            nodeIdxAndScoreList.add(i);
        }

        //// parallel to compute branchScore for every branch
        nodeIdxAndScoreList.parallelStream()
                .forEach(nodeIdx -> {
                    double score = computeBranchScore(nodeIdx, nodeQseq, qname);
                });

        ArrayList<FlexibleNode> selectedNodeList = minGlobalMemoryBranchScoreNodeList;

        //// try to remove ambiguous results
        selectedNodeList = reduceAmbiguousBranch(selectedNodeList, isMultiplePlacements);

        if (DEBUG)
            System.out.println(
                    "minQScore: " + minGlobalMemoryBranchScore + ", selectedNodeList:" + selectedNodeList.size());

        MyFlexibleTree best_mynewTree = null;
        String placementStrings = "";
        /// the first element of selectedNodeList is the best node to be inserted
        for (int k = 0; k < selectedNodeList.size(); ++k) {
            /// get sequences
            FlexibleNode selectedNode = selectedNodeList.get(k);
            int selectedNodeBIndex = selectedNode.getNumber();
            ConcurrentHashMap<Integer, Byte> nodeAseq = getVariantSequenceByNode(
                    (FlexibleNode) selectedNode.getParent());
            ConcurrentHashMap<Integer, Byte> nodeBseq = getVariantSequenceByNode(selectedNode);

            Double[] selectedScores = new Double[3];
            selectedScores[0] = (double) get_num_leaves(selectedNode.getParent());
            selectedScores[1] = (double) get_num_leaves(selectedNode);
            ConcurrentHashMap<Integer, Byte> nodePseq = getStringAndScoreFromNodeABQSeq(nodeAseq, nodeBseq, nodeQseq,
                    selectedScores);

            String selectedBnid = (String) selectedNode.getAttribute(internalnode_nidname);

            if (selectedBnid == null) {
                selectedBnid = selectedNode.getTaxon().getId();
            }

            String selectAName = (String) selectedNode.getParent().getAttribute(internalnode_nidname);

            // Q-P pendent length local estimation
            double pqlen = 0.0;
            if (selectedScores[2] <= MinDoubleNumLimit)
                pqlen = 0.0;
            else /// selectedScores[2] > Double.MIN_VALUE
            {
                double scoreAB = computeNodeScore(nodeAseq, nodeBseq);
                FlexibleNode myNodeB = selectedNode;
                FlexibleNode myNodeA = myNodeB.getParent();
                /// iteratively consider upper branch of A’s parent to A for scaling if
                /// selectedScores[2] > Double.MIN_VALUE and scoreAB <= MinDoubleNumLimit.
                while ((scoreAB <= MinDoubleNumLimit || myNodeB.getLength() <= MinDoubleNumLimit)
                        && !myNodeA.isRoot()) {
                    myNodeB = myNodeA;
                    myNodeA = myNodeB.getParent();
                    scoreAB = computeNodeScore(getVariantSequenceByNode(myNodeA), getVariantSequenceByNode(myNodeB));
                }
                if (scoreAB > MinDoubleNumLimit && myNodeB.getLength() > MinDoubleNumLimit)
                    pqlen = localEstimation(selectedScores[2], scoreAB, myNodeB.getLength());
                else {
                    double p = selectedScores[2] / ((double) getAlignmentLength());
                    pqlen = JC69(p);
                }
            }

            double original_branchAB = selectedNode.getLength();
            if (selectedScores[1] <= MinDoubleNumLimit) {
                ABQ_brlen[0] = original_branchAB;
                ABQ_brlen[1] = 0.0;
                ABQ_brlen[2] = pqlen;
            } else if (selectedScores[0] <= MinDoubleNumLimit) {
                ABQ_brlen[0] = 0.0;
                ABQ_brlen[1] = original_branchAB;
                ABQ_brlen[2] = pqlen;
            } else {
                double Pratio = selectedScores[0] / ((double) (selectedScores[0] + selectedScores[1]));
                double newNodePLength = original_branchAB * Pratio;
                double newNodeBLength = original_branchAB * (1.0 - Pratio);

                ABQ_brlen[0] = newNodePLength;
                ABQ_brlen[1] = newNodeBLength;
                ABQ_brlen[2] = pqlen;
            }

            if (printDisInfoOnScreen) {
                if (k == 0) // the best inserted branch
                    return qname + "\t" + "*" + selectAName + "-" + selectedBnid + "\t" + "ABQ_brlen: " + ABQ_brlen[0]
                            + "\t" + ABQ_brlen[1] + "\t" + ABQ_brlen[2];
                // else
                // System.out.println(qname + "\t" + selectAName + "-" + selectedBnid + "\t" +
                // "ABQ_brlen: " + ABQ_brlen[0] + "\t" + ABQ_brlen[1] + "\t" + ABQ_brlen[2]);
            }
        }
        return "error 02";
    }

    /// compute the branchscore that nodeQseq difference from both nodeA and nodeB

    //// for fasta file
    public double computeBranchScore(int nodeIdx, String nodeQseq, String qname) {
        FlexibleNode nodeB = (FlexibleNode) mytree.getNode(nodeIdx);
        return computeBranchScore(nodeB, nodeQseq, qname);
    }

    public double computeBranchScore(FlexibleNode nodeB, String nodeQseq, String qname) {
        // We use this annotation:
        // nodeA
        // |
        // nodeP
        // | \
        // | \
        // nodeB nodeQ

        if (nodeB.isRoot())
            return Double.MAX_VALUE;

        // nodeA-nodeB pairs
        FlexibleNode nodeA = (FlexibleNode) nodeB.getParent();
        if (nodeA == null) {
            System.out.println("nodeA is null, and its nodeB is " + node2seqName.get(nodeB));
        }

        String nodeAseq = getStringSequenceByNode(nodeA);
        String nodeBseq = getStringSequenceByNode(nodeB);

        if (nodeAseq == null) {
            if (DEBUG)
                System.out.println("nodeAseq cannot find seq - " + node2seqName.get(nodeA));

        } else if (nodeBseq == null) {
            if (DEBUG)
                System.out.println("nodeBseq cannot find seq - " + node2seqName.get(nodeB));
        }

        double score = 0;
        for (int i = 0; i < getAlignmentLength(); i++) {
            char a_i = nodeAseq.charAt(i);
            char b_i = nodeBseq.charAt(i);
            char c_i = nodeQseq.charAt(i);

            // get the substitution scores
            double score_ab = _used_scoreTable[(byte) a_i][(byte) b_i];
            double score_ac = _used_scoreTable[(byte) a_i][(byte) c_i];
            double score_bc = _used_scoreTable[(byte) b_i][(byte) c_i];

            if (a_i != b_i && a_i != c_i && b_i != c_i) { // ATC
                score += (score_ac + score_bc) / 2;
            } else if (a_i == b_i && b_i != c_i) { // AAT
                score += score_bc;
            }

            /// stop when score larger than current minGlobalMemoryBranchScore
            if (score - minGlobalMemoryBranchScore > MinDoubleNumLimit)
                return score;
        }

        lock.lock();
        if (minGlobalMemoryBranchScore - score > MinDoubleNumLimit) {
            minGlobalMemoryBranchScore = score;
            minGlobalMemoryBranchScoreNodeList.clear();
            minGlobalMemoryBranchScoreNodeList.add(nodeB);
        } else if (Math.abs(minGlobalMemoryBranchScore - score) <= MinDoubleNumLimit) {
            minGlobalMemoryBranchScoreNodeList.add(nodeB);
        }
        lock.unlock();

        return score;
    }

    public double computeBranchScore(String nodeAseq, String nodeBseq, String nodeQseq) {
        double score = 0;
        for (int i = 0; i < getAlignmentLength(); i++) {
            char a_i = nodeAseq.charAt(i);
            char b_i = nodeBseq.charAt(i);
            char c_i = nodeQseq.charAt(i);

            double score_ab = _used_scoreTable[(byte) a_i][(byte) b_i];
            double score_ac = _used_scoreTable[(byte) a_i][(byte) c_i];
            double score_bc = _used_scoreTable[(byte) b_i][(byte) c_i];

            if (a_i != b_i && a_i != c_i && b_i != c_i) { // ATC
                score += (score_ac + score_bc) / 2;
            } else if (a_i == b_i && b_i != c_i) { // AAT
                score += score_bc;
            }
        }

        return score;
    }

    //// compute the branch score and generate the P node string.
    public String getStringAndScoreFromNodeABQSeq(String aSeq, String bSeq, String cSeq, Double[] scores) {
        scores[0] = 0.0;
        scores[1] = 0.0;
        scores[2] = 0.0;
        StringBuilder pSeq = new StringBuilder(aSeq);

        // scores[3] is the difference between A and B
        for (int i = 0; i < getAlignmentLength(); i++) {
            char a_i = aSeq.charAt(i);
            char b_i = bSeq.charAt(i);
            char c_i = cSeq.charAt(i);

            double score_ac = _used_scoreTable[(byte) a_i][(byte) c_i];
            double score_bc = _used_scoreTable[(byte) b_i][(byte) c_i];
            double score_ab = _used_scoreTable[(byte) a_i][(byte) b_i];

            if (a_i == b_i && a_i == c_i) {
                continue;
                // do nothing
            } else if (a_i != b_i && a_i != c_i && b_i != c_i) { // ATC
                scores[2] += (score_ac + score_bc) / 2.0;
                if (score_ac > score_bc)
                    pSeq.setCharAt(i, b_i);
            } else if (a_i == b_i && b_i != c_i) { // AAT
                scores[2] += score_bc;
                pSeq.setCharAt(i, a_i);
            } else if (a_i != b_i && b_i == c_i) { // ATT
                scores[0] += score_ab;
                pSeq.setCharAt(i, b_i);
            } else if (a_i == c_i && b_i != c_i) { // TAT
                scores[1] += score_ab;
                pSeq.setCharAt(i, a_i);
            } else if (DEBUG) {
                System.out.println("Unmatched ABQ type");
            }
        }
        return pSeq.toString();
    }

    //// for vcf file
    public double computeBranchScore(int nodeIdx, ConcurrentHashMap<Integer, Byte> nodeQseq, String qname) {
        FlexibleNode nodeB = (FlexibleNode) mytree.getNode(nodeIdx);
        return computeBranchScore(nodeB, nodeQseq, qname);
    }

    public double computeBranchScore(FlexibleNode nodeB, ConcurrentHashMap<Integer, Byte> nodeQseq, String qname) {
        // We use this annotation:
        // nodeA
        // |
        // nodeP
        // | \
        // | \
        // nodeB nodeQ

        if (nodeB.isRoot())
            return Double.MAX_VALUE;

        // nodeA-nodeB pairs
        FlexibleNode nodeA = (FlexibleNode) nodeB.getParent();
        ConcurrentHashMap<Integer, Byte> nodeAseq = getVariantSequenceByNode(nodeA);
        ConcurrentHashMap<Integer, Byte> nodeBseq = getVariantSequenceByNode(nodeB);

        double score = 0;
        HashSet<Integer> mergeIdxSet = new HashSet<Integer>(nodeQseq.keySet());
        mergeIdxSet.addAll(nodeAseq.keySet());
        mergeIdxSet.addAll(nodeBseq.keySet());

        for (Integer key : mergeIdxSet) {
            byte a = nodeAseq.containsKey(key) ? nodeAseq.get(key) : ref_sequence[key];
            byte b = nodeBseq.containsKey(key) ? nodeBseq.get(key) : ref_sequence[key];
            byte c = nodeQseq.containsKey(key) ? nodeQseq.get(key) : ref_sequence[key];

            if (a != b && a != c && b != c) { // ATC
                score += (_used_scoreTable[a][c] + _used_scoreTable[b][c]) / 2;
            } else if (a == b && b != c) { // AAT
                score += _used_scoreTable[b][c];
            }

            if (score - minGlobalMemoryBranchScore > MinDoubleNumLimit)
                return score;
        }

        lock.lock();
        if (minGlobalMemoryBranchScore - score > MinDoubleNumLimit) {
            minGlobalMemoryBranchScore = score;
            minGlobalMemoryBranchScoreNodeList.clear();
            minGlobalMemoryBranchScoreNodeList.add(nodeB);
        } else if (Math.abs(minGlobalMemoryBranchScore - score) <= MinDoubleNumLimit) {
            minGlobalMemoryBranchScoreNodeList.add(nodeB);
        }
        lock.unlock();

        return score;
    }

    public double computeBranchScore(ConcurrentHashMap<Integer, Byte> nodeAseq,
            ConcurrentHashMap<Integer, Byte> nodeBseq, ConcurrentHashMap<Integer, Byte> nodeQseq) {
        double score = 0;
        HashSet<Integer> mergeIdxSet = new HashSet<Integer>(nodeQseq.keySet());
        mergeIdxSet.addAll(nodeAseq.keySet());
        mergeIdxSet.addAll(nodeBseq.keySet());

        for (Integer key : mergeIdxSet) {
            byte a = nodeAseq.containsKey(key) ? nodeAseq.get(key) : ref_sequence[key];
            byte b = nodeBseq.containsKey(key) ? nodeBseq.get(key) : ref_sequence[key];
            byte c = nodeQseq.containsKey(key) ? nodeQseq.get(key) : ref_sequence[key];

            double score_ab = _used_scoreTable[a][b];
            double score_ac = _used_scoreTable[a][c];
            double score_bc = _used_scoreTable[b][c];

            if (a != b && a != c && b != c) { // ATC
                score += (score_ac + score_bc) / 2;
            } else if (a == b && b != c) { // AAT
                score += score_bc;
            }
        }
        return score;
    }

    public ConcurrentHashMap<Integer, Byte> getStringAndScoreFromNodeABQSeq(ConcurrentHashMap<Integer, Byte> nodeAseq,
            ConcurrentHashMap<Integer, Byte> nodeBseq, ConcurrentHashMap<Integer, Byte> nodeQseq, Double[] scores) {
        scores[0] = 0.0;
        scores[1] = 0.0;
        scores[2] = 0.0;
        ConcurrentHashMap<Integer, Byte> pSeq = new ConcurrentHashMap<Integer, Byte>();

        HashSet<Integer> mergeIdxSet = new HashSet<Integer>(nodeQseq.keySet());
        mergeIdxSet.addAll(nodeAseq.keySet());
        mergeIdxSet.addAll(nodeBseq.keySet());

        for (Integer key : mergeIdxSet) {
            byte a = nodeAseq.containsKey(key) ? nodeAseq.get(key) : ref_sequence[key];
            byte b = nodeBseq.containsKey(key) ? nodeBseq.get(key) : ref_sequence[key];
            byte c = nodeQseq.containsKey(key) ? nodeQseq.get(key) : ref_sequence[key];

            double score_ab = _used_scoreTable[(byte) a][(byte) b];
            double score_ac = _used_scoreTable[(byte) a][(byte) c];
            double score_bc = _used_scoreTable[(byte) b][(byte) c];

            byte placeCharacter = a; /// default place a
            if (a == b && a == c) {
                continue;
            } else if (a != b && a != c && b != c) { // ATC
                scores[2] += (score_ac + score_bc) / 2.0;
                if (score_ac > score_bc)
                    placeCharacter = b;
            } else if (a == b && b != c) { // AAT
                scores[2] += score_bc;
            } else if (a != b && b == c) { // ATT
                scores[0] += score_ab;
                placeCharacter = b;
            } else if (a == c && b != c) { // TAT
                scores[1] += score_ab;
            }

            if (placeCharacter != ref_sequence[key]) {
                pSeq.put(key, placeCharacter);
            }
        }

        return pSeq;
    }

    public double computeNodeScore(FlexibleNode nodeB, String nodeQseq) {
        if (nodeB.isRoot())
            return Double.MAX_VALUE;

        String nodeBseq = getStringSequenceByNode(nodeB);
        double scoreNode = 0;
        for (int i = 0; i < getAlignmentLength(); i++) {
            char a_i = nodeBseq.charAt(i);
            char b_i = nodeQseq.charAt(i);
            scoreNode += _used_scoreTable[(byte) a_i][(byte) b_i];
        }

        return scoreNode;
    }

    public double computeNodeScore(String nodeBseq, String nodeQseq) {
        double score = 0;
        for (int i = 0; i < getAlignmentLength(); i++) {
            char a_i = nodeBseq.charAt(i);
            char b_i = nodeQseq.charAt(i);
            score += _used_scoreTable[(byte) a_i][(byte) b_i];
        }
        return score;
    }

    public double computeNodeScore(ConcurrentHashMap<Integer, Byte> nodeBseq,
            ConcurrentHashMap<Integer, Byte> nodeQseq) {
        double score = 0;
        HashSet<Integer> mergeIdxSet = new HashSet<Integer>(nodeQseq.keySet());
        mergeIdxSet.addAll(nodeBseq.keySet());

        for (Integer key : mergeIdxSet) {
            byte b = nodeBseq.containsKey(key) ? nodeBseq.get(key) : ref_sequence[key];
            byte c = nodeQseq.containsKey(key) ? nodeQseq.get(key) : ref_sequence[key];

            score += _used_scoreTable[b][c];
        }
        return score;
    }

    /// filter rules for multiple placements
    public ArrayList<FlexibleNode> reduceAmbiguousBranch(ArrayList<FlexibleNode> selectedNodeList) {
        return reduceAmbiguousBranch(selectedNodeList, false);
    }

    @SuppressWarnings("unchecked")
    public ArrayList<FlexibleNode> reduceAmbiguousBranch(ArrayList<FlexibleNode> selectedNodeList,
            Boolean returnMulti) {
        if (selectedNodeList.size() < 2)
            return selectedNodeList;
        ArrayList<FlexibleNode> finalSelectedNode = new ArrayList<FlexibleNode>();

        /// mapping slected nodes' parents to themselves
        HashMap<FlexibleNode, ArrayList<FlexibleNode>> parentCluster = new HashMap<FlexibleNode, ArrayList<FlexibleNode>>();
        for (int i = 0; i < selectedNodeList.size(); ++i) {
            FlexibleNode nodeB = selectedNodeList.get(i);
            FlexibleNode nodeA = (FlexibleNode) mytree.getParent(nodeB);
            if (!parentCluster.containsKey(nodeA)) {
                ArrayList<FlexibleNode> temp = new ArrayList<FlexibleNode>();
                temp.add(nodeB);
                parentCluster.put(nodeA, temp);
            } else {
                parentCluster.get(nodeA).add(nodeB);
            }
        }

        /// remove duplicate if potential placements are relationships of parent and
        /// child
        ArrayList<FlexibleNode> parentCluste_keySet = new ArrayList<FlexibleNode>(parentCluster.keySet());
        /// post transversal order
        parentCluste_keySet.sort((d1, d2) -> (d2.getNumber() - d1.getNumber()));
        for (int k = 0; k < parentCluste_keySet.size(); ++k) {
            FlexibleNode nodeA = parentCluste_keySet.get(k);
            if (parentCluster.containsKey(nodeA)) {
                ArrayList<FlexibleNode> temp = parentCluster.get(nodeA);
                for (int i = 0; i < temp.size(); ++i) {
                    /// both parent node and its child node are potential placements
                    /// select the node that contains more unique leaf nodes
                    FlexibleNode nodeB = temp.get(i);
                    if (parentCluster.containsKey(nodeB)) {
                        int numLeavesA = get_num_leaves(nodeA);
                        int numLeavesB = get_num_leaves(nodeB);
                        if (numLeavesA < 2 * numLeavesB) {
                            parentCluster.get(nodeA).remove(i);
                            if (parentCluster.get(nodeA).size() < 1)
                                parentCluster.remove(nodeA);
                        } else
                            parentCluster.remove(nodeB);
                    }
                }
            }
        }

        /// select the node contain maximum children, if the same choose minimum height
        /// branch length
        /// for the case that several potential placements share the same parent
        for (FlexibleNode nodeA : parentCluster.keySet()) {
            ArrayList<FlexibleNode> temp = parentCluster.get(nodeA);
            int maxChildrenB = 0;
            ArrayList<Integer> selectIdxs = new ArrayList<Integer>();
            for (int i = 0; i < temp.size(); ++i) {
                if (temp.get(i).getChildCount() > maxChildrenB) {
                    maxChildrenB = temp.get(i).getChildCount();
                    selectIdxs.clear();
                    selectIdxs.add(i);
                } else if (temp.get(i).getChildCount() == maxChildrenB) {
                    selectIdxs.add(i);
                }
            }
            if (selectIdxs.size() < 2)
                finalSelectedNode.add(temp.get(selectIdxs.get(0)));
            else {
                double minHeight = Double.MAX_VALUE;
                for (int i = 0; i < selectIdxs.size(); ++i) {
                    if (temp.get(selectIdxs.get(i)).getHeight() < minHeight) {
                        minHeight = temp.get(selectIdxs.get(i)).getHeight();
                    }
                }
                for (int i = 0; i < selectIdxs.size(); ++i) {
                    if (Math.abs(minHeight - temp.get(selectIdxs.get(i)).getHeight()) <= MinDoubleNumLimit)
                        finalSelectedNode.add(temp.get(selectIdxs.get(i)));
                }
            }
        }

        if (finalSelectedNode.size() < 2)
            return finalSelectedNode;
        selectedNodeList = (ArrayList<FlexibleNode>) finalSelectedNode.clone();
        finalSelectedNode.clear();

        /// return the best placement
        ArrayList<FlexibleNode> bestBranch = reduceBestBranch(selectedNodeList, parentCluster);
        if (returnMulti) {
            finalSelectedNode.add(bestBranch.get(0));
            // if multiple placements are allowed to print out, add all others
            for (int i = 0; i < selectedNodeList.size(); ++i) {
                if (selectedNodeList.get(i).getNumber() != bestBranch.get(0).getNumber()) {
                    finalSelectedNode.add(selectedNodeList.get(i));
                }
            }
            return finalSelectedNode;
        } else
            return bestBranch;
    }

    /// filter different real multiple placements
    public ArrayList<FlexibleNode> reduceBestBranch(ArrayList<FlexibleNode> selectedNodeList,
            HashMap<FlexibleNode, ArrayList<FlexibleNode>> parentCluster) {
        if (selectedNodeList.size() < 2)
            return selectedNodeList;
        ArrayList<FlexibleNode> finalSelectedNode = new ArrayList<FlexibleNode>();

        //// select best answer ammong multiple placements
        /// select the maximum children then minimum height branch length
        int maxChildren = 0;
        for (int i = 0; i < selectedNodeList.size(); ++i) {
            FlexibleNode nodeB = selectedNodeList.get(i);
            FlexibleNode nodeA = (FlexibleNode) mytree.getParent(nodeB);

            int total = nodeA.getChildCount();
            if (parentCluster.containsKey(nodeB))
                total = total - 1; // unique children

            if (total > maxChildren) {
                maxChildren = total;
                finalSelectedNode.clear();
                finalSelectedNode.add(nodeB);
            } else if (maxChildren == total) {
                finalSelectedNode.add(nodeB);
            }
        }

        if (finalSelectedNode.size() < 2)
            return finalSelectedNode;
        selectedNodeList = new ArrayList<FlexibleNode>(finalSelectedNode);
        finalSelectedNode.clear();

        double minHeight = Double.MAX_VALUE;
        for (int i = 0; i < selectedNodeList.size(); ++i) {
            FlexibleNode nodeB = selectedNodeList.get(i);
            FlexibleNode nodeA = (FlexibleNode) mytree.getParent(nodeB);

            double total = nodeA.getHeight();

            if (minHeight - total > MinDoubleNumLimit) {
                minHeight = total;
                finalSelectedNode.clear();
                finalSelectedNode.add(nodeB);
            } else if (Math.abs(minHeight - total) <= MinDoubleNumLimit) {
                finalSelectedNode.add(nodeB);
            }
        }

        if (finalSelectedNode.size() > 1) {
            int randomIdx = (int) (Math.random() * finalSelectedNode.size());
            selectedNodeList = (ArrayList<FlexibleNode>) finalSelectedNode.clone();
            finalSelectedNode.clear();
            finalSelectedNode.add(selectedNodeList.get(randomIdx)); // random select
        }

        return finalSelectedNode;
    }

    public int get_num_leaves(FlexibleNode node) {
        if (!node.hasChildren())
            return 1;
        int num_leaves = 0;
        for (int i = 0; i < node.getChildCount(); ++i) {
            num_leaves += get_num_leaves(node.getChild(i));
        }
        return num_leaves;
    }

    class MyFlexibleTree extends FlexibleTree {
        MyFlexibleTree(Tree t) {
            super(t);
        }

        MyFlexibleTree(Tree t, boolean keepAttribute) {
            super(t, keepAttribute);
        }

        /**
         * Adopt a node hierarchy as its own. Only called by the
         * FlexibleTree(FlexibleNode, TaxonList).
         * This creates the node list and stores the nodes in post-traversal order.
         */
        public void toAdoptNodes(FlexibleNode n) {
            super.adoptNodes(n);
        }

    }

    // JC69 model
    public static double JC69(double p) {
        double d = -3.0 * Math.log(1 - 4.0 * p / 3.0) / 4.0;
        return d;
    }

    // K2P model
    public static double K2P(String nodePseq, String nodeQseq) {
        int S = 0;
        int V = 0;
        int n = nodePseq.length();

        boolean precedingGapQ = true;

        for (int i = 0; i < n; i++) {
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

        double s = 1.0 * S / (1.0 * n);
        double v = 1.0 * V / (1.0 * n);
        double d = -1.0 * Math.log(1.0 - 2.0 * s - v) / 2.0 - 1.0 * Math.log(1.0 - 2.0 * v) / 4.0;
        // System.out.println(s + "\t" + v + "\t" + d);
        return d;
    }

    // local scale model
    public static double localEstimation(Integer[] selectedScores_int, double ABbranch) {
        if (selectedScores_int[0] + selectedScores_int[1] == 0) {
            return 0;
        }

        double d = ABbranch * ((double) selectedScores_int[2])
                / ((double) (selectedScores_int[0] + selectedScores_int[1]));

        return d;
    }

    public static double localEstimation(Double[] selectedScores, double ABbranch) {
        if (selectedScores[0] + selectedScores[1] == 0) {
            return 0;
        }

        double d = ABbranch * selectedScores[2] / (selectedScores[0] + selectedScores[1]);

        return d;
    }

    public static double localEstimation(double score_c2ab, double score_ab, double ABbranch) {
        if (score_c2ab <= MinDoubleNumLimit) {
            return 0;
        }

        double d = ABbranch * score_c2ab / score_ab;

        return d;
    }

    // return common ancestor
    public static FlexibleNodeBranch returnCommonAncestor(FlexibleNode a, FlexibleNode b, Tree tree) {
        FlexibleNodeBranch nodeAndBranch = new TIPars().new FlexibleNodeBranch();
        if (a == b) {
            nodeAndBranch.a = a;
            nodeAndBranch.b = 0.0;
            return nodeAndBranch;
        } else {
            HashMap<FlexibleNode, Double> bAncestors = new HashMap<FlexibleNode, Double>();
            FlexibleNode parent = b;
            double dist = 0.0;
            bAncestors.put(b, dist);
            while (parent != null && !parent.isRoot()) {
                dist = dist + parent.getLength();
                bAncestors.put(parent.getParent(), dist);
                parent = parent.getParent();
            }
            parent = a;
            dist = 0.0;
            while (!bAncestors.containsKey(parent) && !tree.isRoot(parent)) {
                dist = dist + parent.getLength();
                parent = parent.getParent();
            }
            nodeAndBranch.a = parent;
            nodeAndBranch.b = dist;
            if (bAncestors.containsKey(parent))
                nodeAndBranch.b += bAncestors.get(parent);
            return nodeAndBranch;
        }
    }

    // remove single taxon
    public static FlexibleNodeBranch removeTaxon(FlexibleTree tree, FlexibleNode node) {
        try {
            MyFlexibleTree t = new TIPars().new MyFlexibleTree(tree, true);
            copyAttributeFromOneNodeToAnother((FlexibleNode) tree.getRoot(), (FlexibleNode) t.getRoot());

            t.beginTreeEdit();
            FlexibleNode n = (FlexibleNode) t.getNode(node.getNumber());
            FlexibleNode p = n.getParent();

            FlexibleNode pReturn = (FlexibleNode) tree.getNode(p.getNumber());
            FlexibleNode aReturn = (FlexibleNode) (pReturn.getParent());
            double true_dist1 = p.getLength() + n.getLength();
            double true_dist2 = n.getLength();

            p.removeChild(n);
            n.setParent(null);
            if (p.getChildCount() == 1) { // Remove this p node if if has less than 2 child nodes
                if (!p.isRoot()) {
                    FlexibleNode b = (FlexibleNode) (p.getChild(0));
                    FlexibleNode a = (FlexibleNode) (p.getParent());

                    a.addChild(b);
                    a.removeChild(p);
                    b.setParent(a);
                    double oldb2p = b.getLength();
                    b.setLength(b.getLength() + p.getLength());
                    p.setParent(null);
                    t.endTreeEdit();
                    FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(aReturn, true_dist1);
                    result.tree = t;

                    t.toAdoptNodes((FlexibleNode) t.getRoot());
                    return result;
                } else {
                    FlexibleNode b = (FlexibleNode) (p.getChild(0));
                    for (int i = 0; i < b.getChildCount(); ++i) {
                        FlexibleNode c = (FlexibleNode) (b.getChild(i));
                        c.setLength(c.getLength() + b.getLength());
                        p.addChild(c);
                        c.setParent(p);
                    }
                    p.removeChild(b);
                    b.setParent(null);
                    t.endTreeEdit();
                    FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(pReturn, true_dist1);
                    result.tree = t;

                    t.toAdoptNodes((FlexibleNode) t.getRoot());
                    return result;
                }
            } else {
                t.endTreeEdit();
                FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(pReturn, true_dist2);
                result.tree = t;

                t.toAdoptNodes((FlexibleNode) t.getRoot());
                return result;// ? p, p.getLentth ?
            }
        } catch (Exception e) {
            e.printStackTrace();

            return null;
        }
    }

    // remove a list of taxa
    public static FlexibleNodeBranchList removeTaxon(FlexibleTree tree, ArrayList<FlexibleNode> nodes) {
        try {
            MyFlexibleTree t = new TIPars().new MyFlexibleTree(tree, true);
            copyAttributeFromOneNodeToAnother((FlexibleNode) tree.getRoot(), (FlexibleNode) t.getRoot());

            ArrayList<FlexibleNode> removeNodes = new ArrayList<FlexibleNode>();
            for (int i = 0; i < nodes.size(); ++i)
                removeNodes.add((FlexibleNode) t.getNode(nodes.get(i).getNumber()));

            FlexibleNodeBranchList flexibleNodeBranchList = new TIPars().new FlexibleNodeBranchList();
            flexibleNodeBranchList.nodeBranchList = new ArrayList<FlexibleNodeBranch>();

            t.beginTreeEdit();

            for (FlexibleNode node : removeNodes) {

                FlexibleNode n = node; // (FlexibleNode) t.getNode(node.getNumber());
                FlexibleNode p = n.getParent();
                System.out.println(p.getNumber());

                FlexibleNode pReturn = (FlexibleNode) tree.getNode(p.getNumber());
                String pName = (String) p.getAttribute(internalnode_nidname);
                String pReturnName = (String) pReturn.getAttribute(internalnode_nidname);
                if (!pName.equals(pReturnName)) {
                    System.out.println("removeTaxon: " + pName + "|" + pReturnName);
                    return null;
                }

                FlexibleNode aReturn = (FlexibleNode) (pReturn.getParent());
                double true_dist1 = p.getLength() + n.getLength();
                double true_dist2 = n.getLength();

                p.removeChild(n);
                n.setParent(null);
                if (p.getChildCount() == 1) { // Remove this p node if if has less than 2 child nodes
                    FlexibleNode b = (FlexibleNode) (p.getChild(0));
                    FlexibleNode a = (FlexibleNode) (p.getParent());

                    a.addChild(b);
                    b.setParent(a);
                    double oldb2p = b.getLength();
                    b.setLength(b.getLength() + p.getLength());

                    a.removeChild(p);
                    p.setParent(null);

                    FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(aReturn, true_dist1);
                    result.tree = null;
                    flexibleNodeBranchList.nodeBranchList.add(result);
                } else {

                    FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(pReturn, true_dist2);
                    result.tree = null;
                    flexibleNodeBranchList.nodeBranchList.add(result);
                }
            }
            t.endTreeEdit();
            t.toAdoptNodes((FlexibleNode) t.getRoot());
            flexibleNodeBranchList.tree = t;
            return flexibleNodeBranchList;
        } catch (Exception e) {
            e.printStackTrace();

            return null;
        }
    }

    // remove a list of taxa and return a resulting tree
    public static MyFlexibleTree removeTaxonReturnTree(FlexibleTree tree, ArrayList<FlexibleNode> nodes) {
        try {
            MyFlexibleTree t = new TIPars().new MyFlexibleTree(tree, true);
            copyAttributeFromOneNodeToAnother((FlexibleNode) tree.getRoot(), (FlexibleNode) t.getRoot());

            ArrayList<FlexibleNode> removeNodes = new ArrayList<FlexibleNode>();
            for (int i = 0; i < nodes.size(); ++i) {
                removeNodes.add((FlexibleNode) t.getNode(nodes.get(i).getNumber()));
                String sequenceName1 = nodes.get(i).getTaxon().getId();
                String sequenceName2 = ((FlexibleNode) (t.getNode(nodes.get(i).getNumber()))).getTaxon().getId();
                if (!sequenceName1.equals(sequenceName2)) {
                    System.out.println(sequenceName1 + " != " + sequenceName2);
                    return null;
                }
                if (((FlexibleNode) (t.getNode(nodes.get(i).getNumber()))).isRoot()) {
                    System.out.println(sequenceName2 + " original is root");
                    return null;
                }
            }

            FlexibleNodeBranchList flexibleNodeBranchList = new TIPars().new FlexibleNodeBranchList();
            flexibleNodeBranchList.nodeBranchList = new ArrayList<FlexibleNodeBranch>();

            t.beginTreeEdit();
            for (FlexibleNode n : removeNodes) {
                if (n.isRoot()) {
                    System.out.println(n.getTaxon().getId() + " is root");
                    if (n.getParent() != null)
                        System.out.println("root but not null");
                    return null;
                }
                FlexibleNode p = n.getParent();
                if (p == null) {
                    System.out.println(n.getTaxon().getId() + " no parent");
                    return null;
                }
                if (p.getChildCount() > 2) {
                    p.removeChild(n);
                } else if (p.getChildCount() == 2) {
                    for (int j = 0; j < 2; ++j) {
                        if (!p.getChild(j).equals(n)) {
                            if (p.isRoot()) {
                                System.out
                                        .println(n.getTaxon().getId() + "'s parent is root and has only two children");
                                return null;
                            }
                            FlexibleNode a = (FlexibleNode) (p.getParent());
                            FlexibleNode b = (FlexibleNode) (p.getChild(j));
                            double newbl = b.getLength() + p.getLength();
                            b.setLength(newbl);
                            a.addChild(b);
                            a.removeChild(p);
                        }
                    }
                } else {
                    System.out.println(n.getTaxon().getId() + " has no sister");
                    return null;
                }
            }
            t.endTreeEdit();
            t.toAdoptNodes((FlexibleNode) t.getRoot());
            return t;
        } catch (Exception e) {
            e.printStackTrace();

            return null;
        }
    }

    // copy a node
    private static void copyAttributeFromOneNodeToAnother(FlexibleNode n, FlexibleNode m) {
        Iterator attnames = n.getAttributeNames();
        while (attnames != null && attnames.hasNext()) {
            String ahname = (String) attnames.next();
            m.setAttribute(ahname, n.getAttribute(ahname));
        }
    }

    // get string sequence
    private static String getStringSequenceByNode(FlexibleNode a) {
        if (node2seqName.containsKey(a)) {
            if (!seqIdxMap.containsKey((node2seqName.get(a)))) {
                System.out.println("seqIdxMap can not access " + (node2seqName.get(a)));
            }
            String seq = stringSequencesList.get(seqIdxMap.get((node2seqName.get(a))));
            return seq;
        } else {
            String selectedBnid = (String) a.getAttribute(internalnode_nidname);
            if (selectedBnid == null) {
                selectedBnid = a.getTaxon().getId();
            }
            System.out.println("node2seqName can not access " + a.getNumber() + "/" + selectedBnid);
            System.exit(-1);
            return null;
        }
    }

    // get variants sequence
    private static ConcurrentHashMap<Integer, Byte> getVariantSequenceByNode(FlexibleNode a) {
        ConcurrentHashMap<Integer, Byte> seq = multationSequencesMap.get(seqIdxMap.get((node2seqName.get(a))));
        return seq;
    }

    public class FlexibleNodeBranch {
        public FlexibleNode a;
        public Double b;
        public Tree tree;
        public String nodeAName;

        public FlexibleNodeBranch(FlexibleNode a, Double b) {
            this.a = a;
            this.b = b;
        }

        public FlexibleNodeBranch() {
            this.a = null;
            this.b = null;
        }

        public void setNode(FlexibleNode a2) {
            this.a = a2;
        }

        public void setBrlen(Double b2) {
            this.b = b2;
        }

        public FlexibleNode getNode() {
            return this.a;
        }

        public Double getBrlen() {
            return this.b;
        }
    }

    public static int getAlignmentLength() {
        return sequence_character_length;
    }

    // read fasta file to a hashmap
    public static HashMap<Integer, String> readFastaFile2Alignment(String fn) {
        HashMap<Integer, String> alignmnetIdxList = new HashMap<Integer, String>();
        int startIndex = stringSequencesList.size();
        try {
            BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
            String fasline = br2.readLine().trim();
            String desc;
            while (fasline != null) {
                StringBuffer sequence = new StringBuffer();
                if (fasline.matches("^>.+")) {
                    desc = fasline;
                    fasline = br2.readLine().trim();
                    while (fasline != null && !fasline.matches("^>.+")) {
                        sequence.append(fasline);
                        fasline = br2.readLine();
                    }
                    desc = desc.replaceAll(">", "");
                    String seqseq = sequence.toString();
                    seqIdxMap.put(desc, startIndex);
                    stringSequencesList.add(seqseq.toUpperCase());
                    alignmnetIdxList.put(startIndex, desc);
                    startIndex++;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.gc();
        return alignmnetIdxList;
    }

    // read fasta file
    public static void readFastaAlignmentFile(String fn) {
        int startIndex = stringSequencesList.size();
        try {
            BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
            String fasline = br2.readLine().trim();
            String desc;
            while (fasline != null) {
                StringBuffer sequence = new StringBuffer();
                if (fasline.matches("^>.+")) {
                    desc = fasline;
                    fasline = br2.readLine().trim();
                    while (fasline != null && !fasline.matches("^>.+")) {
                        sequence.append(fasline);
                        fasline = br2.readLine();
                    }
                    desc = desc.replaceAll(">", "");
                    String seqseq = sequence.toString();
                    if (sequence_character_length < 0)
                        sequence_character_length = seqseq.length();
                    seqIdxMap.put(desc, startIndex++);
                    stringSequencesList.add(seqseq.toUpperCase());
                }
            }
        } catch (Exception e) {
            e.printStackTrace();

        }
        System.gc();
    }

    // read vcf file
    public static void readVCFAlignmentFile(String fn) {
        int startIndex = multationSequencesMap.size();
        boolean header_found = false;
        try {
            BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
            //// read header
            while (!header_found) {
                String vcfline = br2.readLine().trim();
                if (vcfline == null || vcfline.isEmpty())
                    continue;
                String[] words = vcfline.trim().split("\\s+");
                if (words.length > 2) {
                    if (words[1].contains("POS")) {
                        // Sample names start from the 10th word in the header
                        for (int j = 9; j < words.length; j++) {
                            seqIdxMap.put(words[j], startIndex + j - 9);
                            multationSequencesMap.add(new ConcurrentHashMap<Integer, Byte>());
                        }
                        header_found = true;
                    }
                }
            }

            if (!header_found) {
                System.out.println("Error! Incorrect VCF file. No Header");
                System.exit(-1);
            }
            int numberofThreads = Runtime.getRuntime().availableProcessors();
            // System.out.println("numberofThreads = " + numberofThreads);
            ArrayList<Integer> nodeIdxList = new ArrayList<Integer>(numberofThreads);
            for (int i = 0; i < numberofThreads; ++i)
                nodeIdxList.add(i);

            nodeIdxList.parallelStream()
                    .forEach(nodeIdx -> {
                        String vcfline = null;
                        try {
                            while ((vcfline = br2.readLine()) != null) {
                                if (vcfline.isEmpty())
                                    continue;
                                String[] words = vcfline.trim().split("\\s+");
                                if (words.length != 9 + multationSequencesMap.size() - startIndex) {
                                    System.out.println("Error! Incorrect VCF file.");
                                    System.exit(-1);
                                }
                                int variant_pos = Integer.parseInt(words[1]);
                                byte ref_nuc = (byte) words[3].charAt(0);
                                if (ref_sequence[variant_pos] == (byte) '\0') {
                                    ref_sequence[variant_pos] = ref_nuc;
                                } else if (ref_nuc != ref_sequence[variant_pos]) {
                                    System.out.println("Error! Different referen sequence.");
                                    System.exit(-1);
                                }

                                String[] alleles = words[4].split(",");
                                for (int j = 9; j < words.length; j++) {
                                    if (Character.isDigit(words[j].charAt(0))) {
                                        int allele_id = Integer.parseInt(words[j]);
                                        if (allele_id > 0) {
                                            byte allele = (byte) alleles[allele_id - 1].charAt(0);
                                            multationSequencesMap.get(startIndex + j - 9).put(variant_pos, allele);
                                        }
                                    }
                                }
                            }
                        } catch (NumberFormatException | IOException e) {
                            // TODO Auto-generated catch block
                            e.printStackTrace();
                        }
                    });
            br2.close();
            for (int i = 0; i < ref_sequence.length; ++i)
                if (ref_sequence[i] != (byte) '\0')
                    sequence_character_length = i;
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.gc();
    }

    //// read vcf file to a hashmap
    public static HashMap<Integer, String> readVCFFile2Alignment(InputStream input) {
        HashMap<Integer, String> alignmnetIdxList = new HashMap<Integer, String>();
        int startIndex = multationSequencesMap.size();
        boolean header_found = false;
        try {
            BufferedReader br2 = new BufferedReader(new InputStreamReader(input));

            //// read header
            while (!header_found) {
                String vcfline = br2.readLine().trim();
                if (vcfline == null || vcfline.isEmpty())
                    continue;
                String[] words = vcfline.trim().split("\\s+");
                if (words.length > 2) {
                    if (words[1].contains("POS")) {
                        // Sample names start from the 10th word in the header
                        for (int j = 9; j < words.length; j++) {
                            seqIdxMap.put(words[j], startIndex + j - 9);
                            multationSequencesMap.add(new ConcurrentHashMap<Integer, Byte>());
                            alignmnetIdxList.put(startIndex + j - 9, words[j]);
                        }
                        header_found = true;
                    }
                }
            }

            if (!header_found) {
                System.out.println("Error! Incorrect VCF file. No Header");
                System.exit(-1);
            }
            int numberofThreads = Runtime.getRuntime().availableProcessors();
            ArrayList<Integer> nodeIdxList = new ArrayList<Integer>(numberofThreads);
            for (int i = 0; i < numberofThreads; ++i)
                nodeIdxList.add(i);

            nodeIdxList.parallelStream()
                    .forEach(nodeIdx -> {
                        String vcfline = null;
                        try {
                            while ((vcfline = br2.readLine()) != null) {
                                if (vcfline == null)
                                    break;
                                if (vcfline.isEmpty())
                                    continue;
                                String[] words = vcfline.trim().split("\\s+");
                                if (words.length > 1) {
                                    if (words.length != 9 + multationSequencesMap.size() - startIndex) {
                                        System.out.println("Error! Incorrect VCF file.");
                                        System.exit(-1);
                                    }
                                    int variant_pos = Integer.parseInt(words[1]);
                                    byte ref_nuc = (byte) words[3].charAt(0);
                                    if (ref_sequence[variant_pos] == (byte) '\0') {
                                        ref_sequence[variant_pos] = ref_nuc;
                                    } else if (ref_nuc != ref_sequence[variant_pos]) {
                                        System.out.println("Error! Different referen sequence.");
                                        System.exit(-1);
                                    }

                                    String[] alleles = words[4].split(",");
                                    for (int j = 9; j < words.length; j++) {
                                        if (Character.isDigit(words[j].charAt(0))) {
                                            int allele_id = Integer.parseInt(words[j]);
                                            if (allele_id > 0) {
                                                byte allele = (byte) alleles[allele_id - 1].charAt(0);
                                                multationSequencesMap.get(startIndex + j - 9).put(variant_pos, allele);
                                            }
                                        }
                                    }
                                }
                            }
                        } catch (NumberFormatException e) {
                            // TODO Auto-generated catch block
                            e.printStackTrace();
                        } catch (IOException e) {
                            // TODO Auto-generated catch block
                            e.printStackTrace();
                        }
                    });
            br2.close();
            for (int i = 0; i < ref_sequence.length; ++i)
                if (ref_sequence[i] != (byte) '\0')
                    sequence_character_length = i;
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.gc();
        return alignmnetIdxList;
    }

    public static String getFolder(String filename) {
        String dir = filename.substring(0, filename.lastIndexOf("/") + 1);
        return dir;
    }

    // write nwk tree
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
                for (int i = 0; i < placements.length; i++) {
                    buffer.append(placements[i]);
                    if (i == placements.length - 1) {
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

        } catch (Exception e) {
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
            for (int i = 0; i < n; i++) {
                toNewick(tree, (FlexibleNode) tree.getChild(node, i), buffer, otype);
                if (i == (n - 1)) {
                    buffer.append(')');
                    String label = (String) node.getAttribute(internalnode_nidname);
                    if (node.isExternal())
                        label = node.getTaxon().getId();
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
        }
        ;
    }

    private static void appendLength(Tree tree, FlexibleNode node, StringBuilder buffer) {
        if (tree.hasBranchLengths()) {
            buffer.append(':');
            buffer.append(tree.getBranchLength(node));
        }
    }

    // write nexus tree
    public static void writeNexusTree(Tree t, String fn) {
        try {
            PrintStream fw = new PrintStream(new FileOutputStream(new File(fn)));
            NexusExporter kne = new NexusExporter(fw); // export the tree *with attributes* to the output nexus file.
            Tree[] ts = new Tree[1];
            ts[0] = t;
            kne.exportTrees(ts, true);
            fw.close();
        } catch (Exception e) {
            e.printStackTrace();

        }
    }

    // write fasta
    public static void writeFASTA(String desc, String seq, String output_folder) {
        String output = output_folder + desc + ".fas";
        StringBuilder buffer = new StringBuilder();
        buffer.append(">");
        buffer.append(desc);
        buffer.append("\n");
        buffer.append(seq);
        try {
            PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
            out.println(buffer.toString());
            out.close();

        } catch (Exception e) {
            e.printStackTrace();

        }
    }

    public static void writeFASTAwithFullPath(String desc, String seq, String output_folder) {
        String output = output_folder;
        StringBuilder buffer = new StringBuilder();
        buffer.append(">");
        buffer.append(desc);
        buffer.append("\n");
        buffer.append(seq);
        try {
            PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
            out.println(buffer.toString());
            out.close();

        } catch (Exception e) {
            e.printStackTrace();

        }
    }

    public static void writeFASTA(ArrayList<String> desc, ArrayList<String> seq, String output_folder) {
        String output = output_folder;
        StringBuilder buffer = new StringBuilder();
        for (int i = 0; i < desc.size(); ++i) {
            buffer.append(">");
            buffer.append(desc.get(i));
            buffer.append("\n");
            buffer.append(seq.get(i));
            buffer.append("\n");
        }
        try {
            PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
            out.print(buffer.toString());
            out.close();

        } catch (Exception e) {
            e.printStackTrace();

        }
    }

    public static void writeFASTA(String[] desc, String[] seq, String output_folder) {
        String output = output_folder;
        StringBuilder buffer = new StringBuilder();
        for (int i = 0; i < desc.length; ++i) {
            buffer.append(">");
            buffer.append(desc[i]);
            buffer.append("\n");
            buffer.append(seq[i]);
            buffer.append("\n");
        }
        try {
            PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
            out.print(buffer.toString());
            out.close();

        } catch (Exception e) {
            e.printStackTrace();

        }
    }

    // write variants
    public static void writeMutation(String desc, ConcurrentHashMap<Integer, Byte> seqMutaion, String output_folder) {
        String output = output_folder + desc + ".fas";
        StringBuilder buffer = new StringBuilder();
        buffer.append(desc);
        buffer.append("\n");
        for (ConcurrentHashMap.Entry<Integer, Byte> entry : seqMutaion.entrySet()) {
            System.out.println(entry.getKey());
            System.out.println(entry.getValue());
            buffer.append(ref_sequence[entry.getKey()] + "|" + entry.getKey() + "|" + entry.getValue() + "\n");
        }
        try {
            PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
            out.println(buffer.toString());
            out.close();

        } catch (Exception e) {
            e.printStackTrace();

        }
    }

    public class FlexibleNodeBranchList {
        public ArrayList<FlexibleNodeBranch> nodeBranchList = new ArrayList<FlexibleNodeBranch>();
        public Tree tree;
    }

    public static boolean isNumeric(String str) {
        for (int i = 0; i < str.length(); i++) {
            // System.out.println(str.charAt(i));
            if (!Character.isDigit(str.charAt(i))) {
                return false;
            }
        }
        return true;
    }

    // generate IUPAC nucleotide codes substitution table
    public static double[][] generate_nucleotide_nomenclature_scoreTable() {
        int size = 0;
        for (int i = 0; i < alphabet_nt.length; ++i)
            if ((byte) alphabet_nt[i] > size)
                size = (byte) alphabet_nt[i];
        size = size + 1;
        double[][] NodeScoreBigTable = new double[size][size];
        for (int i = 0; i < size; ++i)
            Arrays.fill(NodeScoreBigTable[0], 0);

        // {'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N',
        // '-'};

        for (int i = 0; i < alphabet_nt.length - 1; ++i) {
            for (int j = i + 1; j < alphabet_nt.length - 1; ++j) {
                int counts = _nucleotide_nomenclature.get((byte) alphabet_nt[i]).size()
                        * _nucleotide_nomenclature.get((byte) alphabet_nt[j]).size();
                int commoncounts = 0;
                // if(counts == 1)
                // {
                for (byte key : _nucleotide_nomenclature.get((byte) alphabet_nt[i])) {
                    if (_nucleotide_nomenclature.get((byte) alphabet_nt[j]).contains(key))
                        commoncounts++;
                }
                if (commoncounts < 1)
                    NodeScoreBigTable[alphabet_nt[i]][alphabet_nt[j]] = 1.0; // - commoncounts / (double)counts;
                // }
            }
        }

        for (int i = 0; i < alphabet_nt.length - 1; ++i) {
            for (int j = 0; j < i; ++j) {
                NodeScoreBigTable[alphabet_nt[i]][alphabet_nt[j]] = NodeScoreBigTable[alphabet_nt[j]][alphabet_nt[i]];
            }
        }

        // for(int i=0; i<alphabet_nt.length; ++i)
        // {
        // for(int j=0; j<alphabet_nt.length; ++j)
        // {
        // System.out.print(NodeScoreBigTable[alphabet_nt[i]][alphabet_nt[j]] + "\t");
        // }
        // System.out.print("\n");
        // }

        return NodeScoreBigTable;
    }

    //// A: 0001
    //// T: 0010
    //// C: 0100
    //// G: 1000
    public static BitSet encodeNucleotide(byte b_i, BitSet nucleotides) {
        for (byte key : _nucleotide_nomenclature.get((byte) b_i)) {
            switch (key) {
                case 'A':
                    nucleotides.set(3);
                    break;
                case 'T':
                    nucleotides.set(2);
                    break;
                case 'C':
                    nucleotides.set(1);
                    break;
                case 'G':
                    nucleotides.set(0);
                    break;
                default:
            }
        }
        return nucleotides;
    }

    public static BitSet encodeNucleotide(byte nucleotide) {
        BitSet nucleotides = new BitSet(4);
        for (byte key : _nucleotide_nomenclature.get(nucleotide)) {
            switch (key) {
                case 'A':
                    nucleotides.set(3);
                    break;
                case 'T':
                    nucleotides.set(2);
                    break;
                case 'C':
                    nucleotides.set(1);
                    break;
                case 'G':
                    nucleotides.set(0);
                    break;
                default:
            }
        }
        return nucleotides;
    }

    public static int NucleotideBitSet2Int(BitSet nucleotides) {
        int nucleotides_int = 0;
        for (int i = 0; i < 4; ++i) {
            if (nucleotides.get(i))
                nucleotides_int += Math.pow(2, i);
        }
        return nucleotides_int;
    }

    public static HashMap<BitSet, Byte> generate_nucleotide_nomenclature_characterTable() {
        HashMap<BitSet, Byte> nucleotide_nomenclature_map2char = new HashMap<BitSet, Byte>();
        Byte[] nucleotide_nomenclature_array = { 'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D',
                'N' };
        for (int i = 0; i < nucleotide_nomenclature_array.length; ++i) {
            Byte nucleotide = nucleotide_nomenclature_array[i];
            BitSet bitSet = encodeNucleotide(nucleotide);
            nucleotide_nomenclature_map2char.put(bitSet, nucleotide);
        }
        return nucleotide_nomenclature_map2char;
    }

    // IUPAC nucleotide codes mapping
    public static HashMap<Byte, HashSet<Byte>> generate_nucleotide_nomenclature() {
        HashMap<Byte, HashSet<Byte>> nucleotide_nomenclature = new HashMap<Byte, HashSet<Byte>>();
        // A
        HashSet<Byte> tempA = new HashSet<Byte>();
        tempA.add((byte) 'A');
        nucleotide_nomenclature.put((byte) 'A', tempA);
        // C
        HashSet<Byte> tempC = new HashSet<Byte>();
        tempC.add((byte) 'C');
        nucleotide_nomenclature.put((byte) 'C', tempC);
        // G
        HashSet<Byte> tempG = new HashSet<Byte>();
        tempG.add((byte) 'G');
        nucleotide_nomenclature.put((byte) 'G', tempG);
        // T
        HashSet<Byte> tempT = new HashSet<Byte>();
        tempT.add((byte) 'T');
        nucleotide_nomenclature.put((byte) 'T', tempT);
        // R
        HashSet<Byte> tempR = new HashSet<Byte>();
        tempR.add((byte) 'A');
        tempR.add((byte) 'G');
        nucleotide_nomenclature.put((byte) 'R', tempR);
        // R
        HashSet<Byte> tempY = new HashSet<Byte>();
        tempY.add((byte) 'C');
        tempY.add((byte) 'T');
        nucleotide_nomenclature.put((byte) 'Y', tempY);
        // M
        HashSet<Byte> tempM = new HashSet<Byte>();
        tempM.add((byte) 'A');
        tempM.add((byte) 'C');
        nucleotide_nomenclature.put((byte) 'M', tempM);
        // K
        HashSet<Byte> tempK = new HashSet<Byte>();
        tempK.add((byte) 'G');
        tempK.add((byte) 'T');
        nucleotide_nomenclature.put((byte) 'K', tempK);
        // S
        HashSet<Byte> tempS = new HashSet<Byte>();
        tempS.add((byte) 'C');
        tempS.add((byte) 'G');
        nucleotide_nomenclature.put((byte) 'S', tempS);
        // W
        HashSet<Byte> tempW = new HashSet<Byte>();
        tempW.add((byte) 'A');
        tempW.add((byte) 'T');
        nucleotide_nomenclature.put((byte) 'W', tempW);
        // H
        HashSet<Byte> tempH = new HashSet<Byte>();
        tempH.add((byte) 'A');
        tempH.add((byte) 'C');
        tempH.add((byte) 'T');
        nucleotide_nomenclature.put((byte) 'H', tempH);
        // B
        HashSet<Byte> tempB = new HashSet<Byte>();
        tempB.add((byte) 'C');
        tempB.add((byte) 'G');
        tempB.add((byte) 'T');
        nucleotide_nomenclature.put((byte) 'B', tempB);
        // V
        HashSet<Byte> tempV = new HashSet<Byte>();
        tempV.add((byte) 'A');
        tempV.add((byte) 'C');
        tempV.add((byte) 'G');
        nucleotide_nomenclature.put((byte) 'V', tempV);
        // D
        HashSet<Byte> tempD = new HashSet<Byte>();
        tempD.add((byte) 'A');
        tempD.add((byte) 'G');
        tempD.add((byte) 'T');
        nucleotide_nomenclature.put((byte) 'D', tempD);
        // N
        HashSet<Byte> tempN = new HashSet<Byte>();
        tempN.add((byte) 'A');
        tempN.add((byte) 'C');
        tempN.add((byte) 'G');
        tempN.add((byte) 'T');
        nucleotide_nomenclature.put((byte) 'N', tempN);

        return nucleotide_nomenclature;
    }

    // generate blosum62 substitution table
    public static double[][] generate_aminoacid_scoreTable() {
        int size = 0;
        for (int i = 0; i < alphabet_aa.length; ++i)
            if ((byte) alphabet_aa[i] > size)
                size = (byte) alphabet_aa[i]; // change char to integer based on ascii code
        size = size + 1;
        double[][] NodeScoreBigTable = new double[size][size];
        for (int i = 0; i < size; ++i)
            Arrays.fill(NodeScoreBigTable[0], 0);

        double[][] blosum62 = {
                // A R N D C Q E G H I L K M F P S T W Y V B Z X -
                { 4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4 },
                { -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4 },
                { -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4 },
                { -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4 },
                { 0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4 },
                { -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4 },
                { -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4 },
                { 0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4 },
                { -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4 },
                { -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4 },
                { -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4 },
                { -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4 },
                { -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4 },
                { -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4 },
                { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4 },
                { 1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4 },
                { 0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4 },
                { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4 },
                { -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4 },
                { 0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4 },
                { -2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4 },
                { -1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4 },
                { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4 },
                { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1 },
        };

        for (int i = 0; i < alphabet_aa.length; ++i) {
            for (int j = 0; j < alphabet_aa.length; ++j) {
                NodeScoreBigTable[alphabet_aa[i]][alphabet_aa[j]] = blosum62[i][j];
            }
        }

        // for(int i=0; i<alphabet_aa.length; ++i)
        // {
        // for(int j=0; j<alphabet_aa.length; ++j)
        // {
        // System.out.print(NodeScoreBigTable[alphabet_aa[i]][alphabet_aa[j]] + "\t");
        // }
        // System.out.print("\n");
        // }

        return NodeScoreBigTable;
    }

    public static String[] runMainVCF(InputStream queryVCFInputStream, boolean isMultiplePlacements,
            InputStream multationSequenceMapInputStream, InputStream seqIdxMapInputStream,
            InputStream refSequenceInputStream, BufferedReader treeBufferedReader) {
        String insfn = "";
        // String intfn = "/tipars/ser_obj/input.tree";
        String inafn = "";
        // String inqfn = "";
        String informat = "vcf";
        String outfn = "";
        String otype = "insertion";
        Boolean aa_flag = false;
        TIPars.isMultiplePlacements = isMultiplePlacements;

        String nidname = "label";
        String attname = "GenName";
        boolean printDisInfoOnScreen = true;

        _nucleotide_nomenclature = generate_nucleotide_nomenclature();
        _nucleotide_nomenclature_scoreTable = generate_nucleotide_nomenclature_scoreTable(); // IUPAC nucleotide codes
                                                                                             // substitution table
        _nucleotide_nomenclature_map2char = generate_nucleotide_nomenclature_characterTable();

        _aminoacid_scoreTable = generate_aminoacid_scoreTable(); // Amino acid substitution table using blosum62 matrix

        try {
            Runtime run = Runtime.getRuntime();
            long startTime = System.currentTimeMillis();
            String output_folder = "/tmp";

            _used_scoreTable = _nucleotide_nomenclature_scoreTable;
            if (aa_flag)
                _used_scoreTable = _aminoacid_scoreTable;

            HashMap<Integer, String> queryList = null;
            long parseStartTime = System.currentTimeMillis();

            try {
                ObjectInputStream ois = new ObjectInputStream(multationSequenceMapInputStream);
                multationSequencesMap = (ArrayList) ois.readObject();
                ois.close();
                multationSequenceMapInputStream.close();
                System.out.println("multationSequencesMap loaded");

                ois = new ObjectInputStream(seqIdxMapInputStream);
                seqIdxMap = (HashMap) ois.readObject();
                ois.close();
                seqIdxMapInputStream.close();
                System.out.println("seqIdxMap loaded");

                ois = new ObjectInputStream(refSequenceInputStream);
                ref_sequence = (byte[]) ois.readObject();
                ois.close();
                refSequenceInputStream.close();
                System.out.println("ref_sequence loaded");
            } catch (IOException ioe) {
                ioe.printStackTrace();
            } catch (ClassNotFoundException c) {
                System.out.println("Class not found");
                c.printStackTrace();
            }

            long serializeLoadTime = System.currentTimeMillis() - parseStartTime;
            System.out.println("load serilization items time: " + (double) serializeLoadTime / 1000);

            queryList = readVCFFile2Alignment(queryVCFInputStream);
            long parseTotalTime = System.currentTimeMillis() - parseStartTime;
            System.out.println("vcf/fasta total parse time: " + (double) parseTotalTime / 1000);

            long parseTreeStartTime = System.currentTimeMillis();
            NewickImporter tni = new NewickImporter(treeBufferedReader);
            Tree tree = tni.importTree(null);
            long parseTreeTotalTime = System.currentTimeMillis() - parseTreeStartTime;
            System.out.println("tree parse time: " + (double) parseTreeTotalTime / 1000);

            // init TIPars
            TIPars myAdd = new TIPars(tree, otype, output_folder);
            Tree outtree = null;

            long startTime2 = System.currentTimeMillis();

            if (!printDisInfoOnScreen)
                System.out.print("Progress: ");

            String progressInfo = "";
            // mutiple taxa insertion for vcf input
            ArrayList<Integer> queryIdxs = new ArrayList<Integer>(queryList.keySet());
            String[] tiparsOutput = new String[queryIdxs.size()];
            queryIdxs.sort(Comparator.naturalOrder());
            for (int i = 0; i < queryIdxs.size(); i++) {
                int idx = queryIdxs.get(i);
                String name = queryList.get(idx);
                ConcurrentHashMap<Integer, Byte> query = multationSequencesMap.get(idx);
                String qid = "q" + (i + 1);
                String pid = "p" + (i + 1);
                tiparsOutput[i] = myAdd.addQuerySequence(name, query, qid, pid, printDisInfoOnScreen, new double[3], otype,
                        0);
            }
            long endTime2 = System.currentTimeMillis();
            long totalTime2 = endTime2 - startTime2;

            System.out.println("Insertion time: " + (double) totalTime2 / 1000);
            return tiparsOutput;
        } catch (Exception e) {
            e.printStackTrace();
            return new String[1];
        }
    }

    // private int getNumberOfQueries(InputStream queryVCFInputStream) {
    //     BufferedReader reader = new BufferedReader(new InputStreamReader(queryVCFInputStream));
    //     int i = 1;
    //     String line;
    //     while(reader.ready() && i <= 4) {
    //         line = reader.readLine();
    //         i++;
    //     }
    //     int tabCount = 0;
    //     for (i = 0; i < line.length(); i++) {
    //         if (line.charAt(i) == '\t') {
    //             tabCount++;
    //         }
    //     }
    //     final int vcfDefaultTabs = 8;
    //     return tabCount - vcfDefaultTabs; 
    // }

    public static void main(String[] args) {
        
        try {
            FileInputStream ifsq = new FileInputStream(args[0]);
            FileInputStream ifsm = new FileInputStream(args[2]);
            FileInputStream ifss = new FileInputStream(args[3]);
            FileInputStream ifsr = new FileInputStream(args[4]);
            BufferedReader ifst = new BufferedReader(new FileReader(args[1]));
            String[] outputs = runMainVCF(ifsq, false, ifsm, ifss, ifsr, ifst);
            System.out.println(Arrays.toString(outputs));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
