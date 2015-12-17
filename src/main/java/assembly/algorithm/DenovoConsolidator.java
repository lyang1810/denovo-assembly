package assembly.algorithm;

import assembly.data.peptide.DenovoPeptide;
import assembly.data.residue.AminoAcid;
import assembly.data.residue.AminoAcidFactory;

import java.util.*;

/**
 * Core algorithm for consolidating de novo sequencing results
 * Created by Lian on 2015-12-15.
 */
public class DenovoConsolidator {

    private float totalResidueMass;
    private float fragmentMzTolerance;
    private TreeMap<Float, Node> network;

    /**
     * Initiate a DenovoConsolidator
     *
     * @param fragmentMzTolerance fragment m/z tolerance in Dalton. Suggest 0.1 Da
     */
    public DenovoConsolidator(float totalResidueMass, float fragmentMzTolerance) {
        this.totalResidueMass = totalResidueMass;
        this.fragmentMzTolerance = fragmentMzTolerance;
        this.network = new TreeMap<>();
        this.network.put(totalResidueMass, new Node(totalResidueMass, null));
    }

    /**
     * Add tags from a de novo peptide to the network
     *
     * @param peptide A de novo peptide
     * @param shift   Add prefix mass to the peptide.
     *                This is for aligning the peptide to the tags already in the network.
     *                This prefix can be either positive (shift right) or negative (shift left)
     */
    public void addDenovoTags(DenovoPeptide peptide, float shift) {
        float nodeMass = shift;
        for (int i = 0; i < peptide.getSequence().length; i++) {
            AminoAcid residue = peptide.getSequence()[i];
            float confidence = peptide.getConfidence()[i];
            if (nodeMass >= 0) {
                if (network.get(nodeMass) == null) {
                    network.put(nodeMass, new Node(nodeMass, null));
                }
                Node node = network.get(nodeMass);
                if (node.tags == null) {
                    node.tags = new ArrayList<>();
                }
                node.tags.add(new DenovoTag(residue, confidence));
            }
            nodeMass += residue.getMass();
        }
    }

    /**
     * Batched version of {@linkplain #addDenovoTags(DenovoPeptide, float)}
     *
     * @param peptides a list of de novo peptides
     * @param shift    Add prefix mass to each of the peptides.
     *                 This is for aligning the peptide to the tags already in the network.
     *                 This prefix can be either positive (shift right) or negative (shift left)
     */
    public void addDenovoTags(List<DenovoPeptide> peptides, float shift) {
        for (DenovoPeptide peptide : peptides) {
            addDenovoTags(peptide, shift);
        }
    }

    public DenovoConsolidatorResult run() {
        // dynamic programming to calculate the best path
        Set<AminoAcid> aaSet = AminoAcidFactory.getInstance().getAminoAcidSet();

        NavigableSet<Float> nodeMasses = network.navigableKeySet();
        for (Float nodeMass : nodeMasses) {
            float maxNodeScore = Float.NaN;
            Node maxNodeScorePrevNode = null;
            DenovoTag maxNodeScorePrevTag = null;

            for (AminoAcid aa : aaSet) {
                float prevMass = nodeMass - aa.getMass();
                if (prevMass < -1e-6) {
                    continue;
                }
                float lowerBound = prevMass - fragmentMzTolerance;
                float upperBound = prevMass + fragmentMzTolerance;
                Collection<Node> prevNodes = network.subMap(lowerBound, true, upperBound, true).values();

                if (prevNodes.size() > 0) {
                    // get the prev node with best score
                    Node bestPrevNode = null;
                    for (Node prevNode : prevNodes) {
                        if (bestPrevNode == null || bestPrevNode.score < prevNode.score) {
                            bestPrevNode = prevNode;
                        }
                    }
                    // get the best tag from a prev node to current node
                    DenovoTag bestPrevTag = null;
                    for (Node prevNode : prevNodes) {
                        for (DenovoTag tag : prevNode.tags) {
                            if (tag.residue == aa &&
                                    (bestPrevTag == null || bestPrevTag.confidence < tag.confidence)) {
                                bestPrevTag = tag;
                            }
                        }
                    }
                    // calculate score for the current node. Only the maximum score is kept.
                    if (bestPrevNode != null && bestPrevTag != null) {
                        float nodeScore = (bestPrevNode.score * bestPrevNode.length + bestPrevTag.confidence) / (bestPrevNode.length + 1);
                        if (Float.isNaN(maxNodeScore) || nodeScore > maxNodeScore) {
                            maxNodeScore = nodeScore;
                            maxNodeScorePrevNode = bestPrevNode;
                            maxNodeScorePrevTag = bestPrevTag;
                        }
                    }
                }
            }
            if (!Float.isNaN(maxNodeScore)) {
                network.get(nodeMass).score = maxNodeScore;
                network.get(nodeMass).length = maxNodeScorePrevNode.length + 1;
                network.get(nodeMass).bestPrevNode = maxNodeScorePrevNode;
                network.get(nodeMass).bestPrevTag = maxNodeScorePrevTag;
            }
        }

        // retrieve the best sequence path and confidence scores
        LinkedList<AminoAcid> sequence = new LinkedList<>();
        LinkedList<Float> confidence = new LinkedList<>();

        Node node = network.get(totalResidueMass);
        while (node.mass > 0) {
            DenovoTag tag = node.bestPrevTag;
            sequence.addFirst(tag.residue);
            confidence.addFirst(tag.confidence);

            node = node.bestPrevNode;
        }

        // pack result
        AminoAcid[] seqArray = sequence.toArray(new AminoAcid[sequence.size()]);
        float[] confArray = new float[confidence.size()];
        for (int i = 0; i < confArray.length; i++) {
            confArray[i] = confidence.get(i);
        }
        return new DenovoConsolidatorResult(seqArray, confArray);
    }

    /**
     * Represent a de novo residue tag in network
     */
    class DenovoTag {
        AminoAcid residue;
        float confidence;

        public DenovoTag(AminoAcid residue, float confidence) {
            this.residue = residue;
            this.confidence = confidence;
        }
    }

    /**
     * Represent a mass node in network
     */
    class Node {
        float mass;     // node mass
        List<DenovoTag> tags;   // de novo tags starting from this node

        float score;    // average residue confidence in the best path to this node
        float length;   // number of residues in the best path to this node

        Node bestPrevNode; // previous node in the best path to this node
        DenovoTag bestPrevTag; // last tag in the best path to this node

        public Node(float mass, List<DenovoTag> tags) {
            this.mass = mass;
            this.tags = tags;
            this.score = 0;
            this.length = 0;
            this.bestPrevNode = null;
            this.bestPrevTag = null;
        }
    }
}

