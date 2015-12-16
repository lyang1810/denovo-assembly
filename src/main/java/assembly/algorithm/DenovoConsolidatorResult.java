package assembly.algorithm;

import assembly.data.residue.AminoAcid;

/**
 * Result of DenovoConsolidator
 * Created by yangl on 2015-12-15.
 */
public class DenovoConsolidatorResult {
    AminoAcid[] sequence;
    float[] confidence;

    public DenovoConsolidatorResult(AminoAcid[] sequence, float[] confidence) {
        this.sequence = sequence;
        this.confidence = confidence;
    }

    public AminoAcid[] getSequence() {
        return sequence;
    }

    public float[] getConfidence() {
        return confidence;
    }
}
