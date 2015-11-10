package assembly.data;

import java.util.Arrays;

/**
 * De novo peptide
 */
public class DenovoPeptide implements Peptide {
    private int fractionIdx;
    private int spectrumIdx;
    private float retentionTime;

    private AminoAcid[] sequence;
    private float[] confidence;

    public DenovoPeptide(AminoAcid[] sequence, float[] confidence) {
        this.sequence = sequence;
        this.confidence = confidence;
    }

    public void setFractionIdx(int fractionIdx) {
        this.fractionIdx = fractionIdx;
    }

    public void setSpectrumIdx(int spectrumIdx) {
        this.spectrumIdx = spectrumIdx;
    }

    public void setRetentionTime(float retentionTime) {
        this.retentionTime = retentionTime;
    }

    public int getFractionIdx() {
        return fractionIdx;
    }

    public int getSpectrumIdx() {
        return spectrumIdx;
    }

    public float getRetentionTime() {
        return retentionTime;
    }

    public AminoAcid[] getSequence() {
        return sequence;
    }

    public float[] getConfidence() {
        return confidence;
    }

    @Override
    public String toString() {
        return Arrays.toString(sequence);
    }
}
