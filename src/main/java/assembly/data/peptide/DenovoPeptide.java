package assembly.data.peptide;

import assembly.data.residue.AminoAcid;

import java.util.Arrays;

/**
 * De novo peptide
 */
public class DenovoPeptide implements Peptide {
    private int fractionIdx;
    private int spectrumIdx;

    private float precursorMz;
    private float retentionTime;

    private AminoAcid[] sequence;
    private float[] confidence;
    private float score;

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

    public void setPrecursorMz(float precursorMz) {
        this.precursorMz = precursorMz;
    }

    public void setRetentionTime(float retentionTime) {
        this.retentionTime = retentionTime;
    }

    public void setScore(float score) {
        this.score = score;
    }

    @Override
    public int getFractionIdx() {
        return fractionIdx;
    }

    @Override
    public int getSpectrumIdx() {
        return spectrumIdx;
    }

    @Override
    public float getPrecursorMz() {
        return precursorMz;
    }

    @Override
    public float getRetentionTime() {
        return retentionTime;
    }

    @Override
    public AminoAcid[] getSequence() {
        return sequence;
    }

    @Override
    public float getScore() {
        return score;
    }

    public float[] getConfidence() {
        return confidence;
    }

    @Override
    public String toString() {
        return Arrays.toString(sequence);
    }

}
