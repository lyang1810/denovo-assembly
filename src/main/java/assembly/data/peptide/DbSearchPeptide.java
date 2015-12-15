package assembly.data.peptide;

import assembly.data.residue.AminoAcid;

import java.util.Arrays;

/**
 * database search peptide
 */
public class DbSearchPeptide implements Peptide {
    private int fractionIdx;
    private int spectrumIdx;
    private float retentionTime;

    private AminoAcid[] sequence;
    private float score;

    public DbSearchPeptide(AminoAcid[] sequence) {
        this.sequence = sequence;
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

    @Override
    public String toString() {
        return Arrays.toString(sequence);
    }
}
