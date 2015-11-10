package assembly.data;

public interface Peptide {

    int getFractionIdx();

    int getSpectrumIdx();

    float getRetentionTime();

    AminoAcid[] getSequence();

    @Override
    String toString();
}
