package assembly.data.peptide;

import assembly.data.residue.AminoAcid;

public interface Peptide {

    int getFractionIdx();

    int getSpectrumIdx();

    float getPrecursorMz();

    float getRetentionTime();

    AminoAcid[] getSequence();

    float getScore();

    @Override
    String toString();
}
