package assembly.util;

import assembly.data.residue.AminoAcid;

/**
 * Utility for mass calculation
 * Created by yangl on 2015-12-15.
 */
public class MassUtil {
    /**
     * Get total residue mass in a peptide
     *
     * @param sequence amino acid sequence of the peptide
     * @return total residue mass
     */
    public static float getTotalResidueMass(AminoAcid[] sequence) {
        float mass = 0;
        for (AminoAcid aa : sequence) {
            mass += aa.getMass();
        }
        return mass;
    }
}
