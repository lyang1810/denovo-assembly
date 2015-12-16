package assembly.experiment;

import assembly.data.residue.AminoAcid;
import assembly.data.residue.AminoAcidFactory;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.*;

/**
 * Test for Evaluator
 * Created by Lian on 2015-12-16.
 */
public class EvaluatorTest {

    @Test
    public void testGetResidueCorrectness() throws Exception {
        AminoAcidFactory factory = AminoAcidFactory.getInstance();

        AminoAcid[] dnSeq = new AminoAcid[]{
                factory.getAminoAcidByCode("E"),
                factory.getAminoAcidByCode("P"),
                factory.getAminoAcidByCode("P"),
                factory.getAminoAcidByCode("T"),
                factory.getAminoAcidByCode("L"),
                factory.getAminoAcidByCode("D"),
                factory.getAminoAcidByCode("E")
        };

        AminoAcid[] dbSeq = new AminoAcid[]{
                factory.getAminoAcidByCode("P"),
                factory.getAminoAcidByCode("E"),
                factory.getAminoAcidByCode("P"),
                factory.getAminoAcidByCode("T"),
                factory.getAminoAcidByCode("I"),
                factory.getAminoAcidByCode("E"),
                factory.getAminoAcidByCode("D")
        };

        boolean[] correctness = Evaluator.checkResidueCorrectness(dnSeq, dbSeq);
        boolean[] expected = new boolean[]{false, false, true, true, true, false, false};
        assertArrayEquals(expected, correctness);
    }
}