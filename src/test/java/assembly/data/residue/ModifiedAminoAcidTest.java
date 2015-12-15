package assembly.data.residue;

import assembly.data.residue.AminoAcid;
import assembly.data.residue.ModifiedAminoAcid;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class ModifiedAminoAcidTest {
    private AminoAcid aa;

    @Before
    public void setUp() throws Exception {
        aa = new ModifiedAminoAcid('X', 1.234f, 2.345f);
    }

    @Test
    public void testGetCode() throws Exception {
        assertEquals("X(+2.35)", aa.getCode());
    }

    @Test
    public void testGetMass() throws Exception {
        assertEquals(1.234f + 2.345f, aa.getMass(), 1e-6);
    }

    @Test
    public void testToString() throws Exception {
        assertEquals("X(+2.35)", aa.toString());
    }
}