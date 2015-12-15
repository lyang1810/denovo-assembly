package assembly.data.residue;

import assembly.data.residue.AminoAcid;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class AminoAcidTest {
    private AminoAcid aa;

    @Before
    public void setUp() throws Exception {
        aa = new AminoAcid('X', 1.234f);
    }

    @Test
    public void testGetCode() throws Exception {
        assertEquals("X", aa.getCode());
    }

    @Test
    public void testGetMass() throws Exception {
        assertEquals(1.234f, aa.getMass(), 1e-6);
    }

    @Test
    public void testToString() throws Exception {
        assertEquals("X", aa.toString());
    }
}