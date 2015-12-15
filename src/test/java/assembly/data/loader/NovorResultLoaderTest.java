package assembly.data.loader;

import assembly.data.peptide.DenovoPeptide;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class NovorResultLoaderTest {

    @Test
    public void testLoad() throws Exception {
        ResultLoader loader = new NovorResultLoader();
        List<DenovoPeptide> peptides = loader.load(
                new File[]{
                        new File("data/novor Light-Chain-Trypsin-1.mgf.csv")
                });
        assertEquals(9160, peptides.size());
    }
}