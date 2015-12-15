package assembly.data.loader;

import assembly.data.peptide.DbSearchPeptide;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class PeaksDbResultLoaderTest {

    @Test
    public void testLoad() throws Exception {
        PeaksDbResultLoader loader = new PeaksDbResultLoader();
        List<DbSearchPeptide> peptides = loader.load(
                new File[]{
                        new File("data/peaks DB search psm.csv")
                });
        assertEquals(29948, peptides.size());
    }
}