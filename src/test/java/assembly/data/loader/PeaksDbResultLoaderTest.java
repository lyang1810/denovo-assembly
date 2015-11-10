package assembly.data.loader;

import assembly.data.DbSearchPeptide;

import java.io.File;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class PeaksDbResultLoaderTest {

    @org.junit.Test
    public void testLoad() throws Exception {
        PeaksDbResultLoader loader = new PeaksDbResultLoader();
        List<DbSearchPeptide> peptides = loader.load(
                new File[]{
                        new File("data/peaks DB search psm.csv")
                });
        assertEquals(29948, peptides.size());
    }
}