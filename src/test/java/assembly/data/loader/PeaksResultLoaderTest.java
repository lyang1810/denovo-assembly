package assembly.data.loader;

import assembly.data.DenovoPeptide;

import java.io.File;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class PeaksResultLoaderTest {

    @org.junit.Test
    public void testLoad() throws Exception {
        ResultLoader loader = new PeaksResultLoader();
        List<DenovoPeptide> peptides = loader.load(
                new File[]{
                        new File("data/peaks de novo peptides.csv"),
                        new File("data/peaks all de novo candidates.csv")
                });
        assertEquals(1098 + 16520, peptides.size());
    }
}