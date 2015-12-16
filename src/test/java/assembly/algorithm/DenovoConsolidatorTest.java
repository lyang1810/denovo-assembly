package assembly.algorithm;

import assembly.data.loader.PeaksResultLoader;
import assembly.data.loader.ResultLoader;
import assembly.data.peptide.DenovoPeptide;
import assembly.util.MassUtil;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Test for DenovoConsolidator
 * Created by yangl on 2015-12-15.
 */
public class DenovoConsolidatorTest {

    @Test
    public void testRun() throws Exception {
        ResultLoader loader = new PeaksResultLoader();
        List<DenovoPeptide> peptides = loader.load(new File[]{new File("data/peaks all de novo candidates.csv")});

        float totalResidueMass = MassUtil.getTotalResidueMass(peptides.get(5).getSequence());
        float fragmentMzTolerance = 0.1f;
        DenovoConsolidator consolidator = new DenovoConsolidator(totalResidueMass, fragmentMzTolerance);

        consolidator.addDenovoTags(peptides.get(5), 0);
        consolidator.addDenovoTags(peptides.get(6), 0);
        consolidator.addDenovoTags(peptides.get(7), 0);
        consolidator.addDenovoTags(peptides.get(8), 0);
        consolidator.addDenovoTags(peptides.get(9), 0);

        DenovoConsolidatorResult result = consolidator.run();
        System.out.println(Arrays.toString(result.getSequence()));
        System.out.println(Arrays.toString(result.getConfidence()));
    }
}