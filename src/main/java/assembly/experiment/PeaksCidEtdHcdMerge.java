package assembly.experiment;

import assembly.data.peptide.DbSearchPeptide;
import assembly.data.peptide.DenovoPeptide;
import assembly.data.loader.PeaksDbResultLoader;
import assembly.data.loader.PeaksResultLoader;

import java.io.File;
import java.util.List;

public class PeaksCidEtdHcdMerge {

    public static void main(String[] args) {
        PeaksResultLoader dnLoader = new PeaksResultLoader();
        List<DenovoPeptide> dnPeptides = dnLoader.load(new File[]{new File("CID_ETD_HCD/dataset-1/all de novo candidates.csv")});

        PeaksDbResultLoader dbLoader = new PeaksDbResultLoader();
        List<DbSearchPeptide> dbPeptides = dbLoader.load(new File[]{new File("CID_ETD_HCD/dataset-1/DB search psm.csv")});

        System.out.println("Done");
    }
}
