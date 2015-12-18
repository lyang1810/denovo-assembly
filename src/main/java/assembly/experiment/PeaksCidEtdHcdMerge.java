package assembly.experiment;

import assembly.algorithm.DenovoConsolidator;
import assembly.algorithm.DenovoConsolidatorResult;
import assembly.data.loader.PeaksDbResultLoader;
import assembly.data.loader.PeaksResultLoader;
import assembly.data.peptide.DbSearchPeptide;
import assembly.data.peptide.DenovoPeptide;
import assembly.data.residue.AminoAcid;
import assembly.util.MassUtil;

import java.io.File;
import java.util.*;

public class PeaksCidEtdHcdMerge {
    private static float fragmentMzTolerance = 0.1f;

    public static void main(String[] args) {
        run("CID_ETD_HCD/dataset-1");
        run("CID_ETD_HCD/dataset-2");
        run("CID_ETD_HCD/dataset-3");
    }

    private static void run(String dataset) {
        System.out.println("Running dataset " + dataset);

        // Load PEAKS de novo results
        PeaksResultLoader dnLoader = new PeaksResultLoader();
        List<DenovoPeptide> dnMatches = dnLoader.load(new File[]{new File(dataset + "/all de novo candidates.csv")});

        // Load PEAKS DB reference result (1% FDR)
        PeaksDbResultLoader dbLoader = new PeaksDbResultLoader();
        List<DbSearchPeptide> dbMatches = dbLoader.load(new File[]{new File(dataset + "/DB search psm.csv")});

        // Output baseline 1 metrics: For each scan, use the best PEAKS de novo candidate from the scan
        {
            System.out.println("Baseline 1");
            List<DenovoPeptide> processedDnMatches = prepareBaseline1(dnMatches);
            evaluate(dbMatches, processedDnMatches);
        }

        // Output baseline 2 metrics: For each scan, use the best PEAKS de novo candidate from the CID/HCD/ETD scan group
        {
            System.out.println("Baseline 2");
            List<DenovoPeptide> processedDnMatches = prepareBaseline2(dnMatches);
            evaluate(dbMatches, processedDnMatches);
        }

        // Output experimental metrics: For each scan, use the sequence consolidated from all PEAKS de novo candidates in the CID/HCD/ETD group
        {
            System.out.println("Experimental");
            List<DenovoPeptide> processedDnMatches = prepareExperimental(dnMatches);
            evaluate(dbMatches, processedDnMatches);
        }
    }

    private static void evaluate(List<DbSearchPeptide> dbMatches, List<DenovoPeptide> dnMatchesConsolidated) {
        Evaluator eval = new Evaluator(dnMatchesConsolidated, dbMatches, fragmentMzTolerance);
        for (int t = 0; t <= 100; t++) {
            float precision = eval.getPrecisionRate(t);
            float recall = eval.getRecallRate(t);
            System.out.println(t + "\t" + precision + "\t" + recall);
        }
    }

    private static List<DenovoPeptide> prepareBaseline1(List<DenovoPeptide> dnMatches) {
        LinkedHashMap<String, DenovoPeptide> dnMatchesFiltered = new LinkedHashMap<>();

        // keep the top candidate for every scan
        for (DenovoPeptide match : dnMatches) {
            String scanTag = match.getFractionIdx() + ":" + match.getSpectrumIdx();
            if (dnMatchesFiltered.containsKey(scanTag)) {
                if (match.getScore() > dnMatchesFiltered.get(scanTag).getScore()) {
                    dnMatchesFiltered.put(scanTag, match);
                }
            } else {
                dnMatchesFiltered.put(scanTag, match);
            }
        }
        return new ArrayList<>(dnMatchesFiltered.values());
    }

    private static List<DenovoPeptide> prepareBaseline2(List<DenovoPeptide> dnMatches) {
        List<DenovoPeptide> dnMatchesConsolidated = new ArrayList<>();

        // find the best de novo candidate of scans in a CID/ETD/HCD group
        List<DenovoPeptide> group = new ArrayList<>();
        float groupMz = -1;

        for (DenovoPeptide dnMatch : dnMatches) {
            if (Math.abs(dnMatch.getPrecursorMz() - groupMz) > 1e-6) {
                dnMatchesConsolidated.addAll(processBaseline2Group(group));
                group.clear();
                groupMz = dnMatch.getPrecursorMz();
                group.add(dnMatch);
            } else {
                group.add(dnMatch);
            }
        }

        // last group
        dnMatchesConsolidated.addAll(processBaseline2Group(group));
        return dnMatchesConsolidated;
    }

    private static List<DenovoPeptide> processBaseline2Group(List<DenovoPeptide> group) {
        List<DenovoPeptide> dnMatchedMerged = new ArrayList<>();
        if (group == null || group.size() == 0) {
            return dnMatchedMerged;
        }

        DenovoPeptide bestMatch = null;
        for (DenovoPeptide match : group) {
            if (bestMatch == null || match.getScore() > bestMatch.getScore()) {
                bestMatch = match;
            }
        }

        // assemble a merged de novo match for each scan
        Set<String> encountered = new HashSet<>();
        for (DenovoPeptide match : group) {
            String scanTag = match.getFractionIdx() + ":" + match.getSpectrumIdx();
            if (!encountered.contains(scanTag)) {
                DenovoPeptide updated = new DenovoPeptide(bestMatch.getSequence(), bestMatch.getConfidence());
                updated.setScore(bestMatch.getScore());
                updated.setPrecursorMz(match.getPrecursorMz());
                updated.setRetentionTime(match.getRetentionTime());
                updated.setFractionIdx(match.getFractionIdx());
                updated.setSpectrumIdx(match.getSpectrumIdx());
                dnMatchedMerged.add(updated);
                encountered.add(scanTag);
            }
        }
        return dnMatchedMerged;
    }

    private static List<DenovoPeptide> prepareExperimental(List<DenovoPeptide> dnMatches) {
        List<DenovoPeptide> dnMatchesConsolidated = new ArrayList<>();

        // consolidate de novo candidates of scans in a CID/ETD/HCD group
        List<DenovoPeptide> group = new ArrayList<>();
        float groupMz = -1;

        for (DenovoPeptide dnMatch : dnMatches) {
            if (Math.abs(dnMatch.getPrecursorMz() - groupMz) > 1e-6) {
                dnMatchesConsolidated.addAll(processExperimentalGroup(group));
                group.clear();
                groupMz = dnMatch.getPrecursorMz();
                group.add(dnMatch);
            } else {
                group.add(dnMatch);
            }
        }

        // last group
        dnMatchesConsolidated.addAll(processExperimentalGroup(group));
        return dnMatchesConsolidated;
    }

    private static List<DenovoPeptide> processExperimentalGroup(List<DenovoPeptide> group) {
        List<DenovoPeptide> dnMatchedMerged = new ArrayList<>();
        if (group == null || group.size() == 0) {
            return dnMatchedMerged;
        }

        // run the consolidating algorithm
        float totalResidueMass = MassUtil.getTotalResidueMass(group.get(0).getSequence());
        DenovoConsolidator consolidator = new DenovoConsolidator(totalResidueMass, fragmentMzTolerance);
        consolidator.addDenovoTags(group, 0);
        DenovoConsolidatorResult result = consolidator.run();

        // assemble a merged de novo match for each scan
        Set<String> encountered = new HashSet<>();
        for (DenovoPeptide match : group) {
            String scanTag = match.getFractionIdx() + ":" + match.getSpectrumIdx();
            if (!encountered.contains(scanTag)) {
                DenovoPeptide updated = new DenovoPeptide(result.getSequence(), result.getConfidence());
                updated.setPrecursorMz(match.getPrecursorMz());
                updated.setRetentionTime(match.getRetentionTime());
                updated.setFractionIdx(match.getFractionIdx());
                updated.setSpectrumIdx(match.getSpectrumIdx());
                float score = 0;
                for (float residueScore : result.getConfidence()) {
                    score += residueScore;
                }
                updated.setScore(score / result.getConfidence().length);
                dnMatchedMerged.add(updated);
                encountered.add(scanTag);
            }
        }
        return dnMatchedMerged;
    }
}
