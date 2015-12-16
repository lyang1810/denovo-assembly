package assembly.experiment;

import assembly.data.peptide.DbSearchPeptide;
import assembly.data.peptide.DenovoPeptide;
import assembly.data.peptide.Peptide;
import assembly.data.residue.AminoAcid;
import assembly.util.MassUtil;

import java.util.*;
import java.util.logging.Logger;

/**
 * Generate metrics for measuring the quality of de novo sequencing result, using database search result as reference.
 * Created by Lian on 2015-12-16.
 */
public class Evaluator {
    private static final Logger log = Logger.getLogger(Evaluator.class.getName());

    private Map<DenovoPeptide, boolean[]> dnMatches; // de novo PSM -> residue correctness
    private Map<String, DbSearchPeptide> dbMatches; // scan tag -> db search PSM

    private int totalReferenceResidueCount;

    /**
     * de novo sequencing PSMs and reference database search PSMs. Assume only one entry for each scan
     *
     * @param dnMatches           de novo sequences PSMs
     * @param dbMatches           datbase search PSMs as a reference. Suggest 1% FDR.
     * @param fragmentMzTolerance fragment m/z tolerance in Dalton. Suggest 0.1 Da
     */
    public Evaluator(List<DenovoPeptide> dnMatches, List<DbSearchPeptide> dbMatches, float fragmentMzTolerance) {
        this.totalReferenceResidueCount = 0;
        // Select de novo matches with reference db search result on the same scan
        this.dnMatches = new HashMap<>();
        this.dbMatches = new HashMap<>();
        for (DbSearchPeptide dbMatch : dbMatches) {
            this.dbMatches.put(getScanTag(dbMatch), dbMatch);
            this.totalReferenceResidueCount += dbMatch.getSequence().length;
        }
        for (DenovoPeptide dnMatch : dnMatches) {
            String scanTag = getScanTag(dnMatch);
            if (this.dbMatches.containsKey(scanTag)) {
                AminoAcid[] dnSeq = dnMatch.getSequence();
                AminoAcid[] dbSeq = this.dbMatches.get(scanTag).getSequence();

                float massDiff = Math.abs(MassUtil.getTotalResidueMass(dnSeq) - MassUtil.getTotalResidueMass(dbSeq));
                if (massDiff > fragmentMzTolerance) {
                    log.warning("Unmatched total residue mass. Diff = " + massDiff);
                }
                boolean[] residueCorrectness = checkResidueCorrectness(dnSeq, dbSeq, fragmentMzTolerance);
                this.dnMatches.put(dnMatch, residueCorrectness);
            }
        }
    }

    public float getPrecisionRate(float residueScoreThreshold) {
        int correctResidueCount = 0;
        int totalResidueCount = 0;

        for (DenovoPeptide dnMatch : dnMatches.keySet()) {
            float[] residueScore = dnMatch.getConfidence();
            boolean[] correctness = dnMatches.get(dnMatch);
            for (int i = 0; i < residueScore.length; i++) {
                if (residueScore[i] >= residueScoreThreshold) {
                    totalResidueCount++;
                    if (correctness[i]) {
                        correctResidueCount++;
                    }
                }
            }
        }
        return (float) correctResidueCount / totalResidueCount;
    }

    public float getRecallRate(float residueScoreThreshold) {
        int correctResidueCount = 0;

        for (DenovoPeptide dnMatch : dnMatches.keySet()) {
            float[] residueScore = dnMatch.getConfidence();
            boolean[] correctness = dnMatches.get(dnMatch);
            for (int i = 0; i < residueScore.length; i++) {
                if (residueScore[i] >= residueScoreThreshold && correctness[i]) {
                    correctResidueCount++;
                }
            }
        }
        return (float) correctResidueCount / totalReferenceResidueCount;
    }

    static boolean[] checkResidueCorrectness(AminoAcid[] dnSeq, AminoAcid[] dbSeq, float fragmentMzTolerance) {
        boolean[] correctness = new boolean[dnSeq.length];
        Arrays.fill(correctness, false);

        float[] dnPrefix = new float[dnSeq.length + 1];
        float[] dbPrefix = new float[dbSeq.length + 1];

        dnPrefix[0] = 0;
        for (int i = 0; i < dnSeq.length; i++) {
            dnPrefix[i + 1] = dnPrefix[i] + dnSeq[i].getMass();
        }

        dbPrefix[0] = 0;
        for (int i = 0; i < dbSeq.length; i++) {
            dbPrefix[i + 1] = dbPrefix[i] + dbSeq[i].getMass();
        }

        int p = 1;
        int q = 1;
        while (p < dnPrefix.length && q < dbPrefix.length) {
            if (Math.abs(dnPrefix[p] - dbPrefix[q]) < fragmentMzTolerance) {
                if (Math.abs(dnPrefix[p - 1] - dbPrefix[q - 1]) < fragmentMzTolerance) {
                    correctness[p - 1] = true;
                }
                p++;
                q++;
            } else if (dnPrefix[p] < dbPrefix[q]) {
                p++;
            } else if (dnPrefix[p] > dbPrefix[q]) {
                q++;
            }
        }
        return correctness;
    }

    static String getScanTag(Peptide match) {
        return match.getFractionIdx() + ":" + match.getSpectrumIdx();
    }
}
