package assembly.data.loader;

import assembly.data.residue.AminoAcid;
import assembly.data.residue.AminoAcidFactory;
import assembly.data.peptide.DenovoPeptide;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings("Duplicates")
public class PeaksResultLoader implements ResultLoader {
    private static final Logger log = Logger.getLogger(PeaksResultLoader.class.getName());

    @Override
    public List<DenovoPeptide> load(File[] files) {
        List<DenovoPeptide> peptides = new ArrayList<>();

        for (int i = 0; i < files.length; ++i) {
            File file = files[i];
            try (BufferedReader br = new BufferedReader(new FileReader(file))) {
                // Acquire column names and index
                Map<String, Integer> headerIdxMap = new HashMap<>();
                String header = br.readLine();
                String[] headers = header.split(",");
                for (int idx = 0; idx < headers.length; ++idx) {
                    headerIdxMap.put(headers[idx], idx);
                }
                int seqColIdx = headerIdxMap.get("Peptide");
                int scoreColIdx = headerIdxMap.get("ALC (%)");
                int confColIdx = headerIdxMap.get("local confidence (%)");
                int scanColIdx = headerIdxMap.get("Scan");
                int precursorMzColIdx = headerIdxMap.get("m/z");
                int rtColIdx = headerIdxMap.get("RT");

                String line;
                while ((line = br.readLine()) != null) {
                    if (line.length() == 0) {
                        continue;
                    }
                    String[] splits = line.split(",");

                    try {
                        AminoAcid[] seq = parseSequence(splits[seqColIdx]);
                        float[] conf = parseConfidence(splits[confColIdx]);
                        float score = parseScore(splits[scoreColIdx]);
                        int[] fractionSpectrumIdx = parseFractionSpectrumIdx(splits[scanColIdx]);
                        int fractionIdx = (fractionSpectrumIdx[0] == -1) ? i : fractionSpectrumIdx[0];
                        int spectrumIdx = fractionSpectrumIdx[1];
                        float precursorMz = parsePrecursorMz(splits[precursorMzColIdx]);
                        float rt = parseRetentionTime(splits[rtColIdx]);

                        DenovoPeptide pep = new DenovoPeptide(seq, conf);
                        pep.setFractionIdx(fractionIdx);
                        pep.setSpectrumIdx(spectrumIdx);
                        pep.setPrecursorMz(precursorMz);
                        pep.setRetentionTime(rt);
                        pep.setScore(score);

                        peptides.add(pep);
                    } catch (NumberFormatException e) {
                        if (!splits[12].contains("/")) {
                            // write log only if the parsing error is not caused by PEAKS merged results
                            log.log(Level.INFO, "Failed to parse line " + line);
                        }

                    }
                }
            } catch (IOException e) {
                log.log(Level.SEVERE, "Failed to load file " + file, e);
            }
        }
        return peptides;
    }

    private AminoAcid[] parseSequence(String seqStr) {
        List<AminoAcid> seq = new ArrayList<>();
        AminoAcidFactory factory = AminoAcidFactory.getInstance();

        for (int i = 0; i < seqStr.length(); ++i) {
            String code = seqStr.substring(i, i + 1);
            if (i + 1 < seqStr.length() && seqStr.charAt(i + 1) == '(') {
                int j = seqStr.indexOf(')', i + 1);
                code = seqStr.substring(i, j + 1);
                i = j;
            }
            AminoAcid aa = factory.getAminoAcidByCode(code);
            if (aa == null) {
                log.log(Level.SEVERE, "Unknown amino acid " + code);
            }
            seq.add(aa);
        }
        return seq.toArray(new AminoAcid[seq.size()]);
    }

    private float[] parseConfidence(String confStr) {
        String[] splits = confStr.split(" ");
        float[] conf = new float[splits.length];
        for (int i = 0; i < splits.length; i++) {
            conf[i] = Float.parseFloat(splits[i]);
        }
        return conf;
    }

    private float parseScore(String scoreStr) {
        return Float.parseFloat(scoreStr);
    }

    private int[] parseFractionSpectrumIdx(String scanStr) {
        if (scanStr.contains(":")) {
            String[] splits = scanStr.split(":");
            int fractionIdx = Integer.parseInt(splits[0].substring(1));
            int spectrumIdx = Integer.parseInt(splits[1]);
            return new int[]{fractionIdx, spectrumIdx};
        } else {
            int spectrumIdx = Integer.parseInt(scanStr);
            return new int[]{-1, spectrumIdx};
        }
    }

    private float parsePrecursorMz(String precursorMzStr) {
        return Float.parseFloat(precursorMzStr);
    }

    private float parseRetentionTime(String rtStr) {
        return Float.parseFloat(rtStr);
    }
}
