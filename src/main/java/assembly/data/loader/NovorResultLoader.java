package assembly.data.loader;

import assembly.data.AminoAcid;
import assembly.data.AminoAcidFactory;
import assembly.data.DenovoPeptide;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

public class NovorResultLoader implements ResultLoader {
    private static final Logger log = Logger.getLogger(NovorResultLoader.class.getName());

    @Override
    public List<DenovoPeptide> load(File[] files) {
        List<DenovoPeptide> peptides = new ArrayList<>();

        for (int i = 0; i < files.length; ++i) {
            File file = files[i];
            try (BufferedReader br = new BufferedReader(new FileReader(file))) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.length() == 0 || line.startsWith("#")) {
                        continue;
                    }
                    String[] splits = line.split(", ");

                    try {
                        AminoAcid[] seq = parseSequence(splits[9]);
                        float[] conf = parseConfidence(splits[10]);
                        int spectrumIdx = parseFractionSpectrumIdx(splits[1]);
                        float rt = parseRetentionTime(splits[2]);

                        DenovoPeptide pep = new DenovoPeptide(seq, conf);
                        pep.setFractionIdx(i);
                        pep.setSpectrumIdx(spectrumIdx);
                        pep.setRetentionTime(rt);

                        peptides.add(pep);
                    } catch (NumberFormatException e) {
                        log.log(Level.SEVERE, "Failed to parse line " + line, e);
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
                String mod = seqStr.substring(i + 1, j + 1);
                switch (mod) {
                    case "(Cam)":
                        code += "(+57.02)";
                        break;
                    case "(O)":
                        code += "(+15.99)";
                        break;
                    default:
                        log.log(Level.SEVERE, "Unknown modification " + mod);
                }
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
        String[] splits = confStr.split("-");
        float[] conf = new float[splits.length];
        for (int i = 0; i < splits.length; i++) {
            conf[i] = Float.parseFloat(splits[i]);
        }
        return conf;
    }

    private int parseFractionSpectrumIdx(String scanStr) {
        return Integer.parseInt(scanStr);
    }

    private float parseRetentionTime(String rtStr) {
        return Float.parseFloat(rtStr);
    }
}
