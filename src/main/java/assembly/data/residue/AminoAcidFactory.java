package assembly.data.residue;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class AminoAcidFactory {
    private static AminoAcidFactory instance = null;

    public static synchronized AminoAcidFactory getInstance() {
        if (instance == null) {
            instance = new AminoAcidFactory();
        }
        return instance;
    }

    private Map<String, AminoAcid> aaMap = null;

    private AminoAcidFactory() {
        aaMap = new HashMap<>();
        aaMap.put("A", new AminoAcid('A', 71.03711f));
        aaMap.put("R", new AminoAcid('R', 156.10111f));
        aaMap.put("N", new AminoAcid('N', 114.04293f));
        aaMap.put("D", new AminoAcid('D', 115.02694f));
        aaMap.put("C", new AminoAcid('C', 103.00919f));
        aaMap.put("E", new AminoAcid('E', 129.04259f));
        aaMap.put("Q", new AminoAcid('Q', 128.05858f));
        aaMap.put("G", new AminoAcid('G', 57.02146f));
        aaMap.put("H", new AminoAcid('H', 137.05891f));
        aaMap.put("I", new AminoAcid('I', 113.08406f));
        aaMap.put("L", new AminoAcid('L', 113.08406f));
        aaMap.put("K", new AminoAcid('K', 128.09496f));
        aaMap.put("M", new AminoAcid('M', 131.04049f));
        aaMap.put("F", new AminoAcid('F', 147.06841f));
        aaMap.put("P", new AminoAcid('P', 97.05276f));
        aaMap.put("S", new AminoAcid('S', 87.03203f));
        aaMap.put("T", new AminoAcid('T', 101.04768f));
        aaMap.put("W", new AminoAcid('W', 186.07931f));
        aaMap.put("Y", new AminoAcid('Y', 163.06333f));
        aaMap.put("V", new AminoAcid('V', 99.06841f));

        aaMap.put("C(+57.02)", new ModifiedAminoAcid('C', 103.00919f, 57.021464f));
        aaMap.put("C(+125.05)", new ModifiedAminoAcid('C', 103.00919f, 125.047676f));
        aaMap.put("M(+15.99)", new ModifiedAminoAcid('M', 131.04049f, 15.994915f));
        aaMap.put("N(+.98)", new ModifiedAminoAcid('N', 114.04293f, 0.984016f));
        aaMap.put("Q(+.98)", new ModifiedAminoAcid('Q', 128.05858f, 0.984016f));
    }

    public AminoAcid getAminoAcidByCode(String code) {
        return aaMap.get(code);
    }

    public Set<AminoAcid> getAminoAcidSet() {
        return new HashSet<>(aaMap.values());
    }
}
