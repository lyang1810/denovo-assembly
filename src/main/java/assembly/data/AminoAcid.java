package assembly.data;

/**
 * Amino Acid
 */
public class AminoAcid {
    protected char code;
    protected float mass; // monoisotopic mass

    public AminoAcid(char code, float mass) {
        this.code = code;
        this.mass = mass;
    }

    public String getCode() {
        return String.valueOf(code);
    }

    public float getMass() {
        return mass;
    }

    @Override
    public String toString() {
        return getCode();
    }
}
