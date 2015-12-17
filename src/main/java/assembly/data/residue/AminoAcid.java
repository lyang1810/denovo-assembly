package assembly.data.residue;

import java.util.Objects;

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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        AminoAcid aminoAcid = (AminoAcid) o;
        return code == aminoAcid.code &&
                Float.compare(aminoAcid.mass, mass) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(code, mass);
    }
}
