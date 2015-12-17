package assembly.data.residue;

import java.util.Objects;

public class ModifiedAminoAcid extends AminoAcid {
    protected float modMass;

    public ModifiedAminoAcid(char code, float mass, float modMass) {
        super(code, mass);
        this.modMass = modMass;
    }

    @Override
    public String getCode() {
        return String.format("%c(%+.2f)", code, modMass);
    }

    @Override
    public float getMass() {
        return mass + modMass;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;
        ModifiedAminoAcid that = (ModifiedAminoAcid) o;
        return Float.compare(that.modMass, modMass) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(super.hashCode(), modMass);
    }
}
