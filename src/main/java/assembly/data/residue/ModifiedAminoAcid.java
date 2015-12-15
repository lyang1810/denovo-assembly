package assembly.data.residue;

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
}
