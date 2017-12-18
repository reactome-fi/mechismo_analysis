package org.reactome.mechismo;

/**
 * Describe a specific residue in a structure.
 * @author wug
 *
 */
public class Residue {
    private int position;
    private String residue;
    private double mechismoScore;
    private ResidueStructure structure;
    private Entity protein;
    
    public Residue() {
    }

    public ResidueStructure getStructure() {
        return structure;
    }

    public void setStructure(ResidueStructure structure) {
        this.structure = structure;
    }

    public Entity getProtein() {
        return protein;
    }

    public void setProtein(Entity protein) {
        this.protein = protein;
    }

    public int getPosition() {
        return position;
    }

    public void setPosition(int position) {
        this.position = position;
    }

    public String getResidue() {
        return residue;
    }

    public void setResidue(String residue) {
        this.residue = residue;
    }

    public double getMechismoScore() {
        return mechismoScore;
    }

    public void setMechismoScore(double mechismoScore) {
        this.mechismoScore = mechismoScore;
    }

}
