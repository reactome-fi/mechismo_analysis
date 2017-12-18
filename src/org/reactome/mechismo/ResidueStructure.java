package org.reactome.mechismo;

/**
 * Describe the structure for one specific {@link Residue Residue}.
 * @author wug
 *
 */
public class ResidueStructure {
    private String pdb;
    private int pdbPosition;
    private String pdbChain;
    private String pdbResidue;
    
    public ResidueStructure() {
    }

    public String getPdb() {
        return pdb;
    }

    public void setPdb(String pdb) {
        this.pdb = pdb;
    }

    public int getPdbPosition() {
        return pdbPosition;
    }

    public void setPdbPosition(int pdbPosition) {
        this.pdbPosition = pdbPosition;
    }

    public String getPdbChain() {
        return pdbChain;
    }

    public void setPdbChain(String pdbChain) {
        this.pdbChain = pdbChain;
    }

    public String getPdbResidue() {
        return pdbResidue;
    }

    public void setPdbResidue(String pdbResidue) {
        this.pdbResidue = pdbResidue;
    }

}
