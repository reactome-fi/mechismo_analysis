package org.reactome.structure.model;

public class PDBUniProtMatch {

    private String chainID;
    private String chainSequence;
    private String uniprot;
    private String gene;
    private int pdbStart;
    private int uniprotStart;
    private int offset; // Add this number to pdb cooridnate to get UniProt coordiate
    
    public PDBUniProtMatch() {
    }

    public void setChainID(String chainID) {
        this.chainID = chainID;
    }

    public void setChainSequence(String chainSequence) {
        this.chainSequence = chainSequence;
    }

    public void setUniprot(String uniprot) {
        this.uniprot = uniprot;
    }

    public void setGene(String gene) {
        this.gene = gene;
    }

    public void setPdbStart(int pdbStart) {
        this.pdbStart = pdbStart;
    }

    public void setUniprotStart(int uniprotStart) {
        this.uniprotStart = uniprotStart;
    }

    public void setOffset(int offset) {
        this.offset = offset;
    }

    public String getChainID() {
        return chainID;
    }

    public String getChainSequence(){
        return chainSequence;
    }

    public String getUniprot() {
        return uniprot;
    }

    public String getGene() {
        return gene;
    }

    public int getPdbStart() {
        return pdbStart;
    }

    public int getUniprotStart() {
        return uniprotStart;
    }

    public int getOffset() {
        return offset;
    }
    
}
