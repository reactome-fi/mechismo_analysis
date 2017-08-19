package org.reactome.structure.model;

import org.biojava.nbio.structure.Structure;

/**
 * A 3D structure to describe a protein/protein interaction.
 * @author wug
 *
 */
public class InteractionStructure {
    private Structure structure;
    private String pdbFileName;
    private ProteinChainInfo firstChain;
    private ProteinChainInfo secondChain;
    
    public InteractionStructure() {
    }

    public Structure getStructure() {
        return structure;
    }

    public void setStructure(Structure structure) {
        this.structure = structure;
    }

    public String getPdbFileName() {
        return pdbFileName;
    }

    public void setPdbFileName(String pdbFileName) {
        this.pdbFileName = pdbFileName;
    }

    public ProteinChainInfo getFirstChain() {
        return firstChain;
    }

    public void setFirstChain(ProteinChainInfo firstChain) {
        this.firstChain = firstChain;
    }

    public ProteinChainInfo getSecondChain() {
        return secondChain;
    }

    public void setSecondChain(ProteinChainInfo secondChain) {
        this.secondChain = secondChain;
    }

}
