package org.reactome.structure.model;

import java.util.Set;

public class ProteinChainSummary {
    
    private Set<Integer> interfaceCoordinates;
    private int chainLength;
    private String proteinName;
    private String uniprotId;
    
    public ProteinChainSummary() {
    }
    
    public ProteinChainSummary(Set<Integer> interfaceCoordinates,
                               int chainLength) {
        if (chainLength <= 0)
            throw new IllegalArgumentException("Invalid chain length detected!");
        this.interfaceCoordinates = interfaceCoordinates;
        this.chainLength = chainLength;
    }

    public String getProteinName() {
        return proteinName;
    }

    public String getUniprotId() {
        return uniprotId;
    }

    public void setUniprotId(String uniprotId) {
        this.uniprotId = uniprotId;
    }

    public void setProteinName(String proteinName) {
        this.proteinName = proteinName;
    }

    public Set<Integer> getInterfaceCoordinates(){
        return this.interfaceCoordinates;
    }
    
    public int getChainLength(){
        return this.chainLength;
    }

    public void setInterfaceCoordinates(Set<Integer> interfaceCoordinates) {
        this.interfaceCoordinates = interfaceCoordinates;
    }

    public void setChainLength(int chainLength) {
        this.chainLength = chainLength;
    }

}
