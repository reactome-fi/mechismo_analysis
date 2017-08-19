package org.reactome.structure.model;

import java.util.Collection;
import java.util.Set;

import org.reactome.funcInt.Protein;

/**
 * Structure information for a protein involved in a specific interaction.
 * @author wug
 *
 */
public class ProteinChainInfo {
    
    private Collection<Integer> interfaceCoordinates;
    private int chainLength;
    private Protein protein;
    
    public ProteinChainInfo() {
    }
    
    public ProteinChainInfo(Set<Integer> interfaceCoordinates,
                               int chainLength) {
        this.interfaceCoordinates = interfaceCoordinates;
        this.chainLength = chainLength;
    }

    public Protein getProtein() {
        return protein;
    }

    public void setProtein(Protein protein) {
        this.protein = protein;
    }

    public String getProteinName() {
        if (protein != null)
            return protein.getShortName();
        return null;
    }

    public String getUniprotId() {
        if (protein != null)
            return protein.getPrimaryAccession();
        return null;
    }

    public void setUniprotId(String uniprotId) {
        if (protein == null)
            protein = new Protein();
        protein.setPrimaryAccession(uniprotId);
    }

    public void setProteinName(String proteinName) {
        if (protein == null)
            protein = new Protein();
        protein.setShortName(proteinName);
    }

    public Collection<Integer> getInterfaceCoordinates(){
        return this.interfaceCoordinates;
    }
    
    public int getChainLength(){
        return this.chainLength;
    }

    public void setInterfaceCoordinates(Collection<Integer> interfaceCoordinates) {
        this.interfaceCoordinates = interfaceCoordinates;
    }

    public void setChainLength(int chainLength) {
        this.chainLength = chainLength;
    }

}
