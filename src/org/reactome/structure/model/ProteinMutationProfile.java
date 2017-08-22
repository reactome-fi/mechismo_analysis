package org.reactome.structure.model;

import java.util.Collection;

import org.reactome.funcInt.Protein;

public class ProteinMutationProfile {
    
    private int interfaceMutationCount;
    private Protein protein;
    private ProteinChainInfo chainInfo;
    private Collection<ProteinMutation> mutations;
    // These two p-values are related to interface enrichment results
    private Double enrichmentPValue = -1.0d;
    private ResiduleMutationProfile minAAProfile;
    
    public ProteinMutationProfile() {
    }
    
    @Deprecated
    public ProteinMutationProfile(String geneName,
                                  int interfaceLength,
                                  int geneLength,
                                  int interfaceMutationCount,
                                  int geneMutationCount){
        this.interfaceMutationCount = interfaceMutationCount;
        setGeneName(geneName);
//        setInterfaceLength(interfaceLength);
//        setInteractionMuationCount(interactionMutationCount);
    }
    
    public Double getEnrichmentPValue() {
        return enrichmentPValue;
    }

    public void setEnrichmentPValue(Double enrichmentPValue) {
        this.enrichmentPValue = enrichmentPValue;
    }

    public ResiduleMutationProfile getMinAAProfile() {
        return minAAProfile;
    }
    
    public boolean isEmpty() {
        return enrichmentPValue < 0.0d;
    }

    /**
     * Set the mutation profile with the minimum p-value in the interaction interface.
     * @param profile
     */
    public void setMinAAProfile(ResiduleMutationProfile profile) {
        this.minAAProfile = profile;
    }

    public void setGeneName(String geneName) {
        if (protein != null)
            return;
        protein = new Protein();
        protein.setShortName(geneName);
    }
    
    public Protein getProtein() {
        return protein;
    }

    public void setProtein(Protein protein) {
        this.protein = protein;
    }

    public ProteinChainInfo getChainInfo() {
        return chainInfo;
    }

    public void setChainInfo(ProteinChainInfo chainInfo) {
        this.chainInfo = chainInfo;
    }

    public Collection<ProteinMutation> getMutations() {
        return mutations;
    }

    public void setMutations(Collection<ProteinMutation> mutations) {
        this.mutations = mutations;
    }

    public void setInterfaceMutationCount(int interfaceMutationCount) {
        this.interfaceMutationCount = interfaceMutationCount;
    }

    public String getGeneName(){
        if (protein != null)
            return protein.getShortName();
        return null;
    }
    
    public int getInterfaceLength() {
        if (chainInfo != null)
            return chainInfo.getInterfaceCoordinates().size();
        return -1;
    }
    
    public int getGeneLength(){
        if (protein != null)
            return protein.getSequence().length();
        return -1;
    }
    
    public int getInterfaceMutationCount(){
        return this.interfaceMutationCount;
    }
    
    public int getGeneMutationCount(){
        if (mutations != null)
            return mutations.size();
        return -1;
    }

    public boolean isValid(){
        return getInterfaceLength() > -1 && getGeneLength() > -1 &&
               getGeneLength() > getInterfaceLength();
    }
    
    public String toString(){
        return String.format("%s,%d,%d,%d,%d",
                                getGeneName(),
                                getInterfaceLength(),
                                getGeneLength(),
                                getInterfaceMutationCount(),
                                getGeneMutationCount());
    }

}
