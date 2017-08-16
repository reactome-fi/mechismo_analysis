package org.reactome.structure.model;

public class InteractionMutationProfile{
    
    private ProteinMutationProfile gene1MutationProfile;
    private ProteinMutationProfile gene2MutationProfile;
    private double p_value;
    private String pdbFileName;
    
    public InteractionMutationProfile(ProteinMutationProfile gene1MutationProfile,
                                      ProteinMutationProfile gene2MutationProfile,
                                      double p_value){
        this.gene1MutationProfile = gene1MutationProfile;
        this.gene2MutationProfile = gene2MutationProfile;
        this.p_value = p_value;
    }
    
    public ProteinMutationProfile getGene1MutationProfile(){
        return this.gene1MutationProfile;
    }
    
    public ProteinMutationProfile getGene2MutationProfile(){
        return this.gene2MutationProfile;
    }
    
    public double getP_value(){
        return this.p_value;
    }

    public String getPdbFileName() {
        return pdbFileName;
    }

    public void setPdbFileName(String pdbFileName) {
        this.pdbFileName = pdbFileName;
    }

    public void setGene1MutationProfile(ProteinMutationProfile gene1MutationProfile) {
        this.gene1MutationProfile = gene1MutationProfile;
    }

    public void setGene2MutationProfile(ProteinMutationProfile gene2MutationProfile) {
        this.gene2MutationProfile = gene2MutationProfile;
    }

    public void setP_value(double p_value) {
        this.p_value = p_value;
    }
}
