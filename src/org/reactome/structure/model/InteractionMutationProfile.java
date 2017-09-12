package org.reactome.structure.model;

public class InteractionMutationProfile{
    
    private ProteinMutationProfile profile1;
    private ProteinMutationProfile profile2;
    private Double p_value = null; // Make as negative to avoid 0.0
    private InteractionStructure structure;
    
    public InteractionMutationProfile() {
    }
    
    @Deprecated
    public InteractionMutationProfile(ProteinMutationProfile gene1MutationProfile,
                                      ProteinMutationProfile gene2MutationProfile,
                                      double p_value){
        this.profile1 = gene1MutationProfile;
        this.profile2 = gene2MutationProfile;
        this.p_value = p_value;
    }
    
    public InteractionStructure getStructure() {
        return structure;
    }

    public void setStructure(InteractionStructure structure) {
        this.structure = structure;
    }

    public ProteinMutationProfile getFirstProteinProfile(){
        return this.profile1;
    }
    
    public ProteinMutationProfile getSecondProteinProfile(){
        return this.profile2;
    }
    
    public double getP_value(){
        return this.p_value;
    }

    public void setFirstProteinProfile(ProteinMutationProfile gene1MutationProfile) {
        this.profile1 = gene1MutationProfile;
    }

    public void setSecondProteinProfile(ProteinMutationProfile gene2MutationProfile) {
        this.profile2 = gene2MutationProfile;
    }

    public void setP_value(double p_value) {
        this.p_value = p_value;
    }
}
