package org.reactome.structure.model;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.reactome.funcInt.Protein;

public class ReactionMutationProfile {
    private Long dbId;
    private String name;
    private Map<String, InteractionMutationProfile> ppiToProfile;
    private Set<String> extractFIs;
    
    public ReactionMutationProfile() {
    }

    public Set<String> getExtractFIs() {
        return extractFIs;
    }

    public void setExtractFIs(Set<String> extractFIs) {
        this.extractFIs = extractFIs;
    }

    public Long getDbId() {
        return dbId;
    }

    public void setDbId(Long dbId) {
        this.dbId = dbId;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Map<String, InteractionMutationProfile> getPpiToProfile() {
        return ppiToProfile;
    }

    public void setPpiToProfile(Map<String, InteractionMutationProfile> ppiToProfile) {
        this.ppiToProfile = ppiToProfile;
    }
    
    public void addPpiToProfile(String ppi,
                                InteractionMutationProfile profile) {
        if (ppiToProfile == null)
            ppiToProfile = new HashMap<>();
        ppiToProfile.put(ppi, profile);
    }
    
    public String exportResults() {
        if (ppiToProfile == null)
            return null;
        StringBuilder builder = new StringBuilder();
        ppiToProfile.forEach((ppi, profile) -> {
            builder.append(ppi).append(": ").append("\n");
            exportResults(profile.getFirstProteinProfile(), builder);
            builder.append("\n");
            exportResults(profile.getSecondProteinProfile(), builder);
            builder.append("\n");
        });
        return builder.toString();
    }
    
    private void exportResults(ProteinMutationProfile profile, 
                               StringBuilder builder) {
        if (profile.isEmpty())
            return;
        builder.append("\t");
        Protein protein = profile.getProtein();
        builder.append(protein.getPrimaryAccession()).append("\t");
        builder.append(protein.getShortName()).append("\t");
        builder.append(profile.getEnrichmentPValue()).append("\t");
        builder.append(profile.getMinAAProfile());
    }
    
}
