package org.reactome.structure.model;

import java.util.Map;
import java.util.Set;

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
    
}
