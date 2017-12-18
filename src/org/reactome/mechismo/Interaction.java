package org.reactome.mechismo;

import java.util.Set;

/**
 * A very simple interaction really.
 * @author wug
 *
 */
public class Interaction {
    private String id; // A key for each search
    private Set<Entity> partners;
    private Set<Mutation> mutations;
    private Set<AnalysisResult> analysisResults;
    
    public Interaction() {
    }

    public Set<AnalysisResult> getAnalysisResults() {
        return analysisResults;
    }

    public void setAnalysisResults(Set<AnalysisResult> analysisResults) {
        this.analysisResults = analysisResults;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public Set<Entity> getPartners() {
        return partners;
    }

    public void setPartners(Set<Entity> partners) {
        this.partners = partners;
    }

    public Set<Mutation> getMutations() {
        return mutations;
    }

    public void setMutations(Set<Mutation> mutations) {
        this.mutations = mutations;
    }

}
