package org.reactome.structure.model;

import java.util.Collection;

/**
 * Mutation profile occured in a single aminio acid residue.
 * @author wug
 *
 */
public class ResiduleMutationProfile {
    private Integer coordinate;
    private Collection<ProteinMutation> mutations;
    private Double pvalue; // Use null to avoid using 0.0d as the default
    
    public ResiduleMutationProfile() {
    }

    public int getCoordinate() {
        return coordinate;
    }

    public void setCoordinate(int coordinate) {
        this.coordinate = coordinate;
    }

    public Collection<ProteinMutation> getMutations() {
        return mutations;
    }

    public void setMutations(Collection<ProteinMutation> mutations) {
        this.mutations = mutations;
    }

    public Double getPvalue() {
        return pvalue;
    }

    public void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }
    
    public String toString() {
        return coordinate + ": " + pvalue;
    }

}
