package org.reactome.mechismo;

import java.util.Set;

/**
 * Describe a mutation in an protein entity.
 * @author wug
 *
 */
public class Mutation {
    private Set<Sample> samples;
    private Residue residue;
    private String variant;
    
    public Mutation() {
    }

    public Set<Sample> getSamples() {
        return samples;
    }

    public void setSamples(Set<Sample> samples) {
        this.samples = samples;
    }

    public Residue getResidue() {
        return residue;
    }

    public void setResidue(Residue residue) {
        this.residue = residue;
    }

    public String getVariant() {
        return variant;
    }

    public void setVariant(String variant) {
        this.variant = variant;
    }

}
