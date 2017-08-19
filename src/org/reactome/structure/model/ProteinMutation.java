package org.reactome.structure.model;

import org.reactome.funcInt.Protein;

/**
 * Created by burkhart on 3/10/17.
 */
public class ProteinMutation {

    private Protein protein;
    
    private Integer coordinate;
    private String mutation;
    private String mutationType;
    private String fathmmType;
    private String mutationZygosity;
    private String sample;

    public ProteinMutation() {
    }

    public ProteinMutation(String gene,
                           Integer coordinate,
                           String mutation,
                           String sample){
        setGene(gene);
        this.coordinate = coordinate;
        this.mutation = mutation;
        this.sample = sample;
    }

    public String getSample() {
        return sample;
    }

    public void setSample(String sample) {
        this.sample = sample;
    }

    public String getMutationZygosity() {
        return mutationZygosity;
    }

    public void setMutationZygosity(String mutationZygosity) {
        this.mutationZygosity = mutationZygosity;
    }

    public String getGene() {
        if (protein != null)
            return protein.getShortName();
        return null;
    }

    public void setGene(String gene) {
        protein = new Protein();
        protein.setShortName(gene);
    }
    
    public void setProtein(Protein protein) {
        this.protein = protein;
    }
    
    public Protein getProtein() {
        return this.protein;
    }

    public Integer getCoordinate() {
        return coordinate;
    }

    public void setCoordinate(Integer coordinate) {
        this.coordinate = coordinate;
    }

    public String getMutation() {
        return mutation;
    }

    public void setMutation(String mutation) {
        this.mutation = mutation;
    }

    public String getMutationType() {
        return mutationType;
    }

    public void setMutationType(String mutationType) {
        this.mutationType = mutationType;
    }

    public String getFathmmType() {
        return fathmmType;
    }

    public void setFathmmType(String fathmmType) {
        this.fathmmType = fathmmType;
    }

    public String toString() {
        return protein.getShortName() + " " + mutation + " " + fathmmType + " " + sample + " " + mutationZygosity;
    }
}
