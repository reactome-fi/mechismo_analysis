package org.reactome.r3.util;

/**
 * Created by burkhart on 3/10/17.
 */
public class MutationObservation {

    private String gene;
    private Integer coordinate;
    private String mutation;
    private String mutationType;
    private String fathmmType;
    private String mutationZygosity;
    private String sample;

    public MutationObservation() {
    }

    public MutationObservation(String gene,Integer coordinate,String mutation,String sample){
        this.gene = gene;
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
        return gene;
    }

    public void setGene(String gene) {
        this.gene = gene;
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
        return gene + " " + mutation + " " + fathmmType + " " + sample + " " + mutationZygosity;
    }
}
