package org.reactome.cancer;

import java.util.Objects;

public class Mutation{
    private Gene gene;
    private Integer position;
    private Character normalResidue;
    private Character mutResidue;
    public Mutation(Gene gene,
                    Integer position,
                    Character normalResidue,
                    Character mutResidue){
     this.gene = gene;
     this.position = position;
     this.normalResidue = normalResidue;
     this.mutResidue = mutResidue;
    }

    public Gene getGene() {
        return gene;
    }

    public int getPosition() {
        return position;
    }

    public char getNormalResidue() {
        return normalResidue;
    }

    public char getMutResidue() {
        return mutResidue;
    }

    @Override
    public int hashCode(){
        return Objects.hash(this.gene,this.position,this.normalResidue,this.mutResidue);
    }

    @Override
    public boolean equals(Object o){
        return o == null
                ? false
                : o.hashCode() == this.hashCode();
    }
}