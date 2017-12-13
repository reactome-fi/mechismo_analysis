package org.reactome.cancer;

import java.util.Objects;

public final class Mutation{
    private final Gene gene;
    private final Integer position;
    private final Character normalResidue;
    private final Character mutResidue;
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
                : o instanceof Mutation && o.hashCode() == this.hashCode();
    }

    @Override
    public String toString(){
        return String.format("gene=%s:position=%d:residue=%c:mutation=%c",
                this.gene,
                this.position,
                this.normalResidue,
                this.mutResidue);
    }
}
