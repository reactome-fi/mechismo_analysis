package org.reactome.cancer;

import java.util.Arrays;
import java.util.Objects;

public class FI{
    private Gene[] genes;
    public FI(Gene gene1, Gene gene2){
        this.genes = new Gene[]{gene1,gene2};
        Arrays.sort(this.genes);
    }

    public Gene[] getGenes() {
        return genes;
    }

    @Override
    public int hashCode(){
        return Objects.hash(this.genes);
    }

    @Override
    public boolean equals(Object o){
        return o == null
                ? false
                : o.hashCode() == this.hashCode();
    }
}
