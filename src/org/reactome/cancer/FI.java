package org.reactome.cancer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

public final class FI implements Comparable<FI>{
    private final List<Gene> genes;
    public FI(Gene gene1, Gene gene2){
        this.genes = new ArrayList<>();
        this.genes.add(gene1);
        this.genes.add(gene2);
        Collections.sort(this.genes);
    }

    public List<Gene> getGenes() {
        return genes;
    }

    @Override
    public int hashCode(){
        return Objects.hash(this.genes.get(0),this.genes.get(1));
    }

    @Override
    public int compareTo(FI o) {
        return o.getGenes().get(0).getHgncName().compareTo(this.genes.get(0).getHgncName());
    }

    @Override
    public boolean equals(Object o){
        return o == null
                ? false
                : o instanceof  FI && o.hashCode() == this.hashCode();
    }

    @Override
    public String toString(){
        return String.format("%s-%s",
                genes.get(0).getHgncName(),
                genes.get(1).getHgncName());
    }

    public String toString(String delim){
        return String.format("%s%s%s",
                genes.get(0).getHgncName(),
                delim,
                genes.get(1).getHgncName());
    }

    public static String convertGeneNamePairToFIName(String delim, String geneName1, String geneName2){
        List<String> fiGenes = new ArrayList<>();
        fiGenes.add(geneName1.trim().toUpperCase());
        fiGenes.add(geneName2.trim().toUpperCase());
        Collections.sort(fiGenes);
        return String.join(delim,fiGenes);
    }
}
