package org.reactome.cancer;

import java.util.Objects;

public final class Gene implements Comparable<Gene> {
    private final String hgncName;
    private final String uniprotID;
    public Gene(String hgncName,
                String uniprotID){
        this.hgncName = hgncName.toUpperCase();
        this.uniprotID = uniprotID.toUpperCase();
    }

    public String getHgncName() {
        return hgncName;
    }

    public String getUniprotID() {
        return uniprotID;
    }

    @Override
    public int compareTo(Gene o) {
        return o.hgncName.compareTo(this.hgncName);
    }

    @Override
    public int hashCode(){
        return Objects.hash(this.hgncName,this.uniprotID);
    }

    @Override
    public boolean equals(Object o){
        return o == null
                ? false
                : o instanceof  Gene && o.hashCode() == this.hashCode();
    }

    @Override
    public String toString(){
        return String.format("%s",
                this.hgncName);
    }
}
