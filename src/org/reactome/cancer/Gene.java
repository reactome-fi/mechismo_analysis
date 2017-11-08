package org.reactome.cancer;

import java.util.Objects;

public class Gene implements Comparable<Gene> {
    private String hgncName;
    private String uniprotID;
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
        return Objects.hash(this.hgncName,
                this.uniprotID);
    }
}
