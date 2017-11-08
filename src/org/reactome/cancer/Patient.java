package org.reactome.cancer;

import java.util.Objects;

public class Patient {
    private String tcgaBarcode;
    private String cancerType;
    public Patient(String tcgaBarcode,String cancerType){
        this.tcgaBarcode = tcgaBarcode;
        this.cancerType = cancerType;
    }

    public String getTcgaBarcode() {
        return tcgaBarcode;
    }

    public String getCancerType() {
        return cancerType;
    }

    @Override
    public int hashCode(){
        return Objects.hash(this.tcgaBarcode,this.cancerType);
    }

    @Override
    public String toString(){
        return this.tcgaBarcode;
    }
}
