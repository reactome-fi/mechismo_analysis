package org.reactome.cancer;

import java.util.Objects;

public final class Patient {
    private final String tcgaBarcode;
    private final String cancerType;
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

    /* I don't think we should include cancertype here bc
    some patients are included more than once (e.g. COAD
    and COADREAD).
    */
    @Override
    public int hashCode(){
        return Objects.hash(this.tcgaBarcode);
    }

    @Override
    public boolean equals(Object o){
        return o == null
                ? false
                : o instanceof Patient && o.hashCode() == this.hashCode();
    }

    @Override
    public String toString(){
        return String.format("%s:%s",
                this.cancerType,
                this.tcgaBarcode);
    }
}
