/*
 * Created on Jul 26, 2016
 *
 */
package org.reactome.genome;

/**
 * Used to model coding regions in a transcript.
 * @author gwu
 *
 */
public class CDS extends Segment {
    private int bpLength;
    
    /**
     * Default constructor.
     */
    public CDS() {
    }

//    public int getPeptideLength() {
//        return peptideLength;
//    }
//
//    public void setPeptideLength(int peptideLength) {
//        this.peptideLength = peptideLength;
//    }
    
    public void setBPLength(int length) {
//        System.out.println("BP Length: " + length);
//        if (length % 3 != 0)
//            throw new IllegalArgumentException("The passed length should be folds of 3! The current length is " + length);
//        peptideLength = length / 3 - 1; // Remove the stop codon
        this.bpLength = length;
    }

    @Override
    public int getLength() {
        return bpLength;
    }
    
}
