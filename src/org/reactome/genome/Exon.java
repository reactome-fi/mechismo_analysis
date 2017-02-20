/*
 * Created on Jul 26, 2016
 *
 */
package org.reactome.genome;

/**
 * Used to describe an Exon.
 * @author gwu
 *
 */
public class Exon extends Segment {
    // coordinate in transcript
    private int txStart;
    private int txEnd;
    
    public Exon() {
    }
    
    public int getTxStart() {
        return txStart;
    }


    public void setTxStart(int txStart) {
        this.txStart = txStart;
    }


    public int getTxEnd() {
        return txEnd;
    }


    public void setTxEnd(int txEnd) {
        this.txEnd = txEnd;
    }

    /**
     * Map from a genomic position to a coordinate of the transcript this exon belong to.
     * @param genomicCoord
     * @return
     */
    public int mapToTranscriptCoordinate(int genomicCoord) {
        return txStart + genomicCoord - start;
    }

    /**
     * @param start
     * @param end
     * @param chromosome
     */
    public Exon(int start, int end, String chromosome) {
        super(start, end, chromosome);
    }
    
}
