/*
 * Created on Jul 26, 2016
 *
 */
package org.reactome.genome;

import java.util.List;

/**
 * @author gwu
 *
 */
public class Transcript extends Segment {
    private List<Exon> exons;
    private CDS cds;
    
    /**
     * Default constructor
     */
    public Transcript() {
    }
    
    public int mapToCDSCoordinate(int genomicCoord) {
        // Find the exon first
        Exon foundExon = null;
        for (Exon exon : exons) {
            if (exon.contains(genomicCoord)) {
                foundExon = exon;
                break;
            }
        }
        if (foundExon == null)
            throw new IllegalArgumentException("Cannot find an exon for " + genomicCoord + " in " + getName());
        int txCoord = foundExon.mapToTranscriptCoordinate(genomicCoord);
        return txCoord - cds.getStart() + 1;
    }
    
    public void setCDSCoordinates(int start, int end) {
        if (exons == null || exons.size() == 0)
            throw new IllegalStateException("Set the exons first before setting the CDS coordinates.");
        // Need to calculate CDS
        cds = new CDS();
        int transcriptLength = getLength();
        // Find the first and the last CDS exons
        int firstIndex = -1;
        int lastIndex = -1;
        for (int i = 0; i < exons.size(); i++) {
            Exon exon = exons.get(i);
            if (exon.contains(start))
                firstIndex = i;
            if (exon.contains(end))
                lastIndex = i;
            if (firstIndex >= 0 && lastIndex >= 0)
                break;
        }
        int cdsBPLength = 0;
        for (int i = firstIndex; i <= lastIndex; i++) {
            Exon exon = exons.get(i);
            if (i == firstIndex)
                cdsBPLength += (exon.getEnd() - start + 1);
            else if (i == lastIndex)
                cdsBPLength += (end - exon.getStart() + 1);
            else
                cdsBPLength += exon.getLength();
        }
        cds.setBPLength(cdsBPLength);
        // CDS' coordinates are based on transcript
        cds.setStart(exons.get(firstIndex).mapToTranscriptCoordinate(start));
        cds.setEnd(exons.get(lastIndex).mapToTranscriptCoordinate(end));
    }
    
    public CDS getCds() {
        return cds;
    }

    public List<Exon> getExons() {
        return exons;
    }

    public void setExons(List<Exon> exons) {
        if (start == 0 && end == 0) {
            throw new IllegalStateException("Set the transcript coordinate first before setting the exons.");
        }
        this.exons = exons;
        calculateExonsTranscriptCoordinates();
    }
    
    private void calculateExonsTranscriptCoordinates() {
        int preLength = 0;
        for (int i = 0; i < exons.size(); i++) {
            Exon exon = exons.get(i);
            exon.setTxStart(preLength + 1); // Start coordinate is 1
            preLength += exon.getLength();
            exon.setTxEnd(preLength - 1); // Make sure to remove 1 so that this exon's length is right
        }
    }
    
    @Override
    public int getLength() {
        int length = 0;
        if (exons != null) {
            for (Exon exon : exons)
                length += exon.getLength();
        }
        return length;
    }
    
}
