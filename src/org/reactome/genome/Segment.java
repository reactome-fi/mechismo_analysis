/*
 * Created on Jun 5, 2007
 *
 */
package org.reactome.genome;


/**
 * This class is used to describe a segment in a chromosome.
 * @author guanming
 *
 */
public class Segment {
    protected String name;
    protected String chrom;
    protected String strand;
    protected int start;
    protected int end;
    
    public Segment() {
    }
    
    public Segment(int start,
                   int end,
                   String chromosome) {
        this.start = start;
        this.end = end;
        this.chrom = chromosome;
    }
    
    public int getLength() {
        return end - start + 1;
    }
    
    public String getChromosome() {
        return chrom;
    }

    public void setChromosome(String chrom) {
        this.chrom = chrom;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }
    
    /** 
     * Check if this Seqment object is overlapped with a chromsome segment.
     * @param sStart
     * @param sEnd
     * @param chrom
     * @return
     */
    public boolean overlap(int sStart, int sEnd, String chrom) {
        // Should be in the same chromosome
        if (!chrom.equals(this.chrom))
            return false;
        // this object is inside of the provided segment
        if (sStart < this.start && sEnd > this.start) {
            return true;
        }
        // this object is at the left side
        else if (sStart >= this.start && sStart <= this.end) {
            return true;
        }
        // this object is at the right side
        else if (sEnd <= this.end && sEnd >= this.start) {
            return true;
        }
        return false;
    }
    
    public boolean overlap(Segment segment) {
        return overlap(segment.getStart(),
                       segment.getEnd(),
                       segment.getChromosome());
    }
    
    /**
     * Check if this Segment object is contained by a chromosome segment.
     * @param sStart
     * @param sEnd
     * @param chromosome
     * @return
     */
    public boolean containedBy(int sStart,
                               int sEnd,
                               String chromosome) {
        if (!chrom.equals(chromosome))
            return false;
        if (sStart <= this.start && sEnd >= this.end)
            return true;
        return false;
    }
    
    /**
     * Check the passed genomic coordinate is in this Segment. It is assumed
     * that this coordinate should be in the same chromosome.
     * @param genomicCoord
     * @return
     */
    public boolean contains(int genomicCoord) {
        return genomicCoord >= this.start && genomicCoord <= this.end;
    }
    
    public boolean containedBy(Segment segment) {
        return containedBy(segment.getStart(),
                           segment.getEnd(),
                           segment.getChromosome());
    }
    
    public String toString() {
        return name + " " + chrom + " " + start + " " + end;
    }
}
