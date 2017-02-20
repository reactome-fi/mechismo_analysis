/*
 * Created on Jun 5, 2007
 *
 */
package org.reactome.genome;

/**
 * This class is used to describe RefSeq information.
 * @author guanming
 *
 */
public class RefSeqInfo extends Segment {
    private String geneName;
    private Transcript transcript;
    
    public Transcript getTranscript() {
        return transcript;
    }

    public void setTranscript(Transcript transcript) {
        this.transcript = transcript;
    }

    public RefSeqInfo() {
    }
    
    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }    
    
}
