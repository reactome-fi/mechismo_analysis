package org.reactome.mechismo;

/**
 * Model cancer samples.
 * @author wug
 *
 */
public class Sample {
    private String source; // e.g. TCGA or COSMIC
    private String name; // For TCGA, this should be the barcode
    private CancerType cancerType;
    
    public Sample() {
    }
    
    public String getSource() {
        return source;
    }
    
    public void setSource(String source) {
        this.source = source;
    }
    
    public String getName() {
        return name;
    }
    
    public void setName(String name) {
        this.name = name;
    }
    
    public CancerType getCancerType() {
        return cancerType;
    }
    
    public void setCancerType(CancerType cancerType) {
        this.cancerType = cancerType;
    }
    
}
