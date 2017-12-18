package org.reactome.mechismo;

public class AnalysisResult {
    private double pvalue;
    private CancerType cancerType;
    
    public AnalysisResult() {
        
    }

    public double getPvalue() {
        return pvalue;
    }

    public void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }

    public CancerType getCancerType() {
        return cancerType;
    }

    public void setCancerType(CancerType cancerType) {
        this.cancerType = cancerType;
    }

}
