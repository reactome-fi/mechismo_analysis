/*
 * Created on Jun 3, 2009
 *
 */
package org.reactome.cancer;

/**
 * This method is used to model gene expression value for one gene.
 * @author wgm
 *
 */
public class GeneExpressionValue {
    
    private String gene;
    private double tvalue;
    private double zvalue;
    private double pvalue;
    
    public GeneExpressionValue() {
    }

    public String getGene() {
        return gene;
    }

    public void setGene(String gene) {
        this.gene = gene;
    }

    public double getTvalue() {
        return tvalue;
    }

    public void setTvalue(double tvalue) {
        this.tvalue = tvalue;
    }

    public double getZvalue() {
        return zvalue;
    }

    public void setZvalue(double zvalue) {
        this.zvalue = zvalue;
    }

    public double getPvalue() {
        return pvalue;
    }

    public void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }
    
}
