package org.reactome.mechismo;

/**
 * Type of cancers (e.g. BRCA, COADREAD, OV). This list most likely will be pulled
 * from TCGA or COSMIC.
 * @author wug
 *
 */
public class CancerType {

    private String abbreviation;
    private String name;
    
    public CancerType() {
    }

    public String getAbbreviation() {
        return abbreviation;
    }

    public void setAbbreviation(String abbreviation) {
        this.abbreviation = abbreviation;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
}
