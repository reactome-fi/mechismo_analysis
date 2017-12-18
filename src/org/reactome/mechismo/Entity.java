package org.reactome.mechismo;

/**
 * Describe targets invovled in interactions. They are usually proteins
 * but may be chemicals.
 * @author wug
 *
 */
public class Entity {
    
    private String name; // For protein, this is gene symbol
    private String identifier; // For protein, this is UniProt id
    
    public Entity() {
    }
    
    public String getName() {
        return name;
    }
    public void setName(String name) {
        this.name = name;
    }
    public String getIdentifier() {
        return identifier;
    }
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }
    
    

}
