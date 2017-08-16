package org.reactome.structure.model;

public class ProteinMutationProfile {
    
    private String proteinName;
    private int interfaceLength;
    private int geneLength;
    private int interfaceMutationCount;
    private int geneMutationCount;
    
    public ProteinMutationProfile() {
    }
    
    public ProteinMutationProfile(String geneName,
                                  int interfaceLength,
                                  int geneLength,
                                  int interfaceMutationCount,
                                  int geneMutationCount){
        this.proteinName = geneName;
        this.interfaceLength = interfaceLength;
        this.geneLength = geneLength;
        this.interfaceMutationCount = interfaceMutationCount;
        this.geneMutationCount = geneMutationCount;
    }
    
    public String getGeneName(){
        return this.proteinName;
    }
    
    public int getInterfaceLength(){
        return this.interfaceLength;
    }
    
    public int getGeneLength(){
        return this.geneLength;
    }
    
    public int getInterfaceMutationCount(){
        return this.interfaceMutationCount;
    }
    
    public int getGeneMutationCount(){
        return this.geneMutationCount;
    }

    public boolean isValid(){
        return this.proteinName != null
                && !this.proteinName.equals("")
                && this.geneLength > 0
                && this.interfaceLength <= this.geneLength;
    }
    
    public String toString(){
        return String.format("%s,%d,%d,%d,%d",
                                this.proteinName,
                                this.interfaceLength,
                                this.geneLength,
                                this.interfaceMutationCount,
                                this.geneMutationCount);
    }

}
