package org.reactome.cancer;

import org.apache.log4j.Logger;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.Pseudograph;
import org.reactome.cancer.driver.CancerDriverReactomeAnalyzer;
import org.reactome.r3.ReactomeAnalyzer;
import org.reactome.r3.util.FileUtility;

import java.io.IOException;
import java.util.*;

public class FIGraphLoader {
    private final static Logger logger = Logger.getLogger(MechismoOutputLoader.class);
    private String fiNetworkFilePath = null;
    private Pseudograph<Gene, DefaultEdge> fiGraph = null;
    private Pseudograph<FI, DefaultEdge> edgeticFIGraph = null;
    private Map<String, String> geneNameToUniprot = null;

    public FIGraphLoader(String fiNetworkFilePath,
                         CancerDriverReactomeAnalyzer cancerDriverReactomeAnalyzer) throws Exception {
        this.fiNetworkFilePath = fiNetworkFilePath;
        this.geneNameToUniprot = new ReactomeAnalyzer().getGeneToUniprotMap(cancerDriverReactomeAnalyzer.getDBA());
    }

    private void ParseFIFile() {
        FileUtility fileUtility = new FileUtility();
        this.fiGraph = new Pseudograph<>(DefaultEdge.class);
        int gene1NameIdx = 0;
        int gene2NameIdx = 1;
        int annotationIdx = 2;
        try {
            fileUtility.setInput(this.fiNetworkFilePath);
            String line;
            String[] tokens;
            while ((line = fileUtility.readLine()) != null) {
                tokens = line.split("\t");
                if(tokens[annotationIdx].contains("complex")) {
                    String gene1Name = tokens[gene1NameIdx];
                    String gene2Name = tokens[gene2NameIdx];

                    if (gene1Name != null && !gene1Name.isEmpty() &&
                            gene2Name != null && !gene2Name.isEmpty()) {

                        String gene1Uniprot = this.geneNameToUniprot.get(gene1Name);
                        String gene2Uniprot = this.geneNameToUniprot.get(gene2Name);

                        if (gene1Uniprot != null && !gene1Uniprot.isEmpty() &&
                                gene2Uniprot != null && !gene2Uniprot.isEmpty()) {

                            Gene gene1 = new Gene(gene1Name, gene1Uniprot);
                            Gene gene2 = new Gene(gene2Name, gene2Uniprot);

                            if (!this.fiGraph.containsVertex(gene1)) {
                                this.fiGraph.addVertex(gene1);
                            }
                            if (!this.fiGraph.containsVertex(gene2)) {
                                this.fiGraph.addVertex(gene2);
                            }
                            if (!this.fiGraph.containsEdge(gene1, gene2)) {
                                this.fiGraph.addEdge(gene1, gene2);
                            }
                        }
                    }
                }
            }
            fileUtility.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    this.fiNetworkFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    private void ConvertToEdgeticGraph() {
        this.edgeticFIGraph = new Pseudograph<>(DefaultEdge.class);
        int counter = 0;
        for (Gene gene : this.fiGraph.vertexSet()) {
            Set<FI> fis = new HashSet<>();
            Set<Gene> genesInComplex = new HashSet<>();
            //find connected genes
            for (DefaultEdge geneConnection : this.fiGraph.edgesOf(gene)) {
                genesInComplex.add(this.fiGraph.getEdgeSource(geneConnection));
                genesInComplex.add(this.fiGraph.getEdgeTarget(geneConnection));
            }
            //create FIs for each connected gene
            for (Gene connectedGene : genesInComplex) {
                if (!gene.equals(connectedGene)) {
                    FI fi = new FI(gene, connectedGene);
                    fis.add(fi);
                    if (!this.edgeticFIGraph.containsVertex(fi)) {
                        this.edgeticFIGraph.addVertex(fi);
                    }
                }
            }
            //create FI connections
            List<FI> fisList = new ArrayList<>(fis);
            for (int i = 0; i < fisList.size() - 1; i++) {
                for (int j = i+1; j < fisList.size(); j++) {
                    FI fi1 = fisList.get(i);
                    FI fi2 = fisList.get(j);
                    if(!fi1.equals(fi2)) {
                        if (!this.edgeticFIGraph.containsEdge(fi1, fi2)) {
                            this.edgeticFIGraph.addEdge(fi1, fi2);
                        }
                    }
                }
            }
            counter++;
                System.out.println(String.format("Processed %d of %d genes...",
                        counter,this.fiGraph.vertexSet().size()));
        }
    }

    public Pseudograph<Gene, DefaultEdge> getFiGraph() {
        if (this.fiGraph == null) {
            this.ParseFIFile();
        }
        return this.fiGraph;
    }

    public Pseudograph<FI, DefaultEdge> getEdgeticFIGraph() {
        if (this.fiGraph == null) {
            this.ParseFIFile();
            this.ConvertToEdgeticGraph();
        }
        return this.edgeticFIGraph;
    }
}
