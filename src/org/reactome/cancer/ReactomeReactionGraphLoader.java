package org.reactome.cancer;

import org.apache.log4j.Logger;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.reactome.r3.util.FileUtility;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class ReactomeReactionGraphLoader {
    private final static Logger logger = Logger.getLogger(MechismoOutputLoader.class);
    private String rxnFilter = null;
    private String reactomeReactionNetworkFilePath = null;
    private String mechismoSamples2SigReactionsFilePath = null;
    private DefaultDirectedGraph<Long, DefaultEdge> reactionGraph = null;
    private Set<Long> reactionSet = null;
    private Set<Long> sigReactionSet = null;

    public ReactomeReactionGraphLoader(String reactomeReactionNetworkFilePath){
        this(reactomeReactionNetworkFilePath,null,null);
    }

    public ReactomeReactionGraphLoader(String reactomeReactionNetworkFilePath,
                                       String mechismoSamples2SigReactionsFilePath,
                                       String rxnFilter) {
        this.reactomeReactionNetworkFilePath = reactomeReactionNetworkFilePath;
        this.mechismoSamples2SigReactionsFilePath = mechismoSamples2SigReactionsFilePath;
        this.rxnFilter = rxnFilter;
    }

    private void ParseMechismoSamples2SigReactionsFile() {
        if(this.mechismoSamples2SigReactionsFilePath != null) {
            FileUtility fileUtility = new FileUtility();
            this.sigReactionSet = new HashSet<>();
            try {
                fileUtility.setInput(this.mechismoSamples2SigReactionsFilePath);
                String line;
                String[] tokens;
                line = fileUtility.readLine();
                tokens = line.split("\t");
                for (int i = 2; i < tokens.length; i++) {
                    this.sigReactionSet.add(Long.parseLong(tokens[i]));
                }
                fileUtility.close();
            } catch (IOException ioe) {
                logger.error(String.format("Couldn't use %s, %s: %s",
                        this.reactomeReactionNetworkFilePath,
                        ioe.getMessage(),
                        Arrays.toString(ioe.getStackTrace())));
            }
        }
    }

    private void ParseReactomeReactionFile() {
        ParseMechismoSamples2SigReactionsFile();
        FileUtility fileUtility = new FileUtility();
        this.reactionGraph = new DefaultDirectedGraph<>(DefaultEdge.class);
        this.reactionSet = new HashSet<>();
        try {
            fileUtility.setInput(this.reactomeReactionNetworkFilePath);
            String line;
            String[] tokens;
            boolean includeRxns;
            while ((line = fileUtility.readLine()) != null) {
                tokens = line.split(" ");
                int rxn1DbIdIdx = 0;
                Long rxn1DbId = Long.parseLong(tokens[rxn1DbIdIdx]);
                int rxn2DbIdIdx = 2;
                Long rxn2DbId = Long.parseLong(tokens[rxn2DbIdIdx]);
                includeRxns = false;

                if(this.rxnFilter != null) {
                    if (rxnFilter.equals("Yes") &&
                            (this.sigReactionSet.contains(rxn1DbId) &&
                                    this.sigReactionSet.contains(rxn2DbId))) {
                        includeRxns = true;
                    } else if (rxnFilter.equals("Linker") &&
                            (this.sigReactionSet.contains(rxn1DbId) ||
                                    this.sigReactionSet.contains(rxn2DbId))) {
                        includeRxns = true;
                    } else if (rxnFilter.equals("No")) {
                        includeRxns = true;
                    }
                }else{
                    includeRxns = true;
                }

                if (includeRxns) {
                    this.reactionSet.add(rxn1DbId);
                    this.reactionSet.add(rxn2DbId);
                    if (!this.reactionGraph.containsVertex(rxn1DbId)) {
                        this.reactionGraph.addVertex(rxn1DbId);
                    }
                    if (!this.reactionGraph.containsVertex(rxn2DbId)) {
                        this.reactionGraph.addVertex(rxn2DbId);
                    }
                    if (!this.reactionGraph.containsEdge(rxn1DbId, rxn2DbId)) {
                        this.reactionGraph.addEdge(rxn1DbId, rxn2DbId);
                    }
                }
            }
            fileUtility.close();
        } catch (IOException ioe) {
            logger.error(String.format("Couldn't use %s, %s: %s",
                    this.reactomeReactionNetworkFilePath,
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    public DefaultDirectedGraph<Long, DefaultEdge> getReactionGraph() {
        if (this.reactionGraph == null) {
            this.ParseReactomeReactionFile();
        }
        return this.reactionGraph;
    }

    Set<Long> getReactionSet() {
        if (this.reactionSet == null) {
            this.ParseReactomeReactionFile();
        }
        return this.reactionSet;
    }
}
