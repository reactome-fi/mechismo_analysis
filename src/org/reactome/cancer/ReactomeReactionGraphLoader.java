package org.reactome.cancer;
import java.io.IOException;
import java.net.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.apache.log4j.Logger;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.reactome.r3.util.FileUtility;

public class ReactomeReactionGraphLoader {
    private final static Logger logger = Logger.getLogger(MechismoOutputLoader.class);
    private final int rxn1DbIdIdx = 0;
    private final int rxn2DbIdIdx = 2;
    private String reactomeReactionNetworkFilePath = null;
    private DefaultDirectedGraph<Long,DefaultEdge> reactionGraph = null;
    private Set<Long> reactionSet = null;

    public ReactomeReactionGraphLoader(String reactomeReactionNetworkFilePath){
       this.reactomeReactionNetworkFilePath = reactomeReactionNetworkFilePath;
    }

    public void ParseReactomeReactionFile(){
        FileUtility fileUtility = new FileUtility();
        this.reactionGraph = new DefaultDirectedGraph<>(DefaultEdge.class);
        this.reactionSet = new HashSet<>();
        try {
            fileUtility.setInput(this.reactomeReactionNetworkFilePath);
            String line = null;
            String[] tokens = null;
            while ((line = fileUtility.readLine()) != null) {
                tokens = line.split(" ");
                Long rxn1DbId = Long.parseLong(tokens[this.rxn1DbIdIdx]);
                Long rxn2DbId = Long.parseLong(tokens[this.rxn2DbIdIdx]);
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
            fileUtility.close();
        }catch(IOException ioe){
            logger.error(String.format("Couldn't use %s, %s: %s",
                    this.reactomeReactionNetworkFilePath.toString(),
                    ioe.getMessage(),
                    Arrays.toString(ioe.getStackTrace())));
        }
    }

    public DefaultDirectedGraph<Long,DefaultEdge> getReactionGraph(){
        if(this.reactionGraph == null){
            this.ParseReactomeReactionFile();
        }
        return this.reactionGraph;
    }

    public Set<Long> getReactionSet(){
        if(this.reactionSet == null){
            this.ParseReactomeReactionFile();
        }
        return this.reactionSet;
    }
}
