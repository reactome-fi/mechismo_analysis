/*
 * Created on Apr 15, 2016
 *
 */
package org.reactome.r3;

import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.util.FileUtilities;
import org.junit.Test;
import org.reactome.fi.util.FileUtility;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.graph.NetworkBuilderForGeneSet;
import org.reactome.r3.util.Configuration;

/**
 * This test class is used to generate a network of reactions. Reaction1 and Reaction2 will be linked
 * together if reaction1 is annotated as a precedingEvent of Reaction2, or one of outputs, which is not
 * a SimpleEntity in the list of ATP, ADP, Pi, H2O, GTP, GDP, CO2, H+, , is an input, catalyst, regulator of Reaction2.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class ReactionMapGenerator {
    private static final Logger logger = Logger.getLogger(ReactionMapGenerator.class);
    
    private final String DIR_NAME = "resources/";
    private final String REACTION_NETWORK_NAME = DIR_NAME + "ReactionNetwork_070517.txt";
//    private final String REACTION_NETWORK_NAME = DIR_NAME + "ReactionNetwork_090120.txt";
//    private final String REACTION_NETWORK_NAME = DIR_NAME + "ReactionNetwork_Rel_71_122820.txt";
    private Set<String> entityEscapeNames;
    
    /**
     * Default constructor.
     */
    public ReactionMapGenerator() {
    }
    
    private MySQLAdaptor getDBA() throws Exception {
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "gk_central_041416",
//                                            "root",
//                                            "macmysql01");
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "reactome_55_plus_i",
//                                            "root",
//                                            "macmysql01");
        
        MySQLAdaptor dba = Configuration.getConfiguration().getReactomeDBA();
        
        return dba;
    }
    
    @Test
    public void checkMCLClusterResults() throws Exception {
        String file = DIR_NAME + "ReactionNewtork_MCL.txt";
        String outFile = DIR_NAME + "ReactionNetwork_MCL_Reaction.txt";
        FileUtilities fu = new FileUtilities();
        fu.setInput(file);
        fu.setOutput(outFile);
        fu.printLine("Reactome\tMCL_Cluster");
        String line = null;
        int cluster = 1;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            System.out.println(cluster + "\t" + tokens.length);
            for (String token : tokens)
                fu.printLine(token + "\t" + cluster);
            cluster ++;
        }
        fu.close();
    }
    
    /**
     * Load a network node feature for reaction nodes. The file should be generated by Cytoscape's
     * Analyze Network... feature and exported into a simple text file.
     * @param fileName A csv file exported from Cytoscape.
     * @param featureName
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadNodeFeature(String fileName,
                                               String featureName) throws IOException {
        Map<String, Double> reactionToFeature = new HashMap<String, Double>();
        Reader reader = new FileReader(fileName);
        Iterable<CSVRecord> records = CSVFormat.DEFAULT.parse(reader);
        int featureIndex = -1;
        int nameIndex = -1;
        boolean isFirst = true;
        for (CSVRecord record : records) {
            if (isFirst) {
                for (int i = 0; i < record.size(); i++) {
                    if (record.get(i).equalsIgnoreCase(featureName))
                        featureIndex = i;
                    else if (record.get(i).equalsIgnoreCase("name"))
                        nameIndex = i;
                }
                if (featureIndex == -1)
                    throw new IllegalArgumentException("Cannot find a column for " + featureName);
                if (nameIndex == -1)
                    throw new IllegalArgumentException("Cannot find a column for name!");
                isFirst = false;
            }
            else {
                reactionToFeature.put(record.get(nameIndex),
                                      new Double(record.get(featureIndex)));
            }
        }
        return reactionToFeature;
    }
    
    @Test
    public void simplifyNetwork() throws Exception {
        String source = DIR_NAME + "ReactionNetwork.txt";
        String target = DIR_NAME + "ReactionNetworkPair.txt";
        FileUtilities fu = new FileUtilities();
        fu.setInput(source);
        fu.setOutput(target);
        String line = null;
        fu.printLine("Rxt1\tRxt2");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(" ");
            fu.printLine(tokens[0] + "\t" + tokens[2]);
        }
        fu.close();
    }
    
    public void generateSubNetwork(Set<String> reactionIds,
                                   PrintStream ps) throws IOException {
        try (Stream<String> stream = Files.lines(Paths.get(REACTION_NETWORK_NAME))){
            stream.forEach(line -> {
                String[] tokens = line.split(" ");
                if (reactionIds.contains(tokens[0]) && reactionIds.contains(tokens[2])) {
                      ps.println(line);
                  }
            });
        }
    }
    
    public void generateSubNetwork(Set<String> reactionIds) throws IOException {
        generateSubNetwork(reactionIds, System.out);
    }
    
    /**
     * Use this method to connect all reactions together. Extract reactions that may be
     * used to connect all reactions together.
     * @param reactionIds
     * @throws IOException
     */
    public void generateSubNetworkForAll(Set<String> reactionIds,
                                         Map<String, Double> reactionToScore,
                                         PrintStream ps) throws IOException {
        try (Stream<String> stream = Files.lines(Paths.get(REACTION_NETWORK_NAME))) {
            // Convert the reaction fis as in the format for the FI network
            Set<String> edges = stream.map(line -> line.split(" "))
                                      .map(tokens -> tokens[0] + "\t" + tokens[2])
                                      .collect(Collectors.toSet());
            NetworkBuilderForGeneSet builder = new NetworkBuilderForGeneSet();
            builder.setAllFIs(edges);
            Set<String> reactionEdges = builder.constructFINetworkForGeneSet(reactionIds, reactionToScore);
            reactionEdges.forEach(edge -> {
                String[] tokens = edge.split("\t");
                if (edges.contains(edge))
                    ps.println(tokens[0] + " preceding " + tokens[1]);
                else
                    ps.println(tokens[1] + " preceding " + tokens[0]);
            });
        }
    }
    
    public void generateSubNetworkForAll(Set<String> reactionIds,
                                         Map<String, Double> reactionToScore) throws IOException {
        generateSubNetworkForAll(reactionIds, reactionToScore, System.out);
    }
    
    /**
     * Load the network without directions.
     * @return
     * @throws IOException
     */
    public Set<String> loadSimpleNetwork() throws IOException {
       return loadNetwork(REACTION_NETWORK_NAME," ", 2);
    }

    /**
     * Load the network without directions.
     * @return
     * @throws IOException
     */
    public Set<String> loadNetwork(String filePath, String delim, int pairIdx) throws IOException {
        try (Stream<String> stream = Files.lines(Paths.get(filePath))){
            Set<String> pairs = stream.map(line -> {
                String[] tokens = line.split(delim);
                // Do a sort
                int compare = tokens[0].compareTo(tokens[pairIdx]);
                if (compare < 0)
                    return tokens[0] + "\t" + tokens[pairIdx];
                else
                    return tokens[pairIdx] + "\t" + tokens[0];
            })
                    .collect(Collectors.toSet());
            return pairs;
        }
    }

    @Test
    public void analyzeGraphComponents() throws Exception {
        Set<String> network = loadSimpleNetwork();
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        List<Set<String>> components = graphAnalyzer.calculateGraphComponents(network);
        System.out.println("Total components: " + components.size());
        components.forEach(comp -> System.out.println(comp.size()));
    }
    
    public Set<String> loadLargestComponents() throws IOException {
        Set<String> network = loadSimpleNetwork();
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        List<Set<String>> components = graphAnalyzer.calculateGraphComponents(network);
        // Filter network so that both reactions are in the largest component
        Set<String> largest = components.get(0);
        Set<String> rtn = network.stream()
                                 .filter(pair -> largest.contains(pair.split("\t")[0]))
                                 .collect(Collectors.toSet());
        return rtn;
    }
    
    @Test
    public void generateReactionToPathwayMap() throws Exception {
        MySQLAdaptor dba = getDBA();
        // Get the top level pathways
        Collection<GKInstance> frontPages = dba.fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
        // There should be only one
        if (frontPages.size() != 1)
            throw new IllegalStateException("Too many FrontPage instances.");
        GKInstance frontPage = frontPages.stream().findAny().get();
        List<GKInstance> frontPageItems = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
        logger.info("Total frontPageItems: " + frontPageItems.size());
        // Exclude all pathways under the disease topic
        Map<GKInstance, Set<GKInstance>> pathwayToReactions = new HashMap<>();
        for (GKInstance topic : frontPageItems) {
            if (topic.getDisplayName().equals("Disease")) {
                logger.info("The Disease topic is excluded.");
                continue;
            }
            generatePathwayToReactionMap(topic, pathwayToReactions);
        }
        // Output
        String outFileName = DIR_NAME + "ReactionToPathway_Rel_71_122820.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(outFileName);
        fu.printLine("ReactionlikeEvent\tPathway");
        for (GKInstance pathway : pathwayToReactions.keySet()) {
            Set<GKInstance> reactions = pathwayToReactions.get(pathway);
            String pathwayId = getStableId(pathway);
            for (GKInstance reaction : reactions) {
                fu.printLine(getStableId(reaction) + "\t" + pathwayId);
            }
        }
        fu.close();
    }
    
    private void generatePathwayToReactionMap(GKInstance pathway,
                                             Map<GKInstance, Set<GKInstance>> pathwayToReactions) throws Exception {
        if (pathwayToReactions.containsKey(pathway))
            return;
        Set<GKInstance> contained = InstanceUtilities.getContainedEvents(pathway);
        pathwayToReactions.put(pathway,
                               contained.stream().filter(e -> e.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)).collect(Collectors.toSet()));
        for (GKInstance event : contained) {
            if (event.getSchemClass().isa(ReactomeJavaConstants.Pathway))
                generatePathwayToReactionMap(event, pathwayToReactions);
        }
    }
    
    @Test
    public void generate() throws Exception {
        MySQLAdaptor dba = getDBA();
        
        Collection<GKInstance> reactions = new ReactomeAnalyzer().loadHumanReactions(dba);
        Collection<GKInstance> regulations = dba.fetchInstancesByClass(ReactomeJavaConstants.Regulation);
        // Load attributes
        dba.loadInstanceAttributeValues(reactions, 
                                        new String[] {ReactomeJavaConstants.input,
                                                      ReactomeJavaConstants.output,
                                                      ReactomeJavaConstants.catalystActivity,
                                                      ReactomeJavaConstants.precedingEvent});
        dba.loadInstanceAttributeValues(regulations,
                                        new String[] {ReactomeJavaConstants.physicalEntity,
                                                      ReactomeJavaConstants.regulatedEntity});
        List<GKInstance> reactionList = new ArrayList<GKInstance>(reactions);
        
        logger.info("Total reactions: " + reactionList.size());
//        if (true)
//            return;
        
        FileUtilities fu = new FileUtilities();
//        fu.setOutput("tmp/ReactionNetwork.txt");
//        fu.setOutput(DIR_NAME + "ReactionNetwork_082916.txt");
//        fu.setOutput(DIR_NAME + "ReactionNetwork_101316.txt");
        fu.setOutput(REACTION_NETWORK_NAME);
        for (int i = 0; i < reactionList.size(); i++) {
            GKInstance reaction1 = reactionList.get(i);
            logger.info(i + ": " + reaction1);
            for (int j = 0; j < reactionList.size(); j++) { // Since it is possible two reactions may point to each other
                                                            // We need to check all
                GKInstance reaction2 = reactionList.get(j);
                if (i == j)
                    continue; // Escape itself
                boolean isPreceding = isPrecedingTo(reaction1, reaction2);
                if (isPreceding) {
//                    fu.printLine(reaction1.getDBID() + " preceding " + reaction2.getDBID());
                    fu.printLine(getStableId(reaction1) + " preceding " + getStableId(reaction2));
                }
            }
        }
        fu.close();
    }
    
    private String getStableId(GKInstance inst) throws Exception {
        GKInstance stableId = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.stableIdentifier);
        if (stableId == null)
            throw new IllegalStateException(inst + " doesn't have a stable id!");
        String rtn = (String) stableId.getAttributeValue(ReactomeJavaConstants.identifier);
        if (rtn == null)
            throw new IllegalStateException(stableId + " doesn't have the id value!");
        return rtn;
    }
    
    private boolean shouldEscape(GKInstance output) throws Exception {
//        if (output.getSchemClass().isa(ReactomeJavaConstants.SimpleEntity)) {
            String name = (String) output.getAttributeValue(ReactomeJavaConstants.name);
            return (getEntityEscapeNames().contains(name));
//        }
//        if (output.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
//            List<GKInstance> hasMember = output.getAttributeValuesList(ReactomeJavaConstants.hasMember);
//            for (GKInstance hasMember1 : hasMember) {
//                if (hasMember1.getSchemClass().isa(ReactomeJavaConstants.SimpleEntity))
//                    return true;
//            }
//        }
//        return false;
    }
    
    private Set<String> getEntityEscapeNames() {
        if (entityEscapeNames != null)
            return entityEscapeNames;
        entityEscapeNames = new HashSet<String>();
        String names = "ATP, ADP, Pi, H2O, GTP, GDP, CO2, H+";
//        // Ub-based protein degradation is applied to all proteins
//        // For this analysis, we escape it
//        String names = "ATP, ADP, Pi, H2O, GTP, GDP, CO2, H+, Ub";
        String[] tokens = names.split(", ");
        for (String token : tokens)
            entityEscapeNames.add(token);
        return entityEscapeNames;
    }
    
    /**
     * Check if reaction1 is preceding to reaction2.
     * @param rxt1
     * @param rxt2
     * @return
     * @throws Exception
     */
    private boolean isPrecedingTo(GKInstance rxt1,
                                  GKInstance rxt2) throws Exception {
        // If rxt1 is in the rxt2's precedingEvent list
        List<GKInstance> precedingEvent2 = rxt2.getAttributeValuesList(ReactomeJavaConstants.precedingEvent);
        if (precedingEvent2.contains(rxt1))
            return true;
        // Check if rxt1's output is in rxt2 input, catalyst, or regulator. 
        // For this test and simplicity, only non-simple molecule entity is tested
        List<GKInstance> output1 = rxt1.getAttributeValuesList(ReactomeJavaConstants.output);
        if (output1.size() == 0)
            return false;
        Set<GKInstance> lfhEntities = getLeftHandEntities(rxt2);
        for (GKInstance output : output1) {
            if (shouldEscape(output))
                continue;
            for (GKInstance lfhEntity : lfhEntities) {
                if (shouldEscape(lfhEntity))
                    continue;
                if (isEquivalent(lfhEntity, output))
                    return true;
            }
        }
        return false;
    }
    
    @Test
    public void testIsPrecedingTo() throws Exception {
        Long preId = 450592L;
        Long folId = 450499L;
        // These two should not
        preId = 5684273L;
        folId = 451634L;
        MySQLAdaptor dba = getDBA();
        GKInstance preRxt = dba.fetchInstance(preId);
        GKInstance folRxt = dba.fetchInstance(folId);
        System.out.println(isPrecedingTo(preRxt, folRxt));
    }
    
    private boolean isEquivalent(GKInstance lfhEntity,
                                 GKInstance output) throws Exception {
        // If they are the same, return true
        if (lfhEntity == output)
            return true;
        Set<GKInstance> lfhMembers = null;
        if (lfhEntity.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
            // If both are EntitySets having shared member, return true
            lfhMembers = InstanceUtilities.getContainedInstances(lfhEntity,
                                                                 ReactomeJavaConstants.hasMember);
            if (lfhMembers.contains(output))
                return true; // If the first reaction output is a member of the second reaction input set
        }
        if (output.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
            // If output is an EntitySet having lfhEntity as a member, return true
            Set<GKInstance> members = InstanceUtilities.getContainedInstances(output,
                                                                              ReactomeJavaConstants.hasMember);
            if (members.contains(lfhEntity))
                return true; // If some member of the output in the first reaction is input to the second reaction
            if (lfhMembers != null) {
                // If both are EntitySets having shared member, return true
                lfhMembers.retainAll(members);
                if (lfhMembers.size() > 0)
                    return true; // There is at least one shared member
            }
        }
        return false;
    }
    
    private Set<GKInstance> getLeftHandEntities(GKInstance reaction) throws Exception {
        Set<GKInstance> rtn = new HashSet<GKInstance>();
        List<GKInstance> input = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        if (input != null)
            rtn.addAll(input);
        GKInstance cas = (GKInstance) reaction.getAttributeValue(ReactomeJavaConstants.catalystActivity);
        if (cas != null) {
            GKInstance ca = (GKInstance) cas.getAttributeValue(ReactomeJavaConstants.physicalEntity);
            if (ca != null)
                rtn.add(ca);
        }
        Collection<GKInstance> regulations = InstanceUtilities.getRegulations(reaction);
        if (regulations != null && regulations.size() > 0) {
            for (GKInstance regulation : regulations) {
                GKInstance regulator = (GKInstance) regulation.getAttributeValue(ReactomeJavaConstants.regulator);
                if (regulator != null)
                    rtn.add(regulator);
            }
        }
        return rtn;
    }
}
