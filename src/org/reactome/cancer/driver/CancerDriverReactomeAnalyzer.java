/*
 * Created on Aug 18, 2016
 *
 */
package org.reactome.cancer.driver;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.ReactomeReactionExpander;
import org.reactome.r3.ReactionMapGenerator;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.util.Configuration;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.FisherExact;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;

/**
 * Analyze cancer drivers distribution in Reactome.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class CancerDriverReactomeAnalyzer {
    // FDR cutoff to choose cancer driver reactions
    // The following value is based on results generated from method
    // searchForBestFDRCutoffForEnrichedReactions() (output was plotted
    // in R and then manually chosen)
    private final static double REACTION_FDR_CUTOFF = 0.017d;
    private FileUtility fu = new FileUtility();
    private MySQLAdaptor dba;

    /**
     * Default constructor.
     */
    public CancerDriverReactomeAnalyzer() {
    }
    
    public List<GKInstance> loadHumanReactions() throws Exception {
        MySQLAdaptor dba = getDBA();
        return new org.reactome.r3.ReactomeAnalyzer().loadHumanReactions(dba);
    }
    
    @Test
    public void dumpReactionIdToName() throws Exception {
        List<GKInstance> humanReactions = loadHumanReactions();
        System.out.println("DB_ID\tDisplay_Name");
        humanReactions.stream().forEach(reaction -> System.out.println(reaction.getDBID() + "\t" + reaction.getDisplayName()));
    }

    @Test
    public void testLoadReactionIdToFIsWithFeatures() throws Exception {
        Map<Long, Set<String>> idsToLines = loadReactionIdToFIsWithFeatures();
        System.out.println("Total reaction ids: " + idsToLines.size());
    }
    
    public Set<String> loadFIsWithPPIFeature(String fileName) throws IOException {
        Set<String> fis = new HashSet<String>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            // Column 4 is for PositiveFeature
            String[] tokens = line.split("\t");
            if (tokens[4].equals("true")) {
                fis.add(tokens[0] + "\t" + tokens[1]);
            }
        }
        fu.close();
        return fis;
    }
    
    /**
     * Load FIs with features for reactions from a pre-generated files.
     * @return
     * @throws Exception
     */
    public Map<Long, Set<String>> loadReactionIdToFIsWithFeatures() throws Exception {
        Map<Long, Set<String>> rxtIdToFIsWithFeatures = new HashMap<Long, Set<String>>();
        String fileName = "results/DriverGenes/Drivers_0816/FIsInSelectedReactions_FDR_05_092516.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // Reaction ids are at the end of the line
            String[] ids = tokens[tokens.length - 1].split(", ");
            for (String id : ids) {
                InteractionUtilities.addElementToSet(rxtIdToFIsWithFeatures,
                                                     new Long(id),
                                                     line);
            }
        }
        fu.close();
        return rxtIdToFIsWithFeatures;
    }
    
    @Test
    public void checkGenesInFIFile() throws IOException {
        String dirName = "results/DriverGenes/Drivers_0816/";
        String fileName = dirName + "FIsInSelectedForInteractome3d_091416.txt";
        Set<String> totalGenes = loadGenesInFIFile(fileName);
        fileName = dirName + "FIsInSelectedForInteractome3d_FDR_05_01_filtered_091416.txt";
        totalGenes.addAll(loadGenesInFIFile(fileName));
        Set<String> cancerGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        Set<String> shared = InteractionUtilities.getShared(cancerGenes, totalGenes);
        System.out.println("Total genes: " + totalGenes.size());
        System.out.println("Cancer genes: " + cancerGenes.size());
        System.out.println("Shared: " + shared.size());
    }

    private Set<String> loadGenesInFIFile(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> totalGenes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            totalGenes.add(tokens[2]);
            totalGenes.add(tokens[3]);
        }
        fu.close();
        return totalGenes;
    }
    
    /**
     * Choose reactions based on enrichment analyses saved in file, MergedReactionEnrichmentAnalysis.
     * @throws IOException
     */
    @Test
    public void chooseHitReactions() throws IOException {
        String dirName = "datasets/ICGC/2016_04/Drivers/";
        String fileName = dirName + "MergedReactionEnrichmentAnalysis_090116.txt";
        String outFileName = dirName + "SelectedHitReactions_090116.txt";
        
        double threshold = 0.01;
        
        fu.setInput(fileName);
        fu.setOutput(outFileName);
        fu.printLine("DB_ID\tReaction\tGeneEnrichment\tFIEnrichment");
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Double geneFDR = new Double(tokens[7]);
            Double fiFDR = 1.0d;
            if (tokens[13].length() > 0)
                fiFDR = new Double(tokens[13]);
            fu.printLine(tokens[0] + "\t" + 
                         tokens[1] + "\t" + 
                         (geneFDR > threshold ? "false" : "true") + "\t" + 
                         (fiFDR > threshold ? "false" : "true"));
        }
        fu.close();
    }
    
    @Test
    public void chooseReactionsForStructureAnalysis() throws IOException {
        String dirName = "results/DriverGenes/Drivers_0816/";
        String fileName = dirName + "MergedReactionEnrichmentAnalysis_090116.txt";
        
        double threshold = 0.05;
        fu.setInput(fileName);
        String line = fu.readLine();
        System.out.println(line);
        int total = 0;
        Set<String> totalGenes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Double geneFDR = new Double(tokens[7]);
            Double fiFDR = 1.0d;
            if (tokens[13].length() > 0)
                fiFDR = new Double(tokens[13]);
            Integer totalFIs = 0;
            if (tokens[8].length() > 0)
                totalFIs = new Integer(tokens[8]);
            // We have to have FIs for analysis
            if (totalFIs == 0)
                continue;
            // Choose either fiFDR or geneFDR for structural analysis
            if (fiFDR > threshold && geneFDR > threshold)
                continue;
            totalGenes.add(tokens[2]);
            totalGenes.add(tokens[3]);
            System.out.println(line);
            total ++;
        }
        fu.close();
        System.out.println("Total selected lines: " + total);
    }
    
    /**
     * Check a reaction subnetwork for hit reactions.
     * @throws IOException
     */
    @Test
    public void checkReactionComponents() throws IOException {
        String dirName = "datasets/ICGC/2016_04/Drivers/";
//        // This file has a bug
//        String fileName = dirName + "MergedReactionEnrichmentAnalysis_083016.txt";
        // Updated file after the bug was fixed
        String fileName = dirName + "MergedReactionEnrichmentAnalysis_090116.txt";
        double fdrCutoff = 0.01d;
        
        Set<String> reactionIds = new HashSet<String>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            double geneFDR = new Double(tokens[7]);
            
//            // Use this code to get FI FDR enriched reactions only
//            geneFDR = 1.0d;
            
            double fiFDR = 1.0d;
            if (tokens[13].length() > 0)
                fiFDR = new Double(tokens[13]);
            if (geneFDR <= fdrCutoff || fiFDR <= fdrCutoff)
                reactionIds.add(tokens[0]);
        }
        fu.close();
        System.out.println("Total selected reactions: " + reactionIds.size());
        
        // Fetch a sub-network containing the above reactionIds
        String reactionMapFile = "/Users/gwu/Documents/wgm/work/reactome/ReactionNetwork/ReactionNetworkWithoutUb_082916.txt";
        fu.setInput(reactionMapFile);
        int count = 0;
        Set<String> relatedReactions = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(" ");
            if (reactionIds.contains(tokens[0]) && reactionIds.contains(tokens[2])) {
//            if (reactionIds.contains(tokens[0]) || reactionIds.contains(tokens[2])) {
                System.out.println(line);
                count ++;
                relatedReactions.add(tokens[0]);
                relatedReactions.add(tokens[2]);
            }
        }
        fu.close();
        System.out.println("Total rection edges: " + count);
        relatedReactions.retainAll(reactionIds);
        System.out.println("Hit reactions: " + relatedReactions.size());
    }
    
    @Test
    public void computeCorrForEnrichmentAndNetworkFeature() throws Exception {
        String reactionMapFeatureFile = "/Users/gwu/Documents/wgm/work/reactome/ReactionNetwork/ReactionNetworkWithoutUb_082916_1_Node.csv";
        reactionMapFeatureFile = "/Users/gwu/Documents/wgm/work/reactome/ReactionNetwork/ReactionNetworkWithUb_082916_1_Node.csv";
        
        String enrichmentName = "datasets/ICGC/2016_04/Drivers/DriverFIReactionEnrichmentAnalysis_082316.txt";
        enrichmentName = "datasets/ICGC/2016_04/Drivers/DriverReactionEnrichmentAnalysis_083016.txt";
        Map<String, Double> rxtToEnrichment = new HashMap<String, Double>();
        FileUtility fu = new FileUtility();
        fu.setInput(enrichmentName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String token = tokens[tokens.length - 1];
            Double enrichment = new Double(token);
            if (enrichment > 0.05)
                continue;
            rxtToEnrichment.put(tokens[0], -Math.log10(enrichment));
        }
        fu.close();
        
        Map<String, String> dbIdToName = loadReactionDBIDToName();
        
        String[] featureNames = new String[] {
                "AverageShortestPathLength",
                "BetweennessCentrality",
                "ClosenessCentrality",
                "ClusteringCoefficient",  
                "Eccentricity",
                "EdgeCount",
                "Indegree", 
                "NeighborhoodConnectivity",
                "Outdegree",
//                "PartnerOfMultiEdgedNodePairs",
                "Stress"
        };
        System.out.println("NetworkNodeFeature\tCorrelation\tP-value\tData");
        for (String featureName : featureNames) {
            Map<String, Double> rxtToFeature = new ReactionMapGenerator().loadNodeFeature(reactionMapFeatureFile, 
                                                                                          featureName);
            
            List<Double> features = new ArrayList<Double>();
            List<Double> enrichments = new ArrayList<Double>();
            for (String dbId : dbIdToName.keySet()) {
                Double feature = rxtToFeature.get(dbId);
                if (feature == null)
                    continue;
                Double enrichment = rxtToEnrichment.get(dbIdToName.get(dbId));
                if (enrichment == null)
                    continue;
                features.add(feature);
                enrichments.add(enrichment);
            }
            PearsonsCorrelation correlation = MathUtilities.constructPearsonCorrelation(features,
                                                                                        enrichments);
            System.out.println(featureName + "\t" + 
                               correlation.getCorrelationMatrix().getEntry(0, 1) + "\t" + 
                               correlation.getCorrelationPValues().getEntry(0, 1) + "\t" +
                               features.size());
        }
    }
    
    @Test
    public void mergeTwoReactionEnrichmentFiles() throws Exception {
        Map<String, String> dbIdToName = loadReactionDBIDToName();
        Map<String, String> nameToDBID = new HashMap<String, String>();
        for (String dbId : dbIdToName.keySet())
            nameToDBID.put(dbIdToName.get(dbId), dbId);
        
        String reactionMapFeatureFile = "/Users/gwu/Documents/wgm/work/reactome/ReactionNetwork/ReactionNetworkWithoutUb_082916_1_Node.csv";
        Map<String, Double> dbIdToFeature = new ReactionMapGenerator().loadNodeFeature(reactionMapFeatureFile, 
                                                                                       "BetweennessCentrality");
        
        String dir = "datasets/ICGC/2016_04/Drivers/";
        String fileName = dir + "DriverFIReactionEnrichmentAnalysis_082316.txt";
        String fileName1 = dir + "DriverReactionEnrichmentAnalysis_083016.txt";
//        String output = dir + "MergedReactionEnrichmentAnalysis_083016.txt";
        String output = dir + "MergedReactionEnrichmentAnalysis_090116.txt";
        Map<String, String> reactionToLine = new HashMap<String, String>();
        fu.setInput(fileName);
        fu.setOutput(output);
        FileUtility fu1 = new FileUtility();
        fu1.setInput(fileName1);
        String line = fu.readLine();
        String line1 = fu1.readLine();
        String[] tokens = line.split("\t");
        int fiHeaderSize = tokens.length - 1; // Don't use reaction
        StringBuilder builder = new StringBuilder();
        builder.append("DB_ID\t").append(line1).append("\t");
        for (int i = 1; i < tokens.length; i++)
            builder.append(tokens[i]).append("\t");
        builder.append("EdgeBetweenness");
        fu.printLine(builder.toString());
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            reactionToLine.put(tokens[0], line);
        }
        while ((line1 = fu1.readLine()) != null) {
            tokens = line1.split("\t");
            line = reactionToLine.get(tokens[0]);
            String dbId = nameToDBID.get(tokens[0]);
            Double feature = dbIdToFeature.get(dbId);
            builder.setLength(0);
            builder.append(dbId).append("\t").append(line1).append("\t");
            if (line == null) {
                for (int i = 0; i < fiHeaderSize; i++)
                    builder.append("\t");
            }
            else {
                tokens = line.split("\t");
                for (int i = 1; i < tokens.length; i++)
                    builder.append(tokens[i]).append("\t");
            }
            builder.append(feature);
            fu.printLine(builder.toString());
        }
        fu.close();
        fu1.close();
    }
    
    @Test
    public void generateReactionNetworkForEnrichedReactions() throws IOException {
        Set<String> enrichedReactionIds = loadEnrichedReactionIds();
        ReactionMapGenerator networkGenerator = new ReactionMapGenerator();
        networkGenerator.generateSubNetwork(enrichedReactionIds);
    }
    
    @Test
    public void searchForBestFDRCutoffForEnrichedReactions() throws IOException {
        Set<String> network = new ReactionMapGenerator().loadSimpleNetwork();
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(network);
        System.out.println("Total reaction ids in network: " + ids.size());
        
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        
        StringBuilder output = new StringBuilder();
        output.append("FDR\tTotalSelected\tLargetsComponentSize\tPercent\n");
        for (double fdr = 0.005d; fdr <= 0.05d; fdr += 0.001d) {
            System.out.println("FDR: " + fdr);
            Set<String> selectedIds = loadEnrichedReactionIds(fdr);
            System.out.println("Total selected ids: " + selectedIds.size());
            selectedIds.retainAll(ids);
            System.out.println("\tIn the network: " + selectedIds.size());
            Set<String> selectedNetwork = InteractionUtilities.getFIs(selectedIds, network);
            List<Set<String>> components = graphAnalyzer.calculateGraphComponents(selectedNetwork);
            components.forEach(comp -> System.out.println(comp.size()));
            output.append(fdr + "\t" + 
                          selectedIds.size() + "\t" + 
                          components.get(0).size() + "\t" + 
                          (double)components.get(0).size() / selectedIds.size() + "\n");
        }
        System.out.println("\n" + output.toString());
    }
    
    @Test
    public void performNetworkComponentSizeTest() throws IOException {
        Set<String> network = new ReactionMapGenerator().loadSimpleNetwork();
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(network);
        System.out.println("Total reaction ids in network: " + ids.size());
        
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        
        Set<String> selectedIds = loadEnrichedReactionIds();
        System.out.println("Total selected ids: " + selectedIds.size());
        selectedIds.retainAll(ids);
        System.out.println("\tIn the network: " + selectedIds.size());
        
        Set<String> selectedNetwork = InteractionUtilities.getFIs(selectedIds, network);
        List<Set<String>> components = graphAnalyzer.calculateGraphComponents(selectedNetwork);
        components.forEach(comp -> System.out.println(comp.size()));
        
        System.out.println("\nRandom permutations:");
        int permutation = 10000;
        List<Integer> randomFirstSizes = new ArrayList<>();
        for (int i = 0; i < permutation; i++) {
            Set<String> randomIds = MathUtilities.randomSampling(ids, selectedIds.size());
            Set<String> randomSelectedNetwork = InteractionUtilities.getFIs(randomIds, network);
            List<Set<String>> randomComps = graphAnalyzer.calculateGraphComponents(randomSelectedNetwork);
            // Just use the first component
//            System.out.println(i + "\t" + randomComps.get(0).size());
            randomFirstSizes.add(randomComps.get(0).size());
        }
        randomFirstSizes.sort(Comparator.reverseOrder());
        for (int i = 0; i < randomFirstSizes.size(); i++)
            System.out.println(i + "\t" + randomFirstSizes.get(i));
    }
    
    private Set<String> loadEnrichedReactionIds() throws IOException {
        return loadEnrichedReactionIds(REACTION_FDR_CUTOFF);
    }
    
    private Set<String> loadEnrichedReactionIds(Double fdrCutoff) throws IOException {
        // Choose FDRs with expansion
        int index = 7;
        Map<String, Double> reactionIdToFDRViaExp = loadReactionIdToCancerGeneEnrichment(index);
        Set<String> selectedReactionIdsViaExp = new HashSet<>();
        reactionIdToFDRViaExp.forEach((id, fdr) -> {
            if (fdr <= fdrCutoff)
                selectedReactionIdsViaExp.add(id);
        });
        System.out.println("Selected with expansion: " + selectedReactionIdsViaExp.size());
        
        index = 11; // FDRs without expanding reactions
        Map<String, Double> reactionIdToFDR = loadReactionIdToCancerGeneEnrichment(index);
        Set<String> selectedReactionIds = new HashSet<>();
        reactionIdToFDR.forEach((id, fdr) -> {
            if (fdr <= fdrCutoff)
                selectedReactionIds.add(id);
        });
        System.out.println("Selected without expansion: " + selectedReactionIds.size());
        
        Set<String> shared = InteractionUtilities.getShared(selectedReactionIdsViaExp, 
                                                            selectedReactionIds);
        System.out.println("\tShared: " + shared.size());
        
        // Use reactions having single genes that are cancer driver genes too
        Set<String> oneHitReactions = loadReactionIdForOneHitGene();
        System.out.println("One hit gene reactions: " + oneHitReactions.size());
        
        Set<String> totalSelected = new HashSet<>();
        totalSelected.addAll(selectedReactionIds);
        totalSelected.addAll(selectedReactionIdsViaExp);
        totalSelected.addAll(oneHitReactions);
        System.out.println("Total selected reactions: " + totalSelected.size());
        
        System.out.println("All checked reactions: " + reactionIdToFDRViaExp.size());
        double percentage = (double) totalSelected.size() / reactionIdToFDRViaExp.size();
        System.out.println("\tPercentage: " + percentage * 100.0d);
        
        return totalSelected;
    }
    
    @Test
    public void checkEnrichedReactions() throws IOException {
        loadEnrichedReactionIds();
    }
    
    @Test
    public void checkPathwaysForEnrichedReactions() throws Exception {
        // Reaction to enrichment scores
//        Map<String, Double> reactionNameToScore = loadReactionToCancerGeneEnrichment();
//        double cutoff = 1.30d; // fdr <= 0.05
        
        // Choose FDRs with expansion
//        int index = 7;
//        index = 11; // FDRs without expanding reactions
//        Map<String, Double> reactionIdToFDR = loadReactionIdToCancerGeneEnrichment(index);
//        Set<String> selectedReactionIds = new HashSet<>();
//        
//        double fdrCutff = 0.05d;
//        
//        reactionIdToFDR.forEach((id, fdr) -> {
//            if (fdr <= fdrCutff)
//                selectedReactionIds.add(id);
//        });
//        
//        // Use to map from name to id since ids are used for reaction to pathway map
////        Map<String, String> reactionIdToName = loadReactionDBIDToName();
////        Map<String, String> reactionNameToId = new HashMap<>();
////        reactionIdToName.forEach((id, name) -> reactionNameToId.put(name, id));
////        Set<String> selectedReactionIds =  reactionNameToScore.keySet().stream().map(name -> reactionNameToId.get(name)).collect(Collectors.toSet());
//        
//        // Use reactions having single genes that are cancer driver genes too
//        Set<String> oneHitReactions = loadReactionIdForOneHitGene();
//        System.out.println("One hit gene reactions: " + oneHitReactions.size());
//        
//        System.out.println("Total selected reactions: " + selectedReactionIds.size());
//        selectedReactionIds.addAll(oneHitReactions);
        
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        AnnotationHelper helper = new AnnotationHelper();
        helper.setReactionIdToPathwayFile("resources/ReactomeReactionsToPathways_051017.txt");
        annotator.setAnnotationHelper(helper);
        annotator.setUseBenjaminiHochbergForFDR(true);
        
        Set<String> selectedReactionIds = loadEnrichedReactionIds();
        
        List<GeneSetAnnotation> results = annotator.annotateReactionsWithReactomePathways(selectedReactionIds);
        System.out.println("Pathway\tNumberInPathway\tRatioOfPathway\tHitNumber\tpValue\tFDR\tHitIds");
        for (GeneSetAnnotation annotation : results) {
            System.out.println(annotation.getTopic() + "\t" + 
                               annotation.getNumberInTopic() + "\t" + 
                               annotation.getRatioOfTopic() + "\t" + 
                               annotation.getHitNumber() + "\t" + 
                               annotation.getPValue() + "\t" + 
                               annotation.getFdr() + "\t" + 
                               annotation.getHitIds());
        }
        
    }
    
    /**
     * In this analysis, a reaction is expanded to multiple implementations if EntitySet is involved
     * in the reaction participants or complexes contained by participants (recursively).
     * @throws Exception
     */
    @Test
    public void checkCancerDrivesInReactionsViaExpand() throws Exception {
        // Load all reaction genes for use
        Set<String> allGenes = loadAllReactionGenes();
        System.out.println("Total reaction genes: " + allGenes.size());
        FisherExact fisher = new FisherExact(allGenes.size());
        // Load known cancer driver genes
        Set<String> driverGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        System.out.println("Total cancer driver genes: " + driverGenes.size());
        driverGenes.retainAll(allGenes);
        System.out.println("\tFilter to reaction genes: " + driverGenes.size());
        // Load all human reactions
        List<GKInstance> reactions = loadHumanReactions();
        System.out.println("Total human reactions: " + reactions.size());
        
        if (true)
            return;
        
        // helper
        ReactomeReactionExpander expander = new ReactomeReactionExpander();
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
         int c = 0;
        SummaryStatistics genesStat = new SummaryStatistics();
        SummaryStatistics sharedStat = new SummaryStatistics();
        List<Double> pvalues = new ArrayList<>();
        // Hold results to calculate FDRs
        List<String> lines = new ArrayList<>();
//        // The following reactions cannot be handled because the generated graph is too big
        // There is no need after introducing a limit to expanding
//        Set<Long> toBeExcluded = new HashSet<>();
//        toBeExcluded.add(983709L); // A BlackboxEvent that can be expanded to 586 PEs and 4664 Relations, failing to get all paths!!!
//        toBeExcluded.add(8868659L); // Too many proteins in one input
//        toBeExcluded.add(983707L); // [Reaction:983707] SYK autophosphorylates at the activated BCR...
//        toBeExcluded.add(8871193L); // [BlackBoxEvent:8871193] Dissociation of AAK1 and dephosphorylation of AP-2 mu2
//        toBeExcluded.add(156912L); // [Reaction:156912] Peptide transfer from P-site tRNA to the A-site tRNA...
//        toBeExcluded.add(2029476L); // [BlackBoxEvent:2029476] Role of myosins in phagosome formation
//        toBeExcluded.add(8868658L); // [Reaction:8868658] HSPA8-mediated ATP hydrolysis promotes vesicle uncoating
        for (GKInstance reaction : reactions) {
            System.out.println(c + ": " + reaction);
//            if (toBeExcluded.contains(reaction.getDBID()))
//                continue; 
            
            // Enrichment analysis using expanding reactions
            Set<Set<String>> setOfGenes = expander.extractGenesFromReaction(reaction);
            if (setOfGenes == null || setOfGenes.size() == 0)
                continue; // In case there is no input! (e.g. [BlackBoxEvent:191072] Synthesis of Cx43)
            pvalues.clear();
            genesStat.clear();
            sharedStat.clear();
            for (Set<String> genes : setOfGenes) {
                Set<String> shared = InteractionUtilities.getShared(genes, driverGenes);
                double pvalue = fisher.getRightTailedP(shared.size(),
                                                       genes.size() - shared.size(), 
                                                       driverGenes.size() - shared.size(), 
                                                       allGenes.size() - genes.size() - driverGenes.size() + shared.size());
                genesStat.addValue(genes.size());
                sharedStat.addValue(shared.size());
                pvalues.add(pvalue);
            }
            double combinedPValue = MathUtilities.combinePValuesWithFisherMethod(pvalues);
            double minPValue = pvalues.stream().min(Comparator.naturalOrder()).get();
            
            // Enrichment analysis using all invovled genes in the reaction.
            Set<String> allReactionGenes = reactomeAnalyzer.grepGenesFromReaction(reaction);
            Set<String> allShared = InteractionUtilities.getShared(allReactionGenes, driverGenes);
            double allPValue = fisher.getRightTailedP(allShared.size(),
                                                      allReactionGenes.size() - allShared.size(),
                                                      driverGenes.size() - allShared.size(), 
                                                      allGenes.size() - allReactionGenes.size() - driverGenes.size() + allShared.size());
                    
            String line = (reaction.getDBID() + "\t" +
                    reaction.getDisplayName() + "\t" + 
                    setOfGenes.size() + "\t" + 
                    (int)(genesStat.getMean() + 0.5d) + "\t" + // Force to be an integer after rounding
                    (int)(sharedStat.getMean() + 0.5d) + "\t" + 
                    minPValue + "\t" + 
                    combinedPValue + "\t" + 
                    allReactionGenes.size() + "\t" + 
                    allShared.size() + "\t" + 
                    allPValue);
            lines.add(line);
            c ++;
//            if (c == 10)
//                break;
        }
        
        Map<String, Double> lineToExpFDR = calculateFDRsForLines(lines, 5);
        Map<String, Double> lineToAllFDR = calculateFDRsForLines(lines, 9);
        
        // Output into a file
        String resultFileName = "results/CancerDriversReactionEnrichmentWithExpand_060917.txt";
        fu.setOutput(resultFileName);
        fu.printLine("DB_ID\tName\tTotalSets\tExp_TotalGenes\tExp_Shared\tMin_PValue\tFisher_PValue\tExp_FDR\tTotalGenes\tShared\tPValue\tFDR");
        StringBuilder builder = new StringBuilder();
        lineToExpFDR.forEach((line, expFDR) -> {
            String[] tokens = line.split("\t");
            builder.setLength(0);
            for (int i = 0; i < 7; i++)
                builder.append(tokens[i]).append("\t");
            builder.append(expFDR);
            for (int i = 7; i < tokens.length; i++)
                builder.append("\t").append(tokens[i]);
            builder.append("\t").append(lineToAllFDR.get(line));
            try {
                fu.printLine(builder.toString());
            }
            catch(IOException e) {
                e.printStackTrace();
            }
        });
        fu.close();
    }

    private Map<String, Double> calculateFDRsForLines(List<String> lines,
                                                      int pvalueIndex) {
        // Calculate FDRs for min-P 
        lines.sort((line1, line2) -> {
            String[] tokens1 = line1.split("\t");
            String[] tokens2 = line2.split("\t");
            // Use min-pvalues to calculate FDRs
            Double minP1 = new Double(tokens1[pvalueIndex]);
            Double minP2 = new Double(tokens2[pvalueIndex]);
            return minP1.compareTo(minP2);
        });
        // P values in the list should have been sorted already. Don't sort them again!
        List<Double> sortedPValue = lines.stream()
                                         .map(line -> new Double(line.split("\t")[pvalueIndex]))
                                         .collect(Collectors.toList());
        List<Double> fdrs = MathUtilities.calculateFDRWithBenjaminiHochberg(sortedPValue);   
        // Lines -> expanded fdrs
        Map<String, Double> lineToFDR = new HashMap<>();
        for (int i = 0; i < lines.size(); i++) {
            lineToFDR.put(lines.get(i), fdrs.get(i));
        }
        return lineToFDR;
    }
    
    private Set<String> loadAllReactionGenes() throws IOException {
        String fiFile = "resources/ReactomeGenesToReactions022717.txt";
        try (Stream<String> stream = Files.lines(Paths.get(fiFile))) {
            Set<String> genes = stream.map(line -> line.split("\t")[0]).collect(Collectors.toSet());
            return genes;
        }
    }
    
    @Test
    public void checkCancerDriversInReactions() throws Exception {
//        String dir = "../FINetworkBuild/results/2015/";
        Set<String> driverGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        System.out.println("Total cancer driver genes: " + driverGenes.size());
        
        // For IL1 pathway
//        String ilGenes = "/Users/gwu/git/PGM-IL1/PGM_IL1_workspace/results/GeneListInIL1.txt";
//        Set<String> driverGenes = fu.loadInteractions(ilGenes);
//        System.out.println("Total driver genes: " + driverGenes.size());
        
        //String fiFile = dir + "ReactomeFIsToReactions_082216.txt";
//        String fiFile = dir + "ReactomeGenesToReactions_082316.txt";
        
        Map<String, Integer> reactionToGeneCount = new HashMap<String, Integer>();
        Set<String> allGenes = new HashSet<String>();
        Map<String, Integer> reactionToCancerGeneCount = new HashMap<String, Integer>();
        Set<String> cancerGenes = new HashSet<String>();

        String fiFile = "resources/ReactomeGenesToReactions022717.txt";
        Files.lines(Paths.get(fiFile))
             .map(line -> line.split("\t"))
             .forEach(tokens -> {
                 if (driverGenes.contains(tokens[0])) {
                     addCount(reactionToCancerGeneCount, tokens[1]);
                     cancerGenes.add(tokens[0]);
                 }
                 addCount(reactionToGeneCount, tokens[1]);
                 allGenes.add(tokens[0]);
             });

        performEnrichmentAnalysis(reactionToGeneCount,
                                  allGenes,
                                  reactionToCancerGeneCount,
                                  cancerGenes,
                                  "Genes");
    }
    
    @Test
    public void checkCancerDriverFIReactions() throws Exception {
        String dir = "../FINetworkBuild/results/2015/";
        Set<String> driverGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        System.out.println("Total driver genes: " + driverGenes.size());
        
        String fiFile = dir + "ReactomeFIsToReactions_082216.txt";
//        String fiFile = dir + "ReactomeFIsToReactionsWithComplexes_082516.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fiFile);
        String line = null;
        Map<String, Integer> reactionToFICount = new HashMap<String, Integer>();
        Set<String> fis = new HashSet<String>();
        Map<String, Integer> reactionToCancerFICount = new HashMap<String, Integer>();
        Set<String> cancerFIs = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (driverGenes.contains(tokens[0]) && driverGenes.contains(tokens[1])) {
//                System.out.println(line);
                addCount(reactionToCancerFICount, tokens[3]);
                cancerFIs.add(tokens[0] + "\t" + tokens[1]);
            }
            addCount(reactionToFICount, tokens[3]);
            fis.add(tokens[0] + "\t" + tokens[1]);
        }
        fu.close();
        
        performEnrichmentAnalysis(reactionToFICount,
                                  fis,
                                  reactionToCancerFICount,
                                  cancerFIs,
                                  "FIs");
    }
    
    private Map<String, Double> loadReactionIdToCancerGeneEnrichment(int index) throws IOException {
        String fileName = "results/CancerDriversReactionEnrichmentWithExpand_060917.txt";
        try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
            Map<String, Double> idToFDR = stream.skip(1) // Skip the first header line
                                                .map(line -> line.split("\t"))
                                                .collect(Collectors.toMap(tokens -> tokens[0],
                                                                          tokens -> new Double(tokens[index])));
            return idToFDR;
        }
    }
    
    /**
     * Select reactions that have one gene only and this one gene is a cancer driver gene. Such a reaction
     * should be regarded as important cancer genes.
     * @return
     * @throws IOException
     */
    private Set<String> loadReactionIdForOneHitGene() throws IOException {
        String fileName = "results/CancerDriversReactionEnrichmentWithExpand_060917.txt";
        try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
            Set<String> selected = stream.skip(1)
                                         .map(line -> line.split("\t"))
                                         .filter(tokens -> tokens[8].equals("1") && tokens[9].equals("1"))
                                         .map(tokens -> tokens[0])
                                         .collect(Collectors.toSet());
            return selected;
        }
    }
    
    public Map<String, Double> loadReactionToCancerGeneEnrichment() throws IOException {
        String driverEnrichmentFile = "results/CancerDriversReactionEnrichment_052317.txt";
        Map<String, Double> reactionToDriverEnrichment = Files.lines(Paths.get(driverEnrichmentFile))
                .skip(1) // Skip the first header line
                .map(line -> line.split("\t"))
                .collect(Collectors.toMap(tokens -> tokens[0],
                                          tokens -> -Math.log10(new Double(tokens[6]))));
        System.out.println("Total reactions in cancer driver enrichment analysis: " + reactionToDriverEnrichment.size());
        return reactionToDriverEnrichment;
    }

    private void performEnrichmentAnalysis(Map<String, Integer> reactionToCount,
                                           Set<String> allEntities,
                                           Map<String, Integer> reactionToCancerCount,
                                           Set<String> cancerEntities,
                                           String type) {
        System.out.println("Total " + type + ": " + allEntities.size());
        System.out.println("Cancer " + type + ": " + cancerEntities.size());
        double allRatio = (double) cancerEntities.size() / allEntities.size();
        System.out.println("Ratio: " + allRatio);
        FisherExact fisher = new FisherExact(allEntities.size());
        System.out.println("\nReaction\tTotal" + type + "\tDriver" + type + "\tRatio\tpValue(binomial)\tpValue(Fisher)\tFDR");
        List<String> lines = new ArrayList<String>();
        for (String reaction : reactionToCount.keySet()) {
            Integer totalFIs = reactionToCount.get(reaction);
            Integer driverFIs = reactionToCancerCount.get(reaction);
            if (driverFIs == null)
                driverFIs = 0;
            double ratio = (double) driverFIs / totalFIs;
            double pvalue = MathUtilities.calculateBinomialPValue(allRatio, totalFIs, driverFIs);
            double fisherPvalue = fisher.getRightTailedP(driverFIs,
                                                         totalFIs - driverFIs,
                                                         cancerEntities.size() - driverFIs,
                                                         allEntities.size() - totalFIs - cancerEntities.size() + driverFIs);
            lines.add(reaction + "\t" + totalFIs + "\t" + 
                      driverFIs + "\t" + ratio + "\t" + 
                      pvalue + "\t" + fisherPvalue);
        }
        Collections.sort(lines, new Comparator<String>() {
            public int compare(String line1, String line2) {
                int index = line1.lastIndexOf("\t");
                Double pvalue1 = new Double(line1.substring(index + 1));
                index = line2.lastIndexOf("\t");
                Double pvalue2 = new Double(line2.substring(index + 1));
                return pvalue1.compareTo(pvalue2);
            }
        });
        List<Double> pvalues = new ArrayList<Double>();
        for (String line1 : lines) {
            int index = line1.lastIndexOf("\t");
            Double pvalue = new Double(line1.substring(index + 1));
            pvalues.add(pvalue);
        }
        List<Double> fdrs = MathUtilities.calculateFDRWithBenjaminiHochberg(pvalues);
        for (int i = 0; i < lines.size(); i++) {
            System.out.println(lines.get(i) + "\t" + fdrs.get(i));
        }
    }
    
    private void addCount(Map<String, Integer> keyToCount,
                          String key) {
        Integer count = keyToCount.get(key);
        if (count == null)
            keyToCount.put(key, 1);
        else
            keyToCount.put(key, ++count);
    }
    
    /**
     * Use another method checkCancerDriversInReactions().
     * @throws Exception
     */
    @Deprecated
    @Test
    public void performCancerGenesReactionEnrichment() throws Exception {
        // Set up annotator
        String dir = "../FINetworkBuild/results/2015/";
        AnnotationHelper helper = new AnnotationHelper();
        helper.setProteinNameToPathwayFile(dir + "ReactomeGenesToReactions_082316.txt");
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setAnnotationHelper(helper);
        annotator.setUseBenjaminiHochbergForFDR(true);
        
        Set<String> driverGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        System.out.println("Total driver genes: " + driverGenes.size());
        
        // Perform enrichment analysis
        List<GeneSetAnnotation> annotations = annotator.annotateGenesWithFDR(driverGenes, AnnotationType.Pathway);
        System.out.println("\nReaction\tHitGenes\tTotalReactionGenes\tRatio\tP-value\tFDR");
        for (GeneSetAnnotation annotation : annotations) {
            System.out.println(annotation.getTopic() + "\t" + 
                               annotation.getHitNumber() + "\t" + 
                               annotation.getNumberInTopic() + "\t" + 
                               annotation.getRatioOfTopic() + "\t" + 
                               annotation.getPValue() + "\t" + 
                               annotation.getFdr());
        }
    }
    
    public MySQLAdaptor getDBA() throws Exception {
        if (dba == null)
            dba = Configuration.getConfiguration().getReactomeDBA();
        return dba;
    }
    
    public void setDBA(MySQLAdaptor dba) {
        this.dba = dba;
    }
    
    @Test
    public void checkDriverDistributionInReactome() throws Exception {
        MySQLAdaptor dba = getDBA();
        Collection<GKInstance> ewases = dba.fetchInstanceByAttribute(ReactomeJavaConstants.EntityWithAccessionedSequence,
                                                                     ReactomeJavaConstants.dataSource,
                                                                     "IS NULL",
                                                                     null);
        System.out.println("Total Reactome EWASes: " + ewases.size());
        dba.loadInstanceAttributeValues(ewases, new String[]{ReactomeJavaConstants.referenceEntity,
                                                             ReactomeJavaConstants.species});
        Set<String> ewasGenes = new HashSet<String>();
        for (GKInstance ewas : ewases) {
            GKInstance species = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.species);
            if (!species.getDisplayName().equals("Homo sapiens"))
                continue;
            GKInstance refEntity = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            String gene = (String) refEntity.getAttributeValue(ReactomeJavaConstants.geneName);
            ewasGenes.add(gene);
        }
        System.out.println("Total genes: " + ewasGenes.size());
        
        // Check drivers in the whole EWAS gene set
//        Set<String> driverGenes = new FICancerDriverPredictor().loadCancerCensusGenes();
        Set<String> driverGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        
        System.out.println("Total driver genes: " + driverGenes.size());
        driverGenes.retainAll(ewasGenes);
        System.out.println("In Reactome: " + driverGenes.size());
        double ratio = (double) driverGenes.size() / ewasGenes.size();
        System.out.println("ratio: " + ratio);
        
        // Check Genes having Catalyst and Regulator roles
        System.out.println("\nChecking CatalystActivities:");
        Set<String> caGenes = fetchRegulationGenes(ReactomeJavaConstants.CatalystActivity,
                                                     ReactomeJavaConstants.physicalEntity,
                                                     dba);
        performBinomialTest(caGenes, driverGenes, ratio);

        // Check Genes having Regulation roles
        System.out.println("\nChecking Regulations:");
        Set<String> regulationGenes = fetchRegulationGenes(ReactomeJavaConstants.Regulation,
                                                           ReactomeJavaConstants.regulator,
                                                           dba);
        performBinomialTest(regulationGenes, driverGenes, ratio);
        
        // Shared genes
        System.out.println("\nGenes having both CA and Regulation roles:");
        Set<String> shared = InteractionUtilities.getShared(regulationGenes, caGenes);
        performBinomialTest(shared, driverGenes, ratio);
        
        // Genes having either CA or Regulation roles
        System.out.println("\nGenes having either CA or Regulation roles:");
        Set<String> bothGenes = new HashSet<String>(caGenes);
        bothGenes.addAll(regulationGenes);
        performBinomialTest(bothGenes, driverGenes, ratio);
        
        // No regulation role genes
        System.out.println("\nGenes having no regulation/ca roles:");
        ewasGenes.removeAll(caGenes);
        ewasGenes.removeAll(regulationGenes);
        performBinomialTest(ewasGenes, driverGenes, ratio);
    }
    
    private void performBinomialTest(Set<String> testGenes,
                                     Set<String> driverGenes,
                                     double ratio) {
        System.out.println("Total test genes: " + testGenes.size());
        Set<String> shared = InteractionUtilities.getShared(testGenes, driverGenes);
        System.out.println("\tare drivers: " + shared.size());
        double caRatio = (double) shared.size() / testGenes.size();
        System.out.println("\tratio: " + caRatio);
        double pvalue = MathUtilities.calculateBinomialPValue(ratio,
                                                              testGenes.size(),
                                                              shared.size());
        System.out.println("\tp-value from binomial test: " + pvalue);
    }

    private Set<String> fetchRegulationGenes(String clsName, String attName, MySQLAdaptor dba)
            throws Exception, InvalidAttributeException {
        // Check Genes having Catalyst and Regulator roles
        Collection<GKInstance> cas = dba.fetchInstanceByAttribute(clsName,
                                                                  ReactomeJavaConstants.dataSource,
                                                                  "IS NULL",
                                                                  null);
        dba.loadInstanceAttributeValues(cas, new String[]{attName});
        Set<String> caGenes = new HashSet<String>();
        for (GKInstance ca : cas) {
            GKInstance pe = (GKInstance) ca.getAttributeValue(attName);
            if (pe == null || !pe.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                continue; // Just in case
            Set<GKInstance> refEntities = InstanceUtilities.grepReferenceEntitiesForPE(pe);
            for (GKInstance refEntity : refEntities) {
                if (refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct)) {
                    String gene = (String) refEntity.getAttributeValue(ReactomeJavaConstants.geneName);
                    caGenes.add(gene);
                }
            }
        }
        return caGenes;
    }
    
    private Map<String, String> loadReactionDBIDToName() throws Exception {
        MySQLAdaptor dba = getDBA();
        // Load instances
        Collection<GKInstance> reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReactionlikeEvent,
                                                                        ReactomeJavaConstants.species,
                                                                        "=",
                                                                        48887L);
        Map<String, String> dbIdToName = new HashMap<String, String>();
        for (GKInstance rxt : reactions)
            dbIdToName.put(rxt.getDBID().toString(),
                           rxt.getDisplayName());
        return dbIdToName;
    }
    
}
