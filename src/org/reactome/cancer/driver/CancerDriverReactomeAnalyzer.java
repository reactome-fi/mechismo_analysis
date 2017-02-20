/*
 * Created on Aug 18, 2016
 *
 */
package org.reactome.cancer.driver;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
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
import org.reactome.r3.ReactionMapGenerator;
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
    private FileUtility fu = new FileUtility();
    
    /**
     * Default constructor.
     */
    public CancerDriverReactomeAnalyzer() {
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
    public void checkCancerDriversInReactions() throws Exception {
        String dir = "../FINetworkBuild/results/2015/";
//        Set<String> driverGenes = new CancerDriverAnalyzer().getDriverGenes(null);
        
        // For IL1 pathway
        String ilGenes = "/Users/gwu/git/PGM-IL1/PGM_IL1_workspace/results/GeneListInIL1.txt";
        Set<String> driverGenes = fu.loadInteractions(ilGenes);
        System.out.println("Total driver genes: " + driverGenes.size());
        
        //String fiFile = dir + "ReactomeFIsToReactions_082216.txt";
        String fiFile = dir + "ReactomeGenesToReactions_082316.txt";
        fu.setInput(fiFile);
        Map<String, Integer> reactionToGeneCount = new HashMap<String, Integer>();
        Set<String> allGenes = new HashSet<String>();
        Map<String, Integer> reactionToCancerGeneCount = new HashMap<String, Integer>();
        Set<String> cancerGenes = new HashSet<String>();
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (driverGenes.contains(tokens[0])) {
                addCount(reactionToCancerGeneCount, tokens[1]);
                cancerGenes.add(tokens[0]);
            }
            addCount(reactionToGeneCount, tokens[1]);
            allGenes.add(tokens[0]);
        }
        fu.close();
        
        performEnrichmentAnalysis(reactionToGeneCount,
                                  allGenes,
                                  reactionToCancerGeneCount,
                                  cancerGenes);
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
                                  cancerFIs);
    }

    private void performEnrichmentAnalysis(Map<String, Integer> reactionToCount,
                                           Set<String> allEntities,
                                           Map<String, Integer> reactionToCancerCount,
                                           Set<String> cancerEntities) {
        System.out.println("Total FIs: " + allEntities.size());
        System.out.println("Cancer FIs: " + cancerEntities.size());
        double allRatio = (double) cancerEntities.size() / allEntities.size();
        System.out.println("Ratio: " + allRatio);
        FisherExact fisher = new FisherExact(allEntities.size());
        System.out.println("\nReaction\tTotalFIs\tDriverFIs\tRatio\tpValue(binomial)\tpValue(Fisher)\tFDR");
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
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_55_plus_i",
                                            "root", 
                                            "macmysql01");
        return dba;
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
